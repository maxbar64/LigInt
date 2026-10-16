#!/usr/bin/env python3
"""Test suite for ligint (chemical rules, geometry, robustness to malformed input).

Usage:
    python3 tests/run_tests.py                      build the sources and run everything
    python3 tests/run_tests.py --bin build/src/ligint   use an already built binary (e.g. from CMake)
    python3 tests/run_tests.py --sanitize           build with ASan+UBSan (detects crashes and undefined behavior)
    python3 tests/run_tests.py -k pi_              only the tests whose name contains "pi_"
    python3 tests/run_tests.py --list               list the tests
    python3 tests/run_tests.py --real-dir DIR       folder with the real PDB files for the "golden" tests (default: tests/data;
                                                    tests with missing files are skipped)

Groups:
    chem   chemical and geometric rules on synthetic fragments (known geometry, known expected outcome)
    tmpl   SDF / MOL2 templates (matching by coordinates, explicit H, fallback)
    pdb    PDB parsing (altlocs, insertion codes, models, ligand selection, CSV)
    fp     interaction fingerprints (CSV/JSON/SVG output, pocket definition files)
    cif    PDBx/mmCIF input (same results as the PDB, long chain ids, malformed files)
    manual the HTML manual (up to date, links resolve, agrees with the program)
    net    KLIFS pocket definitions and RCSB downloads (no real network traffic)
    matrix the interaction matrix of several ligands (--matrix)
    gpcr   pocket definitions of GPCRs from generic residue numbers (--gpcrdb-*)
    robust malformed / anomalous input: the program must not crash nor print nan
    golden real complexes: list of interactions already verified by hand against diagrams (REGRESSION
           tests, not validation tests)

Exit code: 0 if all the tests pass, 1 otherwise.
"""
import argparse
import glob
import html
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import re
import time
import xml.etree.ElementTree as ET
from html.parser import HTMLParser

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _build  # noqa: E402
TESTS = []          # (group, name, function)
CFG = {}            # binary, tmp, ...


# ----------------------------------------------------------------------------------------------
# minimal infrastructure
# ----------------------------------------------------------------------------------------------
class Fail(Exception):
    pass


class Skip(Exception):
    pass


def test(group):
    def deco(fn):
        TESTS.append((group, fn.__name__, fn))
        return fn
    return deco


def check(cond, msg):
    if not cond:
        raise Fail(msg)


def v_add(a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def v_sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def v_mul(a, s): return (a[0] * s, a[1] * s, a[2] * s)
def v_dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def v_norm(a): return math.sqrt(v_dot(a, a))
def v_unit(a):
    n = v_norm(a)
    return (a[0] / n, a[1] / n, a[2] / n)
def v_cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def at(origin, direction, dist):
    """point at 'dist' from 'origin' along 'direction'"""
    return v_add(origin, v_mul(v_unit(direction), dist))


# ----------------------------------------------------------------------------------------------
# synthetic molecules (coordinates in Angstrom, realistic geometry)
# ----------------------------------------------------------------------------------------------
class Mol:
    def __init__(self):
        self.atoms = []   # (name, element, xyz, mol2 type, formal charge)
        self.bonds = []   # (i, j, order)  0-based indices
        self.hs = []      # (heavy atom index, xyz)

    def add(self, name, el, xyz, typ=None, charge=0):
        self.atoms.append((name, el, tuple(xyz), typ or el, charge))
        return len(self.atoms) - 1

    def bond(self, i, j, order=1):
        self.bonds.append((i, j, order))

    def h(self, parent, xyz):
        self.hs.append((parent, tuple(xyz)))

    def pdb_lines(self, sysobj, res="LIG", chain="L", resi=1):
        for name, el, xyz, _, _ in self.atoms:
            sysobj.atom("HETATM", name, res, chain, resi, xyz, el)

    def sdf(self, with_h=True, shift=(0, 0, 0), name="mol"):
        heavy = self.atoms
        hs = self.hs if with_h else []
        lines = [name, "  ligint-test", "", "%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (len(heavy) + len(hs), len(self.bonds) + len(hs))]
        for _, el, xyz, _, q in heavy:
            p = v_add(xyz, shift)
            code = {0: 0, 1: 3, 2: 2, 3: 1, -1: 5, -2: 6, -3: 7}.get(q, 0)
            lines.append("%10.4f%10.4f%10.4f %-3s 0%3d  0  0  0  0  0  0  0  0" % (p[0], p[1], p[2], el, code))
        for _, xyz in hs:
            p = v_add(xyz, shift)
            lines.append("%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0" % (p[0], p[1], p[2], "H"))
        for i, j, o in self.bonds:
            lines.append("%3d%3d%3d  0" % (i + 1, j + 1, o))
        for k, (par, _) in enumerate(hs):
            lines.append("%3d%3d%3d  0" % (par + 1, len(heavy) + k + 1, 1))
        lines.append("M  END")
        lines.append("$$$$")
        return "\n".join(lines) + "\n"

    def mol2(self, name="mol"):
        n = len(self.atoms)
        out = ["@<TRIPOS>MOLECULE", name, "%d %d 0 0 0" % (n, len(self.bonds)), "SMALL", "NO_CHARGES", "",
               "@<TRIPOS>ATOM"]
        for k, (nm, el, xyz, typ, _) in enumerate(self.atoms, 1):
            out.append("%d %s%d %.4f %.4f %.4f %s 1 LIG 0.0000" % (k, el, k, xyz[0], xyz[1], xyz[2], typ))
        out.append("@<TRIPOS>BOND")
        for k, (i, j, o) in enumerate(self.bonds, 1):
            t = "ar" if o == 4 else str(o)
            out.append("%d %d %d %s" % (k, i + 1, j + 1, t))
        return "\n".join(out) + "\n"


def carboxylate():
    m = Mol()
    c1 = m.add("C1", "C", (-1.52, 0, 0), "C.3")
    c2 = m.add("C2", "C", (0, 0, 0), "C.2")
    o1 = m.add("O1", "O", (0.62, 1.07, 0), "O.co2")
    o2 = m.add("O2", "O", (0.62, -1.07, 0), "O.co2")
    m.bond(c1, c2); m.bond(c2, o1, 2); m.bond(c2, o2, 1)
    m.h(o2, (1.12, -1.90, 0))        # present only in the "acid" version
    return m


def amidine():
    m = Mol()
    c1 = m.add("C1", "C", (-1.52, 0, 0))
    c2 = m.add("C2", "C", (0, 0, 0))
    n1 = m.add("N1", "N", (0.67, 1.15, 0))
    n2 = m.add("N2", "N", (0.67, -1.15, 0))
    m.bond(c1, c2); m.bond(c2, n1, 2); m.bond(c2, n2, 1)
    return m


def guanidinium(types=False):
    m = Mol()
    c = m.add("C1", "C", (0, 0, 0), "C.cat" if types else "C")
    n1 = m.add("N1", "N", (0, 1.33, 0), "N.pl3" if types else "N")
    n2 = m.add("N2", "N", (1.152, -0.665, 0), "N.pl3" if types else "N")
    n3 = m.add("N3", "N", (-1.152, -0.665, 0), "N.pl3" if types else "N")
    for n in (n1, n2, n3):
        m.bond(c, n, 1)
    return m


def aminopyrimidine():
    """2-aminopyrimidine: C2 with three N neighbors but IN A RING (not a guanidinium)."""
    m = Mol()
    r = 1.35
    ang = {"C2": 90, "N3": 150, "C4": 210, "C5": 270, "C6": 330, "N1": 30}
    idx = {}
    for nm, a in ang.items():
        idx[nm] = m.add(nm, nm[0], (r * math.cos(math.radians(a)), r * math.sin(math.radians(a)), 0))
    ring = ["C2", "N3", "C4", "C5", "C6", "N1"]
    for k in range(6):
        m.bond(idx[ring[k]], idx[ring[(k + 1) % 6]], 1 + (k % 2))
    n7 = m.add("N7", "N", (0, 2 * r, 0))
    m.bond(idx["C2"], n7)
    return m


def trimethylamine():
    m = Mol()
    n = m.add("N1", "N", (0, 0, 0))
    a = math.radians(70.5)
    for k in range(3):
        phi = math.radians(120 * k)
        p = (1.47 * math.sin(a) * math.cos(phi), 1.47 * math.sin(a) * math.sin(phi), -1.47 * math.cos(a))
        c = m.add("C%d" % (k + 1), "C", p)
        m.bond(n, c)
    return m


def acetamide():
    m = Mol()
    c1 = m.add("C1", "C", (-1.52, 0, 0))
    c2 = m.add("C2", "C", (0, 0, 0))
    o = m.add("O1", "O", (0.60, 1.07, 0))
    n = m.add("N1", "N", (0.70, -1.14, 0))
    m.bond(c1, c2); m.bond(c2, o, 2); m.bond(c2, n, 1)
    return m


HDIR = (0.317, 0.948, 0.0)   # O-H direction of methanol (C-O-H angle = 108.5 degrees)


def methanol():
    m = Mol()
    c = m.add("C1", "C", (0, 0, 0))
    o = m.add("O1", "O", (1.43, 0, 0))
    m.bond(c, o)
    m.h(o, v_add((1.43, 0, 0), v_mul(HDIR, 0.96)))
    return m


def ethane():
    m = Mol()
    a = m.add("C1", "C", (0, 0, 0))
    b = m.add("C2", "C", (1.53, 0, 0))
    m.bond(a, b)
    return m


def chloromethane():
    m = Mol()
    c = m.add("C1", "C", (0, 0, 0))
    cl = m.add("CL1", "Cl", (1.78, 0, 0))
    m.bond(c, cl)
    return m


def ring_points(center, normal, radius, n=6, phase=0.0):
    nz = v_unit(normal)
    ref = (1, 0, 0) if abs(nz[0]) < 0.9 else (0, 1, 0)
    u = v_unit(v_cross(nz, ref))
    w = v_cross(nz, u)
    pts = []
    for k in range(n):
        a = math.radians(phase + 360.0 * k / n)
        pts.append(v_add(center, v_add(v_mul(u, radius * math.cos(a)), v_mul(w, radius * math.sin(a)))))
    return pts


def benzene(types=False):
    m = Mol()
    ids = []
    for k, p in enumerate(ring_points((0, 0, 0), (0, 0, 1), 1.40)):
        ids.append(m.add("C%d" % (k + 1), "C", p, "C.ar" if types else "C"))
    for k in range(6):
        m.bond(ids[k], ids[(k + 1) % 6], 4 if types else 1 + (k % 2))
    return m


PHE_RING = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]


# ----------------------------------------------------------------------------------------------
# construction of the PDB files and execution
# ----------------------------------------------------------------------------------------------
class Sys:
    def __init__(self):
        self.lines = []
        self.n = 0

    def atom(self, rec, name, res, chain, resi, xyz, el, alt=" ", icode=" ", occ=1.0, b=20.0):
        self.n += 1
        nm = name if (len(name) == 4 or len(el) == 2) else " " + name
        self.lines.append("%-6s%5d %-4s%s%3s %s%4d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s" %
                          (rec, self.n % 100000, nm, alt, res, chain, resi, icode, xyz[0], xyz[1], xyz[2], occ, b, el))

    def ligand(self, mol, res="LIG", chain="L", resi=1):
        mol.pdb_lines(self, res, chain, resi)
        return self

    def prot(self, res, name, xyz, resi=10, chain="A", el=None, **kw):
        self.atom("ATOM", name, res, chain, resi, xyz, el or name[0], **kw)
        return self

    def ring(self, res, center, normal, resi=20, chain="A", names=PHE_RING, radius=1.40, el=None):
        for nm, p in zip(names, ring_points(center, normal, radius, len(names))):
            self.atom("ATOM", nm, res, chain, resi, p, el or nm[0])
        return self

    def water(self, xyz, resi=500, chain="A"):
        self.atom("HETATM", "O", "HOH", chain, resi, xyz, "O")
        return self

    def metal(self, el, xyz, resi=300, chain="A"):
        self.atom("HETATM", el.upper(), el.upper(), chain, resi, xyz, el.capitalize())
        return self

    def text(self):
        return "\n".join(self.lines) + "\nEND\n"


class Res:
    def __init__(self, rc, out, err):
        self.rc, self.out, self.err = rc, out, err
        self.rows = []
        lines = [l for l in out.splitlines() if l.strip()]
        if lines and lines[0].startswith("ligand,type"):
            for l in lines[1:]:
                f = l.split(",", 6)
                if len(f) == 7:
                    self.rows.append({"lig": f[0], "kind": f[1], "fav": f[2], "la": f[3], "pa": f[4], "d": f[5], "note": f[6]})

    def find(self, kind=None, la=None, pa=None):
        return [r for r in self.rows if (kind is None or r["kind"] == kind) and (la is None or la in r["la"])
                and (pa is None or pa in r["pa"])]

    def has(self, kind, la=None, pa=None):
        return bool(self.find(kind, la, pa))

    def types(self):
        return sorted({r["kind"] for r in self.rows})


def run_raw(args, timeout=120, env=None):
    t0 = time.time()
    try:
        p = subprocess.run([CFG["bin"]] + args, capture_output=True, timeout=timeout, cwd=CFG["tmp"],
                           env=dict(os.environ, **env) if env else None)
    except subprocess.TimeoutExpired:
        raise Fail("timeout after %d s" % timeout)
    out = p.stdout.decode("utf-8", "replace")
    err = p.stderr.decode("utf-8", "replace")
    CFG["last_time"] = time.time() - t0
    return Res(p.returncode, out, err)


_counter = [0]


def tmpfile(content, suffix=".pdb", binary=False):
    _counter[0] += 1
    path = os.path.join(CFG["tmp"], "f%d%s" % (_counter[0], suffix))
    if binary:
        with open(path, "wb") as f:
            f.write(content)
    else:
        with open(path, "w", newline="", encoding="utf-8") as f:      # explicit: the default is the code page of the machine (cp1252 on Windows)
            f.write(content)
    return path


def go(sysobj, args=(), lig=("LIG",), sdf=None, mol2=None, csv=True, text=None, timeout=120):
    pdb = tmpfile(text if text is not None else sysobj.text())
    a = [pdb] + list(lig)
    if csv:
        a.append("--csv")
    if sdf:
        a += ["--sdf", tmpfile(sdf, ".sdf")]
    if mol2:
        a += ["--mol2", tmpfile(mol2, ".mol2")]
    a += list(args)
    r = run_raw(a, timeout)
    check(r.rc == 0, "exit code %d; stderr: %s" % (r.rc, r.err[:300]))
    return r


def no_crash(r):
    check("nan" not in r.out.lower().replace("annot", ""), "output contains 'nan'")
    check("AddressSanitizer" not in r.err and "runtime error" not in r.err, "sanitizer: " + r.err[:400])


# ----------------------------------------------------------------------------------------------
# GROUP chem: chemical rules
# ----------------------------------------------------------------------------------------------
@test("chem")
def carboxylate_forms_salt_bridge_with_lysine():
    s = Sys().ligand(carboxylate()).prot("LYS", "NZ", at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8))
    r = go(s)
    check(r.has("Salt bridge", la=":O1", pa="LYS"), "expected salt bridge carboxylate-lysine: " + str(r.types()))


@test("chem")
def carboxylate_repels_glutamate():
    s = Sys().ligand(carboxylate()).prot("GLU", "OE1", at((0.62, 1.07, 0), (0.5, 0.866, 0), 3.4))
    r = go(s)
    check(r.has("Electrostatic repulsion"), "expected carboxylate-carboxylate repulsion")
    check(not r.has("Salt bridge"), "salt bridge between like charges")


@test("chem")
def amidine_is_cationic():
    s = Sys().ligand(amidine()).prot("ASP", "OD1", at((0.67, 1.15, 0), (0.5, 0.866, 0), 2.9))
    r = go(s)
    check(r.has("Salt bridge", la=":N1"), "expected cationic amidine: " + str(r.types()))


@test("chem")
def guanidinium_is_cationic():
    s = Sys().ligand(guanidinium()).prot("GLU", "OE1", at((0, 1.33, 0), (0, 1, 0), 2.9))
    r = go(s)
    check(r.has("Salt bridge", la=":N1"), "expected cationic guanidinium: " + str(r.types()))


@test("chem")
def aminopyrimidine_is_not_cationic():
    """Regression: the 'N-C(N)-N' guanidinium in a ring (2-aminopyrimidine) is NOT a cation."""
    m = aminopyrimidine()
    n1 = m.atoms[5][2]                       # N1
    d = v_unit(n1)
    s = Sys().ligand(m).prot("LYS", "NZ", at(n1, d, 3.4))
    r = go(s)
    check(not r.has("Electrostatic repulsion"), "pyrimidine N treated as a cation (false guanidinium)")
    check(not r.has("Salt bridge"), "salt bridge with a neutral ring")


@test("chem")
def ring_nitrogen_is_acceptor():
    m = aminopyrimidine()
    n1 = m.atoms[5][2]
    s = Sys().ligand(m).prot("LYS", "NZ", at(n1, v_unit(n1), 2.9))
    r = go(s)
    check(r.has("Hydrogen bond", la=":N1"), "ring N must be able to accept an H-bond: " + str(r.types()))


@test("chem")
def aliphatic_amine_is_protonated():
    s = Sys().ligand(trimethylamine()).prot("ASP", "OD1", (0, 0, 2.9))
    r = go(s)
    check(r.has("Salt bridge", la=":N1"), "aliphatic amine assumed protonated: " + str(r.types()))


@test("chem")
def amide_nitrogen_is_not_acceptor():
    m = acetamide()
    n = m.atoms[3][2]
    dirn = (0.52, -0.85, 0)
    # donor (LYS NZ) near the amide N: the N is not an acceptor -> no H-bond
    r1 = go(Sys().ligand(m).prot("LYS", "NZ", at(n, dirn, 2.9)))
    check(not r1.has("Hydrogen bond", la=":N1"), "amide N treated as an acceptor")
    # acceptor (ASN OD1) near the N: the amide N-H is a donor -> H-bond
    r2 = go(Sys().ligand(m).prot("ASN", "OD1", at(n, dirn, 2.9)))
    check(r2.has("Hydrogen bond", la=":N1"), "amide N-H must be able to donate")


@test("chem")
def hydroxyl_donates_unless_metal_coordinated():
    m = methanol()
    o1 = m.atoms[1][2]
    glu = at(o1, (1, 0, 0), 2.7)
    r0 = go(Sys().ligand(m).prot("GLU", "OE1", glu))
    check(r0.has("Hydrogen bond", la=":O1"), "free hydroxyl: expected H-bond")
    mg = (o1[0], o1[1] + 2.1, 0)
    r1 = go(Sys().ligand(m).prot("GLU", "OE1", glu).metal("Mg", mg))
    check(r1.has("Metal coordination", la=":O1"), "expected O1-Mg")
    check(not r1.has("Hydrogen bond", la=":O1"), "O coordinated to Mg treated as a donor")
    r2 = go(Sys().ligand(m).prot("GLU", "OE1", glu).metal("Mg", mg), args=["--metal-donors"])
    check(r2.has("Hydrogen bond", la=":O1"), "--metal-donors should restore the donor")


@test("chem")
def metal_distance_limits_per_element():
    m = acetamide()
    o = m.atoms[2][2]
    d = (0.49, 0.87, 0)
    for el, dist, expect in [("Mg", 2.4, True), ("Mg", 2.7, False), ("Zn", 2.55, True), ("Zn", 2.75, False)]:
        r = go(Sys().ligand(m).metal(el, at(o, d, dist)))
        check(r.has("Metal coordination", la=":O1") == expect, "%s at %.2f A: expected %s" % (el, dist, expect))


@test("chem")
def metal_needs_lone_pair():
    m = acetamide()
    n = m.atoms[3][2]
    r = go(Sys().ligand(m).metal("Zn", at(n, (0.52, -0.85, 0), 2.3)))
    check(not r.has("Metal coordination", la=":N1"), "amide N (no lone pair) coordinated to Zn")


@test("chem")
def histidine_bound_to_zinc_does_not_donate():
    m = aminopyrimidine()
    n1 = m.atoms[5][2]
    d = v_unit(n1)
    ne2 = at(n1, d, 3.0)
    r0 = go(Sys().ligand(m).prot("HIS", "NE2", ne2))
    check(r0.has("Hydrogen bond", la=":N1", pa="HIS"), "control: without Zn the H-bond with His must be there")
    r1 = go(Sys().ligand(m).prot("HIS", "NE2", ne2).metal("Zn", at(ne2, d, 2.1)))
    check(not r1.has("Hydrogen bond", la=":N1", pa="HIS"), "His coordinated to Zn treated as a donor")


@test("chem")
def hydrophobic_rules():
    leu = (0.765, 0, 3.72)
    r = go(Sys().ligand(ethane()).prot("LEU", "CD1", leu))
    check(r.has("Hydrophobic contact"), "expected C-C hydrophobic contact")
    r = go(Sys().ligand(methanol()).prot("LEU", "CD1", (-2.0, 0, 3.23)))
    check(not r.has("Hydrophobic contact", la=":C1"), "C bonded to O treated as hydrophobic")
    r = go(Sys().ligand(ethane()).prot("GLU", "CG", leu))
    check(not r.has("Hydrophobic contact"), "GLU CG must not be hydrophobic")
    r = go(Sys().ligand(ethane()).prot("ALA", "CA", leu))
    check(not r.has("Hydrophobic contact") and not r.has("Apolar contact (backbone)"), "CA does not count by default")
    r = go(Sys().ligand(ethane()).prot("ALA", "CA", leu), args=["--backbone-hyd"])
    check(r.has("Apolar contact (backbone)"), "--backbone-hyd must count CA")


@test("chem")
def hydrophobic_thresholds_are_configurable():
    s = lambda: Sys().ligand(ethane()).prot("LEU", "CD1", at((0.765, 0, 0), (0, 0, 1), 4.2))
    check(not go(s()).has("Hydrophobic contact"), "4.2 A outside the default threshold")
    check(go(s(), args=["--hyd-max", "4.5"]).has("Hydrophobic contact"), "--hyd-max 4.5 must include it")


@test("chem")
def pi_stacking_parallel_and_offset_limit():
    for ox, expect in [(1.0, True), (3.0, True), (3.7, False), (6.0, False)]:
        s = Sys().ligand(benzene()).ring("PHE", (ox, 0, 3.6), (0, 0, 1))
        r = go(s)
        found = r.has("pi-stacking (parallel)")
        check(found == expect, "offset %.1f: expected %s, found %s" % (ox, expect, found))


@test("chem")
def pi_stacking_angle_threshold():
    n = (0, -math.sin(math.radians(36)), math.cos(math.radians(36)))
    mk = lambda: Sys().ligand(benzene()).ring("PHE", (0.5, 0, 3.7), n)
    check(not go(mk()).has("pi-stacking (parallel)"), "36 degrees beyond the 30 limit")
    check(go(mk(), args=["--pi-angle", "40"]).has("pi-stacking (parallel)"), "--pi-angle 40 must accept it")


@test("chem")
def pi_stacking_t_shaped():
    s = Sys().ligand(benzene()).ring("PHE", (0, 0, 4.8), (0, 1, 0))
    r = go(s)
    check(r.has("pi-stacking (T-shaped)"), "expected T-shaped stacking: " + str(r.types()))


@test("chem")
def cation_pi_geometry():
    r = go(Sys().ligand(benzene()).prot("LYS", "NZ", (0, 0, 4.0)))
    check(r.has("Cation-pi"), "LYS on the ring axis: expected cation-pi")
    r = go(Sys().ligand(benzene()).prot("LYS", "NZ", (3.46, 0, 2.0)))
    check(not r.has("Cation-pi"), "LYS at 60 degrees from the axis: no cation-pi")
    r = go(Sys().ligand(trimethylamine()).ring("PHE", (0, 0, 4.2), (0, 0, 1)))
    check(r.has("Cation-pi", la=":N1"), "ligand cation above the protein ring: " + str(r.types()))


@test("chem")
def halogen_bond_needs_linear_geometry():
    cl = (1.78, 0, 0)
    r = go(Sys().ligand(chloromethane()).prot("ASN", "OD1", at(cl, (1, 0, 0), 3.2)))
    check(r.has("Halogen bond"), "linear C-Cl...O: expected halogen bond")
    r = go(Sys().ligand(chloromethane()).prot("ASN", "OD1", at(cl, (0, 1, 0), 3.2)))
    check(not r.has("Halogen bond"), "C-Cl...O at 90 degrees: no halogen bond")


@test("chem")
def clash_and_covalent_distances():
    r = go(Sys().ligand(ethane()).prot("LEU", "CD1", (0, 0, 2.5)))
    check(r.has("Steric clash"), "2.5 A between two C: expected clash")
    r = go(Sys().ligand(ethane()).prot("LEU", "CD1", (-1.6, 0, 0)))      # 1.6 A from C1 and 3.1 A from C2
    check(r.has("Possible covalent bond", la=":C1") and not r.has("Steric clash"),
          "1.6 A: expected possible covalent bond: " + str(r.types()))


@test("chem")
def water_bridge():
    m = methanol()
    o1 = m.atoms[1][2]
    w = at(o1, (1, 0, 0), 2.8)
    p = at(w, (0.5, 0.866, 0), 2.8)
    mk = lambda: Sys().ligand(m).water(w).prot("ASN", "ND2", p)
    check(go(mk()).has("Water bridge"), "expected water bridge")
    check(not go(mk(), args=["--no-water"]).has("Water bridge"), "--no-water ignored")
    # direct contact: no bridge
    s = Sys().ligand(m).water(w).prot("ASN", "ND2", at(o1, (1, 0.3, 0), 3.0))
    check(not go(s).has("Water bridge"), "water bridge between atoms already in direct contact")


@test("chem")
def ligand_metal_coordination_removes_hbond_partners():
    """Regression (heme): ligand N and protein N coordinated to the same ligand Fe do not make H-bonds."""
    m = Mol()                                   # control: same ligand WITHOUT Fe
    n1 = m.add("N1", "N", (2.0, 0, 0))
    c1 = m.add("C1", "C", (2.7, 1.1, 0))
    c2 = m.add("C2", "C", (2.7, -1.1, 0))
    m.bond(n1, c1); m.bond(n1, c2)
    nz = (2.0, 0.0, -2.9)
    r = go(Sys().ligand(m).prot("LYS", "NZ", nz))
    check(r.has("Hydrogen bond", la=":N1"), "control: without Fe the ligand N must accept from the lysine")
    # with N coordinated to the ligand Fe the metal must be close: the Fe is at 2.0 A -> N coordinated
    m2 = Mol()
    fe2 = m2.add("FE", "Fe", (0, 0, 0))
    n = m2.add("N1", "N", (2.0, 0, 0))
    a = m2.add("C1", "C", (2.7, 1.1, 0)); b = m2.add("C2", "C", (2.7, -1.1, 0))
    m2.bond(n, a); m2.bond(n, b)
    r = go(Sys().ligand(m2).prot("LYS", "NZ", nz))
    check(not r.has("Hydrogen bond", la=":N1"), "N coordinated to the ligand Fe treated as an acceptor")
    # histidine coordinating the ligand Fe: no H-bond with the ligand N
    r = go(Sys().ligand(m2).prot("HIS", "NE2", (0, 0, 2.1)))
    check(not r.has("Hydrogen bond", pa="HIS"), "His coordinated to the ligand Fe treated as donor/acceptor")


@test("chem")
def coordination_sphere_lists_protein_nitrogens():
    """Regression: a His N bound to a metal must still appear in the coordination sphere of that metal."""
    m = methanol()
    o1 = m.atoms[1][2]
    zn = (o1[0], o1[1] + 2.0, 0)
    ne2 = at(zn, (1, 0.5, 0), 2.1)
    s = Sys().ligand(m).metal("Zn", zn).prot("HIS", "NE2", ne2)
    r = run_raw([tmpfile(s.text()), "LIG"])
    check(r.rc == 0, "exit code %d: %s" % (r.rc, r.err[:200]))
    check("Coordination sphere" in r.out, "coordination sphere missing: " + r.out)
    sphere = r.out.split("Coordination sphere")[1]
    check("HIS 10A:NE2" in sphere and "LIG 1L:O1" in sphere, "sphere without the His N or the ligand O: " + sphere)


@test("chem")
def water_bridge_needs_realistic_angle():
    m = methanol()
    o1 = m.atoms[1][2]
    w = at(o1, (1, 0, 0), 2.8)
    s = lambda p: Sys().ligand(m).water(w).prot("ASN", "ND2", p)
    check(go(s(at(w, (0.5, 0.866, 0), 2.8))).has("Water bridge"), "120 degree angle: bridge expected")
    check(not go(s(at(w, (1, 0, 0), 2.8))).has("Water bridge"), "collinear water (180 degrees): no bridge")
    check(go(s(at(w, (1, 0, 0), 2.8)), args=["--wb-max", "180"]).has("Water bridge"), "--wb-max 180 must readmit the angle")


def rot2(v, deg):
    a = math.radians(deg)
    return (v[0] * math.cos(a) - v[1] * math.sin(a), v[0] * math.sin(a) + v[1] * math.cos(a), 0.0)


@test("chem")
def long_hbond_needs_possible_hydrogen_geometry():
    # hydroxyl (methanol): the H rotates, but the C-O...A angle must be close to 109.5 degrees
    m = methanol()
    o1 = m.atoms[1][2]
    good = at(o1, rot2((-1, 0, 0), -110), 3.9)        # C-O...A = 110 degrees, D...A = 3.9 A
    bad = at(o1, rot2((-1, 0, 0), -40), 3.9)           # 40 degrees: the H cannot point to the acceptor
    s = lambda p: Sys().ligand(m).prot("ASN", "OD1", p)
    check(not go(s(good)).has("Hydrogen bond", la=":O1"), "3.9 A without --hb-long: no H-bond")
    r = go(s(good), args=["--hb-long", "4.1"])
    check(r.has("Hydrogen bond", la=":O1") and "long" in r.find("Hydrogen bond")[0]["note"], "3.9 A with a possible geometry: expected long H-bond")
    check(not go(s(bad), args=["--hb-long", "4.1"]).has("Hydrogen bond", la=":O1"), "impossible geometry: no long H-bond")
    check(not go(s(good), args=["--hb-long", "3.8"]).has("Hydrogen bond", la=":O1"), "beyond the --hb-long limit")
    # amide NH2: the expected C-N...A angle is ~120 degrees
    a = acetamide()
    n = a.atoms[3][2]
    v1 = v_unit(v_sub(a.atoms[1][2], n))
    ok = at(n, rot2(v1, 120), 3.9)
    r = go(Sys().ligand(a).prot("ASN", "OD1", ok), args=["--hb-long", "4.1"])
    check(r.has("Hydrogen bond", la=":N1"), "amide NH2 with a possible geometry: " + str(r.types()))
    no = at(n, v1, 3.9)                               # acceptor on the C side: impossible
    r = go(Sys().ligand(a).prot("ASN", "OD1", no), args=["--hb-long", "4.1"])
    check(not r.has("Hydrogen bond", la=":N1"), "NH2 with the acceptor on the wrong side")


@test("chem")
def salt_bridge_by_group_centroid():
    c = carboxylate()
    near = (0.62 + 4.6, 0, 0)                          # atom-atom 4.7 A (>4.0), centroid-centroid 4.6 A
    s = lambda: Sys().ligand(c).prot("LYS", "NZ", near)
    check(not go(s()).has("Salt bridge"), "default: 4.7 A atom-atom beyond salt-max 4.0")
    r = go(s(), args=["--salt-centroid", "5.0"])
    check(len(r.find("Salt bridge")) == 1, "centroids at 4.6 A: expected a salt bridge: " + str(r.types()))
    check(not go(s(), args=["--salt-centroid", "4.0"]).has("Salt bridge"), "centroids beyond the limit")
    # one charged group = one contact, not one per atom pair
    g = guanidinium()
    mk = lambda: Sys().ligand(g).prot("GLU", "OE1", (-0.5, 4.5, 0)).prot("GLU", "OE2", (0.5, 4.5, 0))
    check(len(go(mk()).find("Salt bridge")) > 1, "default: several atom-atom pairs")
    check(len(go(mk(), args=["--salt-centroid", "5.0"]).find("Salt bridge")) == 1, "per group: a single salt bridge")


@test("chem")
def plip_like_preset_and_explicit_override():
    m = methanol()
    o1 = m.atoms[1][2]
    good = at(o1, rot2((-1, 0, 0), -110), 3.9)
    s = lambda: Sys().ligand(m).prot("ASN", "OD1", good)
    check(not go(s()).has("Hydrogen bond", la=":O1"), "default: no H-bond at 3.9 A")
    check(go(s(), args=["--plip-like"]).has("Hydrogen bond", la=":O1"), "--plip-like must admit the long H-bond")
    check(not go(s(), args=["--plip-like", "--hb-long", "0"]).has("Hydrogen bond", la=":O1"), "explicit option after the preset")
    check(not go(s(), args=["--hb-long", "0", "--plip-like"]).has("Hydrogen bond", la=":O1"), "explicit option before the preset")
    c = carboxylate()
    near = (0.62 + 4.6, 0, 0)
    mk = lambda: Sys().ligand(c).prot("LYS", "NZ", near)
    check(go(mk(), args=["--plip-like"]).has("Salt bridge"), "--plip-like: salt bridge by centroids at 4.6 A")
    check(not go(mk(), args=["--plip-like", "--salt-centroid", "0"]).has("Salt bridge"), "--salt-centroid 0 turns the centroids off")


@test("chem")
def nucleic_acid_stacking_and_phosphate():
    # adenine base (6+5 atoms) at z=0, benzene ligand above it
    c6 = ring_points((0, 0, 0), (0, 0, 1), 1.38)            # N1 C2 N3 C4 C5 C6
    names6 = ["N1", "C2", "N3", "C4", "C5", "C6"]
    s = Sys()
    for nm, p in zip(names6, c6):
        s.atom("ATOM", nm, "DA", "D", 17, p, nm[0])
    c4, c5 = c6[3], c6[4]
    mid = v_mul(v_add(c4, c5), 0.5)
    out = v_unit(mid)
    cen = v_add(mid, v_mul(out, 0.95))
    pent = []
    a0 = math.atan2(c4[1] - cen[1], c4[0] - cen[0])
    sign = 1
    p1 = (cen[0] + 1.174 * math.cos(a0 + math.radians(72)), cen[1] + 1.174 * math.sin(a0 + math.radians(72)))
    if math.hypot(p1[0] - c5[0], p1[1] - c5[1]) > 0.3:
        sign = -1
    for k, nm in [(2, "N7"), (3, "C8"), (4, "N9")]:
        a = a0 + sign * math.radians(72 * k)
        s.atom("ATOM", nm, "DA", "D", 17, (cen[0] + 1.174 * math.cos(a), cen[1] + 1.174 * math.sin(a), 0.0), nm[0])
    lig = Mol()
    for k, p in enumerate(ring_points((0.3, 0, 3.5), (0, 0, 1), 1.40)):
        lig.add("C%d" % (k + 1), "C", p)
    for k in range(6):
        lig.bond(k, (k + 1) % 6, 1 + (k % 2))
    s.ligand(lig)
    r = go(s)
    check(r.has("pi-stacking (parallel)", pa="DA 17D"), "expected stacking on a DNA base: " + str(r.types()))
    # anionic phosphate + ligand amine
    s = Sys().ligand(trimethylamine()).prot("DA", "OP1", (0, 0, 3.0), resi=5, chain="D", el="O")
    r = go(s)
    check(r.has("Salt bridge", pa="DA 5D:OP1"), "expected salt bridge with the DNA phosphate: " + str(r.types()))


# ----------------------------------------------------------------------------------------------
# GROUP tmpl: SDF / MOL2 templates
# ----------------------------------------------------------------------------------------------
@test("tmpl")
def sdf_with_explicit_h_distinguishes_acid_from_carboxylate():
    m = carboxylate()
    probe = at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8)
    s = lambda: Sys().ligand(m).prot("LYS", "NZ", probe)
    check(go(s()).has("Salt bridge"), "without template: anionic carboxylate (physiological pH assumption)")
    r = go(s(), sdf=m.sdf(with_h=True))
    check("template" in r.err, "template not applied: " + r.err)
    check(not r.has("Salt bridge"), "with H on the oxygen the group is an acid: no salt bridge")
    check(r.has("Hydrogen bond", la=":O1"), "the acid O=C must accept the H-bond from the lysine")


@test("tmpl")
def template_hbond_uses_hydrogen_geometry():
    m = methanol()
    o1 = m.atoms[1][2]
    h = v_add(o1, v_mul(HDIR, 0.96))
    aligned = at(o1, HDIR, 2.8)                       # the H points towards the acceptor
    mirror = at(o1, (0.317, -0.948, 0), 2.8)           # the H points elsewhere
    sdf = m.sdf(with_h=True)
    r = go(Sys().ligand(m).prot("ASN", "OD1", aligned), sdf=sdf)
    check(r.has("Hydrogen bond", la=":O1") and "D-H...A" in r.find("Hydrogen bond")[0]["note"],
          "aligned H: expected H-bond with D-H...A angle")
    r = go(Sys().ligand(m).prot("ASN", "OD1", mirror), sdf=sdf)
    check(not r.has("Hydrogen bond", la=":O1"), "H pointing elsewhere: no H-bond with template")
    r = go(Sys().ligand(m).prot("ASN", "OD1", mirror))
    check(r.has("Hydrogen bond", la=":O1"), "without H (heuristics) the contact is accepted")


@test("tmpl")
def template_without_hydrogens_keeps_donors():
    """Regression: an SDF without H must not zero the donors (heuristics are used)."""
    m = methanol()
    o1 = m.atoms[1][2]
    r = go(Sys().ligand(m).prot("ASN", "OD1", at(o1, HDIR, 2.8)), sdf=m.sdf(with_h=False))
    check("template" in r.err, "template not applied")
    check(r.has("Hydrogen bond", la=":O1"), "template without H: the methanol O-H must donate")


@test("tmpl")
def template_mismatch_falls_back_to_heuristics():
    m = methanol()
    s = lambda: Sys().ligand(m).prot("ASN", "OD1", at(m.atoms[1][2], HDIR, 2.8))
    r = go(s(), sdf=m.sdf(shift=(1.0, 0, 0)))
    check("no template matches" in r.err, "SDF shifted by 1 A must not match: " + r.err)
    other = ethane()
    r = go(s(), sdf=other.sdf())
    check("no template matches" in r.err, "different number of atoms: fallback expected")
    r = go(s(), sdf=other.sdf(name="other") + m.sdf(name="right"))
    check("right" in r.err, "multi-molecule file: must choose the molecule that matches: " + r.err)


@test("tmpl")
def mol2_guanidinium_via_cationic_type():
    g = guanidinium(types=True)
    s = Sys().ligand(g).prot("GLU", "OE1", at((0, 1.33, 0), (0, 1, 0), 2.9))
    r = go(s, mol2=g.mol2())
    check("template" in r.err, "mol2 not applied: " + r.err)
    check(r.has("Salt bridge", la=":N1"), "C.cat must make the three N cationic: " + str(r.types()))


@test("tmpl")
def mol2_aromatic_bonds_give_rings():
    b = benzene(types=True)
    s = Sys().ligand(b).ring("PHE", (1.0, 0, 3.6), (0, 0, 1))
    r = go(s, mol2=b.mol2())
    check("template" in r.err and r.has("pi-stacking (parallel)"), "stacking with mol2 template: " + r.err + str(r.types()))


# ----------------------------------------------------------------------------------------------
# GROUP pdb: parsing
# ----------------------------------------------------------------------------------------------
def lys_salt_system(**kw):
    m = carboxylate()
    return Sys().ligand(m), at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8)


@test("pdb")
def altloc_first_encountered_is_used():
    s, p = lys_salt_system()
    far = (30.0, 30.0, 30.0)
    s.prot("LYS", "NZ", p, alt="A").prot("LYS", "NZ", far, alt="B")
    check(go(s).has("Salt bridge"), "altloc A near, B far: expected contact")
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", far, alt="B").prot("LYS", "NZ", p, alt="A")
    check(not go(s).has("Salt bridge"), "the first altloc encountered is used (B, far)")


@test("pdb")
def insertion_codes_are_distinct_residues():
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", p, resi=52, icode=" ")
    s.prot("LYS", "NZ", at((0.62, -1.07, 0), (0.5, -0.866, 0), 2.8), resi=52, icode="A")
    r = go(s)
    check(len(r.find("Salt bridge")) == 2, "expected two salt bridges (52 and 52A): " + str([x["pa"] for x in r.rows]))
    check(any("52AA" in x["pa"] for x in r.rows), "the insertion code must appear in the label")


@test("pdb")
def only_first_model_is_read():
    s, p = lys_salt_system()
    far = (30.0, 30.0, 30.0)
    txt = "MODEL        1\n" + s.prot("LYS", "NZ", far).text().replace("END\n", "") + "ENDMDL\n"
    s2, _ = lys_salt_system()
    txt += "MODEL        2\n" + s2.prot("LYS", "NZ", p).text() + "ENDMDL\nEND\n"
    check(not go(None, text=txt).has("Salt bridge"), "the second model must be ignored")


@test("pdb")
def ligand_selection_and_listing():
    s = Sys()
    s.ligand(carboxylate(), chain="L", resi=1)
    m2 = carboxylate()
    for name, el, xyz, _, _ in m2.atoms:
        s.atom("HETATM", name, "LIG", "M", 2, v_add(xyz, (40, 0, 0)), el)
    s.prot("LYS", "NZ", at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8))
    s.prot("LYS", "NZ", at(v_add((0.62, 1.07, 0), (40, 0, 0)), (0.5, 0.866, 0), 2.8), resi=11)
    check(len({r["lig"] for r in go(s).rows}) == 2, "without filter: two copies of the ligand")
    r = go(s, lig=("LIG", "M", "2"))
    check({r_["lig"] for r_ in r.rows} == {"LIG 2M"}, "chain/residue filter: " + str({r_["lig"] for r_ in r.rows}))
    r = run_raw([tmpfile(s.text())])
    check("LIG(2)" in r.out, "ligand listing with copy count: " + r.out)


@test("pdb")
def linked_hetero_residues_form_one_ligand():
    """Covalently linked HETATM residues (glycans, dipeptides) = a single ligand, as in PLIP."""
    def build(gap):
        s = Sys().ligand(methanol(), res="LIG", chain="L", resi=1)
        x = 1.43 + gap
        s.atom("HETATM", "C1", "MOL", "L", 2, (x, 0, 0), "C")
        s.atom("HETATM", "O1", "MOL", "L", 2, (x + 1.43, 0, 0), "O")
        s.prot("ASN", "OD1", at((x + 1.43, 0, 0), (0.5, 0.866, 0), 2.8))
        return s
    r = go(build(1.42))
    check("Composite" in r.err, "merge not reported: " + r.err)
    check(r.has("Hydrogen bond", la="MOL 2L:O1"), "the atoms of the linked residue must be part of the ligand: " + str(r.types()))
    check(not any("MOL 2L" in x["pa"] for x in r.rows), "the linked residue was also counted as a protein partner")
    r = go(build(1.42), args=["--no-merge"])
    check(not r.has("Hydrogen bond", la="MOL 2L:O1"), "--no-merge must keep the residues separate")
    r = go(build(3.2))
    check("Composite" not in r.err and not any("MOL 2L" in x["la"] for x in r.rows), "unlinked residues merged by mistake")


@test("pdb")
def ligand_selection_with_insertion_code():
    """Residue number with an optional insertion code: 5 (any), 5A (code A), 5- (no code)."""
    s = Sys()
    s.atom("HETATM", "C1", "LIG", "L", 5, (0, 0, 0), "C")
    s.atom("HETATM", "C1", "LIG", "L", 5, (30, 0, 0), "C", icode="A")
    s.atom("HETATM", "C1", "LIG", "L", 6, (60, 0, 0), "C")
    def headers(spec):
        p = run_raw([tmpfile(s.text()), "LIG", "L", spec])
        check(p.rc == 0, "%s: rc=%d %s" % (spec, p.rc, p.err[:150]))
        return re.findall(r"=== Ligand (\S+ \d+\S?) chain", p.out)
    check(sorted(headers("5")) == ["LIG 5", "LIG 5A"], "any code: both residues, got %s" % headers("5"))
    check(headers("5A") == ["LIG 5A"] and headers("5-") == ["LIG 5"], "5A and 5- must select one residue each")
    check(headers("5B") == [], "5B: no such residue")
    for bad in ("5AB", "x", "A5", "", "5-A", "1234567"):
        r = run_raw([tmpfile(s.text()), "LIG", "L", bad])
        check(r.rc == 1 and "Invalid residue number" in r.err, "%r must be rejected: rc=%d %s" % (bad, r.rc, r.err[:100]))


def altloc_system():
    """Methanol ligand; ASN 10 OD1 is far away in conformation A and 2.8 A from O1 in conformation B."""
    m = methanol()
    o1 = m.atoms[1][2]
    s = Sys().ligand(m)
    s.atom("ATOM", "OD1", "ASN", "A", 10, (9, 0, 0), "O", alt="A")
    s.atom("ATOM", "OD1", "ASN", "A", 10, at(o1, (0.5, 0.866, 0), 2.8), "O", alt="B")
    s.atom("ATOM", "OG", "SER", "A", 11, (20, 0, 0), "O", alt="A")
    return s


@test("pdb")
def altloc_option_selects_the_conformation():
    s = altloc_system()
    pdb = tmpfile(s.text())
    default = run_raw([pdb, "LIG", "--csv"])
    check(default.rc == 0 and len(default.out.strip().splitlines()) == 1, "default (first conformation, A): no contact: " + default.out)
    b = run_raw([pdb, "LIG", "--altloc", "B", "--csv"])
    check(b.rc == 0 and "Hydrogen bond" in b.out and "ASN 10A:OD1" in b.out, "conformation B: H-bond expected: " + b.out + b.err)
    check("Alternate conformation B: used for 1 of 2 residues" in b.err, "the choice must be reported on stderr: " + b.err)
    a = run_raw([pdb, "LIG", "--altloc", "A", "--csv"])
    check(a.rc == 0 and "Hydrogen bond" not in a.out, "conformation A: no H-bond")
    c = run_raw([pdb, "LIG", "--altloc", "C", "--csv"])
    check(c.rc == 0 and "alternate conformation C not found" in c.err and "present: A, B" in c.err, "missing conformation must warn: " + c.err)
    # the option combines with the others and with the fingerprint
    fp = run_raw([pdb, "LIG", "--altloc", "B", "--fingerprint", "csv"])
    check(fp.rc == 0 and "ASN 10A" in fp.out, "fingerprint with --altloc: " + fp.out)
    fp0 = run_raw([pdb, "LIG", "--fingerprint", "csv"])
    check(fp0.rc == 0 and "ASN 10A" not in fp0.out, "fingerprint without --altloc: no contact")
    check(run_raw([pdb, "LIG", "--altloc", "B", "--plip-like", "--no-water", "--no-merge"]).rc == 0, "combined options")


@test("pdb")
def altloc_listing_and_argument_errors():
    pdb = tmpfile(altloc_system().text())
    r = run_raw([pdb])
    check(r.rc == 0 and "Alternate conformations present: A, B" in r.out, "listing must show the conformations: " + r.out)
    plain = run_raw([tmpfile(fp_system().text())])
    check("Alternate conformations" not in plain.out, "no line when the structure has none")
    for args, what in [([pdb, "LIG", "--altloc"], "missing value"), ([pdb, "LIG", "--altloc", "AB"], "two characters"),
                       ([pdb, "LIG", "--altloc", ""], "empty"), ([pdb, "LIG", "--altloc", " "], "a blank")]:
        p = run_raw(args)
        check(p.rc == 1 and "altloc" in p.err, "%s: expected code 1, got %d (%s)" % (what, p.rc, p.err[:100]))
        no_crash(p)
    # every conformation can be chosen for a ligand with alternate conformations, and atoms in a second model are ignored
    t = Sys()
    t.atom("HETATM", "C1", "LIG", "L", 1, (0, 0, 0), "C", alt="A")
    t.atom("HETATM", "O1", "LIG", "L", 1, (1.43, 0, 0), "O", alt="A")
    t.atom("HETATM", "C1", "LIG", "L", 1, (50, 0, 0), "C", alt="B")
    t.atom("HETATM", "O1", "LIG", "L", 1, (51.43, 0, 0), "O", alt="B")
    t.prot("ASN", "OD1", at((1.43, 0, 0), (0.5, 0.866, 0), 2.8))
    pdbt = tmpfile(t.text())
    check("Hydrogen bond" in run_raw([pdbt, "LIG", "--altloc", "A", "--csv"]).out, "ligand in conformation A touches the ASN")
    check("Hydrogen bond" not in run_raw([pdbt, "LIG", "--altloc", "B", "--csv"]).out, "ligand in conformation B is far away")
    t.lines.append("ENDMDL")
    t.atom("ATOM", "OD1", "ASN", "A", 10, (0, 0, 0), "O", alt="C")
    check(run_raw([tmpfile(t.text()), "LIG", "--altloc", "C", "--csv"]).rc == 0, "altloc after ENDMDL is ignored without errors")


@test("pdb")
def csv_is_well_formed():
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", p)
    r = go(s)
    lines = r.out.strip().splitlines()
    check(lines[0] == "ligand,type,favorable,ligand_atom,protein_atom,distance,notes", "CSV header: " + lines[0])
    for l in lines[1:]:
        check(len(l.split(",")) == 7, "CSV row with the wrong number of fields: " + l)
    check("annotation" not in r.out.lower() and "B-factor" not in r.out, "service information on stdout in CSV mode")


# ----------------------------------------------------------------------------------------------
# GROUP fp: interaction fingerprints
# ----------------------------------------------------------------------------------------------
SVG_NS = "{http://www.w3.org/2000/svg}"


def fp_system():
    m = methanol()
    return Sys().ligand(m).prot("ASN", "OD1", at(m.atoms[1][2], (0.5, 0.866, 0), 2.8))


def run_fp(system, fmt, extra=(), text=None):
    pdb = tmpfile(text if text is not None else system.text())
    return run_raw([pdb, "LIG", "--fingerprint", fmt] + list(extra))


def pocket_file(lines, newline="\n"):
    return tmpfile(newline.join(lines) + newline, ".tsv")


@test("fp")
def fingerprint_csv_json_svg_are_consistent():
    s = fp_system()
    c = run_fp(s, "csv")
    check(c.rc == 0, "csv rc=%d: %s" % (c.rc, c.err[:200]))
    lines = c.out.strip().splitlines()
    head = lines[0].split(",")
    check(head[:6] == ["ligand", "residue", "chain", "resi", "position", "region"] and "hbond_acceptor" in head, "CSV header: " + lines[0])
    check(len(lines) == 2, "expected a single residue: " + c.out)
    row = dict(zip(head, lines[1].split(",")))
    check(row["residue"] == "ASN 10A" and row["hbond_acceptor"] == "1" and row["hbond_donor"] == "0", "CSV row: " + lines[1])
    j = run_fp(s, "json")
    check(j.rc == 0, "json rc=%d" % j.rc)
    data = json.loads(j.out)
    check(len(data) == 1 and data[0]["ligand"] == "LIG 1L" and len(data[0]["classes"]) == len(head) - 6, "JSON structure")
    res = data[0]["residues"][0]
    check(res["residue"] == "ASN 10A" and res["counts"]["hbond_acceptor"] == 1 and abs(res["min_distance"]["hbond_acceptor"] - 2.8) < 0.01,
          "JSON residue: " + str(res))
    v = run_fp(s, "svg")
    check(v.rc == 0, "svg rc=%d" % v.rc)
    root = ET.fromstring(v.out)                              # well-formed XML
    check(root.tag == SVG_NS + "svg" and root.get("width") and root.get("viewBox"), "svg root")
    check("ASN 10A" in v.out and "LIG 1L" in v.out and "H-bond acceptor" in v.out, "svg content")


@test("fp")
def fingerprint_uses_the_pocket_definition():
    s = fp_system()
    pk = pocket_file(["# example", "@name My pocket", "@color regA #336699",
                      "A\t9\t1\tregA\tALA", "A\t10\t2\tregA\tASN", "A\t11\t3\tregB\tGLY"])
    c = run_fp(s, "csv", ["--pocket", pk])
    check(c.rc == 0, "rc=%d: %s" % (c.rc, c.err[:200]))
    rows = [l.split(",") for l in c.out.strip().splitlines()[1:]]
    check([r[1] for r in rows] == ["ALA 9A", "ASN 10A", "GLY 11A"], "positions in definition order, empty ones included: " + c.out)
    check([r[4] for r in rows] == ["1", "2", "3"] and [r[5] for r in rows] == ["regA", "regA", "regB"], "labels and regions")
    v = run_fp(s, "svg", ["--pocket", pk])
    check(v.rc == 0 and "My pocket" in v.out and "#336699" in v.out, "svg with pocket name and region color")
    for lab in ("ALA 9A", "ASN 10A", "GLY 11A"):
        check(lab in v.out, "cell missing in the svg: " + lab)
    j = json.loads(run_fp(s, "json", ["--pocket", pk]).out)
    check(j[0]["pocket"] == "My pocket" and [r["in_pocket"] for r in j[0]["residues"]] == [True, True, True], "JSON pocket info")
    # a contact residue outside the definition is appended and flagged
    pk2 = pocket_file(["A 9 1 regA ALA"])
    j2 = json.loads(run_fp(s, "json", ["--pocket", pk2]).out)
    check([r["in_pocket"] for r in j2[0]["residues"]] == [True, False], "residue outside the pocket must be flagged")


@test("fp")
def fingerprint_pocket_file_variants_and_warnings():
    s = fp_system()
    variants = {
        "whitespace": ["A 10 p10 reg ASN"],
        "commas": ["A,10,p10,reg,ASN"],
        "any chain": ["* 10 p10 reg ASN"],
        "crlf": None,
    }
    for name, lines in variants.items():
        pk = pocket_file(["A\t10\tp10\treg\tASN"], "\r\n") if lines is None else pocket_file(lines)
        c = run_fp(s, "csv", ["--pocket", pk])
        check(c.rc == 0 and "ASN 10A,A,10,p10,reg" in c.out, "%s: %s | %s" % (name, c.out, c.err[:150]))
    bad = pocket_file(["A 10 p10 reg ASN", "A xyz 1 reg", "A 11", "@bogus 1", "@color reg nocolor"])
    c = run_fp(s, "csv", ["--pocket", bad])
    check(c.rc == 0 and c.err.count("Warning: pocket definition line") == 4, "expected 4 warnings: " + c.err)
    check("ASN 10A" in c.out, "valid lines must still be used")
    # residue name mismatch
    mm = pocket_file(["A 10 p10 reg LEU"])
    c = run_fp(s, "csv", ["--pocket", mm])
    check(c.rc == 0 and "does not match the pocket definition" in c.err, "mismatch warning: " + c.err)


@test("fp")
def fingerprint_command_line_errors():
    pdb = tmpfile(fp_system().text())
    good = pocket_file(["A 10 p10 reg ASN"])
    cases = [([pdb, "LIG", "--fingerprint"], "--fingerprint without a value"),
             ([pdb, "LIG", "--fingerprint", "xml"], "unknown format"),
             ([pdb, "LIG", "--fingerprint", "csv", "--pocket", "/no/such/pocket.tsv"], "nonexistent pocket file"),
             ([pdb, "LIG", "--pocket", good], "--pocket without --fingerprint"),
             ([pdb, "LIG", "--fingerprint", "csv", "--pocket"], "--pocket without a value")]
    for args, what in cases:
        r = run_raw(args)
        check(r.rc == 1, "%s: expected code 1, got %d (%s)" % (what, r.rc, r.err[:120]))
        no_crash(r)
    # no ligand requested: the listing is printed, no fingerprint
    r = run_raw([pdb, "--fingerprint", "csv"])
    check(r.rc == 0 and "Ligands (HETATM) present" in r.out, "listing instead of the fingerprint: " + r.out)


@test("fp")
def fingerprint_robustness():
    s = fp_system()
    random.seed(5)
    garbage = bytes(random.getrandbits(8) for _ in range(6000))
    pk = tmpfile(garbage, ".tsv", binary=True)
    for fmt in ("csv", "json", "svg"):
        r = run_fp(s, fmt, ["--pocket", pk])
        check(r.rc == 0, "%s with a garbage pocket: rc=%d %s" % (fmt, r.rc, r.err[:150]))
        no_crash(r)
    weird = pocket_file(["A 99999999 1 reg ASN", "A 10 " + "x" * 5000 + " reg ASN", "A 10 p<&\"> r&<>\" ASN", "A -5 7 reg ALA"])
    r = run_fp(s, "svg", ["--pocket", weird])
    check(r.rc == 0, "weird pocket: rc=%d" % r.rc)
    ET.fromstring(r.out)                                       # still well-formed XML despite special characters
    check("&amp;" in r.out and "&lt;" in r.out, "special characters must be escaped")
    # empty fingerprint (ligand without contacts)
    alone = Sys().ligand(methanol())
    for fmt in ("csv", "json", "svg"):
        r = run_fp(alone, fmt)
        check(r.rc == 0, "no contacts, %s: rc=%d %s" % (fmt, r.rc, r.err[:100]))
    check("No interactions" in run_fp(alone, "svg").out, "empty figure")
    check(json.loads(run_fp(alone, "json").out)[0]["residues"] == [], "empty JSON residues")
    # many residues in a pocket: figure wraps and stays well-formed
    big = pocket_file(["A %d p%d reg%d ALA" % (i, i, i // 7) for i in range(1, 300)])
    r = run_fp(s, "svg", ["--pocket", big])
    check(r.rc == 0 and len(ET.fromstring(r.out).findall(".//" + SVG_NS + "g")) > 299, "large pocket")


@test("fp")
def fingerprint_several_ligands():
    s = Sys().ligand(methanol(), chain="L", resi=1)
    m2 = methanol()
    for name, el, xyz, _, _ in m2.atoms:
        s.atom("HETATM", name, "LIG", "M", 2, v_add(xyz, (30, 0, 0)), el)
    s.prot("ASN", "OD1", at(m2.atoms[1][2], (0.5, 0.866, 0), 2.8), resi=10)
    s.prot("ASP", "OD1", at(v_add(m2.atoms[1][2], (30, 0, 0)), (0.5, 0.866, 0), 2.8), resi=11)
    c = run_fp(s, "csv")
    check(c.rc == 0 and len({l.split(",")[0] for l in c.out.strip().splitlines()[1:]}) == 2, "CSV must contain both ligands: " + c.out)
    j = json.loads(run_fp(s, "json").out)
    check(len(j) == 2, "JSON must have one object per ligand")
    v = run_fp(s, "svg")
    check(v.rc == 0 and "only the first is drawn" in v.err, "SVG with several ligands: " + v.err)
    ET.fromstring(v.out)


# ----------------------------------------------------------------------------------------------
# GROUP cif: PDBx/mmCIF input
# ----------------------------------------------------------------------------------------------
CIF_COLUMNS = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id", "label_entity_id",
               "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
               "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]


def cif_quote(v):
    return '"%s"' % v if "'" in v else v


def pdb_to_cif(text, chain_map=None, extra_before=""):
    """PDB text -> PDBx/mmCIF text (only the atoms). The label_* columns are made different from the auth_* ones on purpose:
    a reader that used the wrong columns would be caught. chain_map renames chains (e.g. {"A": "AA"})."""
    rows, model, n = [], 1, 0
    for l in text.splitlines():
        if l.startswith("ENDMDL"):
            model += 1
            continue
        if l[:6] not in ("ATOM  ", "HETATM") or len(l) < 54:
            continue
        n += 1
        chain = l[21]
        chain = (chain_map or {}).get(chain, chain)
        el = l[76:78].strip() if len(l) >= 78 else ""
        alt = l[16] if l[16] != " " else "."
        ins = l[26] if l[26] != " " else "?"
        name = l[12:16].strip()
        seq = int(l[22:26])
        rows.append("%s %d %s %s %s %s %sx 1 %d %s %s %s %s 1.00 %s ? %d %s %s %s %d" % (
            l[:6].strip(), n, el or "?", cif_quote(name), alt, l[17:20].strip(), chain if chain != " " else "Z",
            seq + 5000, ins, l[30:38].strip(), l[38:46].strip(), l[46:54].strip(), l[60:66].strip() or "0.00",
            seq, l[17:20].strip(), chain if chain != " " else "?", cif_quote(name), model))
    return ("data_TEST\n#\n_entry.id TEST\n" + extra_before + "loop_\n" + "".join("_atom_site.%s\n" % c for c in CIF_COLUMNS) +
            "\n".join(rows) + "\n#\n")


def both_formats(system_or_text, args, **kw):
    """Runs the same analysis on a PDB text and on its mmCIF conversion; returns both results."""
    text = system_or_text if isinstance(system_or_text, str) else system_or_text.text()
    a = run_raw([tmpfile(text), *args], **kw)
    b = run_raw([tmpfile(pdb_to_cif(text), ".cif"), *args], **kw)
    return a, b


def same_result(a, b, what):
    check(a.rc == b.rc and a.out == b.out and a.err == b.err,
          "%s: PDB and mmCIF differ\n--- PDB\n%s%s\n--- mmCIF\n%s%s" % (what, a.out[:400], a.err[:200], b.out[:400], b.err[:200]))


@test("cif")
def cif_matches_pdb_on_synthetic_systems():
    systems = {"carboxylate+lysine": Sys().ligand(carboxylate()).prot("LYS", "NZ", at(carboxylate().atoms[2][2], (0.5, 0.866, 0), 2.8)),
               "methanol+asn": fp_system(), "alternate conformations": altloc_system()}
    for name, s in systems.items():
        for args in (["LIG", "--csv"], ["LIG"], ["LIG", "--altloc", "B", "--csv"], ["LIG", "--plip-like", "--csv"], ["LIG", "--fingerprint", "csv"],
                     ["LIG", "--fingerprint", "svg"], []):
            a, b = both_formats(s, args)
            same_result(a, b, "%s %s" % (name, " ".join(args)))
    # insertion codes, a second model, a blank chain
    t = Sys()
    t.atom("HETATM", "C1", "LIG", "L", 5, (0, 0, 0), "C")
    t.atom("HETATM", "C1", "LIG", "L", 5, (30, 0, 0), "C", icode="A")
    t.atom("HETATM", "O1", "LIG", "L", 5, (31.43, 0, 0), "O", icode="A")
    t.prot("ASN", "OD1", (33.0, 2.4, 0), resi=10)
    t.lines.append("ENDMDL")
    t.atom("HETATM", "C1", "LIG", "M", 7, (60, 0, 0), "C")
    for spec in ("5", "5A", "5-"):
        a, b = both_formats(t, ["LIG", "L", spec, "--csv"])
        same_result(a, b, "residue " + spec)
    a, b = both_formats(t, [])
    same_result(a, b, "listing with a second model")
    # (a blank chain id does not exist in mmCIF: a null auth_asym_id falls back to label_asym_id, see docs/INPUT_FORMATS.md)


@test("cif")
def cif_matches_pdb_on_real_complexes():
    n = 0
    for f in sorted(glob.glob(os.path.join(CFG["real_dir"], "*.pdb"))):
        text = open(f).read()
        listing = run_raw([f])
        names = [m.group(1) for m in re.finditer(r"(\S+?)\(\d+\)", listing.out.splitlines()[0] if listing.out else "")][:3]
        for nm in names:
            for args in ([nm, "--csv"], [nm, "--plip-like"], [nm, "--fingerprint", "csv"]):
                a, b = both_formats(text, args)
                same_result(a, b, "%s %s" % (os.path.basename(f), " ".join(args)))
                n += 1
    if n == 0:
        raise Skip("no real complexes in " + CFG["real_dir"])


@test("cif")
def cif_real_file_converted_with_gemmi_matches_pdb():
    """tests/data/1AZM.cif was written by gemmi (see tests/data/README.md): a file with the usual mmCIF structure, not made by our converter."""
    pdb, cif = os.path.join(CFG["real_dir"], "1AZM.pdb"), os.path.join(CFG["real_dir"], "1AZM.cif")
    if not (os.path.exists(pdb) and os.path.exists(cif)):
        raise Skip("missing 1AZM.pdb or 1AZM.cif")
    for args in ([], ["AZM", "--csv"], ["AZM"], ["AZM", "--plip-like", "--csv"], ["AZM", "--fingerprint", "csv"], ["ZN", "--csv"]):
        a, b = run_raw([pdb] + args), run_raw([cif] + args)
        same_result(a, b, "1AZM " + " ".join(args))
    r = run_raw([cif, "AZM"])
    check("Coordination sphere" in r.out and "HIS 94A:NE2" in r.out and "Hydrogen bond" in r.out, "the mmCIF analysis must show the Zn sphere: " + r.out[:300])


@test("cif")
def cif_long_chain_ids_and_residue_names():
    m = methanol()
    o1 = m.atoms[1][2]
    # a PDB line cannot hold a five-character component id: convert a normal system and rename the ligand in the mmCIF text
    s = Sys().ligand(m, res="LIG").prot("ASN", "OD1", at(o1, (0.5, 0.866, 0), 2.8), chain="A")
    text = re.sub(r"\bLIG\b", "ABCDE", pdb_to_cif(s.text(), chain_map={"A": "AA"}))
    p = tmpfile(text, ".cif")
    r = run_raw([p, "ABCDE", "--csv"])
    check(r.rc == 0 and "ASN 10 [AA]:OD1" in r.out and r.out.count("ABCDE 1L") >= 1, "long chain label: " + r.out + r.err)
    check("Hydrogen bond" in r.out, "the H-bond must be found")
    f = run_raw([p, "ABCDE", "--fingerprint", "csv"])
    rows = [l.split(",") for l in f.out.strip().splitlines()[1:]]
    check(f.rc == 0 and len(rows) == 1 and rows[0][1] == "ASN 10 [AA]" and rows[0][2] == "AA" and rows[0][3] == "10",
          "fingerprint with a long chain: " + f.out)
    # ligand chain of two characters selected by name, and the listing
    check("ABCDE(1)" in run_raw([p]).out, "listing shows the long name")
    t2 = re.sub(r" L ", " LL ", text).replace(" Lx ", " LLx ")
    sel = run_raw([tmpfile(t2, ".cif"), "ABCDE", "LL", "--csv"])
    check(sel.rc == 0 and "ABCDE 1 [LL]" in sel.out, "ligand with a chain of two characters: " + sel.out + sel.err)
    check("No ligand 'ABCDE' found" in run_raw([tmpfile(t2, ".cif"), "ABCDE", "L", "--csv"]).out, "chain L does not exist")
    # a pocket definition with a long chain
    pk = tmpfile("@name P\nAA\t10\tp10\treg\tASN\nAA\t11\tp11\treg\tGLY\n", ".tsv")
    v = run_raw([p, "ABCDE", "--fingerprint", "csv", "--pocket", pk])
    rows = [l.split(",") for l in v.out.strip().splitlines()[1:]]
    check(v.rc == 0 and [r[1] for r in rows] == ["ASN 10 [AA]", "GLY 11 [AA]"] and rows[0][4] == "p10", "pocket with a long chain: " + v.out + v.err)


@test("cif")
def cif_format_option_and_compressed_files():
    cif = tmpfile(pdb_to_cif(fp_system().text()), ".cif")
    pdb = tmpfile(fp_system().text())
    for args, what in (([cif, "LIG", "--csv"], "auto"), ([cif, "LIG", "--csv", "--format", "cif"], "cif"), ([cif, "LIG", "--csv", "--format", "mmcif"], "mmcif"),
                       ([cif, "LIG", "--csv", "--format", "auto"], "auto explicit")):
        r = run_raw(args)
        check(r.rc == 0 and "Hydrogen bond" in r.out, "%s: %s" % (what, r.out + r.err))
    check("Hydrogen bond" in run_raw([pdb, "LIG", "--csv", "--format", "pdb"]).out, "--format pdb on a PDB")
    forced = run_raw([pdb, "LIG", "--csv", "--format", "cif"])
    check(forced.rc == 0 and "no _atom_site category" in forced.err and "Hydrogen bond" not in forced.out, "PDB read as mmCIF: " + forced.err)
    forced = run_raw([cif, "LIG", "--csv", "--format", "pdb"])
    check(forced.rc == 0 and "Hydrogen bond" not in forced.out, "mmCIF read as PDB finds nothing")
    for args, what in (([cif, "--format"], "missing value"), ([cif, "--format", "xyz"], "unknown format")):
        r = run_raw(args)
        check(r.rc == 1 and "--format" in r.err, "%s: rc=%d %s" % (what, r.rc, r.err[:100]))
        no_crash(r)
    # gzip: refused with a clear message, whatever the extension
    for suffix in (".cif.gz", ".pdb.gz", ".cif"):
        g = tmpfile(b"\x1f\x8b\x08\x00" + bytes(random.getrandbits(8) for _ in range(200)), suffix, binary=True)
        r = run_raw([g, "LIG"])
        check(r.rc == 1 and "compressed" in r.err and "gunzip" in r.err, "%s: rc=%d %s" % (suffix, r.rc, r.err[:150]))
        no_crash(r)
    # the listing shows the alternate conformations of an mmCIF file too
    r = run_raw([tmpfile(pdb_to_cif(altloc_system().text()), ".cif")])
    check("Alternate conformations present: A, B" in r.out, "listing: " + r.out)


@test("cif")
def cif_robustness():
    random.seed(11)
    valid = pdb_to_cif(fp_system().text(), extra_before="_struct.title\n;a title\nwith lines\n;\n_cell.length_a 10\nloop_\n_x.a\n_x.b\n1 'it''s'\n")
    cases = {
        "empty": "", "only data": "data_X\n", "only loop": "data_X\nloop_\n", "tags without values": "data_X\nloop_\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n",
        "unterminated quote": "data_X\nloop_\n_atom_site.group_PDB\n_atom_site.Cartn_x\nATOM 'unterminated\n",
        "unterminated text field": "data_X\n_struct.title\n;never closed\nATOM\n", "text field in the atom loop":
        "data_X\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n;\nHETATM\n;\nC1 1.0\n",
        "huge token": "data_X\n_struct.title " + "x" * 500000 + "\n", "binary": bytes(random.getrandbits(8) for _ in range(5000)).decode("latin-1"),
        "many tags": "data_X\nloop_\n" + "".join("_atom_site.t%d\n" % i for i in range(5000)) + "1 2 3\n",
        "nan coordinates": valid.replace(" 20.00 ", " 20.00 ").replace("0.000 ", "nan ", 1),
        "null coordinates": re.sub(r"(\d+\.\d{3}) (\d+\.\d{3}) (\d+\.\d{3}) 1\.00", r"? ? ? 1.00", valid, count=1),
    }
    for name, text in cases.items():
        r = run_raw([tmpfile(text.encode("latin-1"), ".cif", binary=True) if name == "binary" else tmpfile(text, ".cif"), "LIG", "--csv"])
        check(r.rc == 0, "%s: rc=%d %s" % (name, r.rc, r.err[:150]))
        no_crash(r)
    # random mutations of a valid file
    for i in range(60):
        t = list(valid)
        for _ in range(8):
            k = random.randrange(len(t))
            op = random.choice("dirt")
            if op == "d": del t[k:k + random.randint(1, 10)]
            elif op == "i": t.insert(k, random.choice("';\"#_\n \t?."))
            elif op == "r": t[k] = random.choice("';\"#_\n \t?.0a")
            else: t = t[:k] or ["d"]
        r = run_raw([tmpfile("".join(t), ".cif"), "LIG", "--csv"])
        check(r.rc == 0, "mutation %d: rc=%d %s" % (i, r.rc, r.err[:120]))
        no_crash(r)


@test("cif")
def cif_large_file_runs_in_reasonable_time():
    random.seed(4)
    s = Sys().ligand(carboxylate())
    for i in range(40000):
        base = (random.uniform(-60, 60), random.uniform(-60, 60), random.uniform(-60, 60))
        for nm in ("N", "CA", "C", "O", "CB"):
            s.atom("ATOM", nm, "ALA", "A", i % 9999 + 1, (base[0] + random.uniform(-1.5, 1.5), base[1] + random.uniform(-1.5, 1.5),
                                                          base[2] + random.uniform(-1.5, 1.5)), nm[0])
    text = pdb_to_cif(s.text())
    limit = 120 if CFG["sanitize"] else 40
    t0 = time.time()
    r = run_raw([tmpfile(text, ".cif"), "LIG", "--csv"], timeout=limit * 3)
    check(r.rc == 0 and time.time() - t0 < limit, "200000 atoms in mmCIF: %.1f s (limit %d s)" % (time.time() - t0, limit))


# ----------------------------------------------------------------------------------------------
# GROUP manual: the HTML manual (docs/manual/html), built from Markdown by tools/build_manual.py
# ----------------------------------------------------------------------------------------------
MANUAL_DIR = os.path.normpath(os.path.join(HERE, "..", "docs", "manual", "html"))


class _PageInfo(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids, self.hrefs, self.imgs, self.tags = set(), [], [], set()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self.tags.add(tag)
        if "id" in a: self.ids.add(a["id"])
        if tag == "a" and "name" in a: self.ids.add(a["name"])
        if tag == "a" and "href" in a: self.hrefs.append(a["href"])
        if tag == "img" and "src" in a: self.imgs.append(a["src"])


def manual_pages():
    if not os.path.isdir(MANUAL_DIR):
        raise Skip("no manual in " + MANUAL_DIR)
    info = {}
    for p in sorted(glob.glob(os.path.join(MANUAL_DIR, "*.html"))):
        parser = _PageInfo()
        parser.feed(open(p, encoding="utf-8").read())
        info[os.path.basename(p)] = parser
    return info


@test("manual")
def manual_html_is_up_to_date_with_its_sources():
    script = os.path.normpath(os.path.join(HERE, "..", "tools", "build_manual.py"))
    if not os.path.exists(script) or not os.path.isdir(MANUAL_DIR):
        raise Skip("no manual sources")
    r = subprocess.run([sys.executable, script, "--check"], capture_output=True, text=True)
    if r.returncode == 2 or "needs the Python module" in r.stderr:
        raise Skip(r.stderr.strip()[:200])
    check(r.returncode == 0, "the HTML manual is out of date: run python3 tools/build_manual.py\\n" + r.stderr[:300])


@test("manual")
def manual_links_anchors_and_images_resolve():
    pages = manual_pages()
    check(len(pages) >= 12, "expected at least 12 pages, found %d" % len(pages))
    for name, p in pages.items():
        check("h1" in p.tags and "title" in p.tags, "%s: no title or h1" % name)
        for href in p.hrefs:
            if re.match(r"^(https?:|mailto:)", href):
                continue
            file, _, frag = href.partition("#")
            target = pages.get(file) if file else p
            check(target is not None, "%s: link to a page that does not exist: %s" % (name, href))
            check(not frag or frag in target.ids, "%s: link to an anchor that does not exist: %s" % (name, href))
        for src in p.imgs:
            check(os.path.exists(os.path.join(MANUAL_DIR, src)), "%s: missing image %s" % (name, src))
    toc = json.load(open(os.path.join(MANUAL_DIR, "toc.json"), encoding="utf-8"))
    check([e["file"] for e in toc["pages"]] == sorted([e["file"] for e in toc["pages"]], key=lambda f: [x["file"] for x in toc["pages"]].index(f)) and
          {e["file"] for e in toc["pages"]} == set(pages), "toc.json and the pages on disk differ")
    for e in toc["pages"]:
        for h in e["headings"]:
            check(h["id"] in pages[e["file"]].ids, "toc.json: %s has no anchor %s" % (e["file"], h["id"]))
    qrc = set(re.findall(r'alias="([^"]+)"', open(os.path.join(MANUAL_DIR, "manual.qrc"), encoding="utf-8").read()))
    on_disk = {os.path.relpath(os.path.join(d, f), MANUAL_DIR).replace(os.sep, "/") for d, _, fs in os.walk(MANUAL_DIR) for f in fs} - {"manual.qrc"}
    check(qrc == on_disk, "manual.qrc and the files differ: %s" % sorted(qrc ^ on_disk))


@test("manual")
def manual_command_line_reference_matches_the_program():
    html_text = open(os.path.join(MANUAL_DIR, "cli.html"), encoding="utf-8").read() if os.path.isdir(MANUAL_DIR) else raise_skip("no manual")
    h = run_raw(["--help"])
    text = h.out + h.err
    check("Thresholds (defaults):" in text, "--help has no thresholds section: " + text[:200])
    thr = re.findall(r"--([a-z0-9-]+) (-?[0-9.]+)", text.split("Thresholds (defaults):")[1])
    check(len(thr) >= 20, "too few thresholds in --help: %d" % len(thr))
    rows = dict(re.findall(r"<td><code>--([a-z0-9-]+)</code></td>\s*<td>([^<]*)</td>", html_text))
    for name, default in thr:
        check(name in rows, "cli.html does not document --%s" % name)
        cell = rows[name].split()[0]
        check(float(cell) == float(default), "cli.html says --%s defaults to %s, the program says %s" % (name, rows[name], default))
    options = set(re.findall(r"--([a-z0-9-]+)", text.split("Thresholds (defaults):")[0]))
    check({"csv", "sdf", "mol2", "altloc", "format", "fingerprint", "pocket", "plip-like"} <= options, "unexpected --help options: %s" % sorted(options))
    for o in options:
        check("--" + o in html_text, "cli.html does not mention --%s" % o)
    documented = {n for n in re.findall(r"<td><code>--([a-z0-9-]+)", html_text)}
    known = options | {n for n, _ in thr} | {"help", "version"}
    check(documented <= known | {"hb-long"}, "cli.html documents options that the program does not have: %s" % sorted(documented - known))


def raise_skip(why):
    raise Skip(why)


@test("manual")
def manual_quickstart_example_is_the_real_output():
    pdb = os.path.join(CFG["real_dir"], "1AZM.pdb")
    if not os.path.exists(pdb) or not os.path.isdir(MANUAL_DIR):
        raise Skip("missing 1AZM.pdb or the manual")
    page = html.unescape(open(os.path.join(MANUAL_DIR, "quickstart.html"), encoding="utf-8").read())
    for args in ([], ["AZM"]):
        out = run_raw([pdb] + args).out.replace("\r\n", "\n").rstrip()
        check(out and out in page, "the example of the quick start differs from the output of: ligint 1AZM.pdb %s\n%s" % (" ".join(args), out[:300]))


# ----------------------------------------------------------------------------------------------
# GROUP net: KLIFS pocket definitions and RCSB downloads. Nothing here may touch the network: the transfers that must
# fail are sent to a proxy that does not exist (libcurl honors https_proxy); the cases that succeed use the local cache or a local mapping file.
# ----------------------------------------------------------------------------------------------
NO_NETWORK = {"https_proxy": "http://127.0.0.1:9", "HTTPS_PROXY": "http://127.0.0.1:9", "all_proxy": "http://127.0.0.1:9", "ALL_PROXY": "http://127.0.0.1:9"}
KLIFS_REGIONS = [("I", 3), ("g.l", 6), ("II", 4), ("III", 6), ("aC", 11), ("b.l", 7), ("IV", 4), ("V", 3), ("GK", 1), ("hinge", 3), ("linker", 4),
                 ("aD", 7), ("aE", 5), ("VI", 3), ("c.l", 8), ("VII", 3), ("VIII", 1), ("xDFG", 4), ("a.l", 2)]
KLIFS_HEADER = ("residue.klifs_id", "residue.id", "residue.klifs_region", "residue.aa")


def klifs_rows(first_residue=200):
    rows, p = [], 1
    for region, n in KLIFS_REGIONS:
        for _ in range(n):
            rows.append([p, first_residue + p, region, "ACDEFGHIKLMNPQRSTVWY"[p % 20]])
            p += 1
    rows[10] = [11, "_", "II", "_"]            # a residue that is missing in the crystal
    rows[20] = [21, "221A", "aC", "D"]         # an insertion code
    return rows


def klifs_file(rows, sep="\t", header=KLIFS_HEADER, eol="\n", bom=False, suffix=".tsv"):
    text = eol.join(sep.join(str(x) for x in r) for r in [header] + rows) + eol
    return tmpfile(("\ufeff" if bom else "") + text, suffix)


def klifs_cli(mapping, out, extra=(), chain="A", pdb="3EYG"):
    return run_raw(["--klifs-pdb", pdb, "--klifs-chain", chain, "--klifs-mapping", mapping, "--klifs-output", out] + list(extra))


@test("net")
def klifs_mapping_variants_give_the_same_pocket():
    rows = klifs_rows()
    outs = {}
    for name, mapping in {"tsv": klifs_file(rows), "csv": klifs_file(rows, sep=",", suffix=".csv"), "crlf": klifs_file(rows, eol="\r\n"),
                          "bom": klifs_file(rows, bom=True),
                          "extra columns and order": klifs_file([[r[2], r[1], "x", r[0], r[3]] for r in rows],
                                                                header=("residue.klifs_region", "residue.id", "other", "residue.klifs_id", "residue.aa"))}.items():
        out = os.path.join(CFG["tmp"], "pocket_%s.tsv" % name.replace(" ", "_"))
        r = klifs_cli(mapping, out)
        check(r.rc == 0 and "84 mapped positions" in r.err, "%s: rc=%d %s" % (name, r.rc, r.err[:150]))
        outs[name] = open(out, encoding="utf-8").read()
    check(len(set(outs.values())) == 1, "the variants of the same mapping must give the same pocket: " + str({k: len(v) for k, v in outs.items()}))
    text = outs["tsv"]
    lines = [l for l in text.splitlines() if l and l[0] not in "#@"]
    check(len(lines) == 84 and not any(l.split("\t")[1] == "211" for l in lines), "the missing residue (position 11) must be omitted")
    check("A\t221A\t21\taC\tASP" in lines, "insertion code and residue name: " + str([l for l in lines if "221" in l]))
    check("A\t210\t10\tII\tMET" in lines, "residue name from the one-letter code")
    check(text.count("@color ") == len(KLIFS_REGIONS), "one color per region: %d" % text.count("@color "))
    # every one-letter code gives the right residue name
    names = dict(zip("ACDEFGHIKLMNPQRSTVWY", "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()))
    got = {}
    for l in lines:
        f = l.split("\t")
        got[int(f[2])] = f[4]
    for p in range(1, 86):
        if p in (11,):
            continue
        check(got[p] == names["ACDEFGHIKLMNPQRSTVWY"[p % 20]] or p == 21, "position %d: %s" % (p, got[p]))


@test("net")
def klifs_pocket_works_with_the_fingerprint():
    s = fp_system()                                     # ligand + ASN 10 (chain A) that accepts a hydrogen bond
    mapping = klifs_file([[1, 9, "I", "A"], [2, 10, "hinge", "N"], [3, "_", "hinge", "_"], [4, 11, "hinge", "G"]], suffix=".tsv")
    out = os.path.join(CFG["tmp"], "pocket_fp.tsv")
    check(klifs_cli(mapping, out).rc == 0, "generation")
    pdb = tmpfile(s.text())
    r = run_raw([pdb, "LIG", "--fingerprint", "csv", "--pocket", out])
    rows = [l.split(",") for l in r.out.strip().splitlines()[1:]]
    check(r.rc == 0 and [x[1:6] for x in rows] == [["ALA 9A", "A", "9", "1", "I"], ["ASN 10A", "A", "10", "2", "hinge"], ["GLY 11A", "A", "11", "4", "hinge"]],
          "positions and regions of the fingerprint: " + r.out + r.err)
    check(rows[1][10] == "1", "the hydrogen bond of ASN 10 must be in the position 2: " + str(rows[1]))
    # the pocket can also be generated and used in one command (--klifs-* instead of --pocket)
    r2 = run_raw([pdb, "LIG", "--fingerprint", "csv", "--klifs-pdb", "3EYG", "--klifs-chain", "A", "--klifs-mapping", mapping])
    check(r2.rc == 0 and r2.out == r.out and "3 mapped positions" in r2.err, "one-command form: " + r2.err[:200])
    json.loads(run_raw([pdb, "LIG", "--fingerprint", "json", "--pocket", out]).out)


@test("net")
def klifs_mapping_errors_are_reported_and_nothing_is_written():
    good = [[1, 9, "I", "A"], [2, 10, "hinge", "N"]]
    bad = {"duplicate position": good + [[2, 99, "I", "A"]], "position 0": [[0, 5, "I", "A"]], "position 86": [[86, 5, "I", "A"]],
           "position not a number": [["x", 5, "I", "A"]], "position with a fraction": [["1.5", 5, "I", "A"]],
           "residue not a number": [[1, "abc", "I", "A"]], "residue with a fraction": [[1, "12.5", "I", "A"]], "residue with two letters": [[1, "12AB", "I", "A"]],
           "region missing": [[1, 5, "", "A"]], "only missing residues": [[1, "_", "I", "_"], [2, "_", "I", "_"]]}
    for name, rows in bad.items():
        out = os.path.join(CFG["tmp"], "never_%s.tsv" % name.replace(" ", "_"))
        r = klifs_cli(klifs_file(rows), out)
        check(r.rc == 1 and r.err.startswith("KLIFS: "), "%s: rc=%d %s" % (name, r.rc, r.err[:150]))
        check(not os.path.exists(out), "%s: no file must be written on error" % name)
        no_crash(r)
    texts = {"empty file": "", "header only": "\t".join(KLIFS_HEADER) + "\n", "wrong columns": "a\tb\tc\n1\t2\t3\n",
             "binary": bytes(random.Random(1).getrandbits(8) for _ in range(3000)).decode("latin-1"), "only comments": "# one\n# two\n"}
    for name, text in texts.items():
        out = os.path.join(CFG["tmp"], "never2_%s.tsv" % name.replace(" ", "_"))
        r = klifs_cli(tmpfile(text, ".tsv"), out)
        check(r.rc == 1 and r.err.startswith("KLIFS: "), "%s: rc=%d %s" % (name, r.rc, r.err[:150]))
        no_crash(r)
    # accepted and unusual: a negative residue number, a quoted CSV value with a comma, a very long region name
    for name, rows, suffix, sep in (("negative residue", [[1, "-3", "I", "A"]], ".tsv", "\t"), ("long region", [[1, 5, "x" * 50000, "A"]], ".tsv", "\t")):
        r = klifs_cli(klifs_file(rows, sep=sep, suffix=suffix), os.path.join(CFG["tmp"], "ok_%s.tsv" % name.replace(" ", "_")))
        check(r.rc == 0, "%s: rc=%d %s" % (name, r.rc, r.err[:150]))
    q = tmpfile('residue.klifs_id,residue.id,residue.klifs_region\n1,"5","a,b"\n', ".csv")
    check(klifs_cli(q, os.path.join(CFG["tmp"], "ok_quoted.tsv")).rc == 0, "quoted CSV value")
    # random corruption of a valid mapping never crashes
    rnd = random.Random(8)
    base = open(klifs_file(klifs_rows()), encoding="utf-8").read()
    for i in range(40):
        t = list(base)
        for _ in range(6):
            k = rnd.randrange(len(t))
            t[k] = rnd.choice("\t,\n\"#_ 0A9-.")
        r = klifs_cli(tmpfile("".join(t), ".tsv"), os.path.join(CFG["tmp"], "fuzz.tsv"))
        check(r.rc in (0, 1), "corrupted mapping %d: rc=%d" % (i, r.rc))
        no_crash(r)


@test("net")
def klifs_command_line_rules():
    mapping = klifs_file([[1, 9, "I", "A"], [2, 10, "hinge", "N"]])
    pdb = tmpfile(fp_system().text())
    pocket = os.path.join(CFG["tmp"], "rules_pocket.tsv")
    klifs_cli(mapping, pocket)
    base = ["--klifs-pdb", "3EYG", "--klifs-chain", "A", "--klifs-mapping", mapping]
    r = run_raw(["--klifs-pdb", "3EYG", "--klifs-mapping", mapping, "--klifs-output", os.path.join(CFG["tmp"], "x.tsv")])
    check(r.rc == 1 and "chain is required" in r.err, "chain is required: " + r.err)
    r = run_raw(base)
    check(r.rc == 1 and "--klifs-output" in r.err, "a pocket generated with no structure and no output would be lost: rc=%d %s" % (r.rc, r.err))
    r = run_raw([pdb, "LIG", "--pocket", pocket] + base + ["--fingerprint", "csv"])
    check(r.rc == 1 and "either --pocket or --klifs" in r.err, "--pocket and --klifs-* together: " + r.err)
    r = run_raw([pdb, "LIG"] + base)
    check(r.rc == 0 and "only used with --fingerprint" in r.err, "the pocket is useless without --fingerprint: warn: rc=%d %s" % (r.rc, r.err))
    r = run_raw(base + ["--klifs-output", "/nonexistent/folder/p.tsv"])
    check(r.rc == 1 and "Cannot write" in r.err, "unwritable output: " + r.err)
    for bad in ("abc", "0", "-4", "1.5", ""):
        r = run_raw(["--klifs-id", bad, "--klifs-chain", "A", "--klifs-mapping", mapping, "--klifs-output", os.path.join(CFG["tmp"], "y.tsv")])
        check(r.rc == 1 and "Invalid KLIFS structure ID" in r.err, "--klifs-id %r: rc=%d %s" % (bad, r.rc, r.err[:100]))
    for opt in ("--klifs-pdb", "--klifs-chain", "--klifs-model", "--klifs-mapping", "--klifs-output", "--klifs-id"):
        r = run_raw([opt])
        check(r.rc == 1 and "requires a value" in r.err, "%s without a value: %s" % (opt, r.err[:100]))
    # without a mapping the service is needed: it fails cleanly, offline (a proxy that does not exist) or in a build without it
    r = run_raw(["--klifs-pdb", "3EYG", "--klifs-chain", "A", "--klifs-output", os.path.join(CFG["tmp"], "z.tsv")], env=NO_NETWORK, timeout=90)
    check(r.rc == 1 and r.err.startswith("KLIFS: ") and ("unavailable" in r.err or "request failed" in r.err), "KLIFS service not reachable: rc=%d %s" % (r.rc, r.err[:200]))
    r = run_raw(["--klifs-pdb", "3E/G", "--klifs-chain", "A", "--klifs-output", os.path.join(CFG["tmp"], "z.tsv")], env=NO_NETWORK)
    check(r.rc == 1 and "four-character PDB ID" in r.err, "invalid PDB ID: " + r.err)


# ----------------------------------------------------------------------------------------------
# GROUP matrix: the interaction matrix of several ligands (--matrix LIST)
# ----------------------------------------------------------------------------------------------
MATRIX_POCKET = "@name test\n@color gk #ffcc00\nA\t10\t1\thinge\tASN\nA\t30\t2\tgk\tLYS\nA\t77\t3\tgk\tGLY\n"
CLASSES7 = ["hydrophobic", "face_to_face", "edge_to_face", "hbond_donor", "hbond_acceptor", "ionic_positive", "ionic_negative"]


def matrix_structures():
    """Two structures with one ligand each: ASN 10 accepts a hydrogen bond (A); a carboxylate and LYS 30, which is positive and donates (B)."""
    a = fp_system()
    b = Sys().ligand(carboxylate()).prot("LYS", "NZ", at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8), resi=30)
    return tmpfile(a.text()), tmpfile(b.text())


def matrix_list_file(pa, pb, extra_header="", extra_a="", extra_b="", relative=True):
    na, nb = (os.path.basename(pa), os.path.basename(pb)) if relative else (pa, pb)
    return tmpfile("structure\tligand\tlabel%s\n%s\tLIG\talpha%s\n%s\tLIG\tbeta%s\n" % (extra_header, na, extra_a, nb, extra_b), ".tsv")


def matrix_csv(out):
    lines = out.strip().splitlines()
    head = lines[0].split(",")
    rows = {}
    for l in lines[1:]:
        cells = l.split(",")
        rows[cells[0]] = dict(zip(head, cells))
    return head, rows


def bits_of(row):
    return {k: v for k, v in row.items() if ":" in k and v == "1"}


@test("matrix")
def matrix_csv_json_svg_are_consistent():
    pa, pb = matrix_structures()
    lst = matrix_list_file(pa, pb)                                 # relative paths: resolved against the folder of the list
    pocket = tmpfile(MATRIX_POCKET, ".tsv")
    r = run_raw(["--matrix", lst, "--pocket", pocket])
    check(r.rc == 0 and not r.err.strip(), "csv rc=%d %s" % (r.rc, r.err[:200]))
    head, rows = matrix_csv(r.out)
    check(head == ["label", "structure", "ligand"] + ["%s:%s" % (p, c) for p in ("1", "2", "3") for c in CLASSES7], "header: " + ",".join(head[:8]))
    check(list(rows) == ["alpha", "beta"] and os.path.normpath(rows["alpha"]["structure"]) == os.path.normpath(pa), "rows: %s %s" % (list(rows), rows["alpha"]["structure"]))
    check(set(bits_of(rows["alpha"])) == {"1:hbond_acceptor"}, "alpha: " + str(bits_of(rows["alpha"])))
    check(set(bits_of(rows["beta"])) == {"2:hbond_donor", "2:ionic_positive"}, "beta: " + str(bits_of(rows["beta"])))
    check(rows["alpha"]["ligand"] == "LIG 1L" and rows["beta"]["ligand"] == "LIG 1L", "ligand ids")
    j = run_raw(["--matrix", lst, "--pocket", pocket, "--matrix-format", "json"])
    check(j.rc == 0, "json rc=%d" % j.rc)
    d = json.loads(j.out)
    check(d["classes"] == CLASSES7 and d["pocket"] is True and [p["position"] for p in d["positions"]] == ["1", "2", "3"] and
          [p["region"] for p in d["positions"]] == ["hinge", "gk", "gk"], "json structure")
    check([x["label"] for x in d["rows"]] == ["alpha", "beta"] and all(len(x["bits"]) == 21 and len(x["residues"]) == 3 for x in d["rows"]), "json rows")
    for x in d["rows"]:                                           # the bits of the JSON are the cells of the CSV, in the same order
        csv_bits = [int(rows[x["label"]]["%s:%s" % (p, c)]) for p in ("1", "2", "3") for c in CLASSES7]
        check(x["bits"] == csv_bits, "json bits = csv cells for " + x["label"])
    check(d["rows"][0]["residues"][0] == "ASN 10A" and d["rows"][1]["residues"][1] == "LYS 30A", "residues in the json: " + str([x["residues"] for x in d["rows"]]))
    v = run_raw(["--matrix", lst, "--pocket", pocket, "--matrix-format", "svg", "--matrix-title", "A & B <test>"])
    check(v.rc == 0, "svg rc=%d" % v.rc)
    root = ET.fromstring(v.out.encode("utf-8"))                     # well-formed XML
    check(root.tag == SVG_NS + "svg" and root.get("width") and root.get("viewBox"), "svg root")
    for text in ("alpha", "beta", "hinge", "A &amp; B &lt;test&gt;", "H-bond acceptor", "Ionic positive"):
        check(text in v.out, "svg must contain " + text)
    stripes = [e for e in root.iter(SVG_NS + "rect") if e.find(SVG_NS + "title") is not None]
    check(len(stripes) == 3, "one stripe for each bit that is set (3), got %d" % len(stripes))
    check("#ffcc00" in v.out, "the @color of the pocket definition is used for its region")
    # the same matrix, whatever the way it is asked for
    check(run_raw(["--matrix", lst, "--matrix-format", "csv", "--pocket", pocket]).out == r.out, "the default format is csv")


@test("matrix")
def matrix_without_pocket_uses_the_residues_in_contact():
    pa, pb = matrix_structures()
    lst = matrix_list_file(pa, pb)
    r = run_raw(["--matrix", lst])
    check(r.rc == 0 and not r.err.strip(), "rc=%d %s" % (r.rc, r.err[:200]))
    head, rows = matrix_csv(r.out)
    check(head[3:] == ["%s:%s" % (p, c) for p in ("ASN 10A", "LYS 30A") for c in CLASSES7], "the columns are the residues in contact: " + ",".join(head[3:6]))
    check(set(bits_of(rows["alpha"])) == {"ASN 10A:hbond_acceptor"} and set(bits_of(rows["beta"])) == {"LYS 30A:hbond_donor", "LYS 30A:ionic_positive"}, "bits")
    j = json.loads(run_raw(["--matrix", lst, "--matrix-format", "json"]).out)
    check(j["pocket"] is False and all(p["region"] == "" for p in j["positions"]), "no regions without a pocket")
    svg = run_raw(["--matrix", lst, "--matrix-format", "svg"])
    ET.fromstring(svg.out.encode("utf-8"))
    # the classes
    r = run_raw(["--matrix", lst, "--matrix-classes", "all"])
    head, rows = matrix_csv(r.out)
    check(len(head) == 3 + 2 * 11 and head[3 + 7] == "ASN 10A:cation_pi" and head[3 + 10] == "ASN 10A:water_bridge", "all the classes: " + ",".join(head[3:5]))
    r = run_raw(["--matrix", lst, "--matrix-classes", "hbond_acceptor, hbond_donor"])
    head, rows = matrix_csv(r.out)
    check(head[3:] == ["ASN 10A:hbond_donor", "ASN 10A:hbond_acceptor", "LYS 30A:hbond_donor", "LYS 30A:hbond_acceptor"], "a selection, in the order of the classes: " + str(head[3:]))
    r = run_raw(["--matrix", lst, "--matrix-classes", "klifs"])
    check(matrix_csv(r.out)[0][3:10] == ["ASN 10A:" + c for c in CLASSES7], "klifs = the 7 classes")
    r = run_raw(["--matrix", lst, "--matrix-classes", "hydrophobic,nonsense"])
    check(r.rc == 1 and "--matrix-classes: unknown interaction class 'nonsense'" in r.err and not r.out, "unknown class: " + r.err)
    r = run_raw(["--matrix", lst, "--matrix-classes", ","])
    check(r.rc == 1 and "no interaction class" in r.err, "empty selection")


@test("matrix")
def matrix_pockets_per_entry_and_hiding_empty_positions():
    pa, pb = matrix_structures()
    p1 = tmpfile("A\t10\t1\thinge\tASN\nA\t30\t2\tgk\tLYS\nA\t99\t4\tgk\tGLY\n", ".tsv")
    p2 = tmpfile("A\t10\t1\thinge\tASN\nA\t30\t2\tgk\tLYS\nA\t50\t3\tgk\tALA\nA\t99\t4\tgk\tGLY\n", ".tsv")
    for first, second in ((p1, p2), (p2, p1)):                          # the position that only one pocket has goes between 2 and 4, in either order
        lst = matrix_list_file(pa, pb, extra_header="\tpocket", extra_a="\t" + os.path.basename(first), extra_b="\t" + os.path.basename(second))
        r = run_raw(["--matrix", lst])
        check(r.rc == 0, "rc=%d %s" % (r.rc, r.err[:200]))
        head, rows = matrix_csv(r.out)
        positions = [h.split(":")[0] for h in head[3::7]]
        check(positions == ["1", "2", "3", "4"], "the order of the merged positions: %s (%s first)" % (positions, os.path.basename(first)))
        check(set(bits_of(rows["alpha"])) == {"1:hbond_acceptor"} and set(bits_of(rows["beta"])) == {"2:hbond_donor", "2:ionic_positive"}, "bits")
    # the common pocket of --pocket is used by the entries that do not name one
    lst = matrix_list_file(pa, pb, extra_header="\tpocket", extra_a="\t" + os.path.basename(p1), extra_b="")
    r = run_raw(["--matrix", lst, "--pocket", p2])
    check(r.rc == 0 and [h.split(":")[0] for h in matrix_csv(r.out)[0][3::7]] == ["1", "2", "3", "4"], "mixed: " + r.err[:200])
    # positions without any interaction can be left out
    pocket = tmpfile(MATRIX_POCKET, ".tsv")
    lst = matrix_list_file(pa, pb)
    full = matrix_csv(run_raw(["--matrix", lst, "--pocket", pocket]).out)[0]
    small = run_raw(["--matrix", lst, "--pocket", pocket, "--matrix-hide-empty"])
    check(small.rc == 0, "hide-empty rc")
    head, rows = matrix_csv(small.out)
    check([h.split(":")[0] for h in head[3::7]] == ["1", "2"] and len(full) == 3 + 21 and len(head) == 3 + 14, "position 3 has no interaction and is left out")
    check(set(bits_of(rows["beta"])) == {"2:hbond_donor", "2:ionic_positive"}, "the bits are the same")
    d = json.loads(run_raw(["--matrix", lst, "--pocket", pocket, "--matrix-hide-empty", "--matrix-format", "json"]).out)
    check([p["position"] for p in d["positions"]] == ["1", "2"] and all(len(x["bits"]) == 14 and len(x["residues"]) == 2 for x in d["rows"]), "json with hide-empty")
    ET.fromstring(run_raw(["--matrix", lst, "--pocket", pocket, "--matrix-hide-empty", "--matrix-format", "svg"]).out.encode("utf-8"))


@test("matrix")
def matrix_list_file_variants_and_ligand_selection():
    pa, pb = matrix_structures()
    pocket = tmpfile(MATRIX_POCKET, ".tsv")
    na, nb = os.path.basename(pa), os.path.basename(pb)
    base = run_raw(["--matrix", matrix_list_file(pa, pb), "--pocket", pocket]).out
    # CSV, other column names and order, quotes, comments, blank lines, a byte order mark and Windows line ends: the same matrix
    text = ('# a list of ligands\r\n\r\nName,File,Ligand\r\n"alpha",%s,LIG\r\n  # another comment\r\nbeta,"%s",LIG\r\n' % (na, nb))
    for name, content, suffix in (("csv", text, ".csv"), ("bom", "\ufeff" + text, ".csv")):
        r = run_raw(["--matrix", tmpfile(content, suffix), "--pocket", pocket])
        check(r.rc == 0 and r.out == base, "%s: rc=%d %s" % (name, r.rc, r.err[:200]))
    r = run_raw(["--matrix", tmpfile(text.replace(",", "\t"), ".tsv"), "--pocket", pocket])
    check(r.rc == 0 and r.out == base, "tsv with other columns: " + r.err[:200])
    # the entries can name a label that has a comma (quoted)
    r = run_raw(["--matrix", tmpfile('structure,ligand,label\n%s,LIG,"one, two"\n' % na, ".csv"), "--pocket", pocket])
    check(r.rc == 0 and r.out.splitlines()[1].startswith('"one, two",'), "label with a comma is quoted in the csv: " + r.out[-120:])
    # a default label from the file and the ligand
    r = run_raw(["--matrix", tmpfile("structure\tligand\n%s\tLIG\n" % na, ".tsv"), "--pocket", pocket])
    check(r.rc == 0 and r.out.splitlines()[1].startswith(os.path.splitext(na)[0] + " LIG 1L,"), "default label: " + r.out.splitlines()[1][:60])
    # a structure with two ligands of the same name: the chain and the residue number choose
    m = methanol()
    s2 = Sys().ligand(m, chain="L", resi=1)
    s2.atom("HETATM", "C1", "LIG", "M", 2, (30.0, 0.0, 0.0), "C")              # a second methanol, far from the first (close ligands would be merged into one)
    s2.atom("HETATM", "O1", "LIG", "M", 2, (31.43, 0.0, 0.0), "O")
    two = tmpfile(s2.prot("ASN", "OD1", at(m.atoms[1][2], (0.5, 0.866, 0), 2.8)).text())
    n2 = os.path.basename(two)
    r = run_raw(["--matrix", tmpfile("structure\tligand\n%s\tLIG\n" % n2, ".tsv")])
    check(r.rc == 1 and "matches 2 ligands" in r.err and "LIG 1L" in r.err and "LIG 2M" in r.err and not r.out, "ambiguous: " + r.err)
    for extra, wanted in (("chain\tresseq\nM\t\n", "LIG 2M"), ("chain\tresseq\nL\t1\n", "LIG 1L"), ("chain\tresseq\n\t2\n", "LIG 2M"), ("chain\tresseq\n\t2-\n", "LIG 2M"), ("chain\tresseq\n\t1A\n", None)):
        head, body = extra.split("\n")[0], extra.split("\n")[1]
        r = run_raw(["--matrix", tmpfile("structure\tligand\t%s\n%s\tLIG\t%s\n" % (head, n2, body), ".tsv")])
        if wanted:
            check(r.rc == 0 and r.out.splitlines()[1].split(",")[2] == wanted, "%r: %s %s" % (body, r.out[-60:], r.err[:150]))
        else:
            check(r.rc == 1, "a residue number that does not exist: rc=%d" % r.rc)
    # a column that is not known is ignored
    r = run_raw(["--matrix", tmpfile("structure\tligand\tcomment\n%s\tLIG\twhatever\n" % na, ".tsv")])
    check(r.rc == 0, "unknown column")


@test("matrix")
def matrix_errors_name_the_entry_and_write_nothing():
    pa, pb = matrix_structures()
    na, nb = os.path.basename(pa), os.path.basename(pb)
    pocket = tmpfile(MATRIX_POCKET, ".tsv")
    cases = {
        "missing structure file": ("structure\tligand\n%s\tLIG\nnothing_here.pdb\tLIG\n" % na, ["Entry 2 (line 3, nothing_here.pdb)"]),
        "ligand not found": ("structure\tligand\n%s\tXYZ\n" % na, ["ligand 'XYZ' not found", "LIG(1)"]),
        "ligand missing": ("structure\tligand\n%s\t\n" % na, ["the ligand is required", "LIG(1)"]),
        "duplicate label": ("structure\tligand\tlabel\n%s\tLIG\tx\n%s\tLIG\tx\n" % (na, nb), ["used by another entry"]),
        "no header": ("ligand\tchain\nLIG\tA\n", ["Line 1", "'structure' column"]),
        "empty structure": ("structure\tligand\n\tLIG\n", ["Line 2: the structure is empty"]),
        "bad residue number": ("structure\tligand\tresseq\n%s\tLIG\t12AB\n" % na, ["Line 2: invalid residue number '12AB'"]),
        "bad altloc": ("structure\tligand\taltloc\n%s\tLIG\tAB\n" % na, ["altloc expects a single character"]),
        "header only": ("structure\tligand\n", ["no entries"]),
        "only comments": ("# nothing\n", ["empty"]),
        "pocket does not exist": ("structure\tligand\tpocket\n%s\tLIG\tno_such_pocket.tsv\n" % na, ["cannot open the pocket definition"]),
        "pocket for some entries only": ("structure\tligand\tpocket\n%s\tLIG\t%s\n%s\tLIG\t\n" % (na, os.path.basename(pocket), nb), ["Entry 2", "has no pocket definition while others have one"]),
    }
    for name, (text, wanted) in cases.items():
        r = run_raw(["--matrix", tmpfile(text, ".tsv")])
        check(r.rc == 1 and r.err.startswith("Matrix: ") and not r.out, "%s: rc=%d out=%r err=%s" % (name, r.rc, r.out[:40], r.err[:200]))
        for w in wanted:
            check(w in r.err, "%s: the message must contain %r: %s" % (name, w, r.err[:300]))
    r = run_raw(["--matrix", os.path.join(CFG["tmp"], "no_such_list.tsv")])
    check(r.rc == 1 and "Cannot open the matrix list" in r.err, "missing list: " + r.err)
    # the mixed case is also an error when the common pocket is given for some entries only (through the list)
    r = run_raw(["--matrix", tmpfile("structure\tligand\n%s\tLIG\n" % na, ".tsv"), "--pocket", os.path.join(CFG["tmp"], "no_such_pocket.tsv")])
    check(r.rc == 1 and "cannot open the pocket definition" in r.err, "common pocket that does not exist: " + r.err)


@test("matrix")
def matrix_command_line_rules():
    pa, pb = matrix_structures()
    lst = matrix_list_file(pa, pb)
    pocket = tmpfile(MATRIX_POCKET, ".tsv")
    r = run_raw(["--matrix", lst, pa])
    check(r.rc == 1 and "takes no structure argument" in r.err, "positional argument: " + r.err)
    for other in (["--fingerprint", "csv"], ["--csv"], ["--klifs-pdb", "3EYG", "--klifs-chain", "A", "--klifs-mapping", pocket], ["--download-pdb", "1AZM"]):
        r = run_raw(["--matrix", lst] + other)
        check(r.rc == 1 and "--matrix cannot be combined" in r.err and not r.out, "%s: %s" % (other[0], r.err[:150]))
    for opt, value in (("--matrix-format", "csv"), ("--matrix-classes", "klifs"), ("--matrix-title", "t")):
        r = run_raw([opt, value])
        check(r.rc == 1 and "require --matrix LIST" in r.err, "%s without --matrix: %s" % (opt, r.err[:150]))
    r = run_raw(["--matrix-hide-empty"])
    check(r.rc == 1 and "require --matrix LIST" in r.err, "--matrix-hide-empty without --matrix")
    r = run_raw(["--matrix", lst, "--matrix-format", "xml"])
    check(r.rc == 1 and "expects csv, json or svg" in r.err, "format: " + r.err)
    for opt in ("--matrix", "--matrix-format", "--matrix-classes", "--matrix-title"):
        r = run_raw([opt])
        check(r.rc == 1 and "requires a value" in r.err, "%s without a value" % opt)
    # the thresholds of the analysis apply to every row
    base = matrix_csv(run_raw(["--matrix", lst, "--pocket", pocket]).out)[1]
    check(len(bits_of(base["alpha"])) == 1, "alpha has a hydrogen bond")
    none = run_raw(["--matrix", lst, "--pocket", pocket, "--hb-max", "2.0"])
    rows = matrix_csv(none.out)[1] if none.rc == 0 else {}
    check(none.rc == 0 and len(bits_of(rows["alpha"])) == 0 and set(bits_of(rows["beta"])) == {"2:ionic_positive"},
          "--hb-max 2.0 removes the hydrogen bonds, not the salt bridge: %s %s" % (bits_of(rows.get("alpha", {})), bits_of(rows.get("beta", {}))))
    ok = run_raw(["--matrix", lst, "--pocket", pocket, "--plip-like"])
    check(ok.rc == 0, "--plip-like")
    r = run_raw(["--matrix", lst, "--no-water", "--metal-donors", "--backbone-hyd", "--no-merge", "--altloc", "B", "--format", "pdb"])
    check(r.rc == 0, "the other options: " + r.err[:150])


@test("matrix")
def matrix_reports_contacts_outside_the_pocket():
    pa, pb = matrix_structures()
    small = tmpfile("A\t10\t1\thinge\tASN\n", ".tsv")                     # no position for LYS 30
    r = run_raw(["--matrix", matrix_list_file(pa, pb), "--pocket", small])
    check(r.rc == 0, "rc=%d" % r.rc)
    check("Warning: beta: 1 residue(s) in contact are not in the pocket definition and are not in the matrix: LYS 30A" in r.err, "warning: " + r.err)
    head, rows = matrix_csv(r.out)
    check(len(head) == 3 + 7 and len(bits_of(rows["beta"])) == 0 and len(bits_of(rows["alpha"])) == 1, "the contact is not a column")
    # an entry in the list whose pocket warns about a line
    bad = tmpfile("A\t10\t1\thinge\tASN\nthis line is not valid\n", ".tsv")
    r = run_raw(["--matrix", matrix_list_file(pa, pb), "--pocket", bad])
    check(r.rc == 0 and r.err.count("Warning:") >= 1 and os.path.basename(bad) in r.err, "a malformed pocket line is reported once: " + r.err[:200])


@test("matrix")
def matrix_agrees_with_the_single_fingerprints():
    pa, pb = matrix_structures()
    pocket = tmpfile(MATRIX_POCKET, ".tsv")
    head, rows = matrix_csv(run_raw(["--matrix", matrix_list_file(pa, pb), "--pocket", pocket]).out)
    for label, path in (("alpha", pa), ("beta", pb)):
        fp = run_raw([path, "LIG", "--fingerprint", "csv", "--pocket", pocket])
        lines = fp.out.strip().splitlines()
        fhead = lines[0].split(",")
        for l in lines[1:]:
            c = dict(zip(fhead, l.split(",")))
            if not c["position"]:
                continue
            for cls in CLASSES7:
                check(rows[label]["%s:%s" % (c["position"], cls)] == ("1" if int(c[cls]) > 0 else "0"), "%s %s %s" % (label, c["position"], cls))


@test("matrix")
def matrix_on_real_complexes_agrees_with_the_single_fingerprints():
    names = (("1HSG", "MK1"), ("1HVR", "XK2"), ("1HXB", "ROC"))
    root = os.path.abspath(os.path.join(CFG["real_dir"], "..", ".."))
    pocket = os.path.join(root, "docs", "examples", "hiv-protease-pocket.tsv")
    if not all(os.path.exists(os.path.join(CFG["real_dir"], n + ".pdb")) for n, _ in names) or not os.path.exists(pocket):
        raise Skip("the HIV protease complexes or their pocket definition are not in " + root)
    lst = tmpfile("structure\tligand\tlabel\n" + "".join("%s\t%s\t%s\n" % (os.path.join(CFG["real_dir"], n + ".pdb"), lig, n) for n, lig in names), ".tsv")
    head, rows = matrix_csv(run_raw(["--matrix", lst, "--pocket", pocket]).out)
    check(list(rows) == [n for n, _ in names] and len(head) == 3 + 7 * ((len(head) - 3) // 7) and (len(head) - 3) // 7 > 20, "rows and columns: %d" % len(head))
    total = 0
    for n, lig in names:
        fp = run_raw([os.path.join(CFG["real_dir"], n + ".pdb"), lig, "--fingerprint", "csv", "--pocket", pocket])
        fl = fp.out.strip().splitlines()
        fhead = fl[0].split(",")
        expected = {}
        for l in fl[1:]:
            c = dict(zip(fhead, l.split(",")))
            if c["position"]:
                for cls in CLASSES7:
                    expected["%s:%s" % (c["position"], cls)] = "1" if int(c[cls]) > 0 else "0"
        for col, v in expected.items():
            check(rows[n][col] == v, "%s %s: matrix %s, fingerprint %s" % (n, col, rows[n][col], v))
        total += len(bits_of(rows[n]))
    check(total > 30, "the three inhibitors make many contacts: %d" % total)
    check(any(rows[n]["A25:hbond_acceptor"] == "1" for n, _ in names), "the catalytic Asp25 accepts a hydrogen bond in at least one complex")
    # the example list of the documentation gives the same matrix
    example = os.path.join(root, "docs", "examples", "hiv-protease-matrix.tsv")
    if os.path.exists(example):
        r = run_raw(["--matrix", example, "--pocket", pocket])
        check(r.rc == 0 and len(matrix_csv(r.out)[1]) == 3, "the example of the documentation runs: " + r.err[:200])


@test("matrix")
def matrix_robustness():
    pa, pb = matrix_structures()
    pocket = tmpfile(MATRIX_POCKET, ".tsv")
    good = open(matrix_list_file(pa, pb), encoding="utf-8").read()
    rnd = random.Random(33)
    for i in range(60):
        t = list(good)
        for _ in range(6):
            k = rnd.randrange(len(t))
            t[k] = rnd.choice("\t,\"\n#\r \x00\u00e9a1-")
        r = run_raw(["--matrix", tmpfile("".join(t), ".tsv"), "--pocket", pocket])
        no_crash(r)
        check(r.rc in (0, 1), "corrupted list %d: rc=%d" % (i, r.rc))
        check(r.rc == 1 or r.out.startswith("label,structure,ligand"), "a list that is accepted gives a matrix: %d" % i)
    garbage = bytes(random.Random(4).getrandbits(8) for _ in range(4000))
    r = run_raw(["--matrix", tmpfile(garbage, ".tsv", binary=True)])
    no_crash(r)
    check(r.rc == 1, "binary garbage: rc=%d" % r.rc)
    r = run_raw(["--matrix", tmpfile("structure\tligand\tlabel\n%s\tLIG\t%s\n" % (os.path.basename(pa), "x" * 30000), ".tsv")])
    no_crash(r)
    check(r.rc == 0 and ET.fromstring(run_raw(["--matrix", tmpfile("structure\tligand\tlabel\n%s\tLIG\t%s\n" % (os.path.basename(pa), "x" * 30000), ".tsv"), "--matrix-format", "svg"]).out.encode("utf-8")) is not None,
          "a very long label")
    # many entries
    many = "structure\tligand\tlabel\n" + "".join("%s\tLIG\tr%d\n" % (os.path.basename(pa if i % 2 else pb), i) for i in range(150))
    r = run_raw(["--matrix", tmpfile(many, ".tsv"), "--pocket", pocket], timeout=120)
    check(r.rc == 0 and len(r.out.strip().splitlines()) == 151, "150 rows: rc=%d %.1f s" % (r.rc, CFG["last_time"]))
    check(CFG["last_time"] < 60, "150 rows in a reasonable time: %.1f s" % CFG["last_time"])


# ----------------------------------------------------------------------------------------------
# GROUP gpcr: pocket definitions of GPCRs from the generic residue numbers of GPCRdb (--gpcrdb-*)
# ----------------------------------------------------------------------------------------------
def pdb_atom_line(serial, name, res, chain, resi, xyz, b, el):
    return "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s" % (serial, name, res, chain, resi, xyz[0], xyz[1], xyz[2], 1.0, b, el)


def generic_number(r):
    return 3 + (22 + r) / 100.0                         # residue 1 is 3x23, residue 10 is 3x32, residue 30 is 3x52


def numbered_lines(chain="A", count=40, names=None, annotated=True, bw=False):
    """CA (GPCRdb number) and N (Ballesteros-Weinstein number or 1.00) atoms of `count` residues, as GPCRdb writes them; or B-factors of a structure."""
    out, serial = [], 5000
    for r in range(1, count + 1):
        res = (names or {}).get(r, "GLY")
        b = generic_number(r) if annotated else 40.0 + (r * 37 % 53) / 3.0
        out.append(pdb_atom_line(serial, " N  ", res, chain, r, (r * 1.5, 30.0, 30.0), b if (bw and annotated) else 1.0, " N"))
        out.append(pdb_atom_line(serial + 1, " CA ", res, chain, r, (r * 1.5 + 1.0, 30.0, 30.0), b, " C"))
        serial += 2
    return "\n".join(out) + "\n"


def add_lines(text, lines):
    check(text.rstrip().endswith("END"), "the system must end with END")
    return text.rstrip()[:-3] + lines + "END\n"


def gpcr_structures(annotated=True, bw=False):
    """(A) ASN 10 (3x32) accepts a hydrogen bond from the ligand; (B) LYS 30 (3x52) is positive and donates to a carboxylate."""
    a = fp_system().text()
    b = Sys().ligand(carboxylate()).prot("LYS", "NZ", at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8), resi=30).text()
    return (add_lines(a, numbered_lines(names={10: "ASN"}, annotated=annotated, bw=bw)), add_lines(b, numbered_lines(names={30: "LYS"}, annotated=annotated, bw=bw)))


def gpcr_cli(pdb, positions, out=None, extra=()):
    args = ["--gpcrdb-pdb", pdb, "--gpcrdb-positions", positions] + (["--gpcrdb-output", out] if out else []) + list(extra)
    return run_raw(args)


def pocket_lines(path):
    return [l.split("\t") for l in open(path, encoding="utf-8").read().splitlines() if l and l[0] not in "#@"]


@test("gpcr")
def gpcr_pocket_from_the_annotated_structure():
    a, b = gpcr_structures()
    pa = tmpfile(a)
    out = os.path.join(CFG["tmp"], "gpcr_pocket.tsv")
    r = gpcr_cli(pa, "3x40,3x32, 3X25 ,3x32", out)                       # any order, blanks, upper case, a repeated number
    check(r.rc == 0 and "GPCRdb pocket: 3 positions" in r.err and not r.out, "rc=%d %s" % (r.rc, r.err[:200]))
    text = open(out, encoding="utf-8").read()
    check(text.startswith("@name GPCR pocket (generic residue numbers, chain A)\n@color TM3 #") and text.count("@color") == 1, "the header: " + text[:120])
    check(pocket_lines(out) == [["A", "3", "3x25", "TM3", "GLY"], ["A", "10", "3x32", "TM3", "ASN"], ["A", "18", "3x40", "TM3", "GLY"]], str(pocket_lines(out)))
    # the positions of a file with comments, and the pocket from the structure of the analysis itself
    plist = tmpfile("# the pocket of interest\n3x32\n3x40  # a comment\n", ".txt")
    r2 = run_raw(["--gpcrdb-pdb", pa, "--gpcrdb-positions", plist, "--gpcrdb-output", os.path.join(CFG["tmp"], "gpcr_pocket2.tsv")])
    check(r2.rc == 0 and "2 positions" in r2.err, "positions from a file: " + r2.err)
    # the pocket file is the pocket definition of the fingerprint
    fp = run_raw([pa, "LIG", "--fingerprint", "csv", "--pocket", out])
    rows = {l.split(",")[4]: l.split(",") for l in fp.out.strip().splitlines()[1:]}
    head = fp.out.splitlines()[0].split(",")
    check(fp.rc == 0 and list(rows) == ["3x25", "3x32", "3x40"] and rows["3x32"][head.index("hbond_acceptor")] == "1" and rows["3x32"][5] == "TM3", "fingerprint: " + fp.out[:300] + fp.err[:200])


@test("gpcr")
def gpcr_fingerprint_and_matrix_with_generic_numbers():
    a, b = gpcr_structures()
    pa, pb = tmpfile(a), tmpfile(b)
    pocket = os.path.join(CFG["tmp"], "gpcr_ref.tsv")
    check(gpcr_cli(pa, "3x25,3x32,3x40,3x52", pocket).rc == 0, "generate")
    ref = run_raw([pa, "LIG", "--fingerprint", "csv", "--pocket", pocket]).out
    # the same in one command: the structure of the analysis is the annotated file
    one = run_raw([pa, "LIG", "--fingerprint", "csv", "--gpcrdb-positions", "3x25,3x32,3x40,3x52"])
    check(one.rc == 0 and one.out == ref and "GPCRdb pocket: 4 positions" in one.err, "one command: " + one.err[:200])
    # a structure without numbers (the original of the PDB) and the file annotated by GPCRdb apart
    plain, _ = gpcr_structures(annotated=False)
    original = tmpfile(fp_system().text())
    sep = run_raw([original, "LIG", "--fingerprint", "csv", "--gpcrdb-pdb", pa, "--gpcrdb-positions", "3x25,3x32,3x40,3x52"])
    check(sep.rc == 0 and sep.out == ref, "annotated file apart: " + sep.err[:200])
    for fmt in ("json", "svg"):
        r = run_raw([pa, "LIG", "--fingerprint", fmt, "--gpcrdb-positions", "3x25,3x32"])
        check(r.rc == 0, fmt + ": " + r.err[:150])
        if fmt == "json":
            json.loads(r.out)
        else:
            ET.fromstring(r.out.encode("utf-8"))
    # a matrix of two receptors on the same generic positions
    pocket_a, pocket_b = os.path.join(CFG["tmp"], "gpcr_pa.tsv"), os.path.join(CFG["tmp"], "gpcr_pb.tsv")
    check(gpcr_cli(pa, "3x32,3x52", pocket_a).rc == 0 and gpcr_cli(pb, "3x32,3x52", pocket_b).rc == 0, "pockets")
    lst = tmpfile("structure\tligand\tlabel\tpocket\n%s\tLIG\talpha\t%s\n%s\tLIG\tbeta\t%s\n" % (os.path.basename(pa), os.path.basename(pocket_a), os.path.basename(pb), os.path.basename(pocket_b)), ".tsv")
    m = run_raw(["--matrix", lst])
    check(m.rc == 0, "matrix: " + m.err[:200])
    head, rows = matrix_csv(m.out)
    check([h.split(":")[0] for h in head[3::7]] == ["3x32", "3x52"], "the columns are the generic numbers: " + str(head[3::7]))
    check(set(bits_of(rows["alpha"])) == {"3x32:hbond_acceptor"} and set(bits_of(rows["beta"])) == {"3x52:hbond_donor", "3x52:ionic_positive"}, "bits: %s %s" % (bits_of(rows["alpha"]), bits_of(rows["beta"])))
    svg = run_raw(["--matrix", lst, "--matrix-format", "svg"])
    check(svg.rc == 0 and "TM3" in svg.out, "the segment is the region of the matrix")


@test("gpcr")
def gpcr_table_gives_the_same_pocket_as_the_annotated_structure():
    a, _ = gpcr_structures()
    pa = tmpfile(a)
    out1, out2 = os.path.join(CFG["tmp"], "gpcr_t1.tsv"), os.path.join(CFG["tmp"], "gpcr_t2.tsv")
    check(gpcr_cli(pa, "3x25,3x32,3x40", out1, ["--gpcrdb-chain", "A"]).rc == 0, "from the structure")
    rows = ["chain,resi,gpcrdb,aa"] + ["A,%d,%dx%02d,%s" % (r, 3, 22 + r, "N" if r == 10 else "G") for r in range(1, 41)] + ["A,100,-,G", "A,101,,K"]
    table = tmpfile("\n".join(rows) + "\n", ".csv")
    r = run_raw(["--gpcrdb-mapping", table, "--gpcrdb-chain", "A", "--gpcrdb-positions", "3x25,3x32,3x40", "--gpcrdb-output", out2])
    check(r.rc == 0, "from the table: " + r.err[:200])
    check(open(out1, encoding="utf-8").read() == open(out2, encoding="utf-8").read(), "the pocket is the same from the structure and from the table:\n%s\n%s" % (open(out1).read(), open(out2).read()))
    # the numbers of GPCRdb written together, Ballesteros-Weinstein numbers, a segment and an insertion code
    t2 = tmpfile("resi\tgpcrdb\tbw\tsegment\n113\t3.32x32\t\tTM3\n203\t\t5.42\tTM5\n316A\t7x42\t\t\n", ".tsv")
    out3 = os.path.join(CFG["tmp"], "gpcr_t3.tsv")
    r = run_raw(["--gpcrdb-mapping", t2, "--gpcrdb-positions", "3.32,3x32,5.42,7x42", "--gpcrdb-output", out3])
    check(r.rc == 0, "table with several forms: " + r.err[:200])
    # the same position asked in two notations keeps the order of the request
    check(pocket_lines(out3) == [["", "113", "3.32", "TM3", ""], ["", "113", "3x32", "TM3", ""], ["", "203", "5.42", "TM5", ""], ["", "316A", "7x42", "TM7", ""]], str(pocket_lines(out3)))


@test("gpcr")
def gpcr_errors_and_rules():
    a, b = gpcr_structures()
    pa = tmpfile(a)
    out = os.path.join(CFG["tmp"], "gpcr_never.tsv")
    def fails(args, *texts):
        r = run_raw(args)
        check(r.rc == 1 and r.err.startswith(("GPCRdb: ", "--", "Choose", "Warning")) and not r.out, "%s: rc=%d %s" % (args[:4], r.rc, r.err[:250]))
        for t in texts:
            check(t in r.err, "%s must contain %r: %s" % (args[:4], t, r.err[:300]))
        return r
    # Ballesteros-Weinstein numbers need the N atoms; the numbers with an x do not; they are not the same
    fails(["--gpcrdb-pdb", pa, "--gpcrdb-positions", "3.32", "--gpcrdb-output", out], "Ballesteros-Weinstein numbers", "5.42 and 5x42")
    a_bw, _ = gpcr_structures(bw=True)
    pbw = tmpfile(a_bw)
    r = run_raw(["--gpcrdb-pdb", pbw, "--gpcrdb-positions", "3.32,3x32", "--gpcrdb-output", out])
    check(r.rc == 0 and [l[2] for l in pocket_lines(out)] == ["3.32", "3x32"] and all(l[1] == "10" for l in pocket_lines(out)), "BW and GPCRdb numbers: " + r.err[:200] + str(pocket_lines(out)))
    # a structure that is not annotated: the B-factors of the structure are not numbers
    plain, _ = gpcr_structures(annotated=False)
    fails(["--gpcrdb-pdb", tmpfile(plain), "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "does not have GPCRdb generic numbers", "generic_numbering_index")
    fails(["--gpcrdb-pdb", tmpfile(fp_system().text()), "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "no protein residues with CA atoms")
    # chains
    two = tmpfile(add_lines(a, numbered_lines(chain="B", names={10: "ASN"})))
    fails(["--gpcrdb-pdb", two, "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "Several chains have generic numbers (A, B)", "--gpcrdb-chain")
    r = run_raw(["--gpcrdb-pdb", two, "--gpcrdb-chain", "B", "--gpcrdb-positions", "3x32", "--gpcrdb-output", out])
    check(r.rc == 0 and pocket_lines(out)[0][0] == "B", "chain B: " + r.err[:200])
    fails(["--gpcrdb-pdb", pa, "--gpcrdb-chain", "Z", "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "Chain Z is not in the structure")
    # positions that the structure does not have: a warning; none of them: an error
    r = run_raw(["--gpcrdb-pdb", pa, "--gpcrdb-positions", "3x32,3x99,6x48", "--gpcrdb-output", out])
    check(r.rc == 0 and "Warning: positions that this structure does not have" in r.err and "3x99, 6x48" in r.err and "1 positions" in r.err, "missing: " + r.err)
    fails(["--gpcrdb-pdb", pa, "--gpcrdb-positions", "3x99,6x48", "--gpcrdb-output", out], "None of the requested positions (3x99, 6x48)")
    # the list of positions
    fails(["--gpcrdb-pdb", pa, "--gpcrdb-positions", "3x32,hello", "--gpcrdb-output", out], "Invalid generic residue number 'hello'")
    fails(["--gpcrdb-pdb", pa, "--gpcrdb-positions", "9x50", "--gpcrdb-output", out], "Invalid generic residue number '9x50'")
    fails(["--gpcrdb-pdb", pa, "--gpcrdb-output", out], "--gpcrdb-positions is required")
    fails(["--gpcrdb-pdb", pa, "--gpcrdb-positions", " ", "--gpcrdb-output", out], "No generic residue numbers")
    fails(["--gpcrdb-pdb", os.path.join(CFG["tmp"], "no_such.pdb"), "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "Cannot open")
    fails(["--gpcrdb-mapping", tmpfile("name\tvalue\nx\ty\n", ".tsv"), "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "needs a header with a residue number column")
    fails(["--gpcrdb-mapping", tmpfile("resi\tgpcrdb\n12AB\t3x32\n", ".tsv"), "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "invalid residue number '12AB'")
    fails(["--gpcrdb-mapping", os.path.join(CFG["tmp"], "no_such_table.tsv"), "--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "Cannot open")
    # the rules of the command line
    fails(["--gpcrdb-positions", "3x32", "--gpcrdb-pdb", pa], "needs --gpcrdb-output")
    fails(["--gpcrdb-positions", "3x32", "--gpcrdb-output", out], "--gpcrdb-pdb FILE")
    fails([pa, "LIG", "--pocket", out, "--fingerprint", "csv", "--gpcrdb-positions", "3x32"], "either --pocket or --gpcrdb")
    fails([pa, "LIG", "--fingerprint", "csv", "--gpcrdb-positions", "3x32", "--klifs-pdb", "3EYG", "--klifs-chain", "A"], "either --klifs-* or --gpcrdb-*")
    lst = tmpfile("structure\tligand\n%s\tLIG\n" % os.path.basename(pa), ".tsv")
    fails(["--matrix", lst, "--gpcrdb-positions", "3x32"], "--matrix cannot be combined")
    for opt in ("--gpcrdb-positions", "--gpcrdb-pdb", "--gpcrdb-mapping", "--gpcrdb-chain", "--gpcrdb-name", "--gpcrdb-output"):
        r = run_raw([opt])
        check(r.rc == 1 and "requires a value" in r.err, "%s without a value" % opt)
    # without --fingerprint the pocket is not used: a warning, and the analysis goes on
    r = run_raw([pa, "LIG", "--gpcrdb-positions", "3x32"])
    check(r.rc == 0 and "only used with --fingerprint" in r.err, "warning: " + r.err[:200])
    # a name for the pocket
    r = run_raw(["--gpcrdb-pdb", pa, "--gpcrdb-positions", "3x32", "--gpcrdb-name", "Orthosteric site", "--gpcrdb-output", out])
    check(r.rc == 0 and open(out, encoding="utf-8").read().startswith("@name Orthosteric site\n"), "name")


@test("gpcr")
def gpcr_bulges_loops_and_helix_8_of_a_numbered_structure():
    # numbers of a structure as GPCRdb writes them: a repeated number is a bulge; a negative B-factor is a bulge too; loops are 12, 23, 34, 45, 56, 67; helix 8 is 8
    lines, serial = [], 1
    spec = [(2, 2.53, "ALA"), (3, 2.54, "GLY"), (4, 2.55, "LEU"), (5, 2.55, "ALA"), (6, 2.56, "VAL"), (7, 2.57, "VAL"), (8, 2.58, "PRO"), (9, 3.50, "ARG"), (10, 3.51, "TYR"),
            (11, 34.49, "PRO"), (12, 34.50, "PHE"), (13, 45.50, "CYS"), (14, 56.50, "SER"), (15, 67.50, "GLU"), (16, 23.50, "TRP"), (17, 12.50, "LEU"), (18, 8.47, "SER"), (19, 6.48, "TRP"),
            (20, 6.49, "ALA"), (21, 6.49, "ALA"), (22, 6.50, "PRO"), (23, 6.51, "PHE"), (24, 6.52, "PHE")]
    for i in range(1, 31):
        spec.append((100 + i, 1.0 + (19 + i) / 100.0, "GLY"))                  # a helix of 30 residues so that the structure looks like an annotated one
    spec.append((200, -2.57, "SER"))                                            # a bulge as GPCRdb writes it in its own output: after 2.57
    for resi, b, res in spec:
        lines.append(pdb_atom_line(serial, " N  ", res, "A", resi, (resi * 1.5, 0, 0), 1.0, " N"))
        lines.append(pdb_atom_line(serial + 1, " CA ", res, "A", resi, (resi * 1.5 + 1, 0, 0), b, " C"))
        serial += 2
    pdb = tmpfile("\n".join(lines) + "\nEND\n")
    out = os.path.join(CFG["tmp"], "gpcr_bulge.tsv")
    r = gpcr_cli(pdb, "2x55,2x551,2x56,2x571,3x50,34x50,45x50,56x50,67x50,23x50,12x50,8x47,6x49,6x491", out)
    check(r.rc == 0, "rc=%d %s" % (r.rc, r.err[:300]))
    got = {l[2]: (l[3], l[1]) for l in pocket_lines(out)}
    check(got["2x55"] == ("TM2", "4") and got["2x551"] == ("TM2", "5"), "a repeated number: the second residue is the bulge 2x551: " + str(got))
    check(got["2x571"] == ("TM2", "200"), "a negative B-factor is a bulge: the number of the preceding residue, 2x57, followed by 1: " + str(got))
    check(got["34x50"] == ("ICL2", "12") and got["45x50"] == ("ECL2", "13") and got["56x50"] == ("ICL3", "14") and got["67x50"] == ("ECL3", "15") and got["23x50"] == ("ECL1", "16") and
          got["12x50"] == ("ICL1", "17") and got["8x47"] == ("H8", "18"), "loops and helix 8: " + str(got))
    check(got["6x49"] == ("TM6", "20") and got["6x491"] == ("TM6", "21"), "another repeated number: " + str(got))
    order = [l[2] for l in pocket_lines(out)]
    check(order.index("12x50") < order.index("2x55") < order.index("2x551") < order.index("2x56") < order.index("2x571") < order.index("23x50") < order.index("3x50") < order.index("34x50") <
          order.index("45x50") < order.index("56x50") < order.index("6x49") < order.index("6x491") < order.index("67x50") < order.index("8x47"),
          "in the order of the receptor (ICL1, TM2, ECL1, TM3, ICL2, ECL2, ICL3, TM6, ECL3, H8): " + str(order))


@test("gpcr")
def gpcr_robustness():
    a, _ = gpcr_structures()
    rnd = random.Random(52)
    for i in range(50):
        t = list(a)
        for _ in range(10):
            k = rnd.randrange(len(t))
            t[k] = rnd.choice(" -.0123456789ACN\n\x00\u00e9")
        r = run_raw(["--gpcrdb-pdb", tmpfile("".join(t)), "--gpcrdb-positions", "3x25,3x32,3x40", "--gpcrdb-output", os.path.join(CFG["tmp"], "gpcr_fuzz.tsv")])
        no_crash(r)
        check(r.rc in (0, 1), "corrupted structure %d: rc=%d" % (i, r.rc))
    table = "chain,resi,gpcrdb,aa\n" + "".join("A,%d,3x%02d,G\n" % (r, 22 + r) for r in range(1, 30))
    for i in range(50):
        t = list(table)
        for _ in range(8):
            k = rnd.randrange(len(t))
            t[k] = rnd.choice(",\t\n\"x.0123456789A-#\x00")
        r = run_raw(["--gpcrdb-mapping", tmpfile("".join(t), ".csv"), "--gpcrdb-positions", "3x25,3x32", "--gpcrdb-output", os.path.join(CFG["tmp"], "gpcr_fuzz.tsv")])
        no_crash(r)
        check(r.rc in (0, 1), "corrupted table %d: rc=%d" % (i, r.rc))
    garbage = tmpfile(bytes(random.Random(8).getrandbits(8) for _ in range(5000)), ".pdb", binary=True)
    r = run_raw(["--gpcrdb-pdb", garbage, "--gpcrdb-positions", "3x32", "--gpcrdb-output", os.path.join(CFG["tmp"], "gpcr_fuzz.tsv")])
    no_crash(r)
    check(r.rc == 1, "binary garbage: rc=%d" % r.rc)


def mol2_pocket(residues, eol="\n"):
    """A pocket.mol2 as KLIFS writes it: residues = [(resname, number)]."""
    out = ["@<TRIPOS>MOLECULE", "HUMAN/TEST_1abc_altB_chainA", "%d 0 %d 0 0 " % (2 * len(residues), len(residues)), "BIOPOLYMER", "USER_CHARGES", "", "", "@<TRIPOS>ATOM"]
    for i, (n, num) in enumerate(residues):
        out.append("   %d N        6.2495    17.2789    53.5862 N.3    %d %s%s  -0.6701 BACKBONE" % (2 * i + 1, i + 1, n, num))
    out += ["@<TRIPOS>BOND", "@<TRIPOS>SUBSTRUCTURE"]
    for i, (n, num) in enumerate(residues):
        out.append("  %d %s%s    %d RESIDUE 4 A %s 1" % (i + 1, n, num, 2 * i + 1, n))
    out.append("# MOE 2012.10 (io_trps.svl 2012.10)")
    return eol.join(out) + eol


def pocket_sequence(placed):
    """The 85-character pocket sequence: placed = {position: letter}."""
    return "".join(placed.get(p, "_") for p in range(1, 86))


@test("net")
def klifs_pocket_mol2_gives_the_fingerprint_of_the_kinase_pocket():
    s = fp_system()                                                      # ligand + ALA 9, ASN 10 (accepts a hydrogen bond), GLY 11 of chain A
    mol2 = tmpfile(mol2_pocket([("ALA", 9), ("ASN", 10), ("GLY", 11)]), ".mol2")
    seq = pocket_sequence({1: "A", 20: "N", 46: "G"})                    # ASN 10 sits in the alpha C helix, GLY 11 in the hinge
    out = os.path.join(CFG["tmp"], "mol2_pocket.tsv")
    r = klifs_cli(mol2, out, extra=["--klifs-sequence", seq])
    check(r.rc == 0 and "3 mapped positions" in r.err, "rc=%d %s" % (r.rc, r.err[:200]))
    lines = [l.split("\t") for l in open(out, encoding="utf-8").read().splitlines() if l and l[0] not in "#@"]
    check([(l[1], l[2], l[3], l[4]) for l in lines] == [("9", "1", "I", "ALA"), ("10", "20", "\u03b1C", "ASN"), ("11", "46", "hinge", "GLY")], str(lines))
    check("@color \u03b1C #" in open(out, encoding="utf-8").read(), "the alpha C region must have its own color")
    pdb = tmpfile(s.text())
    fp = run_raw([pdb, "LIG", "--fingerprint", "csv", "--pocket", out])
    rows = [l.split(",") for l in fp.out.strip().splitlines()[1:]]
    check(fp.rc == 0 and [(x[1], x[4], x[5]) for x in rows] == [("ALA 9A", "1", "I"), ("ASN 10A", "20", "\u03b1C"), ("GLY 11A", "46", "hinge")], fp.out + fp.err)
    check(rows[1][10] == "1", "the hydrogen bond of ASN 10 is in the position 20: " + str(rows[1]))
    svg = run_raw([pdb, "LIG", "--fingerprint", "svg", "--pocket", out])
    import xml.etree.ElementTree as ET
    ET.fromstring(svg.out.encode("utf-8"))                               # well-formed
    check("\u03b1C" in svg.out, "the Greek letter must reach the figure")
    js = json.loads(run_raw([pdb, "LIG", "--fingerprint", "json", "--pocket", out]).out)
    check(js, "json")
    # the one-command form and the Windows line ends give the same
    mol2_crlf = tmpfile(mol2_pocket([("ALA", 9), ("ASN", 10), ("GLY", 11)], "\r\n"), ".mol2")
    r2 = run_raw([pdb, "LIG", "--fingerprint", "csv", "--klifs-pdb", "1ABC", "--klifs-chain", "A", "--klifs-mapping", mol2_crlf, "--klifs-sequence", seq])
    check(r2.rc == 0 and r2.out == fp.out, "one command, CRLF: rc=%d %s" % (r2.rc, r2.err[:200]))


@test("net")
def klifs_pocket_mol2_without_sequence_needs_all_85_residues():
    out = os.path.join(CFG["tmp"], "mol2_all.tsv")
    all85 = [("ALA", 100 + p) for p in range(1, 86)]
    r = klifs_cli(tmpfile(mol2_pocket(all85), ".mol2"), out)
    check(r.rc == 0 and "85 mapped positions" in r.err, "85 residues, no sequence: " + r.err)
    lines = [l.split("\t") for l in open(out, encoding="utf-8").read().splitlines() if l and l[0] not in "#@"]
    check([l[3] for l in lines[:4]] == ["I", "I", "I", "g.l"] and lines[44][3] == "GK" and lines[79][3] == "xDFG" and lines[84][3] == "a.l", "regions of the 85 positions")
    regions = [l[3] for l in lines]
    expected = {"I": 3, "g.l": 6, "II": 4, "III": 6, "\u03b1C": 11, "b.l": 7, "IV": 4, "V": 3, "GK": 1, "hinge": 3, "linker": 4, "\u03b1D": 7, "\u03b1E": 5,
                "VI": 3, "c.l": 8, "VII": 3, "VIII": 1, "xDFG": 4, "a.l": 2}
    check({k: regions.count(k) for k in expected} == expected and len(set(regions)) == 19, "the 19 KLIFS regions with their sizes")
    r = klifs_cli(tmpfile(mol2_pocket(all85[:81]), ".mol2"), os.path.join(CFG["tmp"], "never_mol2.tsv"))
    check(r.rc == 1 and "81 residues" in r.err and "--klifs-sequence" in r.err, "81 residues need the sequence: " + r.err)
    check(not os.path.exists(os.path.join(CFG["tmp"], "never_mol2.tsv")), "nothing written")


@test("net")
def klifs_pocket_mol2_errors_and_corruption():
    res = [("ALA", 9), ("ASN", 10), ("GLY", 11)]
    good = tmpfile(mol2_pocket(res), ".mol2")
    seq = pocket_sequence({1: "A", 20: "N", 46: "G"})
    def run(mol2, sequence, name):
        out = os.path.join(CFG["tmp"], "never_%s.tsv" % name)
        if os.path.exists(out):
            os.remove(out)                                   # a file left by an earlier run that succeeded does not count
        r = klifs_cli(mol2, out, extra=["--klifs-sequence", sequence] if sequence is not None else [])
        check(not os.path.exists(out) or r.rc == 0, "%s: no file on error" % name)
        no_crash(r)
        return r
    for name, mol2, sequence, text in (
            ("seq too short", good, seq[1:], "must have 85 characters"), ("seq too long", good, seq + "A", "must have 85 characters"),
            ("seq invalid char", good, seq[:-1] + "1", "Invalid character"), ("seq residues differ in number", good, pocket_sequence({1: "A", 20: "N"}), "has 2 residues but the pocket MOL2 has 3"),
            ("residue differs from its letter", good, pocket_sequence({1: "A", 20: "W", 46: "G"}), "Position 20"),
            ("no sequence, 3 residues", good, None, "pocket sequence"),
            ("no substructures", tmpfile("@<TRIPOS>MOLECULE\nx\n1 0 0 0 0\n@<TRIPOS>ATOM\n", ".mol2"), "", "no residues"),
            ("short substructure line", tmpfile("@<TRIPOS>MOLECULE\n@<TRIPOS>SUBSTRUCTURE\n 1 ALA1 1 RESIDUE\n", ".mol2"), "", "Invalid substructure line"),
            ("bad residue number", tmpfile("@<TRIPOS>MOLECULE\n@<TRIPOS>SUBSTRUCTURE\n 1 ALAxyz 1 RESIDUE 4 A ALA 1\n", ".mol2"), "", "Invalid residue"),
            ("86 residues", tmpfile(mol2_pocket([("ALA", i) for i in range(86)]), ".mol2"), "", "at most 85")):
        r = run(mol2, sequence, name.replace(" ", "_"))
        check(r.rc == 1 and r.err.startswith("KLIFS: ") and text in r.err, "%s: rc=%d %s" % (name, r.rc, r.err[:200]))
    r = run_raw(["--klifs-sequence"])
    check(r.rc == 1 and "requires a value" in r.err, "--klifs-sequence without a value")
    # modified residues, insertion codes and negative numbers are accepted
    odd = tmpfile(mol2_pocket([("ALA", "9A"), ("ASN", -3), ("MSE", 11)]), ".mol2")
    r = klifs_cli(odd, os.path.join(CFG["tmp"], "mol2_odd.tsv"), extra=["--klifs-sequence", pocket_sequence({1: "A", 20: "N", 46: "M"})])
    check(r.rc == 0, "odd residues: " + r.err[:200])
    check("A\t9A\t1\tI\tALA" in open(os.path.join(CFG["tmp"], "mol2_odd.tsv"), encoding="utf-8").read(), "insertion code kept")
    # a random corruption never crashes
    rnd = random.Random(21)
    base = open(good, encoding="utf-8").read()
    for i in range(60):
        t = list(base)
        for _ in range(8):
            k = rnd.randrange(len(t))
            t[k] = rnd.choice("@<>\n\t _#0A9-.\x00\u00e9")
        r = run(tmpfile("".join(t), ".mol2"), seq, "fuzz")
        check(r.rc in (0, 1), "corrupted MOL2 %d: rc=%d" % (i, r.rc))


def rcsb_cache():
    d = os.path.join(CFG["tmp"], "rcsb_cache_%d" % random.randrange(10**9))
    os.makedirs(d)
    shutil.copy(os.path.join(CFG["real_dir"], "1AZM.pdb"), os.path.join(d, "1AZM.pdb")) if os.path.exists(os.path.join(CFG["real_dir"], "1AZM.pdb")) else \
        open(os.path.join(d, "1AZM.pdb"), "w").write(fp_system().text() + "\n" * 5)
    return d


@test("net")
def rcsb_cache_ligand_arguments_and_validation():
    cache = rcsb_cache()
    base = ["--download-pdb", "1azm", "--download-cache", cache]
    r = run_raw(base)
    check(r.rc == 0 and "RCSB structure: " in r.err and "Ligands (HETATM) present:" in r.out, "cached structure, lower-case ID: " + r.out + r.err)
    path = run_raw(base + ["--download-only"]).out.strip()
    check(os.path.basename(path) == "1AZM.pdb" and os.path.exists(path), "--download-only prints the path: " + path)
    names = re.findall(r"(\S+?)\(\d+\)", r.out.splitlines()[0])
    lig = "LIG" if "LIG" in names else "AZM"
    # the ligand can be named in the same command (RESNAME [chain [resseq]])
    direct = run_raw(["--download-pdb", "1AZM", "--download-cache", cache, lig, "--csv"])
    separate = run_raw([path, lig, "--csv"])
    check(direct.rc == 0 and direct.out == separate.out and direct.out.count("\n") > 1, "download + ligand in one command: rc=%d %s" % (direct.rc, direct.err[:200]))
    check(run_raw(["--download-pdb", "1AZM", "--download-cache", cache, lig, "A", "--csv"]).out.count("\n") >= 1, "with a chain")
    r = run_raw(base + [path])
    check(r.rc == 1 and "without a separate structure path" in r.err, "a structure path together with --download-pdb: " + r.err)
    for bad in ("../x", "1AZ", "1AZMM", "1A_M", "1AZ\u00e9", "", "1AZM/../.."):         # (four letters or digits are a valid ID: they would be downloaded)
        r = run_raw(["--download-pdb", bad, "--download-cache", cache])
        check(r.rc == 1 and "four-character alphanumeric PDB ID" in r.err, "ID %r: rc=%d %s" % (bad, r.rc, r.err[:100]))
        check(sorted(os.listdir(cache)) == ["1AZM.pdb"], "%r: nothing may be created: %s" % (bad, os.listdir(cache)))
    r = run_raw(base + ["--download-format", "xyz"])
    check(r.rc == 1 and "format must be pdb or cif" in r.err, "format: " + r.err)
    check(run_raw(["--download-only"]).rc == 1, "--download-only alone")
    for opt in ("--download-pdb", "--download-format", "--download-output", "--download-cache"):
        r = run_raw([opt])
        check(r.rc == 1 and "requires a value" in r.err, "%s without a value" % opt)


@test("net")
def rcsb_failed_transfers_leave_the_cache_untouched():
    cache = rcsb_cache()
    cached = os.path.join(cache, "1AZM.pdb")
    before = open(cached, "rb").read()
    # a refresh, a format that is not cached, and a corrupt cache entry all need the network: here it is unreachable
    for args, what in ((["--download-refresh"], "refresh"), (["--download-format", "cif"], "cif not cached")):
        r = run_raw(["--download-pdb", "1AZM", "--download-cache", cache] + args, env=NO_NETWORK, timeout=90)
        check(r.rc == 1 and r.err.count("RCSB: ") >= 1 and ("unavailable" in r.err or "failed" in r.err), "%s: rc=%d %s" % (what, r.rc, r.err[:200]))
        check(open(cached, "rb").read() == before, "%s: the cached file must not change" % what)
        check(sorted(os.listdir(cache)) == ["1AZM.pdb"], "%s: no partial file may remain: %s" % (what, os.listdir(cache)))
    for name, content in (("2ABC.pdb", b"short"), ("3ABC.pdb", b""), ("4ABC.pdb", b"x" * 500)):       # too short / empty / not a structure
        open(os.path.join(cache, name), "wb").write(content)
        r = run_raw(["--download-pdb", name[:4], "--download-cache", cache], env=NO_NETWORK, timeout=90)
        check(r.rc == 1 and "RCSB: " in r.err, "invalid cache entry %s must not be trusted: rc=%d %s" % (name, r.rc, r.err[:150]))
        check(open(os.path.join(cache, name), "rb").read() == content, "%s: an invalid entry is left as it was" % name)
    # --download-output chooses the file; a valid file there is used as it is
    out = os.path.join(cache, "chosen.pdb")
    shutil.copy(cached, out)
    r = run_raw(["--download-pdb", "1AZM", "--download-output", out, "--download-only"])
    check(r.rc == 0 and r.out.strip() == out, "--download-output: " + r.out + r.err)


# ----------------------------------------------------------------------------------------------
# GROUP robust: malformed input
# ----------------------------------------------------------------------------------------------
def base_system():
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", p)
    return s


@test("robust")
def missing_pdb_file():
    r = run_raw([os.path.join(CFG["tmp"], "non_esiste.pdb")])
    check(r.rc == 1 and "Cannot open" in r.err, "rc=%d err=%s" % (r.rc, r.err))


@test("robust")
def no_arguments_prints_usage():
    r = run_raw([])
    check(r.rc == 1 and "Usage" in r.err, "rc=%d" % r.rc)


@test("robust")
def empty_and_whitespace_files():
    for content in ["", "\n\n\n", "   \n"]:
        r = run_raw([tmpfile(content)])
        check(r.rc == 0, "empty file: rc=%d %s" % (r.rc, r.err[:200]))
        no_crash(r)


@test("robust")
def ligand_not_present():
    r = run_raw([tmpfile(base_system().text()), "XXX"])
    check(r.rc == 0 and "No ligand" in r.out, r.out + r.err)


@test("robust")
def truncated_lines_do_not_crash():
    random.seed(1)
    lines = base_system().text().splitlines()
    cut = [l[:random.randint(0, len(l))] for l in lines]
    for k in range(5):
        cut = [l[:random.randint(0, len(l))] if random.random() < 0.5 else l for l in lines]
        r = run_raw([tmpfile("\n".join(cut)), "LIG", "--csv"])
        check(r.rc == 0, "rc=%d: %s" % (r.rc, r.err[:300]))
        no_crash(r)


@test("robust")
def binary_garbage_does_not_crash():
    random.seed(2)
    data = bytes(random.getrandbits(8) for _ in range(20000))
    r = run_raw([tmpfile(data, binary=True), "LIG"])
    check(r.rc == 0, "rc=%d: %s" % (r.rc, r.err[:300]))
    no_crash(r)


@test("robust")
def non_numeric_coordinates():
    txt = base_system().text()
    bad = "\n".join(l[:30] + "abcdefgh" + l[38:] if l.startswith(("ATOM", "HETATM")) else l for l in txt.splitlines())
    r = run_raw([tmpfile(bad), "LIG", "--csv"])
    check(r.rc == 0, "rc=%d: %s" % (r.rc, r.err[:300]))
    no_crash(r)


@test("robust")
def huge_and_nan_coordinates():
    txt = base_system().text()
    for token in ["1e30    ", "-1e30   ", "nan     ", "inf     ", "99999999"]:
        bad = "\n".join(l[:30] + token + l[38:] if l.startswith("ATOM") else l for l in txt.splitlines())
        r = run_raw([tmpfile(bad), "LIG", "--csv"])
        check(r.rc == 0, "coordinate '%s': rc=%d %s" % (token.strip(), r.rc, r.err[:300]))
        no_crash(r)


@test("robust")
def coincident_atoms_zero_distance():
    s = Sys().ligand(carboxylate())
    s.prot("LYS", "NZ", (0.62, 1.07, 0))                  # exactly above O1
    s.prot("ASP", "OD1", (0.62, -1.07, 0))
    r = go(s)
    no_crash(r)


@test("robust")
def crlf_line_endings_give_same_result():
    txt = base_system().text()
    a = go(None, text=txt)
    b = go(None, text=txt.replace("\n", "\r\n"))
    check(a.rows == b.rows and a.rows, "CRLF changes the result")


@test("robust")
def missing_or_lowercase_element_column():
    txt = base_system().text()
    a = go(None, text=txt)
    nocol = "\n".join(l[:76] if l.startswith(("ATOM", "HETATM")) else l for l in txt.splitlines()) + "\n"
    b = go(None, text=nocol)
    check(sorted(x["kind"] for x in a.rows) == sorted(x["kind"] for x in b.rows) and b.rows,
          "without element column: %s versus %s" % (a.types(), b.types()))
    low = "\n".join(l[:76] + l[76:].lower() if l.startswith(("ATOM", "HETATM")) else l for l in txt.splitlines()) + "\n"
    c = go(None, text=low)
    check(sorted(x["kind"] for x in a.rows) == sorted(x["kind"] for x in c.rows), "lowercase element")


@test("robust")
def single_atom_ligand():
    s = Sys().prot("HIS", "NE2", (2.1, 0, 0), el="N").metal("Zn", (0, 0, 0))
    r = go(s, lig=("ZN",))
    no_crash(r)
    check(r.has("Metal coordination"), "single-atom ligand (Zn): " + str(r.types()))


@test("robust")
def ligand_with_only_hydrogens():
    s = Sys().prot("LYS", "NZ", (5, 5, 5))
    s.atom("HETATM", "H1", "LIG", "L", 1, (0, 0, 0), "H")
    r = run_raw([tmpfile(s.text()), "LIG"])
    check(r.rc == 0 and "No ligand" in r.out, r.out)


def sdf_system():
    m = methanol()
    s = Sys().ligand(m).prot("ASN", "OD1", at(m.atoms[1][2], HDIR, 2.8))
    return m, s


@test("robust")
def malformed_sdf_files():
    m, s = sdf_system()
    good = m.sdf()
    lines = good.splitlines()
    cases = {
        "v3000": good.replace("V2000", "V3000"),
        "troncato_a_meta": "\n".join(lines[:5]) + "\n",
        "senza_blocco_bond": "\n".join(lines[:6]) + "\nM  END\n$$$$\n",
        "legame_indice_fuori_range": good.replace("  1  2  1  0", " 99  2  1  0"),
        "legame_indice_zero": good.replace("  1  2  1  0", "  0  2  1  0"),
        "legame_indice_negativo": good.replace("  1  2  1  0", " -5  2  1  0"),
        "conteggio_atomi_enorme": good.replace("  3  2  0  0", "999999  0  0  0") if False else "\n".join(lines[:3] + ["999999999  1  0  0  0  0  0  0  0  0999 V2000"] + lines[4:]),
        "vuoto": "",
        "solo_testo": "ciao\nmondo\n",
        "riga_atomo_corta": "\n".join(lines[:4] + ["1.0 2.0"] + lines[5:]),
    }
    for name, content in cases.items():
        r = go(s, sdf=content)
        no_crash(r)


@test("robust")
def malformed_mol2_files():
    m, s = sdf_system()
    good = m.mol2()
    lines = good.splitlines()
    cases = {
        "troncato": "\n".join(lines[:6]) + "\n",
        "id_bond_inesistenti": good.replace("1 1 2 1", "1 77 88 1"),
        "senza_molecule": "\n".join(l for l in lines if "MOLECULE" not in l),
        "senza_bond": "\n".join(lines[:lines.index("@<TRIPOS>BOND")]) + "\n",
        "atomo_incompleto": good.replace("C1 0.0000 0.0000 0.0000", "C1"),
        "vuoto": "",
        "tipo_vuoto": good.replace(" C 1 LIG", " . 1 LIG"),
    }
    for name, content in cases.items():
        r = go(s, mol2=content)
        no_crash(r)


@test("robust")
def invalid_command_line_is_rejected():
    pdb = tmpfile(base_system().text())
    for args, what in [([pdb, "--sdf"], "--sdf without a file"), ([pdb, "--foo", "1"], "unknown option"),
                       ([pdb, "LIG", "--hb-max"], "threshold without a value"),
                       ([pdb, "LIG", "--hb-max", "abc"], "non-numeric threshold"),
                       ([pdb, "LIG", "--sdf", "/no/such.sdf"], "nonexistent SDF file"),
                       ([pdb, "LIG", "--mol2", "/no/such.mol2"], "nonexistent MOL2 file")]:
        r = run_raw(args)
        check(r.rc == 1, "%s: expected code 1, got %d (%s)" % (what, r.rc, r.err[:120]))
        no_crash(r)


@test("robust")
def oversized_thresholds_are_clamped_not_silent():
    r = go(base_system(), args=["--hb-max", "9", "--salt-max", "12"])
    check("limited" in r.err.lower() or "clamp" in r.err.lower() or "massimo" in r.err.lower(),
          "threshold beyond the grid limit (6 A) without a warning: " + r.err[:200])


@test("robust")
def large_system_runs_in_reasonable_time():
    random.seed(3)
    s = Sys().ligand(carboxylate())
    n = 40000
    resname = "ALA"
    for i in range(n):
        base = (random.uniform(-60, 60), random.uniform(-60, 60), random.uniform(-60, 60))
        for nm in ("N", "CA", "C", "O", "CB"):
            p = (base[0] + random.uniform(-1.5, 1.5), base[1] + random.uniform(-1.5, 1.5), base[2] + random.uniform(-1.5, 1.5))
            s.atom("ATOM", nm, resname, "A", i % 9999 + 1, p, nm[0])
    limit = 120 if CFG["sanitize"] else 40
    r = go(s, timeout=limit * 3)
    check(CFG["last_time"] < limit, "200000 atoms: %.1f s (limit %d s)" % (CFG["last_time"], limit))
    # choosing an alternate conformation reads the file twice: it must stay in the same order of magnitude
    t0 = time.time()
    r2 = run_raw([tmpfile(s.text()), "LIG", "--altloc", "A", "--csv"], timeout=limit * 3)
    check(r2.rc == 0 and time.time() - t0 < limit, "200000 atoms with --altloc: %.1f s (limit %d s)" % (time.time() - t0, limit))
    no_crash(r)


# ----------------------------------------------------------------------------------------------
# GROUP golden: real complexes (regression of results verified by hand)
# ----------------------------------------------------------------------------------------------
GOLDEN = [
    # (pdb, ligand arguments, options, [expected (kind, ligand contains, protein contains)], [forbidden (kind, ...)])
    ("2OZR.pdb", ["GG1", "A", "2001"], [],
     [("Hydrogen bond", ":O2", "THR 226A:N"), ("Salt bridge", ":O6", "LYS 119A:NZ"),
      ("pi-stacking (parallel)", "", "PHE 231A")], []),
    ("1efy.pdb", ["BZC"], ["--sdf", "1efy.sdf"],
     [("Hydrogen bond", ":N1", "GLY 863A:O"), ("Hydrogen bond", ":O1", "GLY 863A:N"),
      ("pi-stacking (parallel)", "", "TYR 907A")], []),
    ("4IKF.pdb", ["M76"], [],
     [("Metal coordination", ":O14", "MG 405A"), ("Metal coordination", ":O16", "MG 404A"),
      ("Metal coordination", ":O18", "MG 404A"), ("pi-stacking (parallel)", "", "DA 17D"),
      ("pi-stacking (parallel)", "", "DC 16D")],
     [("Hydrogen bond", "", ""), ("Metal coordination", ":N15", "")]),
    ("3PTB.pdb", ["BEN"], ["--no-water"],
     [("Salt bridge", ":N1", "ASP 189A"), ("Salt bridge", ":N2", "ASP 189A"),
      ("Hydrogen bond", ":N1", "GLY 219A:O"), ("Hydrogen bond", ":N2", "SER 190A:OG")], []),
    ("3PTB.pdb", ["BEN"], ["--backbone-hyd"],
     [("Apolar contact (backbone)", "", "TRP 215A:C"), ("Apolar contact (backbone)", "", "CYS 191A:CA")], []),
    ("1STP.pdb", ["BTN"], ["--no-water"],
     [("Hydrogen bond", ":N1", "ASP 128A"), ("Hydrogen bond", ":N2", "SER 45A"),
      ("Hydrogen bond", ":O11", "SER 88A"), ("Hydrogen bond", ":O12", "ASN 49A:N"),
      ("Hydrogen bond", ":O3", "TYR 43A"), ("Hydrogen bond", ":O3", "SER 27A"),
      ("Hydrogen bond", ":O3", "ASN 23A")], []),
    ("1IEP.pdb", ["STI", "A", "201"], ["--no-water", "--pi-angle", "40"],
     [("Hydrogen bond", ":N13", "THR 315A"), ("Hydrogen bond", ":N21", "GLU 286A"),
      ("Hydrogen bond", ":N3", "MET 318A"), ("Hydrogen bond", ":O29", "ASP 381A"),
      ("Hydrogen bond", ":N51", "ILE 360A:O"), ("pi-stacking (parallel)", "", "PHE 317A"),
      ("pi-stacking (parallel)", "", "PHE 382A")],
     [("Cation-pi", "", "")]),
    ("1HSG.pdb", ["MK1"], ["--no-water"],
     [("Hydrogen bond", ":O2", "ASP 25B"), ("Hydrogen bond", ":O2", "ASP 25A"),
      ("Hydrogen bond", ":N4", "GLY 27B:O"), ("Hydrogen bond", ":O4", "ASP 29B:N"),
      ("Hydrogen bond", ":O4", "GLY 27B:O")], []),
    ("1AZM.pdb", ["AZM"], ["--no-water"],
     [("Metal coordination", ":N1", "ZN 261A"), ("Hydrogen bond", ":N1", "THR 199A:OG1"),
      ("Hydrogen bond", ":O2", "THR 199A:N")],
     [("Hydrogen bond", "", "HIS 94A"), ("Hydrogen bond", "", "HIS 96A"), ("Hydrogen bond", "", "HIS 119A")]),
]


def make_golden(idx, entry):
    pdb, lig, opts, expect, forbid = entry

    def fn():
        path = os.path.join(CFG["real_dir"], pdb)
        if not os.path.exists(path):
            raise Skip("missing %s" % path)
        extra = []
        i = 0
        while i < len(opts):
            if opts[i] in ("--sdf", "--mol2"):
                p = os.path.join(CFG["real_dir"], opts[i + 1])
                if not os.path.exists(p):
                    raise Skip("missing %s" % p)
                extra += [opts[i], p]
                i += 2
            else:
                extra.append(opts[i]); i += 1
        r = run_raw([path] + lig + ["--csv"] + extra)
        check(r.rc == 0, "rc=%d %s" % (r.rc, r.err[:200]))
        no_crash(r)
        for kind, la, pa in expect:
            check(r.has(kind, la or None, pa or None), "expected but not found: %s %s %s" % (kind, la, pa))
        for kind, la, pa in forbid:
            check(not r.has(kind, la or None, pa or None), "found a forbidden contact: %s %s %s" % (kind, la, pa))
    fn.__name__ = "golden_%s_%s" % (pdb.split(".")[0], "_".join(lig).replace(" ", "")) + ("_%d" % idx if idx else "")
    return fn


for _i, _e in enumerate(GOLDEN):
    _f = make_golden(_i, _e)
    TESTS.append(("golden", _f.__name__, _f))


# ----------------------------------------------------------------------------------------------
# esecuzione
# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default="", help="already built ligint binary (default: build the sources)")
    ap.add_argument("--sanitize", action="store_true")
    ap.add_argument("-k", default="", help="run only the tests with this substring in the name")
    ap.add_argument("--group", default="", help="chem, tmpl, pdb, robust, golden (comma-separated)")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--real-dir", default=os.path.join(HERE, "data"))
    ap.add_argument("-v", action="store_true", help="show the message also for passed/skipped tests")
    ap.add_argument("--keep", action="store_true", help="keep the temporary files (they are removed when all the tests pass; with failures they are kept)")
    a = ap.parse_args()
    if a.list:
        for g, n, _ in TESTS:
            print("%-7s %s" % (g, n))
        return 0
    tmp = tempfile.mkdtemp(prefix="ligint_test_")
    CFG.update(tmp=tmp, sanitize=a.sanitize, real_dir=os.path.abspath(a.real_dir), bin=os.path.join(tmp, "ligint"))
    warn = ""
    if a.bin:
        CFG["bin"] = os.path.abspath(a.bin)
        if not os.path.exists(CFG["bin"]):
            shutil.rmtree(tmp, ignore_errors=True)
            sys.exit("binary not found: " + CFG["bin"])
        print("binary: %s (prebuilt)" % CFG["bin"])
    else:
        warn = _build.build_ligint(CFG["bin"], a.sanitize)
        print("binary: %s (%s)" % (CFG["bin"], "ASan+UBSan" if a.sanitize else "-O2"))
    if warn:
        print("compiler warnings:\n" + warn)
    groups = {g for g in a.group.split(",") if g}
    ok = fail = skip = 0
    failures = []
    for g, n, fn in TESTS:
        if (groups and g not in groups) or (a.k and a.k not in n):
            continue
        try:
            fn()
            ok += 1
            print("  ok    %-7s %s" % (g, n))
        except Skip as e:
            skip += 1
            print("  skip  %-7s %s  (%s)" % (g, n, e))
        except Fail as e:
            fail += 1
            failures.append((g, n, str(e)))
            print("  FAIL  %-7s %s\n          -> %s" % (g, n, str(e).replace("\n", " ")[:600]))
        except Exception as e:   # error in the test itself
            fail += 1
            failures.append((g, n, "error in the test: %r" % e))
            print("  ERROR %-7s %s\n          -> %r" % (g, n, e))
    print("\nresult: %d passed, %d failed, %d skipped" % (ok, fail, skip))
    if fail or a.keep:
        print("temporary files kept in: " + tmp)
    else:
        shutil.rmtree(tmp, ignore_errors=True)       # the files of the tests are large (hundreds of MB): do not leave them in the temporary folder
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
