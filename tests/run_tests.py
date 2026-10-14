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
    robust malformed / anomalous input: the program must not crash nor print nan
    golden real complexes: list of interactions already verified by hand against diagrams (REGRESSION
           tests, not validation tests)

Exit code: 0 if all the tests pass, 1 otherwise.
"""
import argparse
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

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


def run_raw(args, timeout=120):
    t0 = time.time()
    try:
        p = subprocess.run([CFG["bin"]] + args, capture_output=True, timeout=timeout, cwd=CFG["tmp"])
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
        with open(path, "w", newline="") as f:
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
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
