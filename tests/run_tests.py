#!/usr/bin/env python3
"""Suite di test per ligint (regole chimiche, geometria, robustezza degli input).

Uso:
    python3 tests/run_tests.py                      compila i sorgenti ed esegue tutto
    python3 tests/run_tests.py --bin build/src/ligint   usa un binario gia' compilato (es. da CMake)
    python3 tests/run_tests.py --sanitize           compila con ASan+UBSan (rileva crash e comportamenti indefiniti)
    python3 tests/run_tests.py -k pi_              solo i test il cui nome contiene "pi_"
    python3 tests/run_tests.py --list               elenca i test
    python3 tests/run_tests.py --real-dir DIR       cartella con i PDB reali per i test "golden" (default: tests/data;
                                                    i test con file mancanti vengono saltati)

Gruppi:
    chem   regole chimiche e geometriche su frammenti sintetici (geometria nota, esito atteso noto)
    tmpl   template SDF / MOL2 (abbinamento per coordinate, H espliciti, fallback)
    pdb    parsing del PDB (altloc, codici di inserzione, modelli, selezione del ligando, CSV)
    robust input malformati / anomali: il programma non deve andare in crash ne' stampare nan
    golden complessi reali: elenco di interazioni gia' verificate a mano contro diagrammi (sono test di
           REGRESSIONE, non di validazione)

Codice di uscita: 0 se tutti i test passano, 1 altrimenti.
"""
import argparse
import math
import os
import random
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _build  # noqa: E402
TESTS = []          # (gruppo, nome, funzione)
CFG = {}            # binario, tmp, ...


# ----------------------------------------------------------------------------------------------
# infrastruttura minima
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
    """punto a 'dist' da 'origin' lungo 'direction'"""
    return v_add(origin, v_mul(v_unit(direction), dist))


# ----------------------------------------------------------------------------------------------
# molecole sintetiche (coordinate in Angstrom, geometria realistica)
# ----------------------------------------------------------------------------------------------
class Mol:
    def __init__(self):
        self.atoms = []   # (nome, elemento, xyz, tipo mol2, carica formale)
        self.bonds = []   # (i, j, ordine)  indici 0-based
        self.hs = []      # (indice atomo pesante, xyz)

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
    m.h(o2, (1.12, -1.90, 0))        # presente solo nella versione "acido"
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
    """2-amminopirimidina: C2 con tre N vicini ma IN ANELLO (non e' un guanidinio)."""
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


HDIR = (0.317, 0.948, 0.0)   # direzione O-H del metanolo (angolo C-O-H = 108.5 gradi)


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
# costruzione dei PDB e esecuzione
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
        if lines and lines[0].startswith("ligando,tipo"):
            for l in lines[1:]:
                f = l.split(",", 6)
                if len(f) == 7:
                    self.rows.append({"lig": f[0], "tipo": f[1], "fav": f[2], "la": f[3], "pa": f[4], "d": f[5], "note": f[6]})

    def find(self, tipo=None, la=None, pa=None):
        return [r for r in self.rows if (tipo is None or r["tipo"] == tipo) and (la is None or la in r["la"])
                and (pa is None or pa in r["pa"])]

    def has(self, tipo, la=None, pa=None):
        return bool(self.find(tipo, la, pa))

    def types(self):
        return sorted({r["tipo"] for r in self.rows})


def run_raw(args, timeout=120):
    t0 = time.time()
    try:
        p = subprocess.run([CFG["bin"]] + args, capture_output=True, timeout=timeout, cwd=CFG["tmp"])
    except subprocess.TimeoutExpired:
        raise Fail("timeout dopo %d s" % timeout)
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
    check(r.rc == 0, "codice di uscita %d; stderr: %s" % (r.rc, r.err[:300]))
    return r


def no_crash(r):
    check("nan" not in r.out.lower().replace("annot", ""), "output contiene 'nan'")
    check("AddressSanitizer" not in r.err and "runtime error" not in r.err, "sanitizer: " + r.err[:400])


# ----------------------------------------------------------------------------------------------
# GRUPPO chem: regole chimiche
# ----------------------------------------------------------------------------------------------
@test("chem")
def carboxylate_forms_salt_bridge_with_lysine():
    s = Sys().ligand(carboxylate()).prot("LYS", "NZ", at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8))
    r = go(s)
    check(r.has("Ponte salino", la=":O1", pa="LYS"), "atteso ponte salino carbossilato-lisina: " + str(r.types()))


@test("chem")
def carboxylate_repels_glutamate():
    s = Sys().ligand(carboxylate()).prot("GLU", "OE1", at((0.62, 1.07, 0), (0.5, 0.866, 0), 3.4))
    r = go(s)
    check(r.has("Repulsione elettrostatica"), "attesa repulsione carbossilato-carbossilato")
    check(not r.has("Ponte salino"), "ponte salino tra cariche uguali")


@test("chem")
def amidine_is_cationic():
    s = Sys().ligand(amidine()).prot("ASP", "OD1", at((0.67, 1.15, 0), (0.5, 0.866, 0), 2.9))
    r = go(s)
    check(r.has("Ponte salino", la=":N1"), "attesa amidina cationica: " + str(r.types()))


@test("chem")
def guanidinium_is_cationic():
    s = Sys().ligand(guanidinium()).prot("GLU", "OE1", at((0, 1.33, 0), (0, 1, 0), 2.9))
    r = go(s)
    check(r.has("Ponte salino", la=":N1"), "atteso guanidinio cationico: " + str(r.types()))


@test("chem")
def aminopyrimidine_is_not_cationic():
    """Regressione: il guanidinio 'N-C(N)-N' in anello (2-amminopirimidina) NON e' un catione."""
    m = aminopyrimidine()
    n1 = m.atoms[5][2]                       # N1
    d = v_unit(n1)
    s = Sys().ligand(m).prot("LYS", "NZ", at(n1, d, 3.4))
    r = go(s)
    check(not r.has("Repulsione elettrostatica"), "N di pirimidina trattato come catione (falso guanidinio)")
    check(not r.has("Ponte salino"), "ponte salino con anello neutro")


@test("chem")
def ring_nitrogen_is_acceptor():
    m = aminopyrimidine()
    n1 = m.atoms[5][2]
    s = Sys().ligand(m).prot("LYS", "NZ", at(n1, v_unit(n1), 2.9))
    r = go(s)
    check(r.has("Legame a idrogeno", la=":N1"), "N di anello deve poter accettare un H-bond: " + str(r.types()))


@test("chem")
def aliphatic_amine_is_protonated():
    s = Sys().ligand(trimethylamine()).prot("ASP", "OD1", (0, 0, 2.9))
    r = go(s)
    check(r.has("Ponte salino", la=":N1"), "ammina alifatica assunta protonata: " + str(r.types()))


@test("chem")
def amide_nitrogen_is_not_acceptor():
    m = acetamide()
    n = m.atoms[3][2]
    dirn = (0.52, -0.85, 0)
    # donatore (LYS NZ) vicino all'N ammidico: l'N non e' accettore -> nessun H-bond
    r1 = go(Sys().ligand(m).prot("LYS", "NZ", at(n, dirn, 2.9)))
    check(not r1.has("Legame a idrogeno", la=":N1"), "N ammidico trattato da accettore")
    # accettore (ASN OD1) vicino all'N: l'N-H ammidico e' donatore -> H-bond
    r2 = go(Sys().ligand(m).prot("ASN", "OD1", at(n, dirn, 2.9)))
    check(r2.has("Legame a idrogeno", la=":N1"), "N-H ammidico deve poter donare")


@test("chem")
def hydroxyl_donates_unless_metal_coordinated():
    m = methanol()
    o1 = m.atoms[1][2]
    glu = at(o1, (1, 0, 0), 2.7)
    r0 = go(Sys().ligand(m).prot("GLU", "OE1", glu))
    check(r0.has("Legame a idrogeno", la=":O1"), "idrossile libero: atteso H-bond")
    mg = (o1[0], o1[1] + 2.1, 0)
    r1 = go(Sys().ligand(m).prot("GLU", "OE1", glu).metal("Mg", mg))
    check(r1.has("Coordinazione metallica", la=":O1"), "atteso O1-Mg")
    check(not r1.has("Legame a idrogeno", la=":O1"), "O coordinato a Mg trattato come donatore")
    r2 = go(Sys().ligand(m).prot("GLU", "OE1", glu).metal("Mg", mg), args=["--metal-donors"])
    check(r2.has("Legame a idrogeno", la=":O1"), "--metal-donors dovrebbe ripristinare il donatore")


@test("chem")
def metal_distance_limits_per_element():
    m = acetamide()
    o = m.atoms[2][2]
    d = (0.49, 0.87, 0)
    for el, dist, expect in [("Mg", 2.4, True), ("Mg", 2.7, False), ("Zn", 2.55, True), ("Zn", 2.75, False)]:
        r = go(Sys().ligand(m).metal(el, at(o, d, dist)))
        check(r.has("Coordinazione metallica", la=":O1") == expect, "%s a %.2f A: atteso %s" % (el, dist, expect))


@test("chem")
def metal_needs_lone_pair():
    m = acetamide()
    n = m.atoms[3][2]
    r = go(Sys().ligand(m).metal("Zn", at(n, (0.52, -0.85, 0), 2.3)))
    check(not r.has("Coordinazione metallica", la=":N1"), "N ammidico (senza coppia solitaria) coordinato a Zn")


@test("chem")
def histidine_bound_to_zinc_does_not_donate():
    m = aminopyrimidine()
    n1 = m.atoms[5][2]
    d = v_unit(n1)
    ne2 = at(n1, d, 3.0)
    r0 = go(Sys().ligand(m).prot("HIS", "NE2", ne2))
    check(r0.has("Legame a idrogeno", la=":N1", pa="HIS"), "controllo: senza Zn l'H-bond con His deve esserci")
    r1 = go(Sys().ligand(m).prot("HIS", "NE2", ne2).metal("Zn", at(ne2, d, 2.1)))
    check(not r1.has("Legame a idrogeno", la=":N1", pa="HIS"), "His coordinata a Zn trattata da donatore")


@test("chem")
def hydrophobic_rules():
    leu = (0.765, 0, 3.72)
    r = go(Sys().ligand(ethane()).prot("LEU", "CD1", leu))
    check(r.has("Contatto idrofobico"), "atteso contatto idrofobico C-C")
    r = go(Sys().ligand(methanol()).prot("LEU", "CD1", (-2.0, 0, 3.23)))
    check(not r.has("Contatto idrofobico", la=":C1"), "C legato a O trattato da idrofobico")
    r = go(Sys().ligand(ethane()).prot("GLU", "CG", leu))
    check(not r.has("Contatto idrofobico"), "GLU CG non deve essere idrofobico")
    r = go(Sys().ligand(ethane()).prot("ALA", "CA", leu))
    check(not r.has("Contatto idrofobico") and not r.has("Contatto apolare (backbone)"), "CA non conta di default")
    r = go(Sys().ligand(ethane()).prot("ALA", "CA", leu), args=["--backbone-hyd"])
    check(r.has("Contatto apolare (backbone)"), "--backbone-hyd deve contare CA")


@test("chem")
def hydrophobic_thresholds_are_configurable():
    s = lambda: Sys().ligand(ethane()).prot("LEU", "CD1", at((0.765, 0, 0), (0, 0, 1), 4.2))
    check(not go(s()).has("Contatto idrofobico"), "4,2 A fuori dalla soglia di default")
    check(go(s(), args=["--hyd-max", "4.5"]).has("Contatto idrofobico"), "--hyd-max 4.5 deve includerlo")


@test("chem")
def pi_stacking_parallel_and_offset_limit():
    for ox, expect in [(1.0, True), (3.0, True), (3.7, False), (6.0, False)]:
        s = Sys().ligand(benzene()).ring("PHE", (ox, 0, 3.6), (0, 0, 1))
        r = go(s)
        found = r.has("pi-stacking (parallelo)")
        check(found == expect, "offset %.1f: atteso %s, trovato %s" % (ox, expect, found))


@test("chem")
def pi_stacking_angle_threshold():
    n = (0, -math.sin(math.radians(36)), math.cos(math.radians(36)))
    mk = lambda: Sys().ligand(benzene()).ring("PHE", (0.5, 0, 3.7), n)
    check(not go(mk()).has("pi-stacking (parallelo)"), "36 gradi oltre il limite di 30")
    check(go(mk(), args=["--pi-angle", "40"]).has("pi-stacking (parallelo)"), "--pi-angle 40 deve accettarlo")


@test("chem")
def pi_stacking_t_shaped():
    s = Sys().ligand(benzene()).ring("PHE", (0, 0, 4.8), (0, 1, 0))
    r = go(s)
    check(r.has("pi-stacking (a T)"), "atteso stacking a T: " + str(r.types()))


@test("chem")
def cation_pi_geometry():
    r = go(Sys().ligand(benzene()).prot("LYS", "NZ", (0, 0, 4.0)))
    check(r.has("Cation-pi"), "LYS sull'asse dell'anello: atteso cation-pi")
    r = go(Sys().ligand(benzene()).prot("LYS", "NZ", (3.46, 0, 2.0)))
    check(not r.has("Cation-pi"), "LYS a 60 gradi dall'asse: nessun cation-pi")
    r = go(Sys().ligand(trimethylamine()).ring("PHE", (0, 0, 4.2), (0, 0, 1)))
    check(r.has("Cation-pi", la=":N1"), "catione del ligando sopra l'anello della proteina: " + str(r.types()))


@test("chem")
def halogen_bond_needs_linear_geometry():
    cl = (1.78, 0, 0)
    r = go(Sys().ligand(chloromethane()).prot("ASN", "OD1", at(cl, (1, 0, 0), 3.2)))
    check(r.has("Legame ad alogeno"), "C-Cl...O lineare: atteso legame ad alogeno")
    r = go(Sys().ligand(chloromethane()).prot("ASN", "OD1", at(cl, (0, 1, 0), 3.2)))
    check(not r.has("Legame ad alogeno"), "C-Cl...O a 90 gradi: nessun legame ad alogeno")


@test("chem")
def clash_and_covalent_distances():
    r = go(Sys().ligand(ethane()).prot("LEU", "CD1", (0, 0, 2.5)))
    check(r.has("Clash sterico"), "2,5 A tra due C: atteso clash")
    r = go(Sys().ligand(ethane()).prot("LEU", "CD1", (-1.6, 0, 0)))      # 1,6 A da C1 e 3,1 A da C2
    check(r.has("Possibile legame covalente", la=":C1") and not r.has("Clash sterico"),
          "1,6 A: atteso possibile legame covalente: " + str(r.types()))


@test("chem")
def water_bridge():
    m = methanol()
    o1 = m.atoms[1][2]
    w = at(o1, (1, 0, 0), 2.8)
    p = at(w, (0.5, 0.866, 0), 2.8)
    mk = lambda: Sys().ligand(m).water(w).prot("ASN", "ND2", p)
    check(go(mk()).has("Ponte d'acqua"), "atteso ponte d'acqua")
    check(not go(mk(), args=["--no-water"]).has("Ponte d'acqua"), "--no-water ignorato")
    # contatto diretto: nessun ponte
    s = Sys().ligand(m).water(w).prot("ASN", "ND2", at(o1, (1, 0.3, 0), 3.0))
    check(not go(s).has("Ponte d'acqua"), "ponte d'acqua tra atomi gia' a contatto diretto")


@test("chem")
def ligand_metal_coordination_removes_hbond_partners():
    """Regressione (eme): N del ligando e N della proteina coordinati allo stesso Fe del ligando non fanno H-bond."""
    m = Mol()                                   # controllo: stesso ligando SENZA Fe
    n1 = m.add("N1", "N", (2.0, 0, 0))
    c1 = m.add("C1", "C", (2.7, 1.1, 0))
    c2 = m.add("C2", "C", (2.7, -1.1, 0))
    m.bond(n1, c1); m.bond(n1, c2)
    nz = (2.0, 0.0, -2.9)
    r = go(Sys().ligand(m).prot("LYS", "NZ", nz))
    check(r.has("Legame a idrogeno", la=":N1"), "controllo: senza Fe l'N del ligando deve accettare dalla lisina")
    # con N coordinato al Fe del ligando serve il metallo vicino: il Fe e' a 2,0 A -> N coordinato
    m2 = Mol()
    fe2 = m2.add("FE", "Fe", (0, 0, 0))
    n = m2.add("N1", "N", (2.0, 0, 0))
    a = m2.add("C1", "C", (2.7, 1.1, 0)); b = m2.add("C2", "C", (2.7, -1.1, 0))
    m2.bond(n, a); m2.bond(n, b)
    r = go(Sys().ligand(m2).prot("LYS", "NZ", nz))
    check(not r.has("Legame a idrogeno", la=":N1"), "N coordinato al Fe del ligando trattato da accettore")
    # istidina che coordina il Fe del ligando: nessun H-bond con l'N del ligando
    r = go(Sys().ligand(m2).prot("HIS", "NE2", (0, 0, 2.1)))
    check(not r.has("Legame a idrogeno", pa="HIS"), "His coordinata al Fe del ligando trattata da donatore/accettore")


@test("chem")
def water_bridge_needs_realistic_angle():
    m = methanol()
    o1 = m.atoms[1][2]
    w = at(o1, (1, 0, 0), 2.8)
    s = lambda p: Sys().ligand(m).water(w).prot("ASN", "ND2", p)
    check(go(s(at(w, (0.5, 0.866, 0), 2.8))).has("Ponte d'acqua"), "angolo a 120 gradi: ponte atteso")
    check(not go(s(at(w, (1, 0, 0), 2.8))).has("Ponte d'acqua"), "acqua collineare (180 gradi): nessun ponte")
    check(go(s(at(w, (1, 0, 0), 2.8)), args=["--wb-max", "180"]).has("Ponte d'acqua"), "--wb-max 180 deve riammettere l'angolo")


def rot2(v, deg):
    a = math.radians(deg)
    return (v[0] * math.cos(a) - v[1] * math.sin(a), v[0] * math.sin(a) + v[1] * math.cos(a), 0.0)


@test("chem")
def long_hbond_needs_possible_hydrogen_geometry():
    # idrossile (metanolo): l'H ruota, ma l'angolo C-O...A deve essere vicino a 109,5 gradi
    m = methanol()
    o1 = m.atoms[1][2]
    good = at(o1, rot2((-1, 0, 0), -110), 3.9)        # C-O...A = 110 gradi, D...A = 3,9 A
    bad = at(o1, rot2((-1, 0, 0), -40), 3.9)           # 40 gradi: l'H non puo' puntare all'accettore
    s = lambda p: Sys().ligand(m).prot("ASN", "OD1", p)
    check(not go(s(good)).has("Legame a idrogeno", la=":O1"), "3,9 A senza --hb-long: nessun H-bond")
    r = go(s(good), args=["--hb-long", "4.1"])
    check(r.has("Legame a idrogeno", la=":O1") and "lungo" in r.find("Legame a idrogeno")[0]["note"], "3,9 A con geometria possibile: atteso H-bond lungo")
    check(not go(s(bad), args=["--hb-long", "4.1"]).has("Legame a idrogeno", la=":O1"), "geometria impossibile: nessun H-bond lungo")
    check(not go(s(good), args=["--hb-long", "3.8"]).has("Legame a idrogeno", la=":O1"), "oltre il limite --hb-long")
    # NH2 ammidico: l'angolo C-N...A atteso e' ~120 gradi
    a = acetamide()
    n = a.atoms[3][2]
    v1 = v_unit(v_sub(a.atoms[1][2], n))
    ok = at(n, rot2(v1, 120), 3.9)
    r = go(Sys().ligand(a).prot("ASN", "OD1", ok), args=["--hb-long", "4.1"])
    check(r.has("Legame a idrogeno", la=":N1"), "NH2 ammidico con geometria possibile: " + str(r.types()))
    no = at(n, v1, 3.9)                               # accettore dal lato del C: impossibile
    r = go(Sys().ligand(a).prot("ASN", "OD1", no), args=["--hb-long", "4.1"])
    check(not r.has("Legame a idrogeno", la=":N1"), "NH2 con accettore dal lato sbagliato")


@test("chem")
def salt_bridge_by_group_centroid():
    c = carboxylate()
    near = (0.62 + 4.6, 0, 0)                          # atomo-atomo 4,7 A (>4,0), centroide-centroide 4,6 A
    s = lambda: Sys().ligand(c).prot("LYS", "NZ", near)
    check(not go(s()).has("Ponte salino"), "default: 4,7 A atomo-atomo oltre salt-max 4,0")
    r = go(s(), args=["--salt-centroid", "5.0"])
    check(len(r.find("Ponte salino")) == 1, "centroidi a 4,6 A: atteso un ponte salino: " + str(r.types()))
    check(not go(s(), args=["--salt-centroid", "4.0"]).has("Ponte salino"), "centroidi oltre il limite")
    # un gruppo carico = un solo contatto, non uno per coppia di atomi
    g = guanidinium()
    mk = lambda: Sys().ligand(g).prot("GLU", "OE1", (-0.5, 4.5, 0)).prot("GLU", "OE2", (0.5, 4.5, 0))
    check(len(go(mk()).find("Ponte salino")) > 1, "default: piu' coppie atomo-atomo")
    check(len(go(mk(), args=["--salt-centroid", "5.0"]).find("Ponte salino")) == 1, "per gruppo: un solo ponte salino")


@test("chem")
def plip_like_preset_and_explicit_override():
    m = methanol()
    o1 = m.atoms[1][2]
    good = at(o1, rot2((-1, 0, 0), -110), 3.9)
    s = lambda: Sys().ligand(m).prot("ASN", "OD1", good)
    check(not go(s()).has("Legame a idrogeno", la=":O1"), "default: nessun H-bond a 3,9 A")
    check(go(s(), args=["--plip-like"]).has("Legame a idrogeno", la=":O1"), "--plip-like deve ammettere l'H-bond lungo")
    check(not go(s(), args=["--plip-like", "--hb-long", "0"]).has("Legame a idrogeno", la=":O1"), "opzione esplicita dopo il preset")
    check(not go(s(), args=["--hb-long", "0", "--plip-like"]).has("Legame a idrogeno", la=":O1"), "opzione esplicita prima del preset")
    c = carboxylate()
    near = (0.62 + 4.6, 0, 0)
    mk = lambda: Sys().ligand(c).prot("LYS", "NZ", near)
    check(go(mk(), args=["--plip-like"]).has("Ponte salino"), "--plip-like: ponte salino per centroidi a 4,6 A")
    check(not go(mk(), args=["--plip-like", "--salt-centroid", "0"]).has("Ponte salino"), "--salt-centroid 0 spegne i centroidi")


@test("chem")
def nucleic_acid_stacking_and_phosphate():
    # base adeninica (6+5 atomi) in z=0, ligando benzenico sopra
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
    check(r.has("pi-stacking (parallelo)", pa="DA 17D"), "atteso stacking su base di DNA: " + str(r.types()))
    # fosfato anionico + ammina del ligando
    s = Sys().ligand(trimethylamine()).prot("DA", "OP1", (0, 0, 3.0), resi=5, chain="D", el="O")
    r = go(s)
    check(r.has("Ponte salino", pa="DA 5D:OP1"), "atteso ponte salino con fosfato del DNA: " + str(r.types()))


# ----------------------------------------------------------------------------------------------
# GRUPPO tmpl: template SDF / MOL2
# ----------------------------------------------------------------------------------------------
@test("tmpl")
def sdf_with_explicit_h_distinguishes_acid_from_carboxylate():
    m = carboxylate()
    probe = at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8)
    s = lambda: Sys().ligand(m).prot("LYS", "NZ", probe)
    check(go(s()).has("Ponte salino"), "senza template: carbossilato anionico (assunzione a pH fisiologico)")
    r = go(s(), sdf=m.sdf(with_h=True))
    check("template" in r.err, "template non applicato: " + r.err)
    check(not r.has("Ponte salino"), "con H sull'ossigeno il gruppo e' un acido: nessun ponte salino")
    check(r.has("Legame a idrogeno", la=":O1"), "O=C dell'acido deve accettare l'H-bond dalla lisina")


@test("tmpl")
def template_hbond_uses_hydrogen_geometry():
    m = methanol()
    o1 = m.atoms[1][2]
    h = v_add(o1, v_mul(HDIR, 0.96))
    aligned = at(o1, HDIR, 2.8)                       # H punta verso l'accettore
    mirror = at(o1, (0.317, -0.948, 0), 2.8)           # H punta altrove
    sdf = m.sdf(with_h=True)
    r = go(Sys().ligand(m).prot("ASN", "OD1", aligned), sdf=sdf)
    check(r.has("Legame a idrogeno", la=":O1") and "D-H...A" in r.find("Legame a idrogeno")[0]["note"],
          "H allineato: atteso H-bond con angolo D-H...A")
    r = go(Sys().ligand(m).prot("ASN", "OD1", mirror), sdf=sdf)
    check(not r.has("Legame a idrogeno", la=":O1"), "H che punta altrove: nessun H-bond con template")
    r = go(Sys().ligand(m).prot("ASN", "OD1", mirror))
    check(r.has("Legame a idrogeno", la=":O1"), "senza H (euristica) il contatto e' accettato")


@test("tmpl")
def template_without_hydrogens_keeps_donors():
    """Regressione: un SDF senza H non deve azzerare i donatori (si usa l'euristica)."""
    m = methanol()
    o1 = m.atoms[1][2]
    r = go(Sys().ligand(m).prot("ASN", "OD1", at(o1, HDIR, 2.8)), sdf=m.sdf(with_h=False))
    check("template" in r.err, "template non applicato")
    check(r.has("Legame a idrogeno", la=":O1"), "template senza H: O-H del metanolo deve donare")


@test("tmpl")
def template_mismatch_falls_back_to_heuristics():
    m = methanol()
    s = lambda: Sys().ligand(m).prot("ASN", "OD1", at(m.atoms[1][2], HDIR, 2.8))
    r = go(s(), sdf=m.sdf(shift=(1.0, 0, 0)))
    check("nessun template" in r.err, "SDF traslato di 1 A non deve combaciare: " + r.err)
    other = ethane()
    r = go(s(), sdf=other.sdf())
    check("nessun template" in r.err, "numero di atomi diverso: fallback atteso")
    r = go(s(), sdf=other.sdf(name="altra") + m.sdf(name="giusta"))
    check("giusta" in r.err, "file multi-molecola: deve scegliere la molecola che combacia: " + r.err)


@test("tmpl")
def mol2_guanidinium_via_cationic_type():
    g = guanidinium(types=True)
    s = Sys().ligand(g).prot("GLU", "OE1", at((0, 1.33, 0), (0, 1, 0), 2.9))
    r = go(s, mol2=g.mol2())
    check("template" in r.err, "mol2 non applicato: " + r.err)
    check(r.has("Ponte salino", la=":N1"), "C.cat deve rendere cationici i tre N: " + str(r.types()))


@test("tmpl")
def mol2_aromatic_bonds_give_rings():
    b = benzene(types=True)
    s = Sys().ligand(b).ring("PHE", (1.0, 0, 3.6), (0, 0, 1))
    r = go(s, mol2=b.mol2())
    check("template" in r.err and r.has("pi-stacking (parallelo)"), "stacking con template mol2: " + r.err + str(r.types()))


# ----------------------------------------------------------------------------------------------
# GRUPPO pdb: parsing
# ----------------------------------------------------------------------------------------------
def lys_salt_system(**kw):
    m = carboxylate()
    return Sys().ligand(m), at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8)


@test("pdb")
def altloc_first_encountered_is_used():
    s, p = lys_salt_system()
    far = (30.0, 30.0, 30.0)
    s.prot("LYS", "NZ", p, alt="A").prot("LYS", "NZ", far, alt="B")
    check(go(s).has("Ponte salino"), "altloc A vicino, B lontano: atteso contatto")
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", far, alt="B").prot("LYS", "NZ", p, alt="A")
    check(not go(s).has("Ponte salino"), "si usa la prima altloc incontrata (B, lontana)")


@test("pdb")
def insertion_codes_are_distinct_residues():
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", p, resi=52, icode=" ")
    s.prot("LYS", "NZ", at((0.62, -1.07, 0), (0.5, -0.866, 0), 2.8), resi=52, icode="A")
    r = go(s)
    check(len(r.find("Ponte salino")) == 2, "attesi due ponti salini (52 e 52A): " + str([x["pa"] for x in r.rows]))
    check(any("52AA" in x["pa"] for x in r.rows), "il codice di inserzione deve comparire nell'etichetta")


@test("pdb")
def only_first_model_is_read():
    s, p = lys_salt_system()
    far = (30.0, 30.0, 30.0)
    txt = "MODEL        1\n" + s.prot("LYS", "NZ", far).text().replace("END\n", "") + "ENDMDL\n"
    s2, _ = lys_salt_system()
    txt += "MODEL        2\n" + s2.prot("LYS", "NZ", p).text() + "ENDMDL\nEND\n"
    check(not go(None, text=txt).has("Ponte salino"), "il secondo modello deve essere ignorato")


@test("pdb")
def ligand_selection_and_listing():
    s = Sys()
    s.ligand(carboxylate(), chain="L", resi=1)
    m2 = carboxylate()
    for name, el, xyz, _, _ in m2.atoms:
        s.atom("HETATM", name, "LIG", "M", 2, v_add(xyz, (40, 0, 0)), el)
    s.prot("LYS", "NZ", at((0.62, 1.07, 0), (0.5, 0.866, 0), 2.8))
    s.prot("LYS", "NZ", at(v_add((0.62, 1.07, 0), (40, 0, 0)), (0.5, 0.866, 0), 2.8), resi=11)
    check(len({r["lig"] for r in go(s).rows}) == 2, "senza filtro: due copie del ligando")
    r = go(s, lig=("LIG", "M", "2"))
    check({r_["lig"] for r_ in r.rows} == {"LIG 2M"}, "filtro catena/residuo: " + str({r_["lig"] for r_ in r.rows}))
    r = run_raw([tmpfile(s.text())])
    check("LIG(2)" in r.out, "elenco dei ligandi con conteggio copie: " + r.out)


@test("pdb")
def linked_hetero_residues_form_one_ligand():
    """Residui HETATM legati covalentemente (glicani, dipeptidi) = un solo ligando, come in PLIP."""
    def build(gap):
        s = Sys().ligand(methanol(), res="LIG", chain="L", resi=1)
        x = 1.43 + gap
        s.atom("HETATM", "C1", "MOL", "L", 2, (x, 0, 0), "C")
        s.atom("HETATM", "O1", "MOL", "L", 2, (x + 1.43, 0, 0), "O")
        s.prot("ASN", "OD1", at((x + 1.43, 0, 0), (0.5, 0.866, 0), 2.8))
        return s
    r = go(build(1.42))
    check("composto" in r.err, "unione non segnalata: " + r.err)
    check(r.has("Legame a idrogeno", la="MOL 2L:O1"), "gli atomi del residuo legato devono far parte del ligando: " + str(r.types()))
    check(not any("MOL 2L" in x["pa"] for x in r.rows), "il residuo legato e' stato contato anche come partner proteico")
    r = go(build(1.42), args=["--no-merge"])
    check(not r.has("Legame a idrogeno", la="MOL 2L:O1"), "--no-merge deve tenere separati i residui")
    r = go(build(3.2))
    check("composto" not in r.err and not any("MOL 2L" in x["la"] for x in r.rows), "residui non legati uniti per errore")


@test("pdb")
def csv_is_well_formed():
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", p)
    r = go(s)
    lines = r.out.strip().splitlines()
    check(lines[0] == "ligando,tipo,favorevole,atomo_ligando,atomo_proteina,distanza,note", "intestazione CSV: " + lines[0])
    for l in lines[1:]:
        check(len(l.split(",")) == 7, "riga CSV con numero di campi errato: " + l)
    check("Annotazione" not in r.out and "B-factor" not in r.out, "informazioni di servizio su stdout in modalita' CSV")


# ----------------------------------------------------------------------------------------------
# GRUPPO robust: input malformati
# ----------------------------------------------------------------------------------------------
def base_system():
    s, p = lys_salt_system()
    s.prot("LYS", "NZ", p)
    return s


@test("robust")
def missing_pdb_file():
    r = run_raw([os.path.join(CFG["tmp"], "non_esiste.pdb")])
    check(r.rc == 1 and "Impossibile" in r.err, "rc=%d err=%s" % (r.rc, r.err))


@test("robust")
def no_arguments_prints_usage():
    r = run_raw([])
    check(r.rc == 1 and "Uso" in r.err, "rc=%d" % r.rc)


@test("robust")
def empty_and_whitespace_files():
    for content in ["", "\n\n\n", "   \n"]:
        r = run_raw([tmpfile(content)])
        check(r.rc == 0, "file vuoto: rc=%d %s" % (r.rc, r.err[:200]))
        no_crash(r)


@test("robust")
def ligand_not_present():
    r = run_raw([tmpfile(base_system().text()), "XXX"])
    check(r.rc == 0 and "Nessun ligando" in r.out, r.out + r.err)


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
        check(r.rc == 0, "coordinata '%s': rc=%d %s" % (token.strip(), r.rc, r.err[:300]))
        no_crash(r)


@test("robust")
def coincident_atoms_zero_distance():
    s = Sys().ligand(carboxylate())
    s.prot("LYS", "NZ", (0.62, 1.07, 0))                  # esattamente sopra O1
    s.prot("ASP", "OD1", (0.62, -1.07, 0))
    r = go(s)
    no_crash(r)


@test("robust")
def crlf_line_endings_give_same_result():
    txt = base_system().text()
    a = go(None, text=txt)
    b = go(None, text=txt.replace("\n", "\r\n"))
    check(a.rows == b.rows and a.rows, "CRLF cambia il risultato")


@test("robust")
def missing_or_lowercase_element_column():
    txt = base_system().text()
    a = go(None, text=txt)
    nocol = "\n".join(l[:76] if l.startswith(("ATOM", "HETATM")) else l for l in txt.splitlines()) + "\n"
    b = go(None, text=nocol)
    check(sorted(x["tipo"] for x in a.rows) == sorted(x["tipo"] for x in b.rows) and b.rows,
          "senza colonna elemento: %s contro %s" % (a.types(), b.types()))
    low = "\n".join(l[:76] + l[76:].lower() if l.startswith(("ATOM", "HETATM")) else l for l in txt.splitlines()) + "\n"
    c = go(None, text=low)
    check(sorted(x["tipo"] for x in a.rows) == sorted(x["tipo"] for x in c.rows), "elemento minuscolo")


@test("robust")
def single_atom_ligand():
    s = Sys().prot("HIS", "NE2", (2.1, 0, 0), el="N").metal("Zn", (0, 0, 0))
    r = go(s, lig=("ZN",))
    no_crash(r)
    check(r.has("Coordinazione metallica"), "ligando di un solo atomo (Zn): " + str(r.types()))


@test("robust")
def ligand_with_only_hydrogens():
    s = Sys().prot("LYS", "NZ", (5, 5, 5))
    s.atom("HETATM", "H1", "LIG", "L", 1, (0, 0, 0), "H")
    r = run_raw([tmpfile(s.text()), "LIG"])
    check(r.rc == 0 and "Nessun ligando" in r.out, r.out)


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
    for args, what in [([pdb, "--sdf"], "--sdf senza file"), ([pdb, "--foo", "1"], "opzione sconosciuta"),
                       ([pdb, "LIG", "--hb-max"], "soglia senza valore"),
                       ([pdb, "LIG", "--hb-max", "abc"], "soglia non numerica"),
                       ([pdb, "LIG", "--sdf", "/non/esiste.sdf"], "file SDF inesistente"),
                       ([pdb, "LIG", "--mol2", "/non/esiste.mol2"], "file MOL2 inesistente")]:
        r = run_raw(args)
        check(r.rc == 1, "%s: atteso codice 1, ottenuto %d (%s)" % (what, r.rc, r.err[:120]))
        no_crash(r)


@test("robust")
def oversized_thresholds_are_clamped_not_silent():
    r = go(base_system(), args=["--hb-max", "9", "--salt-max", "12"])
    check("limitat" in r.err.lower() or "clamp" in r.err.lower() or "massimo" in r.err.lower(),
          "soglia oltre il limite della griglia (6 A) senza avviso: " + r.err[:200])


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
    check(CFG["last_time"] < limit, "200000 atomi: %.1f s (limite %d s)" % (CFG["last_time"], limit))
    no_crash(r)


# ----------------------------------------------------------------------------------------------
# GRUPPO golden: complessi reali (regressione di risultati verificati a mano)
# ----------------------------------------------------------------------------------------------
GOLDEN = [
    # (pdb, argomenti ligando, opzioni, [attesi (tipo, ligando contiene, proteina contiene)], [vietati (tipo, ...)])
    ("2OZR.pdb", ["GG1", "A", "2001"], [],
     [("Legame a idrogeno", ":O2", "THR 226A:N"), ("Ponte salino", ":O6", "LYS 119A:NZ"),
      ("pi-stacking (parallelo)", "", "PHE 231A")], []),
    ("1efy.pdb", ["BZC"], ["--sdf", "1efy.sdf"],
     [("Legame a idrogeno", ":N1", "GLY 863A:O"), ("Legame a idrogeno", ":O1", "GLY 863A:N"),
      ("pi-stacking (parallelo)", "", "TYR 907A")], []),
    ("4IKF.pdb", ["M76"], [],
     [("Coordinazione metallica", ":O14", "MG 405A"), ("Coordinazione metallica", ":O16", "MG 404A"),
      ("Coordinazione metallica", ":O18", "MG 404A"), ("pi-stacking (parallelo)", "", "DA 17D"),
      ("pi-stacking (parallelo)", "", "DC 16D")],
     [("Legame a idrogeno", "", ""), ("Coordinazione metallica", ":N15", "")]),
    ("3PTB.pdb", ["BEN"], ["--no-water"],
     [("Ponte salino", ":N1", "ASP 189A"), ("Ponte salino", ":N2", "ASP 189A"),
      ("Legame a idrogeno", ":N1", "GLY 219A:O"), ("Legame a idrogeno", ":N2", "SER 190A:OG")], []),
    ("3PTB.pdb", ["BEN"], ["--backbone-hyd"],
     [("Contatto apolare (backbone)", "", "TRP 215A:C"), ("Contatto apolare (backbone)", "", "CYS 191A:CA")], []),
    ("1STP.pdb", ["BTN"], ["--no-water"],
     [("Legame a idrogeno", ":N1", "ASP 128A"), ("Legame a idrogeno", ":N2", "SER 45A"),
      ("Legame a idrogeno", ":O11", "SER 88A"), ("Legame a idrogeno", ":O12", "ASN 49A:N"),
      ("Legame a idrogeno", ":O3", "TYR 43A"), ("Legame a idrogeno", ":O3", "SER 27A"),
      ("Legame a idrogeno", ":O3", "ASN 23A")], []),
    ("1IEP.pdb", ["STI", "A", "201"], ["--no-water", "--pi-angle", "40"],
     [("Legame a idrogeno", ":N13", "THR 315A"), ("Legame a idrogeno", ":N21", "GLU 286A"),
      ("Legame a idrogeno", ":N3", "MET 318A"), ("Legame a idrogeno", ":O29", "ASP 381A"),
      ("Legame a idrogeno", ":N51", "ILE 360A:O"), ("pi-stacking (parallelo)", "", "PHE 317A"),
      ("pi-stacking (parallelo)", "", "PHE 382A")],
     [("Cation-pi", "", "")]),
    ("1HSG.pdb", ["MK1"], ["--no-water"],
     [("Legame a idrogeno", ":O2", "ASP 25B"), ("Legame a idrogeno", ":O2", "ASP 25A"),
      ("Legame a idrogeno", ":N4", "GLY 27B:O"), ("Legame a idrogeno", ":O4", "ASP 29B:N"),
      ("Legame a idrogeno", ":O4", "GLY 27B:O")], []),
    ("1AZM.pdb", ["AZM"], ["--no-water"],
     [("Coordinazione metallica", ":N1", "ZN 261A"), ("Legame a idrogeno", ":N1", "THR 199A:OG1"),
      ("Legame a idrogeno", ":O2", "THR 199A:N")],
     [("Legame a idrogeno", "", "HIS 94A"), ("Legame a idrogeno", "", "HIS 96A"), ("Legame a idrogeno", "", "HIS 119A")]),
]


def make_golden(idx, entry):
    pdb, lig, opts, expect, forbid = entry

    def fn():
        path = os.path.join(CFG["real_dir"], pdb)
        if not os.path.exists(path):
            raise Skip("manca %s" % path)
        extra = []
        i = 0
        while i < len(opts):
            if opts[i] in ("--sdf", "--mol2"):
                p = os.path.join(CFG["real_dir"], opts[i + 1])
                if not os.path.exists(p):
                    raise Skip("manca %s" % p)
                extra += [opts[i], p]
                i += 2
            else:
                extra.append(opts[i]); i += 1
        r = run_raw([path] + lig + ["--csv"] + extra)
        check(r.rc == 0, "rc=%d %s" % (r.rc, r.err[:200]))
        no_crash(r)
        for tipo, la, pa in expect:
            check(r.has(tipo, la or None, pa or None), "atteso e non trovato: %s %s %s" % (tipo, la, pa))
        for tipo, la, pa in forbid:
            check(not r.has(tipo, la or None, pa or None), "trovato un contatto vietato: %s %s %s" % (tipo, la, pa))
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
    ap.add_argument("--bin", default="", help="binario di ligint gia' compilato (default: compila i sorgenti)")
    ap.add_argument("--sanitize", action="store_true")
    ap.add_argument("-k", default="", help="esegue solo i test con questa sottostringa nel nome")
    ap.add_argument("--group", default="", help="chem, tmpl, pdb, robust, golden (separati da virgola)")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--real-dir", default=os.path.join(HERE, "data"))
    ap.add_argument("-v", action="store_true", help="mostra il messaggio anche dei test riusciti/saltati")
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
            sys.exit("binario non trovato: " + CFG["bin"])
        print("binario: %s (precompilato)" % CFG["bin"])
    else:
        warn = _build.build_ligint(CFG["bin"], a.sanitize)
        print("binario: %s (%s)" % (CFG["bin"], "ASan+UBSan" if a.sanitize else "-O2"))
    if warn:
        print("avvisi del compilatore:\n" + warn)
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
        except Exception as e:   # errore nel test stesso
            fail += 1
            failures.append((g, n, "errore nel test: %r" % e))
            print("  ERROR %-7s %s\n          -> %r" % (g, n, e))
    print("\nrisultato: %d ok, %d falliti, %d saltati" % (ok, fail, skip))
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
