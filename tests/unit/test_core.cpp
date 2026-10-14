// Library tests (no external dependencies): API, thresholds, reports, templates, composite ligands.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ligint/ligint.hpp"

static int failures = 0;
#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            failures++;                                                                \
        }                                                                              \
    } while (0)

// ---- construction of synthetic PDBs -----------------------------------------------------------
static std::string atomLine(const char* rec, int serial, const std::string& name, const std::string& res, char chain, int resi,
                            double x, double y, double z, const std::string& el) {
    char buf[200];
    std::string nm = (name.size() == 4 || el.size() == 2) ? name : " " + name;
    std::snprintf(buf, sizeof buf, "%-6s%5d %-4s %3s %c%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s", rec, serial, nm.c_str(),
                  res.c_str(), chain, resi, x, y, z, el.c_str());
    return buf;
}

// methanol (C1 at 0,0,0; O1 at 1.43,0,0) as ligand LIG + an asparagine with OD1 2.8 A from O1 (C-O...O angle = 120)
static std::string methanolSystem(bool withProbe = true) {
    std::ostringstream s;
    s << atomLine("HETATM", 1, "C1", "LIG", 'L', 1, 0, 0, 0, "C") << "\n";
    s << atomLine("HETATM", 2, "O1", "LIG", 'L', 1, 1.43, 0, 0, "O") << "\n";
    if (withProbe) s << atomLine("ATOM", 3, "OD1", "ASN", 'A', 10, 1.43 + 2.8 * 0.5, 2.8 * 0.866, 0, "O") << "\n";
    s << "END\n";
    return s.str();
}

static ligint::Result run(const std::string& pdb, ligint::Options opt) {
    std::istringstream in(pdb);
    return ligint::analyzeStream(in, opt);
}

static bool hasType(const ligint::LigandResult& l, const std::string& type) {
    for (auto& h : l.hits) if (h.type == type) return true;
    return false;
}

// ---- tests ------------------------------------------------------------------------------------
static void testListing() {
    ligint::Result r = run(methanolSystem(), ligint::Options());          // no ligand requested
    CHECK(r.listingOnly);
    CHECK(r.hetCounts.count("LIG") == 1 && r.hetCounts["LIG"] == 1);
    CHECK(r.hetResidues.size() == 1 && r.hetResidues[0].atoms == 2 && r.hetResidues[0].chain == "L");
    ligint::Options opt;
    opt.ligName = "XXX";
    r = run(methanolSystem(), opt);
    CHECK(r.listingOnly && r.notFound == "XXX");
}

static void testHydrogenBondAndPositions() {
    ligint::Options opt;
    opt.ligName = "LIG";
    ligint::Result r = run(methanolSystem(), opt);
    CHECK(!r.listingOnly && r.ligands.size() == 1);
    const auto& l = r.ligands[0];
    CHECK(l.id == "LIG 1L" && l.atoms.size() == 2 && l.bonds.size() == 1);
    CHECK(hasType(l, ligint::type::HydrogenBond));
    for (auto& h : l.hits) {
        if (h.type != ligint::type::HydrogenBond) continue;
        CHECK(h.ligandAtom == "LIG 1L:O1" && h.proteinAtom == "ASN 10A:OD1" && h.residue == "ASN 10A");
        CHECK(std::fabs(h.distance - 2.8) < 0.01);
        CHECK(std::fabs(h.ligandPos[0] - 1.43) < 2e-3);                   // coordinates for the graphics
        CHECK(std::fabs(h.proteinPos[1] - 2.8 * 0.866) < 2e-3);
    }
    CHECK(!l.contactAtoms.empty());                                       // atoms of the residue in contact
}

static void testParams() {
    ligint::Params p;
    CHECK(p.set("hb-max", 3.8) && p.hbMax == 3.8);
    CHECK(!p.set("no-such-threshold", 1.0));
    ligint::Params q;
    q.applyPlipLike();
    CHECK(q.hbLong == 4.1 && q.saltCentroid == 5.5 && q.halMax == 4.0);
    ligint::Params c;
    c.hbMax = 9;
    auto w = c.clampToGrid();
    CHECK(c.hbMax == ligint::Params::kMaxGridRadius && w.size() == 1);
    CHECK(c.clampToGrid().empty());                                       // idempotent
}

static void testThresholdsChangeResult() {
    // acceptor at 3.9 A: beyond hb-max (3.5), admitted only with hb-long and a possible geometry
    std::ostringstream s;
    s << atomLine("HETATM", 1, "C1", "LIG", 'L', 1, 0, 0, 0, "C") << "\n"
      << atomLine("HETATM", 2, "O1", "LIG", 'L', 1, 1.43, 0, 0, "O") << "\n"
      << atomLine("ATOM", 3, "OD1", "ASN", 'A', 10, 1.43 + 3.9 * 0.342, 3.9 * 0.940, 0, "O") << "\n";   // C-O...A = 110 gradi
    ligint::Options opt;
    opt.ligName = "LIG";
    CHECK(!hasType(run(s.str(), opt).ligands[0], ligint::type::HydrogenBond));
    opt.params.applyPlipLike();
    CHECK(hasType(run(s.str(), opt).ligands[0], ligint::type::HydrogenBond));
}

static void testReports() {
    ligint::Options opt;
    opt.ligName = "LIG";
    ligint::Result r = run(methanolSystem(), opt);
    std::ostringstream text;
    ligint::printText(text, r);
    CHECK(text.str().find("=== Ligand LIG 1 chain L: 2 atoms, 0 rings ===") != std::string::npos);
    CHECK(text.str().find("Summary by type:") != std::string::npos);
    std::ostringstream csv, info;
    ligint::printCsv(csv, info, r);
    std::string out = csv.str();
    CHECK(out.rfind("ligand,type,favorable,ligand_atom,protein_atom,distance,notes\n", 0) == 0);
    size_t rows = 0;
    for (char ch : out) if (ch == '\n') rows++;
    CHECK(rows == 1 + r.ligands[0].hits.size());
    CHECK(info.str().find("Ligand annotation LIG 1L") != std::string::npos);
}

static void testMissingFile() {
    bool thrown = false;
    try {
        ligint::analyzeFile("/no/such.pdb", ligint::Options());
    } catch (const std::runtime_error&) {
        thrown = true;
    }
    CHECK(thrown);
}

static void testTemplateSdf() {
    const char* path = "ligint_test_methanol.sdf";
    {
        std::ofstream f(path);
        f << "methanol\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
             "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
             "    1.4300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
             "    1.7460    0.9100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
             "  1  2  1  0\n  2  3  1  0\nM  END\n$$$$\n";
    }
    ligint::Options opt;
    std::vector<std::string> warnings;
    CHECK(ligint::loadSdf(path, opt.templates, warnings) && opt.templates.size() == 1 && warnings.empty());
    CHECK(!ligint::loadSdf("/no/such.sdf", opt.templates, warnings));
    opt.ligName = "LIG";
    ligint::Result r = run(methanolSystem(), opt);
    CHECK(r.ligands[0].annotation.find("template 'methanol'") != std::string::npos);
    CHECK(!r.info.empty() && r.info[0] == "Templates loaded: 1");
    std::remove(path);
}

static void testLinkedResiduesAreMerged() {
    std::ostringstream s;
    s << atomLine("HETATM", 1, "C1", "LIG", 'L', 1, 0, 0, 0, "C") << "\n"
      << atomLine("HETATM", 2, "O1", "LIG", 'L', 1, 1.43, 0, 0, "O") << "\n"
      << atomLine("HETATM", 3, "C1", "MOL", 'L', 2, 2.85, 0, 0, "C") << "\n"          // bonded to O1 (1.42 A)
      << atomLine("HETATM", 4, "O1", "MOL", 'L', 2, 4.28, 0, 0, "O") << "\n"
      << atomLine("ATOM", 5, "OD1", "ASN", 'A', 10, 4.28 + 1.4, 2.42, 0, "O") << "\n";  // acceptor close to the second residue
    ligint::Options opt;
    opt.ligName = "LIG";
    ligint::Result r = run(s.str(), opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].residues.size() == 2 && r.ligands[0].atoms.size() == 4);
    CHECK(hasType(r.ligands[0], ligint::type::HydrogenBond));
    opt.mergeLinked = false;
    r = run(s.str(), opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].residues.empty() && r.ligands[0].atoms.size() == 2);
}

static void testMalformedInputDoesNotCrash() {
    ligint::Options opt;
    opt.ligName = "LIG";
    for (const std::string& txt : {std::string(), std::string("\n\n"), std::string("ATOM  \nHETATM"),
                                   atomLine("HETATM", 1, "C1", "LIG", 'L', 1, 1e30, 0, 0, "C"),
                                   std::string(5000, 'x')}) {
        ligint::Result r = run(txt, opt);
        CHECK(r.ligands.empty());
    }
}

int main() {
    testListing();
    testHydrogenBondAndPositions();
    testParams();
    testThresholdsChangeResult();
    testReports();
    testMissingFile();
    testTemplateSdf();
    testLinkedResiduesAreMerged();
    testMalformedInputDoesNotCrash();
    if (failures) {
        std::cerr << failures << " checks failed\n";
        return 1;
    }
    std::cout << "ligint_core_tests: all tests passed\n";
    return 0;
}
