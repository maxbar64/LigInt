// Library tests (no external dependencies): API, thresholds, reports, templates, composite ligands.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ligint/fingerprint.hpp"
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

static void testInsertionCodeSelection() {
    // two residues with the same name, chain and number, differing only in the insertion code
    std::ostringstream s;
    s << atomLine("HETATM", 1, "C1", "LIG", 'L', 5, 0, 0, 0, "C") << "\n";
    std::string a = atomLine("HETATM", 2, "C1", "LIG", 'L', 5, 30, 0, 0, "C");
    a[26] = 'A';
    s << a << "\n";
    ligint::Options opt;
    opt.ligName = "LIG"; opt.chain = "L"; opt.resi = 5;
    CHECK(run(s.str(), opt).ligands.size() == 2);          // any insertion code
    opt.icode = 'A';
    ligint::Result r = run(s.str(), opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].id == "LIG 5AL");
    opt.icode = ' ';
    r = run(s.str(), opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].id == "LIG 5L");
    opt.icode = 'B';
    CHECK(run(s.str(), opt).ligands.empty());
}

static std::string withAltloc(std::string line, char alt) {
    line[16] = alt;
    return line;
}

static void testAlternateConformationSelection() {
    // methanol ligand; ASN 10 has an OD1 far away in conformation A and 2.8 A from O1 in conformation B;
    // SER 11 has only conformation A (far away); GLY 12 has no alternate conformations
    std::ostringstream s;
    s << atomLine("HETATM", 1, "C1", "LIG", 'L', 1, 0, 0, 0, "C") << "\n"
      << atomLine("HETATM", 2, "O1", "LIG", 'L', 1, 1.43, 0, 0, "O") << "\n"
      << withAltloc(atomLine("ATOM", 3, "OD1", "ASN", 'A', 10, 9, 0, 0, "O"), 'A') << "\n"
      << withAltloc(atomLine("ATOM", 4, "OD1", "ASN", 'A', 10, 1.43 + 2.8 * 0.5, 2.8 * 0.866, 0, "O"), 'B') << "\n"
      << withAltloc(atomLine("ATOM", 5, "OG", "SER", 'A', 11, 20, 0, 0, "O"), 'A') << "\n"
      << atomLine("ATOM", 6, "O", "GLY", 'A', 12, 30, 0, 0, "O") << "\n";
    ligint::Options opt;
    opt.ligName = "LIG";
    ligint::Result r = run(s.str(), opt);                                     // default: the first conformation (A)
    CHECK(r.altlocs.size() == 2 && r.altlocs[0] == 'A' && r.altlocs[1] == 'B');
    CHECK(r.ligands.size() == 1 && r.ligands[0].hits.empty() && r.warnings.empty());
    opt.altloc = 'A';
    CHECK(run(s.str(), opt).ligands[0].hits.empty());
    opt.altloc = 'B';
    r = run(s.str(), opt);                                                    // conformation B: the H-bond appears
    CHECK(r.ligands[0].hits.size() == 1);
    if (r.ligands[0].hits.size() == 1) {
        CHECK(r.ligands[0].hits[0].type == ligint::type::HydrogenBond);
        CHECK(r.ligands[0].hits[0].proteinAtom == "ASN 10A:OD1" && std::fabs(r.ligands[0].hits[0].distance - 2.8) < 0.01);
    }
    bool info = false;
    for (auto& i : r.info) if (i.find("Alternate conformation B: used for 1 of 2 residues") != std::string::npos) info = true;
    CHECK(info);                                                              // SER 11 has no B and keeps A: reported
    opt.altloc = 'C';                                                         // not present anywhere: warning, same as the default
    r = run(s.str(), opt);
    CHECK(r.warnings.size() == 1 && r.warnings[0].find("alternate conformation C not found") != std::string::npos &&
          r.warnings[0].find("present: A, B") != std::string::npos && r.ligands[0].hits.empty());
    // a structure without alternate conformations
    opt.altloc = 'B';
    r = run(methanolSystem(), opt);
    CHECK(r.warnings.size() == 1 && r.warnings[0].find("no alternate conformations") != std::string::npos && !r.ligands[0].hits.empty());
    // the ligand can have alternate conformations too, and the first model only is read
    std::ostringstream t;
    t << withAltloc(atomLine("HETATM", 1, "C1", "LIG", 'L', 1, 0, 0, 0, "C"), 'A') << "\n"
      << withAltloc(atomLine("HETATM", 2, "O1", "LIG", 'L', 1, 1.43, 0, 0, "O"), 'A') << "\n"
      << withAltloc(atomLine("HETATM", 3, "C1", "LIG", 'L', 1, 50, 0, 0, "C"), 'B') << "\n"
      << withAltloc(atomLine("HETATM", 4, "O1", "LIG", 'L', 1, 51.43, 0, 0, "O"), 'B') << "\n"
      << atomLine("ATOM", 5, "OD1", "ASN", 'A', 10, 1.43 + 2.8 * 0.5, 2.8 * 0.866, 0, "O") << "\n"
      << "ENDMDL\n"
      << withAltloc(atomLine("ATOM", 6, "OD1", "ASN", 'A', 10, 0, 0, 0, "O"), 'C') << "\n";
    opt.altloc = 'A';
    CHECK(run(t.str(), opt).ligands[0].hits.size() == 1);                     // ligand in conformation A, next to the ASN
    opt.altloc = 'B';
    CHECK(run(t.str(), opt).ligands[0].hits.empty());                         // ligand in conformation B, far from it
    opt.altloc = 'C';
    CHECK(!run(t.str(), opt).warnings.empty());                               // C is after ENDMDL: not part of the first model
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


// ---- fingerprints -----------------------------------------------------------------------------
// acetate (carboxylate) as ligand + a lysine NZ 2.8 A from O1, in the direction away from the carboxylate
static std::string carboxylateLysSystem() {
    std::ostringstream s;
    s << atomLine("HETATM", 1, "C1", "LIG", 'L', 1, -1.52, 0, 0, "C") << "\n"
      << atomLine("HETATM", 2, "C2", "LIG", 'L', 1, 0, 0, 0, "C") << "\n"
      << atomLine("HETATM", 3, "O1", "LIG", 'L', 1, 0.62, 1.07, 0, "O") << "\n"
      << atomLine("HETATM", 4, "O2", "LIG", 'L', 1, 0.62, -1.07, 0, "O") << "\n"
      << atomLine("ATOM", 5, "NZ", "LYS", 'A', 30, 0.62 + 2.8 * 0.5, 1.07 + 2.8 * 0.866, 0, "N") << "\n";
    return s.str();
}

static void testFingerprintClasses() {
    ligint::Options opt;
    opt.ligName = "LIG";
    // ligand hydroxyl donates to an asparagine O: the RESIDUE is the acceptor
    ligint::Result r = run(methanolSystem(), opt);
    ligint::Fingerprint fp = ligint::makeFingerprint(r.ligands[0]);
    CHECK(fp.ligandId == "LIG 1L" && fp.entries.size() == 1);
    CHECK(fp.entries[0].residue == "ASN 10A" && fp.entries[0].inPocket && fp.entries[0].position.empty());
    CHECK(fp.entries[0].counts[ligint::FpHBondAcceptor] == 1 && fp.entries[0].counts[ligint::FpHBondDonor] == 0);
    CHECK(std::fabs(fp.entries[0].minDistance[ligint::FpHBondAcceptor] - 2.8) < 0.01);
    // carboxylate + lysine: ionic positive (the residue is positive) and H-bond donor (the lysine donates)
    r = run(carboxylateLysSystem(), opt);
    fp = ligint::makeFingerprint(r.ligands[0]);
    CHECK(fp.entries.size() == 1);
    CHECK(fp.entries[0].counts[ligint::FpIonicPositive] == 1 && fp.entries[0].counts[ligint::FpIonicNegative] == 0);
    CHECK(fp.entries[0].counts[ligint::FpHBondDonor] == 1);
    // the hit records carry the information used by the fingerprint
    bool seen = false;
    for (auto& h : r.ligands[0].hits)
        if (h.type == ligint::type::SaltBridge) { seen = true; CHECK(h.residueCharge == 1 && h.hbResidueRole == 'D'); }
    CHECK(seen);
    // clashes and repulsions do not enter the fingerprint
    ligint::LigandResult fake;
    fake.id = "X";
    ligint::HitRecord bad;
    bad.type = ligint::type::StericClash; bad.residue = "ALA 5A"; bad.favorable = false; bad.distance = 2.0;
    fake.hits.push_back(bad);
    bad.type = ligint::type::ElectrostaticRepulsion;
    fake.hits.push_back(bad);
    CHECK(ligint::makeFingerprint(fake).entries.empty());
    CHECK(std::string(ligint::fpClassName(ligint::FpHBondDonor)) == "hbond_donor" && ligint::oneLetterCode("TRP") == 'W' && ligint::oneLetterCode("XYZ") == '?');
}

static void testPocketDefinition() {
    std::istringstream in(
        "# comment\n"
        "@name Test pocket\n"
        "@color hinge #ff0000\n"
        "@color broken nocolor\n"
        "@bogus x\n"
        "A\t10\t1\thinge\tASN\n"          // tab separated, with the residue name
        "A,11A,2,hinge,gly\n"                 // comma separated, insertion code, lowercase name
        "*  12  3  gk\n"                      // whitespace separated, any chain
        "A abc 4 x\n"                         // malformed residue number
        "A 13\n"                              // too few fields
        "\n");
    ligint::PocketDefinition def;
    std::vector<std::string> warnings;
    CHECK(ligint::parsePocketDefinition(in, def, warnings));
    CHECK(def.name == "Test pocket" && def.colors.size() == 1 && def.colors[0].first == "hinge");
    CHECK(def.positions.size() == 3 && warnings.size() == 4);
    CHECK(def.positions[0].chain == "A" && def.positions[0].resi == 10 && def.positions[0].resName == "ASN" && def.positions[0].region == "hinge");
    CHECK(def.positions[1].icode == 'A' && def.positions[1].resName == "GLY");
    CHECK(def.positions[2].chain.empty() && def.positions[2].label == "3" && def.positions[2].resName.empty());
    CHECK(!ligint::loadPocketDefinition("/no/such/pocket.tsv", def, warnings));
    // binary garbage must not crash
    std::string junk(4000, '\0');
    for (size_t i = 0; i < junk.size(); i++) junk[i] = static_cast<char>((i * 37) % 251);
    std::istringstream g(junk);
    ligint::PocketDefinition d2;
    std::vector<std::string> w2;
    ligint::parsePocketDefinition(g, d2, w2);
}

static void testFingerprintWithPocket() {
    ligint::Options opt;
    opt.ligName = "LIG";
    ligint::Result r = run(methanolSystem(), opt);       // contact: ASN 10 chain A
    ligint::PocketDefinition def;
    def.name = "P";
    def.positions = {{"A", 9, ' ', "p9", "reg1", "ALA"}, {"", 10, ' ', "p10", "reg1", "ASN"}, {"A", 11, ' ', "p11", "reg2", "GLY"}};
    ligint::Fingerprint fp = ligint::makeFingerprint(r.ligands[0], &def);
    CHECK(fp.entries.size() == 3 && fp.pocketName == "P");                     // every position, also empty ones
    CHECK(!fp.entries[0].any() && fp.entries[1].any() && !fp.entries[2].any());
    CHECK(fp.entries[1].position == "p10" && fp.entries[1].region == "reg1" && fp.entries[1].chain == "A");
    CHECK(fp.warnings.empty());
    // a contact residue outside the definition is appended and flagged
    def.positions = {{"A", 9, ' ', "p9", "reg1", "ALA"}};
    fp = ligint::makeFingerprint(r.ligands[0], &def);
    CHECK(fp.entries.size() == 2 && !fp.entries[1].inPocket && fp.entries[1].residue == "ASN 10A");
    // residue name mismatch is reported
    def.positions = {{"A", 10, ' ', "p10", "reg1", "LEU"}};
    fp = ligint::makeFingerprint(r.ligands[0], &def);
    CHECK(fp.entries.size() == 1 && fp.warnings.size() == 1 && fp.warnings[0].find("does not match") != std::string::npos);
}

static void testFingerprintOutputs() {
    ligint::Options opt;
    opt.ligName = "LIG";
    ligint::Result r = run(carboxylateLysSystem(), opt);
    ligint::PocketDefinition def;
    def.name = "a&b <pocket>";
    def.positions = {{"A", 30, ' ', "p<30>", "reg&1", "LYS"}, {"A", 31, ' ', "p31", "reg&1", "GLU"}};
    ligint::Fingerprint fp = ligint::makeFingerprint(r.ligands[0], &def);
    std::ostringstream csv;
    ligint::writeFingerprintCsv(csv, {fp});
    std::string c = csv.str();
    size_t lines = 0, commas = 0;
    for (char ch : c.substr(0, c.find('\n'))) if (ch == ',') commas++;
    for (char ch : c) if (ch == '\n') lines++;
    CHECK(lines == 1 + fp.entries.size() && commas == 5 + ligint::FpClassCount);
    CHECK(c.find("LIG 1L,LYS 30A,A,30,p<30>,reg&1,") != std::string::npos);
    std::ostringstream js;
    ligint::writeFingerprintJson(js, {fp});
    CHECK(js.str().front() == '[' && js.str().find("\"hbond_donor\": 1") != std::string::npos && js.str().find("\"ionic_positive\"") != std::string::npos);
    CHECK(js.str().find("a&b <pocket>") != std::string::npos);
    std::string svg = ligint::fingerprintSvg(fp);
    CHECK(svg.rfind("<?xml", 0) == 0 && svg.find("</svg>") != std::string::npos);
    CHECK(svg.find("a&amp;b &lt;pocket&gt;") != std::string::npos && svg.find("reg&amp;1") != std::string::npos);   // XML escaping
    CHECK(svg.find("a&b") == std::string::npos);
    CHECK(svg.find("Hydrophobic") != std::string::npos && svg.find("Ionic positive") != std::string::npos);        // legend
    ligint::SvgOptions so;
    so.showLegend = false;
    CHECK(ligint::fingerprintSvg(fp, so).find("Ionic positive") == std::string::npos);
    CHECK(ligint::fingerprintSvg(ligint::Fingerprint()).find("No interactions") != std::string::npos);
    // wrapping: many residues in one region fill several rows
    ligint::Fingerprint big;
    for (int i = 0; i < 50; i++) { ligint::FingerprintEntry e; e.residue = "ALA " + std::to_string(i) + "A"; e.resName = "ALA"; e.resi = i; e.region = "r"; big.entries.push_back(e); }
    so.cellsPerRow = 10;
    std::string bs = ligint::fingerprintSvg(big, so);
    size_t headers = 0;
    for (size_t p = 0; (p = bs.find(">r</text>", p)) != std::string::npos; p++) headers++;
    CHECK(headers == 5);
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
    testInsertionCodeSelection();
    testAlternateConformationSelection();
    testFingerprintClasses();
    testPocketDefinition();
    testFingerprintWithPocket();
    testFingerprintOutputs();
    if (failures) {
        std::cerr << failures << " checks failed\n";
        return 1;
    }
    std::cout << "ligint_core_tests: all tests passed\n";
    return 0;
}
