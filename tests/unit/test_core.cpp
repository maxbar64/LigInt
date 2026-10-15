// Library tests (no external dependencies): API, thresholds, reports, templates, composite ligands.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ligint/fingerprint.hpp"
#include "ligint/klifs.hpp"
#include "ligint/rcsb.hpp"
#include <filesystem>
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


// ---- mmCIF input -----------------------------------------------------------------------------
struct CifRow {
    const char* group; const char* el; const char* atom; const char* alt; const char* comp; const char* chain;
    int resi; const char* ins; double x, y, z; int model;
};

// a PDBx/mmCIF file as the PDB writes it: label_* columns different from the auth_* ones on purpose, so that a
// reader that used the wrong columns would be caught
static std::string cifText(const std::vector<CifRow>& rows, const std::string& before = "", const std::string& after = "") {
    std::ostringstream s;
    s << "data_TEST\n#\n_entry.id TEST\n" << before
      << "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n"
         "_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
         "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n"
         "_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n"
         "_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n";
    int id = 0;
    for (const auto& r : rows) {
        char b[512];
        ++id;
        const std::string atom = std::string(r.atom).find('\'') != std::string::npos ? std::string("\"") + r.atom + "\"" : r.atom;
        std::snprintf(b, sizeof b, "%s %d %s %s %s %s %sz 1 %d %s %.3f %.3f %.3f 1.00 20.00 ? %d %s %s %s %d\n", r.group, id, r.el,
                      atom.c_str(), r.alt, r.comp, r.chain, 900 + id, r.ins, r.x, r.y, r.z, r.resi, r.comp, r.chain, atom.c_str(), r.model);
        s << b;
    }
    s << after;
    return s.str();
}

static std::vector<CifRow> cifMethanol(const char* chain = "A", const char* comp = "ASN") {
    const double ox = 1.43 + 2.8 * 0.5, oy = 2.8 * 0.866;
    return {{"HETATM", "C", "C1", ".", "LIG", "L", 1, "?", 0, 0, 0, 1},
            {"HETATM", "O", "O1", ".", "LIG", "L", 1, "?", 1.43, 0, 0, 1},
            {"ATOM", "O", "OD1", ".", comp, chain, 10, "?", ox, oy, 0, 1}};
}

static void testMmcifBasic() {
    ligint::Options opt;
    opt.ligName = "LIG";
    ligint::Result r = run(cifText(cifMethanol()), opt);                        // format recognized from the content
    CHECK(r.format == ligint::StructureFormat::Mmcif && r.warnings.empty());
    CHECK(r.ligands.size() == 1 && r.ligands[0].id == "LIG 1L" && r.ligands[0].hits.size() == 1);
    CHECK(r.ligands[0].hits[0].type == ligint::type::HydrogenBond && r.ligands[0].hits[0].proteinAtom == "ASN 10A:OD1");
    CHECK(std::fabs(r.ligands[0].hits[0].distance - 2.8) < 0.01);               // chain and number are the auth_* ones
    // same system as a PDB file: same result
    ligint::Result p = run(methanolSystem(), opt);
    CHECK(p.format == ligint::StructureFormat::Pdb && p.ligands[0].hits.size() == 1 &&
          p.ligands[0].hits[0].proteinAtom == r.ligands[0].hits[0].proteinAtom);
    // the listing
    ligint::Options none;
    ligint::Result l = run(cifText(cifMethanol()), none);
    CHECK(l.listingOnly && l.hetCounts.count("LIG") == 1 && l.hetCounts.at("LIG") == 1);
    // CRLF line endings and lower/upper case tags and data block keywords
    std::string crlf = cifText(cifMethanol()), up = crlf;
    for (size_t i = 0; (i = crlf.find('\n', i)) != std::string::npos; i += 2) crlf.replace(i, 1, "\r\n");
    CHECK(run(crlf, opt).ligands[0].hits.size() == 1);
    for (size_t i = 0; (i = up.find("_atom_site.", i)) != std::string::npos; i += 11) up.replace(i, 11, "_ATOM_SITE.");
    up.replace(up.find("loop_"), 5, "LOOP_");
    up.replace(up.find("data_"), 5, "DATA_");
    CHECK(run(up, opt).ligands[0].hits.size() == 1);
}

static void testMmcifNamesAndChains() {
    ligint::Options opt;
    opt.ligName = "ABCDE";                                                      // a five-character component id
    std::vector<CifRow> rows = {{"HETATM", "C", "C1", ".", "ABCDE", "L", 1, "?", 0, 0, 0, 1},
                                {"HETATM", "O", "O1", ".", "ABCDE", "L", 1, "?", 1.43, 0, 0, 1},
                                {"ATOM", "O", "OD1", ".", "ASN", "AA", 10, "?", 1.43 + 1.4, 2.42, 0, 1}};   // chain of two characters
    ligint::Result r = run(cifText(rows), opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].id == "ABCDE 1L" && r.ligands[0].resName == "ABCDE");
    CHECK(r.ligands[0].hits.size() == 1 && r.ligands[0].hits[0].residue == "ASN 10 [AA]");   // the label stays unambiguous
    CHECK(r.ligands[0].hits[0].proteinAtom == "ASN 10 [AA]:OD1");
    ligint::Fingerprint fp = ligint::makeFingerprint(r.ligands[0]);
    CHECK(fp.entries.size() == 1 && fp.entries[0].chain == "AA" && fp.entries[0].resi == 10 && fp.entries[0].icode == ' ');
    // insertion code + chain of two characters
    rows[2].ins = "B";
    r = run(cifText(rows), opt);
    CHECK(r.ligands[0].hits[0].residue == "ASN 10B [AA]");
    fp = ligint::makeFingerprint(r.ligands[0]);
    CHECK(fp.entries.size() == 1 && fp.entries[0].chain == "AA" && fp.entries[0].icode == 'B');
    // a ligand with a chain of two characters is selected by that chain
    rows[0].chain = rows[1].chain = "LL";
    opt.chain = "LL";
    CHECK(run(cifText(rows), opt).ligands.size() == 1 && run(cifText(rows), opt).ligands[0].id == "ABCDE 1 [LL]");
    opt.chain = "L";
    CHECK(run(cifText(rows), opt).ligands.empty());
    // atom names with an apostrophe (nucleic acids) are quoted in the file
    std::vector<CifRow> dna = {{"ATOM", "O", "O5'", ".", "DA", "D", 1, "?", 0, 0, 0, 1}, {"ATOM", "C", "C5'", ".", "DA", "D", 1, "?", 1.4, 0, 0, 1},
                               {"HETATM", "C", "C1", ".", "LIG", "L", 1, "?", 30, 0, 0, 1}};
    ligint::Options o2;
    o2.ligName = "LIG";
    ligint::Result d = run(cifText(dna), o2);
    CHECK(d.ligands.size() == 1 && d.warnings.empty());
}

static void testMmcifSyntax() {
    ligint::Options opt;
    opt.ligName = "LIG";
    // text fields, quoted values, comments and other loops before the atoms; another category and a second data block after them
    const std::string before =
        "_struct.title\n;A title: with 'quotes', a # and\nloop_\n_atom_site.group_PDB\n;\n"
        "# a comment between categories\n"
        "loop_\n_citation.id\n_citation.title\n_citation.note\nprimary 'Title with # hash and a quote's' 'x'\n"
        "2 \"double \"quoted\" text\" ?\n3\n;a text\nvalue\n;\n.\n"
        "_cell.length_a 81.85\n";
    const std::string after = "#\nloop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.U11\n1 0.1\n2 0.2\ndata_SECOND\n_entry.id SECOND\n";
    ligint::Result r = run(cifText(cifMethanol(), before, after), opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].hits.size() == 1 && r.warnings.empty());
    // tags in another order, and only the label_* columns for the residue (no auth_*)
    const std::string noAuth =
        "data_X\nloop_\n_atom_site.Cartn_z\n_atom_site.Cartn_y\n_atom_site.Cartn_x\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.group_PDB\n_atom_site.type_symbol\n"
        "0 0 0 C1 LIG L 1 HETATM C\n0 0 1.43 O1 LIG L 1 HETATM O\n0 2.42 2.83 OD1 ASN A 10 ATOM O\n";
    r = run(noAuth, opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].hits.size() == 1 && r.ligands[0].hits[0].proteinAtom == "ASN 10A:OD1");
    // no element column: derived from the atom name; no group column: nothing is HETATM, so no ligand
    const std::string noEl =
        "data_X\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
        "_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "HETATM C1 LIG L 1 0 0 0\nHETATM O1 LIG L 1 1.43 0 0\nATOM OD1 ASN A 10 2.83 2.42 0\n";
    r = run(noEl, opt);
    CHECK(r.ligands.size() == 1 && r.ligands[0].hits.size() == 1);
    // an _atom_site category written as pairs (a one-atom file)
    const std::string single = "data_X\n_atom_site.group_PDB HETATM\n_atom_site.type_symbol ZN\n_atom_site.label_atom_id ZN\n_atom_site.label_comp_id ZN\n"
                               "_atom_site.auth_asym_id A\n_atom_site.auth_seq_id 5\n_atom_site.Cartn_x 1.0\n_atom_site.Cartn_y 2.0\n_atom_site.Cartn_z 3.0\n";
    ligint::Options none;
    r = run(single, none);
    CHECK(r.listingOnly && r.hetCounts.count("ZN") == 1);
}

static void testMmcifModelsAndAltlocs() {
    ligint::Options opt;
    opt.ligName = "LIG";
    // only the first model is read
    std::vector<CifRow> rows = cifMethanol();
    rows.push_back({"HETATM", "C", "C1", ".", "LIG", "L", 2, "?", 0, 0, 0, 2});     // a second ligand, only in model 2
    ligint::Result r = run(cifText(rows), opt);
    CHECK(r.ligands.size() == 1);
    // alternate conformations: ASN 10 has OD1 far away in A and 2.8 A from O1 in B; SER 11 has only A
    std::vector<CifRow> alt = {{"HETATM", "C", "C1", ".", "LIG", "L", 1, "?", 0, 0, 0, 1},
                               {"HETATM", "O", "O1", ".", "LIG", "L", 1, "?", 1.43, 0, 0, 1},
                               {"ATOM", "O", "OD1", "A", "ASN", "A", 10, "?", 9, 0, 0, 1},
                               {"ATOM", "O", "OD1", "B", "ASN", "A", 10, "?", 1.43 + 1.4, 2.42, 0, 1},
                               {"ATOM", "O", "OG", "A", "SER", "A", 11, "?", 20, 0, 0, 1}};
    r = run(cifText(alt), opt);
    CHECK(r.altlocs.size() == 2 && r.altlocs[0] == 'A' && r.altlocs[1] == 'B' && r.ligands[0].hits.empty());
    opt.altloc = 'B';
    r = run(cifText(alt), opt);
    CHECK(r.ligands[0].hits.size() == 1 && r.ligands[0].hits[0].proteinAtom == "ASN 10A:OD1");
    bool info = false;
    for (auto& i : r.info) if (i.find("Alternate conformation B: used for 1 of 2 residues") != std::string::npos) info = true;
    CHECK(info);
    opt.altloc = 'C';
    r = run(cifText(alt), opt);
    CHECK(r.warnings.size() == 1 && r.warnings[0].find("alternate conformation C not found") != std::string::npos);
}

static void testMmcifMalformed() {
    ligint::Options opt;
    opt.ligName = "LIG";
    // incomplete last row
    std::string text = cifText(cifMethanol());
    text.resize(text.size() - 20);
    ligint::Result r = run(text, opt);
    bool partial = false;
    for (auto& w : r.warnings) if (w.find("incomplete last row") != std::string::npos) partial = true;
    CHECK(partial && r.ligands.size() == 1);
    // no atoms at all, a category without coordinates, and text that is not a structure
    r = run("data_X\n_entry.id X\n", opt);
    CHECK(r.listingOnly && !r.warnings.empty() && r.warnings[0].find("no _atom_site") != std::string::npos);
    r = run("data_X\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\nATOM CA\n", opt);
    CHECK(r.listingOnly && !r.warnings.empty());
    r = run("", opt);
    CHECK(r.listingOnly);
    // missing coordinates ("?") count as invalid coordinates, rows without a name are counted
    std::vector<CifRow> rows = cifMethanol();
    std::string t = cifText(rows);
    const size_t p = t.rfind("ASN");
    (void)p;
    rows.push_back({"ATOM", "C", "CA", ".", "ALA", "A", 11, "?", 0, 0, 0, 1});
    t = cifText(rows);
    const size_t q = t.rfind("0.000 0.000 0.000");
    t.replace(q, 5, "?");
    r = run(t, opt);
    bool bad = false;
    for (auto& w : r.warnings) if (w.find("invalid coordinates") != std::string::npos) bad = true;
    CHECK(bad);
    // pseudo-random garbage and truncations of a valid file never crash
    std::string valid = cifText(cifMethanol(), "_struct.title ;x;\n");
    unsigned seed = 12345;
    for (int n = 0; n < 300; n++) {
        std::string g = valid;
        for (int k = 0; k < 6; k++) {
            seed = seed * 1103515245u + 12345u;
            const size_t pos = (seed >> 8) % g.size();
            seed = seed * 1103515245u + 12345u;
            const int mode = (seed >> 8) % 4;
            if (mode == 0) g[pos] = static_cast<char>(seed >> 16);
            else if (mode == 1) g.erase(pos, 1 + (seed >> 12) % 8);
            else if (mode == 2) g.insert(pos, "\n;");
            else g.resize(pos);
            if (g.empty()) g = "d";
        }
        run(g, opt);
    }
    // the format can be forced
    opt.format = ligint::StructureFormat::Pdb;
    CHECK(run(cifText(cifMethanol()), opt).ligands.empty());
    opt.format = ligint::StructureFormat::Mmcif;
    ligint::Result f = run(methanolSystem(), opt);
    CHECK(f.format == ligint::StructureFormat::Mmcif && f.listingOnly && !f.warnings.empty());
    // a gzip file is refused with a clear message
    const std::string gz = "/tmp/ligint_test_compressed.cif.gz";
    {
        std::ofstream out(gz, std::ios::binary);
        out.put(static_cast<char>(0x1f)); out.put(static_cast<char>(0x8b)); out << "garbage";
    }
    bool refused = false;
    try { ligint::analyzeFile(gz, ligint::Options()); } catch (const std::runtime_error& e) { refused = std::string(e.what()).find("compressed") != std::string::npos; }
    CHECK(refused);
    std::remove(gz.c_str());
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

static void testKlifsStructureSelection() {
    ligint::KlifsOptions o;
    o.pdbId = "3EYG";
    o.chain = "A";
    // an explicit structure ID needs no answer at all
    o.structureId = 7;
    CHECK(ligint::selectKlifsStructure("this is not json", o) == 7);
    o.structureId = 0;
    if (!ligint::klifsOnlineAvailable()) {                                   // a build without network support cannot read the answer
        bool thrown = false;
        try { ligint::selectKlifsStructure("[]", o); } catch (const std::runtime_error&) { thrown = true; }
        CHECK(thrown);
        return;
    }
    auto message = [&](const std::string& json, const ligint::KlifsOptions& opt) {
        try { ligint::selectKlifsStructure(json, opt); } catch (const std::exception& e) { return std::string(e.what()); }
        return std::string("(no exception)");
    };
    // the chain chooses
    CHECK(ligint::selectKlifsStructure(R"([{"structure_ID":1445,"chain":"A","alt":"A"},{"structure_ID":1446,"chain":"B","alt":"A"}])", o) == 1445);
    CHECK(ligint::selectKlifsStructure(R"({"structure_ID":9,"chain":"a"})", o) == 9);                          // one object instead of a list, chain in lower case
    CHECK(ligint::selectKlifsStructure(R"([{"structure_ID":"1445","chain":"A"}])", o) == 1445);                // the ID as a string
    // two alternate models of the same chain: ambiguous, unless the model is given; the model can be named in several ways
    const std::string twoModels = R"([{"structure_ID":1445,"chain":"A","alt":"A"},{"structure_ID":1453,"chain":"A","alt":"B"}])";
    try {
        ligint::selectKlifsStructure(twoModels, o);
        CHECK(false);
    } catch (const ligint::KlifsAmbiguousError& e) {
        CHECK(e.candidates.size() == 2 && e.candidates[0].id == 1445 && e.candidates[1].id == 1453 && e.candidates[1].model == "B");
        const std::string m = e.what();
        CHECK(m.find("--klifs-id") != std::string::npos && m.find("--klifs-model") != std::string::npos && m.find("model=B") != std::string::npos);
    }
    o.model = "b";
    CHECK(ligint::selectKlifsStructure(twoModels, o) == 1453);
    CHECK(ligint::selectKlifsStructure(R"([{"structure_ID":1,"chain":"A","alternate_model":"A"},{"structure_ID":2,"chain":"A","alternate_model":"B"}])", o) == 2);
    CHECK(ligint::selectKlifsStructure(R"([{"structure_ID":1,"chain":"A","structure.alternate_model":"A"},{"structure_ID":2,"chain":"A","structure.alternate_model":"B"}])", o) == 2);
    o.model = "C";
    CHECK(message(twoModels, o).find("no structure") != std::string::npos);                                   // that model does not exist
    o.model.clear();
    // the situation of a first real run: two structures of chain A and no field that names the model: the message must say what the answer contains
    const std::string unknownModel = R"([{"structure_ID":1445,"kinase":"X","pdb":"3eyg","chain":"A","rmsd1":1.2},{"structure_ID":1453,"kinase":"X","pdb":"3eyg","chain":"A","rmsd1":1.3}])";
    try {
        ligint::selectKlifsStructure(unknownModel, o);
        CHECK(false);
    } catch (const ligint::KlifsAmbiguousError& e) {
        const std::string m = e.what();
        CHECK(e.candidates.size() == 2 && e.candidates[0].model.empty());
        CHECK(m.find("--klifs-id") != std::string::npos && m.find("or with --klifs-model") == std::string::npos && m.find("--klifs-model cannot choose") != std::string::npos);
        CHECK(m.find("fields of the answer are") != std::string::npos && m.find("rmsd1") != std::string::npos && m.find("structure_ID") != std::string::npos);
    }
    // no match, empty answers and answers that are not understood
    o.chain = "Z";
    CHECK(message(twoModels, o).find("no structure for chain Z") != std::string::npos && message(twoModels, o).find("ID=1453") != std::string::npos);
    CHECK(message("[]", o).find("(none)") != std::string::npos);
    o.chain = "A";
    CHECK(message("<html>not json</html>", o).find("starts with: '<html>not json</html>'") != std::string::npos);
    CHECK(message(R"([{"chain":"A"}])", o).find("no numeric structure ID") != std::string::npos);
    CHECK(message(R"([{"structure_ID":"abc","chain":"A"}])", o).find("no numeric structure ID") != std::string::npos);
    CHECK(message(R"([{"structure_ID":-4,"chain":"A"}])", o).find("no numeric structure ID") != std::string::npos);
    CHECK(message(R"([1, "x", null])", o).find("no structure") != std::string::npos);                          // elements that are not objects are skipped
    CHECK(message(std::string(5000, '{'), o).find("...") != std::string::npos);                               // a long answer is cut in the message
}

static void testRcsbCacheAndValidation() {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "ligint_rcsb_unit";
    fs::remove_all(dir);
    fs::create_directories(dir);
    {
        std::ofstream f(dir / "1ABC.pdb");
        f << "HEADER    TEST ENTRY\n" << std::string(150, ' ') << "\nATOM      1  N   ALA A   1       0.000   0.000   0.000\n";
    }
    ligint::RcsbDownloadOptions o;
    o.pdbId = "1abc";                                              // lower case is accepted and the cached file is found
    o.cacheDir = dir.string();
    CHECK(fs::path(ligint::downloadRcsbStructure(o)).filename() == "1ABC.pdb");
    // identifiers that are not four letters or digits never get to the file system or to the network
    for (const char* bad : {"", "1AB", "1ABCD", "../x", "1A/C", "1A C", "1AB\xC3\xA9", "1AB\n"}) {
        ligint::RcsbDownloadOptions b = o;
        b.pdbId = bad;
        bool thrown = false;
        try { ligint::downloadRcsbStructure(b); } catch (const std::runtime_error&) { thrown = true; }
        CHECK(thrown);
    }
    ligint::RcsbDownloadOptions f = o;
    f.format = "xyz";
    bool thrown = false;
    try { ligint::downloadRcsbStructure(f); } catch (const std::runtime_error&) { thrown = true; }
    CHECK(thrown);
    CHECK(std::distance(fs::directory_iterator(dir), fs::directory_iterator()) == 1);     // nothing was created
    fs::remove_all(dir);
    (void)ligint::rcsbOnlineAvailable();
    (void)ligint::klifsOnlineAvailable();
}

static void testKlifsNative() {
    const std::string filename="ligint_klifs_unit_mapping.tsv";
    {std::ofstream f(filename);f<<"residue.klifs_id\tresidue.id\tresidue.klifs_region\tresidue.aa\n"
       <<"1\t879\tI\tR\n47\t957\thinge\tE\n48\t_\thinge\tF\n81\t1021\tDFG\tD\n";}
    ligint::KlifsOptions o;o.pdbId="3EYG";o.chain="A";o.mappingFile=filename;
    const auto p=ligint::generateKlifsPocket(o);
    CHECK(p.positions.size()==3);
    CHECK(p.positions[1].resi==957);
    CHECK(p.positions[1].label=="47");
    CHECK(p.positions[1].region=="hinge");
    CHECK(p.positions[1].resName=="GLU");
    const std::string out="ligint_klifs_unit_pocket.tsv";
    ligint::saveKlifsPocket(p,out);
    ligint::PocketDefinition read;std::vector<std::string> warnings;
    CHECK(ligint::loadPocketDefinition(out,read,warnings));
    CHECK(read.positions.size()==3);
    CHECK(read.positions[1].region=="hinge");
    {std::ofstream f(filename);f<<"residue.klifs_id\tresidue.id\tresidue.klifs_region\n1\t900\tI\n1\t901\tI\n";}
    bool rejected=false;try{(void)ligint::generateKlifsPocket(o);}catch(const std::exception&){rejected=true;}
    CHECK(rejected);
    std::remove(filename.c_str());std::remove(out.c_str());
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
    testMmcifBasic();
    testMmcifNamesAndChains();
    testMmcifSyntax();
    testMmcifModelsAndAltlocs();
    testMmcifMalformed();
    testAlternateConformationSelection();
    testFingerprintClasses();
    testKlifsNative();
    testRcsbCacheAndValidation();
    testKlifsStructureSelection();
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
