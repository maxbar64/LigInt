// Library tests (no external dependencies): API, thresholds, reports, templates, composite ligands.
#include <cmath>
#include <map>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ligint/fingerprint.hpp"
#include "ligint/gpcrdb.hpp"
#include "ligint/klifs.hpp"
#include "ligint/matrix.hpp"
#include "ligint/rcsb.hpp"
#include "ligint/scene3d.hpp"
#include "ligint/pocket.hpp"
#include <algorithm>
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

// a pocket.mol2 as KLIFS writes it (header, atoms, substructures, MOE comment); residue i is called names[i] and numbered numbers[i]
static std::string pocketMol2(const std::vector<std::pair<std::string, std::string>>& residues, const char* eol = "\n") {
    std::string t;
    auto line = [&](const std::string& s) { t += s + eol; };
    line("@<TRIPOS>MOLECULE");
    line("HUMAN/TEST_1abc_altB_chainA");
    line(std::to_string(residues.size() * 2) + " 0 " + std::to_string(residues.size()) + " 0 0 ");
    line("BIOPOLYMER");
    line("USER_CHARGES");
    line("");
    line("");
    line("@<TRIPOS>ATOM");
    for (size_t i = 0; i < residues.size(); ++i)
        line("   " + std::to_string(2 * i + 1) + " N        6.2495    17.2789    53.5862 N.3    " + std::to_string(i + 1) + " " + residues[i].first + residues[i].second + "  -0.6701 BACKBONE");
    line("@<TRIPOS>BOND");
    line("@<TRIPOS>SUBSTRUCTURE");
    for (size_t i = 0; i < residues.size(); ++i)
        line("  " + std::to_string(i + 1) + " " + residues[i].first + residues[i].second + "    " + std::to_string(2 * i + 1) + " RESIDUE 4 A " + residues[i].first + " 1");
    line("# MOE 2012.10 (io_trps.svl 2012.10)");
    return t;
}

static void writeTextFile(const std::string& path, const std::string& text) {
    std::ofstream f(path, std::ios::binary);
    f << text;
}

// A residue of a receptor as GPCRdb annotates it: CA with the GPCRdb number (or 0), N with the Ballesteros-Weinstein number (or 1.00)
static std::string annotatedAtoms(const std::string& chain, int resi, const std::string& res, double ca, double n, int& serial) {
    auto line = [&](const char* name, double b, double x) {
        char buf[100];
        std::snprintf(buf, sizeof buf, "%-6s%5d %-4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n", "ATOM", serial++, name, res.c_str(), chain[0], resi, x, 0.0, 0.0, 1.0, b,
                      name[1] == 'N' ? " N" : " C");
        return std::string(buf);
    };
    return line(" N  ", n, resi * 1.5) + line(" CA ", ca, resi * 1.5 + 1.0);
}

// 30 residues of TM3 (numbers 3.23-3.52, resi 100-129) with a bulge after 3.35 (the next residue repeats the number), then 5 residues without number,
// ICL2, ECL2 and H8 residues, and a residue with the bulge written as GPCRdb writes it (negative B-factor)
static std::string annotatedReceptor(bool withBw, const std::string& chain = "A", bool numbers = true) {
    std::string t;
    int serial = 1;
    auto add = [&](int resi, const char* res, double ca, double bw) { t += annotatedAtoms(chain, resi, res, numbers ? ca : 20.0 + resi * 0.37, withBw && numbers ? bw : 1.0, serial); };
    int resi = 100;
    for (int k = 23; k <= 52; ++k) {
        const double v = 3.0 + k / 100.0;
        add(resi++, k == 32 ? "ASP" : "ALA", v, v);
        if (k == 35) add(resi++, "SER", v, 3.36);                      // a bulge: the number of the preceding residue again (a file that went through another program)
    }
    for (int i = 0; i < 5; ++i) add(resi++, "GLY", 0.0, 1.0);          // not numbered
    add(resi++, "PRO", 34.49, 0);                                         // loops and helix 8
    add(resi++, "PRO", 34.50, 0);
    add(resi++, "CYS", 45.52, 0);
    add(resi++, "TYR", 8.47, 0);
    add(resi++, "LEU", 2.55, 2.55);
    add(resi++, "ALA", -2.55, 2.56);                                     // a bulge in the raw output: negative, the number of the preceding residue
    t += "END\n";
    return t;
}

static void testGpcrdbNumbers() {
    namespace fs = std::filesystem;
    std::string c;
    char sch = 0;
    CHECK(ligint::parseGenericNumber("3x32", c, sch) && c == "3x32" && sch == 'x');
    CHECK(ligint::parseGenericNumber(" 3X32 ", c, sch) && c == "3x32");
    CHECK(ligint::parseGenericNumber("03x32", c, sch) && c == "3x32");
    CHECK(ligint::parseGenericNumber("3.32", c, sch) && c == "3.32" && sch == '.');
    CHECK(ligint::parseGenericNumber("2x551", c, sch) && c == "2x551");                               // a bulge
    CHECK(ligint::parseGenericNumber("45x50", c, sch) && c == "45x50" && ligint::parseGenericNumber("8x47", c, sch) && ligint::parseGenericNumber("12x50", c, sch));
    for (const char* bad : {"", "3", "x32", "3x", "3x3", "3x32345", "3y32", "9x50", "13x50", "3.321", "45.50", "ab", "3x3a", "3 x32", "-3x32", "100x32"})
        CHECK(!ligint::parseGenericNumber(bad, c, sch));
    CHECK(ligint::genericSegment("3x32") == "TM3" && ligint::genericSegment("8x47") == "H8" && ligint::genericSegment("12x50") == "ICL1" && ligint::genericSegment("23x50") == "ECL1" &&
          ligint::genericSegment("34x50") == "ICL2" && ligint::genericSegment("45x50") == "ECL2" && ligint::genericSegment("56x50") == "ICL3" && ligint::genericSegment("67x50") == "ECL3" &&
          ligint::genericSegment("nonsense").empty());
    // lists of positions
    auto pos = ligint::parseGpcrdbPositions("3x32, 3x33;6x48\n# a comment\n 3X32 7.39 # again\n");
    CHECK((pos == std::vector<std::string>{"3x32", "3x33", "6x48", "7.39"}));
    writeTextFile("gp_list.txt", "# the pocket\n3x32 3x33\n\n6x48\n");
    CHECK((ligint::parseGpcrdbPositions("gp_list.txt") == std::vector<std::string>{"3x32", "3x33", "6x48"}));
    auto listFailure = [](const std::string& text) {
        try { ligint::parseGpcrdbPositions(text); } catch (const std::runtime_error& e) { return std::string(e.what()); }
        return std::string("(no exception)");
    };
    CHECK(listFailure("3x32,xyz").find("Invalid generic residue number 'xyz'") != std::string::npos && listFailure("").find("No generic residue numbers") != std::string::npos &&
          listFailure("# only a comment").find("No generic residue numbers") != std::string::npos);
    std::remove("gp_list.txt");

    // a structure annotated by GPCRdb
    writeTextFile("gp_a.pdb", annotatedReceptor(false));
    auto failure = [](const ligint::GpcrdbOptions& o) {
        try { ligint::generateGpcrdbPocket(o); } catch (const std::runtime_error& e) { return std::string(e.what()); }
        return std::string("(no exception)");
    };
    ligint::GpcrdbOptions o;
    o.numberedStructure = "gp_a.pdb";
    o.positions = {"3x40", "3x32", "3x351", "3x35", "3x36", "34x50", "45x52", "8x47", "2x551", "3x99"};
    auto g = ligint::generateGpcrdbPocket(o);
    std::vector<std::string> labels, regions;
    for (const auto& p : g.pocket.positions) { labels.push_back(p.label); regions.push_back(p.region); }
    // in the order of the receptor: TM2, TM3, ICL2, ECL2, H8; a bulge number comes after the number that it follows
    CHECK((labels == std::vector<std::string>{"2x551", "3x32", "3x35", "3x351", "3x36", "3x40", "34x50", "45x52", "8x47"}));
    CHECK((regions == std::vector<std::string>{"TM2", "TM3", "TM3", "TM3", "TM3", "TM3", "ICL2", "ECL2", "H8"}));
    CHECK(g.pocket.positions[1].resName == "ASP" && g.pocket.positions[1].resi == 109 && g.pocket.positions[1].chain == "A");                  // 3.32 is the 10th numbered residue
    CHECK(g.pocket.positions[3].resName == "SER" && g.pocket.positions[3].resi == g.pocket.positions[2].resi + 1);                             // the repeated number is the bulge: 3x351
    CHECK(g.pocket.positions[0].resName == "ALA" && g.pocket.positions[0].resi == 141);                                  // 2x551 is the residue after Leu 140 (2x55); 100-130 are TM3, 131-135 are not numbered, 136-139 loops and H8
    CHECK((g.missing == std::vector<std::string>{"3x99"}) && g.warnings.size() == 1 && g.warnings[0].find("3x99") != std::string::npos);
    CHECK(g.pocket.colors.size() == 5 && g.pocket.name.find("chain A") != std::string::npos);                    // TM2, TM3, ICL2, ECL2, H8
    // the same with the negative B-factor: the raw output of GPCRdb
    CHECK(g.pocket.positions[0].label == "2x551" && g.pocket.positions[0].chain == "A");
    // the pocket file can be read back
    ligint::saveGpcrdbPocket(g.pocket, "gp_pocket.tsv");
    ligint::PocketDefinition back;
    std::vector<std::string> warnings;
    CHECK(ligint::loadPocketDefinition("gp_pocket.tsv", back, warnings) && warnings.empty() && back.positions.size() == g.pocket.positions.size() && back.name == g.pocket.name &&
          back.positions[3].label == "3x351" && back.positions[3].region == "TM3" && back.positions[6].region == "ICL2" && back.colors.size() == g.pocket.colors.size());
    // Ballesteros-Weinstein numbers need the N atoms
    o.positions = {"3.32", "3.40"};
    CHECK(failure(o).find("Ballesteros-Weinstein numbers (with a point") != std::string::npos && failure(o).find("5.42 and 5x42") != std::string::npos);
    writeTextFile("gp_b.pdb", annotatedReceptor(true));
    o.numberedStructure = "gp_b.pdb";
    g = ligint::generateGpcrdbPocket(o);
    CHECK(g.pocket.positions.size() == 2 && g.pocket.positions[0].label == "3.32" && g.pocket.positions[0].resName == "ASP" && g.pocket.positions[1].label == "3.40");
    o.positions = {"3x32", "3.32"};                                                                    // both schemes in one list
    CHECK(ligint::generateGpcrdbPocket(o).pocket.positions.size() == 2);
    // not annotated: the B-factors of a structure
    writeTextFile("gp_plain.pdb", annotatedReceptor(false, "A", false));
    o.numberedStructure = "gp_plain.pdb";
    o.positions = {"3x32"};
    CHECK(failure(o).find("does not have GPCRdb generic numbers") != std::string::npos && failure(o).find("generic_numbering_index") != std::string::npos);
    // chains
    writeTextFile("gp_two.pdb", annotatedReceptor(false, "A") + annotatedReceptor(false, "B"));
    o.numberedStructure = "gp_two.pdb";
    CHECK(failure(o).find("Several chains have generic numbers (A, B)") != std::string::npos);
    o.chain = "B";
    CHECK(ligint::generateGpcrdbPocket(o).pocket.positions[0].chain == "B");
    o.chain = "Z";
    CHECK(failure(o).find("Chain Z is not in the structure (chains: A: 37 of 42") != std::string::npos);
    o.chain.clear();
    // errors
    o.numberedStructure = "gp_a.pdb";
    o.positions = {"3x98", "3x99"};
    CHECK(failure(o).find("None of the requested positions (3x98, 3x99)") != std::string::npos);
    o.positions.clear();
    CHECK(failure(o).find("No generic residue numbers") != std::string::npos);
    o.positions = {"junk"};
    CHECK(failure(o).find("Invalid generic residue number 'junk'") != std::string::npos);
    o.positions = {"3x32"};
    o.numberedStructure = "no_such_file.pdb";
    CHECK(failure(o).find("Cannot open") != std::string::npos);
    o.numberedStructure.clear();
    CHECK(failure(o).find("annotated with generic numbers, or a table") != std::string::npos);

    // a table: the residue number, the numbers (GPCRdb, Ballesteros-Weinstein or both) and optionally chain, segment and amino acid
    writeTextFile("gp_table.tsv", "chain\tresi\tgpcrdb\tbw\tsegment\taa\nA\t113\t3x32\t3.32\tTM3\tD\nA\t114\t3x33\t3.33\tTM3\tV\nA\t203\t5x43\t5.42\tTM5\tS\nA\t204\t-\t-\t\tS\n"
                                  "A\t312\t7.39x38\t\tTM7\tN\nA\t316A\t7x42\t\t\tY\nB\t113\t3x32\t\tTM3\tD\nA\t118\tjunk\t\t\tT\n");
    ligint::GpcrdbOptions t;
    t.mappingFile = "gp_table.tsv";
    t.chain = "A";
    t.positions = {"3x32", "5.42", "7x38", "7.39", "7x42", "3x33"};
    g = ligint::generateGpcrdbPocket(t);
    std::vector<std::string> tl;
    for (const auto& p : g.pocket.positions) tl.push_back(p.label + "=" + p.chain + std::to_string(p.resi) + (p.icode == ' ' ? "" : std::string(1, p.icode)) + p.resName);
    CHECK((tl == std::vector<std::string>{"3x32=A113ASP", "3x33=A114VAL", "5.42=A203SER", "7x38=A312ASN", "7.39=A312ASN", "7x42=A316ATYR"}));   // in the order of the receptor; "7.39x38" gives both numbers
    CHECK(g.pocket.positions[2].region == "TM5" && g.pocket.positions[4].region == "TM7" && g.missing.empty());
    CHECK(g.warnings.size() == 1 && g.warnings[0].find("'junk' is not a generic residue number") != std::string::npos);
    t.chain = "B";
    t.positions = {"3x32", "3x33"};
    g = ligint::generateGpcrdbPocket(t);
    CHECK(g.pocket.positions.size() == 1 && g.pocket.positions[0].chain == "B" && (g.missing == std::vector<std::string>{"3x33"}));
    writeTextFile("gp_table_bad.tsv", "chain\tname\nA\tx\n");
    t.mappingFile = "gp_table_bad.tsv";
    CHECK(failure(t).find("needs a header with a residue number column") != std::string::npos);
    writeTextFile("gp_table_bad.tsv", "resi\tgpcrdb\n12AB\t3x32\n");
    CHECK(failure(t).find("Line 2: invalid residue number '12AB'") != std::string::npos);
    writeTextFile("gp_table_bad.tsv", "resi,gpcrdb\n");
    CHECK(failure(t).find("has no generic numbers") != std::string::npos);
    writeTextFile("gp_table_bad.tsv", "");
    CHECK(failure(t).find("empty") != std::string::npos);
    for (const char* f : {"gp_a.pdb", "gp_b.pdb", "gp_plain.pdb", "gp_two.pdb", "gp_pocket.tsv", "gp_table.tsv", "gp_table_bad.tsv"}) std::remove(f);
}

static void testMatrixList() {
    // header in any order, aliases, comments, quotes, a BOM and Windows line ends
    std::istringstream in("\xEF\xBB\xBF# a list\r\nlabel,Ligand,FILE,chain,resi,altloc,pocket\r\n"
                          "\"first, one\",STI,\"dir/a.pdb\",A,201,B,p.tsv\r\n\r\n"
                          "second,ATP,/abs/b.cif,,501A,,\r\nthird,ZN,c.pdb,B,12-,,\r\n");
    auto e = ligint::parseMatrixList(in, "base");
    CHECK(e.size() == 3);
    CHECK(e[0].label == "first, one" && e[0].ligand == "STI" && e[0].chain == "A" && e[0].resi == 201 && e[0].icode == 0 && e[0].altloc == 'B');
    CHECK(e[0].structure == (std::filesystem::path("base") / "dir/a.pdb").string() && e[0].pocket == (std::filesystem::path("base") / "p.tsv").string());
    CHECK(e[1].structure == "/abs/b.cif" && e[1].resi == 501 && e[1].icode == 'A' && e[1].pocket.empty() && e[1].altloc == 0 && e[1].line == 5);
    CHECK(e[2].resi == 12 && e[2].icode == ' ' && e[2].chain == "B");
    // tab-separated, with only the required column
    std::istringstream tsv("structure\nx.pdb\n");
    CHECK(ligint::parseMatrixList(tsv, "")[0].structure == "x.pdb");
    // errors name the line
    auto message = [](const std::string& text) {
        std::istringstream i(text);
        try { ligint::parseMatrixList(i, ""); } catch (const std::runtime_error& ex) { return std::string(ex.what()); }
        return std::string("(no exception)");
    };
    CHECK(message("").find("empty") != std::string::npos);
    CHECK(message("# only a comment\n").find("empty") != std::string::npos);
    CHECK(message("structure\tligand\n").find("no entries") != std::string::npos);
    std::istringstream alias("file\tligand\nx.pdb\tSTI\n");                                                       // 'file' is an alias of structure
    CHECK(ligint::parseMatrixList(alias, "")[0].ligand == "STI");
    CHECK(message("ligand\tchain\nSTI\tA\n").find("Line 1") != std::string::npos && message("ligand\tchain\nSTI\tA\n").find("'structure' column") != std::string::npos);
    CHECK(message("structure\tligand\n\tSTI\n").find("Line 2: the structure is empty") != std::string::npos);
    CHECK(message("structure\tresseq\na.pdb\t12AB\n").find("Line 2: invalid residue number '12AB'") != std::string::npos);
    CHECK(message("structure\tresseq\na.pdb\tx\n").find("invalid residue number") != std::string::npos);
    CHECK(message("structure\taltloc\na.pdb\tAB\n").find("altloc expects a single character") != std::string::npos);
    // residue numbers as the command line writes them
    int resi = 0;
    char icode = 9;
    CHECK(ligint::parseLigandResidueNumber("501", resi, icode) && resi == 501 && icode == 0);
    CHECK(ligint::parseLigandResidueNumber("501A", resi, icode) && icode == 'A');
    CHECK(ligint::parseLigandResidueNumber("501-", resi, icode) && icode == ' ');
    CHECK(ligint::parseLigandResidueNumber("-5", resi, icode) && resi == -5);
    CHECK(!ligint::parseLigandResidueNumber("", resi, icode) && !ligint::parseLigandResidueNumber("A", resi, icode) && !ligint::parseLigandResidueNumber("12AB", resi, icode) &&
          !ligint::parseLigandResidueNumber("1234567", resi, icode));
    // missing file
    bool thrown = false;
    try { ligint::readMatrixList("this_list_does_not_exist.tsv"); } catch (const std::runtime_error&) { thrown = true; }
    CHECK(thrown);
    // classes
    std::vector<int> c;
    std::string err;
    CHECK(ligint::parseMatrixClasses("klifs", c, err) && c == ligint::defaultMatrixClasses() && c.size() == 7);
    CHECK(ligint::parseMatrixClasses("ALL", c, err) && c.size() == static_cast<std::size_t>(ligint::FpClassCount));
    CHECK(ligint::parseMatrixClasses("hbond_donor, hydrophobic,hbond_donor", c, err) && c.size() == 2 && c[0] == ligint::FpHydrophobic && c[1] == ligint::FpHBondDonor);
    CHECK(!ligint::parseMatrixClasses("hydrophobic,nonsense", c, err) && err.find("nonsense") != std::string::npos);
    CHECK(!ligint::parseMatrixClasses(" , ", c, err) && !err.empty());
}

static void testMatrixBuild() {
    namespace fs = std::filesystem;
    writeTextFile("lm_a.pdb", methanolSystem(true));              // ligand LIG L1; ASN A10 accepts a hydrogen bond
    writeTextFile("lm_b.pdb", carboxylateLysSystem());           // ligand LIG L1; LYS A30: ionic positive and H-bond donor
    writeTextFile("lm_none.pdb", methanolSystem(false));         // a ligand and no protein
    writeTextFile("lm_two.pdb", methanolSystem(false).substr(0, methanolSystem(false).find("END")) +
                                    atomLine("HETATM", 8, "C1", "LIG", 'M', 2, 20, 0, 0, "C") + "\n" + atomLine("HETATM", 9, "O1", "LIG", 'M', 2, 21.43, 0, 0, "O") + "\nEND\n");
    writeTextFile("lm_pocket.tsv", "@name test\n@color R #ff0000\nA\t10\t1\thinge\tASN\nA\t30\t2\tR\tLYS\nA\t40\t3\tR\tGLY\n");
    writeTextFile("lm_pocket1.tsv", "A\t10\t1\thinge\tASN\nA\t30\t2\tR\tLYS\nA\t99\t4\tR\tGLY\n");
    writeTextFile("lm_pocket2.tsv", "A\t10\t1\thinge\tASN\nA\t30\t2\tR\tLYS\nA\t50\t3\tR\tALA\nA\t99\t4\tR\tGLY\n");
    auto entry = [](const std::string& file, const std::string& label, const std::string& pocket = std::string()) {
        ligint::MatrixEntry e;
        e.structure = file;
        e.ligand = "LIG";
        e.label = label;
        e.pocket = pocket;
        return e;
    };
    auto failure = [](const std::vector<ligint::MatrixEntry>& entries, const ligint::MatrixOptions& o = ligint::MatrixOptions()) {
        try { ligint::makeMatrix(entries, o); } catch (const std::runtime_error& ex) { return std::string(ex.what()); }
        return std::string("(no exception)");
    };
    // without a pocket definition the columns are the residues in contact
    auto m = ligint::makeMatrix({entry("lm_a.pdb", "a"), entry("lm_b.pdb", "b")});
    CHECK(!m.pocketMode && m.rows.size() == 2 && m.positions.size() == 2 && m.classes.size() == 7);
    CHECK(m.positions[0].label == "ASN 10A" && m.positions[1].label == "LYS 30A" && m.positions[0].region.empty());
    CHECK(m.bit(0, 0, 4) && m.rows[0].bitCount() == 1 && !m.bit(0, 1, 3));                                      // a: ASN accepts (class 4)
    CHECK(m.bit(1, 1, 3) && m.bit(1, 1, 5) && m.rows[1].bitCount() == 2 && m.rows[1].residues[0].empty() && m.rows[1].residues[1] == "LYS 30A");
    CHECK(m.rows[0].ligandId == "LIG 1L" && m.rows[0].label == "a" && m.warnings.empty());
    // the same with a pocket definition: every position is a column, also those without contacts
    ligint::MatrixOptions po;
    po.pocket = "lm_pocket.tsv";
    m = ligint::makeMatrix({entry("lm_a.pdb", "a"), entry("lm_b.pdb", "b")}, po);
    CHECK(m.pocketMode && m.positions.size() == 3 && m.positions[0].label == "1" && m.positions[0].region == "hinge" && m.positions[2].region == "R");
    CHECK(m.regionColors.size() == 1 && m.regionColors[0].first == "R");
    CHECK(m.bit(0, 0, 4) && m.bit(1, 1, 3) && m.bit(1, 1, 5) && m.rows[0].bitCount() == 1 && m.rows[1].bitCount() == 2);
    po.hideEmptyPositions = true;                                                                                // position 3 has no contact
    m = ligint::makeMatrix({entry("lm_a.pdb", "a"), entry("lm_b.pdb", "b")}, po);
    CHECK(m.positions.size() == 2 && m.positions[1].label == "2" && m.rows[1].bits.size() == 14 && m.rows[1].residues.size() == 2 && m.bit(1, 1, 5));
    // a selection of classes
    po.hideEmptyPositions = false;
    po.classes = {ligint::FpHBondDonor, ligint::FpIonicPositive};
    m = ligint::makeMatrix({entry("lm_b.pdb", "b")}, po);
    CHECK(m.classes.size() == 2 && m.rows[0].bits.size() == 6 && m.bit(0, 1, 0) && m.bit(0, 1, 1) && m.rows[0].bitCount() == 2);
    // pockets that differ per entry: a position that only one pocket has goes where it belongs
    m = ligint::makeMatrix({entry("lm_a.pdb", "a", "lm_pocket1.tsv"), entry("lm_b.pdb", "b", "lm_pocket2.tsv")});
    CHECK(m.positions.size() == 4 && m.positions[0].label == "1" && m.positions[1].label == "2" && m.positions[2].label == "3" && m.positions[3].label == "4");
    m = ligint::makeMatrix({entry("lm_b.pdb", "b", "lm_pocket2.tsv"), entry("lm_a.pdb", "a", "lm_pocket1.tsv")});                  // the same in the other order
    CHECK(m.positions.size() == 4 && m.positions[2].label == "3" && m.positions[3].label == "4");
    // a contact that is outside the pocket definition is reported and is not a column
    writeTextFile("lm_small.tsv", "A\t30\t2\tR\tLYS\n");
    po = ligint::MatrixOptions();
    po.pocket = "lm_small.tsv";
    m = ligint::makeMatrix({entry("lm_a.pdb", "a")}, po);
    CHECK(m.positions.size() == 1 && m.rows[0].bitCount() == 0 && m.warnings.size() == 1 && m.warnings[0].find("not in the pocket definition") != std::string::npos &&
          m.warnings[0].find("ASN 10A") != std::string::npos);
    // errors name the entry
    auto noLigand = entry("lm_a.pdb", "x");
    noLigand.ligand.clear();
    CHECK(failure({noLigand}).find("the ligand is required; the structure has: LIG(1)") != std::string::npos);
    auto missing = entry("lm_a.pdb", "x");
    missing.ligand = "XYZ";
    CHECK(failure({entry("lm_a.pdb", "a"), missing}).find("Entry 2") != std::string::npos && failure({missing}).find("ligand 'XYZ' not found") != std::string::npos);
    CHECK(failure({entry("lm_two.pdb", "two")}).find("matches 2 ligands") != std::string::npos);
    auto chosen = entry("lm_two.pdb", "two");
    chosen.chain = "M";
    CHECK(ligint::makeMatrix({chosen}).rows[0].ligandId == "LIG 2M");
    CHECK(failure({entry("no_such_file.pdb", "n")}).find("Entry 1") != std::string::npos);
    CHECK(failure({entry("lm_a.pdb", "same"), entry("lm_b.pdb", "same")}).find("used by another entry") != std::string::npos);
    CHECK(failure({entry("lm_a.pdb", "a", "lm_pocket.tsv"), entry("lm_b.pdb", "b")}).find("Entry 2") != std::string::npos &&
          failure({entry("lm_a.pdb", "a", "lm_pocket.tsv"), entry("lm_b.pdb", "b")}).find("has no pocket definition while others have one") != std::string::npos);
    CHECK(failure({entry("lm_a.pdb", "a", "no_such_pocket.tsv")}).find("cannot open the pocket definition") != std::string::npos);
    CHECK(failure({}).find("no entries") != std::string::npos);
    ligint::MatrixOptions bad;
    bad.classes = {99};
    CHECK(failure({entry("lm_a.pdb", "a")}, bad).find("Invalid interaction class") != std::string::npos);
    // a default label from the file and the ligand
    auto unnamed = entry("lm_a.pdb", "");
    CHECK(ligint::makeMatrix({unnamed}).rows[0].label == "lm_a LIG 1L");
    // the outputs
    po = ligint::MatrixOptions();
    po.pocket = "lm_pocket.tsv";
    m = ligint::makeMatrix({entry("lm_a.pdb", "a, with comma"), entry("lm_b.pdb", "b")}, po);
    std::ostringstream csv;
    ligint::writeMatrixCsv(csv, m);
    std::istringstream rows(csv.str());
    std::string header, row1, row2;
    std::getline(rows, header);
    std::getline(rows, row1);
    std::getline(rows, row2);
    CHECK(header.rfind("label,structure,ligand,1:hydrophobic,1:face_to_face", 0) == 0 && header.find("3:ionic_negative") != std::string::npos);
    CHECK(row1.rfind("\"a, with comma\",lm_a.pdb,LIG 1L,", 0) == 0 && row2.rfind("b,lm_b.pdb,LIG 1L,", 0) == 0);
    CHECK(std::count(header.begin(), header.end(), ',') == 3 + 21 - 1 && std::count(row2.begin(), row2.end(), '1') >= 2);
    std::ostringstream json;
    ligint::writeMatrixJson(json, m);
    CHECK(json.str().find("\"classes\": [\"hydrophobic\"") != std::string::npos && json.str().find("\"pocket\": true") != std::string::npos &&
          json.str().find("\"label\": \"a, with comma\"") != std::string::npos && json.str().find("\"region\": \"hinge\"") != std::string::npos);
    const std::string svg = ligint::matrixSvg(m);
    CHECK(svg.rfind("<svg ", 0) == 0 && svg.find("</svg>") != std::string::npos && svg.find("hinge") != std::string::npos && svg.find("Hydrophobic") != std::string::npos &&
          svg.find("a, with comma") != std::string::npos);
    ligint::MatrixSvgOptions so;
    so.title = "T & <x>";
    so.showLegend = false;
    const std::string svg2 = ligint::matrixSvg(m, so);
    CHECK(svg2.find("T &amp; &lt;x&gt;") != std::string::npos && svg2.find("Hydrophobic") == std::string::npos);
    for (const char* f : {"lm_a.pdb", "lm_b.pdb", "lm_none.pdb", "lm_two.pdb", "lm_pocket.tsv", "lm_pocket1.tsv", "lm_pocket2.tsv", "lm_small.tsv"}) std::remove(f);
}

static void testKlifsPocketMol2() {
    namespace fs = std::filesystem;
    const std::string aa3[] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};
    const std::string letters = "ACDEFGHIKLMNPQRSTVWY";
    auto generate = [&](const std::string& mol2, const std::string& sequence, std::string* error = nullptr) {
        const std::string file = "ligint_unit_pocket.mol2";
        { std::ofstream f(file, std::ios::binary); f << mol2; }
        ligint::KlifsOptions o;
        o.pdbId = "1ABC";
        o.chain = "A";
        o.mappingFile = file;
        o.pocketSequence = sequence;
        ligint::PocketDefinition p;
        try { p = ligint::generateKlifsPocket(o); } catch (const std::exception& e) { if (error) *error = e.what(); }
        std::remove(file.c_str());
        return p;
    };
    // all 85 positions: no sequence needed; the positions and their regions are those of KLIFS
    std::vector<std::pair<std::string, std::string>> all;
    std::string seq85;
    for (int p = 1; p <= 85; ++p) {
        all.push_back({aa3[p % 20], std::to_string(100 + p)});
        seq85 += letters[p % 20];
    }
    auto p = generate(pocketMol2(all), "");
    CHECK(p.positions.size() == 85);
    CHECK(p.positions[0].label == "1" && p.positions[0].region == "I" && p.positions[0].resi == 101 && p.positions[0].resName == aa3[1]);
    CHECK(p.positions[2].region == "I" && p.positions[3].region == "g.l" && p.positions[19].region == "\xCE\xB1" "C" && p.positions[44].region == "GK");
    CHECK(p.positions[45].region == "hinge" && p.positions[79].region == "xDFG" && p.positions[84].region == "a.l" && p.positions[84].resi == 185);
    CHECK(generate(pocketMol2(all, "\r\n"), seq85).positions.size() == 85);                      // Windows line ends, and the sequence agrees
    // 81 residues and a sequence with four missing positions: the residues go to the positions that the sequence has
    std::string seq = seq85;
    for (int missing : {6, 7, 8, 9}) seq[missing - 1] = '_';
    std::vector<std::pair<std::string, std::string>> rest;
    for (int q = 1; q <= 85; ++q) if (seq[q - 1] != '_') rest.push_back({aa3[q % 20], std::to_string(100 + q)});
    p = generate(pocketMol2(rest), seq);
    CHECK(p.positions.size() == 81 && p.positions[4].label == "5" && p.positions[5].label == "10" && p.positions[5].resi == 110 && p.positions[5].region == "II");
    {
        std::string e;
        generate(pocketMol2(rest), "", &e);                                                       // fewer than 85 and no sequence: not placeable
        CHECK(e.find("81 residues") != std::string::npos && e.find("pocket sequence") != std::string::npos);
        generate(pocketMol2(rest), seq85, &e);                                                    // the sequence has 85 residues, the file 81
        CHECK(e.find("has 85 residues but the pocket MOL2 has 81") != std::string::npos);
        std::string wrong = seq;
        wrong[10] = wrong[10] == 'A' ? 'C' : 'A';                                                 // a standard residue that disagrees
        generate(pocketMol2(rest), wrong, &e);
        CHECK(e.find("Position 11") != std::string::npos);
        generate(pocketMol2(rest), seq.substr(1), &e);
        CHECK(e.find("must have 85 characters") != std::string::npos);
        generate(pocketMol2(all), std::string(84, 'A') + "1", &e);
        CHECK(e.find("Invalid character") != std::string::npos);
        generate("@<TRIPOS>MOLECULE\nx\n1 0 0 0 0\n@<TRIPOS>ATOM\n", "", &e);                      // no substructures
        CHECK(e.find("no residues") != std::string::npos);
        generate("@<TRIPOS>MOLECULE\n@<TRIPOS>SUBSTRUCTURE\n 1 ALA1 1 RESIDUE\n", "", &e);
        CHECK(e.find("Invalid substructure line") != std::string::npos);
        generate("@<TRIPOS>MOLECULE\n@<TRIPOS>SUBSTRUCTURE\n 1 ALAxyz 1 RESIDUE 4 A ALA 1\n", "", &e);
        CHECK(e.find("Invalid residue") != std::string::npos);
        auto many = all;
        many.push_back({"ALA", "999"});
        generate(pocketMol2(many), "", &e);
        CHECK(e.find("at most 85") != std::string::npos);
    }
    // an insertion code, a negative residue number and a modified residue (its letter in the sequence need not match)
    auto special = rest;
    special[0] = {"CYS", "421A"};                                                              // position 1 is a C in the sequence
    special[1] = {"ASP", "-3"};
    special[2] = {"MSE", "7"};
    p = generate(pocketMol2(special), seq);
    CHECK(p.positions.size() == 81 && p.positions[0].resi == 421 && p.positions[0].icode == 'A' && p.positions[1].resi == -3 && p.positions[2].resName == "MSE");
    // a file that is not a MOL2 still goes through the CSV/TSV reader
    std::string e;
    generate("a\tb\n1\t2\n", "", &e);
    CHECK(e.find("no position") != std::string::npos);
    CHECK(fs::exists("ligint_unit_pocket.mol2") == false);
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


// ---- 3D scene: backbone, secondary structure, cartoon, surfaces ------------------------------------
#ifdef LIGINT_TEST_DATA_DIR
// fraction of the residues whose DSSP-like assignment agrees with the HELIX/SHEET records of the file
static double ssAgreement(const std::string& file, const std::string& lig, int* nHelix = nullptr, int* nStrand = nullptr) {
    ligint::Options o;
    o.ligName = lig;
    o.protein3d = true;
    const ligint::Result r = ligint::analyzeFile(std::string(LIGINT_TEST_DATA_DIR) + "/" + file, o);
    std::map<std::pair<std::string, int>, char> ref;
    std::ifstream in(std::string(LIGINT_TEST_DATA_DIR) + "/" + file);
    std::string l;
    auto range = [&](char kind, const std::string& chain, int a, int b) { for (int i = a; i <= b; i++) ref[{chain, i}] = kind; };
    while (std::getline(in, l)) {
        if (l.compare(0, 5, "HELIX") == 0 && l.size() > 37) range('H', std::string(1, l[19]), std::stoi(l.substr(21, 4)), std::stoi(l.substr(33, 4)));
        if (l.compare(0, 5, "SHEET") == 0 && l.size() > 37) range('E', std::string(1, l[21]), std::stoi(l.substr(22, 4)), std::stoi(l.substr(33, 4)));
    }
    int agree = 0, h = 0, e = 0;
    for (const auto& b : r.backbone) {
        auto it = ref.find({b.chain, b.resi});
        agree += (it == ref.end() ? 'C' : it->second) == b.ss;
        h += b.ss == 'H';
        e += b.ss == 'E';
    }
    if (nHelix) *nHelix = h;
    if (nStrand) *nStrand = e;
    return r.backbone.empty() ? 0.0 : static_cast<double>(agree) / r.backbone.size();
}
#endif

static double radiusOf(const ligint::V3f& p) { return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); }

static void testScene3d() {
    // the protein description is returned only on request (the command line does not need it)
    {
        ligint::Options plain;
        plain.ligName = "LIG";
        const std::string pdb = methanolSystem();
        std::istringstream a(pdb), b(pdb);
        CHECK(ligint::analyzeStream(a, plain).backbone.empty());
        ligint::Options with = plain;
        with.protein3d = true;
        CHECK(!ligint::analyzeStream(b, with).proteinAtoms.empty());
    }
#ifdef LIGINT_TEST_DATA_DIR
    int h = 0, e = 0;
    CHECK(ssAgreement("1ATP.pdb", "ATP", &h, &e) > 0.90);       // measured 97 %: a kinase (helices and strands)
    CHECK(h > 100 && e > 30);
    CHECK(ssAgreement("1AZM.pdb", "AZM") > 0.90);               // carbonic anhydrase: mostly strands
    CHECK(ssAgreement("1HSG.pdb", "MK1", &h, &e) > 0.90);       // HIV protease: strands
    CHECK(e > 60 && h < 40);

    // cartoon of a real protein: finite, consistent mesh with unit normals and tags inside the chain
    {
        ligint::Options o;
        o.ligName = "ATP";
        o.protein3d = true;
        const ligint::Result r = ligint::analyzeFile(std::string(LIGINT_TEST_DATA_DIR) + "/1ATP.pdb", o);
        const ligint::Mesh m = ligint::buildCartoon(r.backbone);
        CHECK(m.tri.size() > 50000 && m.pos.size() == m.nrm.size() && m.pos.size() == m.tag.size());
        bool ok = true;
        for (size_t i = 0; i < m.pos.size(); i++) {
            for (float c : m.pos[i]) ok = ok && std::isfinite(c);
            ok = ok && std::fabs(radiusOf(m.nrm[i]) - 1.0) < 1e-3 && m.tag[i] >= 0 && m.tag[i] <= r.backbone.size() - 1;
        }
        for (const auto& t : m.tri) for (int k : t) ok = ok && k >= 0 && k < static_cast<int>(m.pos.size());
        CHECK(ok);
        ligint::CartoonOptions tube;
        tube.tubeOnly = true;
        CHECK(!ligint::buildCartoon(r.backbone, tube).empty());
        CHECK(ligint::buildCartoon(std::vector<ligint::BackboneResidue>()).empty());
    }

    // surface: converges when the grid is refined, and the pocket clip keeps a part of it
    {
        ligint::Options o;
        o.ligName = "BNZ";
        o.protein3d = true;
        const ligint::Result r = ligint::analyzeFile(std::string(LIGINT_TEST_DATA_DIR) + "/181L.pdb", o);
        std::vector<ligint::ProteinAtom> atoms;
        for (const auto& a : r.proteinAtoms) if (!a.hetero) atoms.push_back(a);
        ligint::SurfaceOptions coarse, fine;
        coarse.spacing = 0.5;
        fine.spacing = 0.3;
        const ligint::Mesh mc = ligint::buildSurface(atoms, coarse), mf = ligint::buildSurface(atoms, fine);
        const double ac = ligint::meshArea(mc), af = ligint::meshArea(mf);
        CHECK(af > 5000 && af < 9000);                           // a 164-residue protein: roughly 7000 A^2
        CHECK(std::fabs(ac - af) / af < 0.06);
        CHECK(mc.tag.size() == mc.pos.size());
        std::vector<ligint::V3> pocket{r.ligands.at(0).atoms.at(0).pos};
        const ligint::Mesh part = ligint::clipMeshNear(mc, pocket, 8.0);
        CHECK(!part.empty() && part.tri.size() < mc.tri.size() / 2);
    }
#endif
    // a single atom: the exact areas are known (SAS: sphere of radius r + probe; SES: sphere of radius r)
    {
        std::vector<ligint::ProteinAtom> one(1);
        one[0].radius = 1.7f;
        ligint::SurfaceOptions sas;
        sas.kind = ligint::SurfaceKind::SolventAccessible;
        sas.spacing = 0.4;
        const ligint::Mesh ms = ligint::buildSurface(one, sas);
        const double exactSas = 4 * M_PI * 3.1 * 3.1;
        CHECK(std::fabs(ligint::meshArea(ms) / exactSas - 1.0) < 0.03);
        int outward = 0;
        for (size_t i = 0; i < ms.pos.size(); i++) outward += ms.nrm[i][0] * ms.pos[i][0] + ms.nrm[i][1] * ms.pos[i][1] + ms.nrm[i][2] * ms.pos[i][2] > 0;
        CHECK(outward == static_cast<int>(ms.pos.size()));       // the normals point to the solvent
        ligint::SurfaceOptions ses;
        ses.spacing = 0.3;
        const double exactSes = 4 * M_PI * 1.7 * 1.7;
        CHECK(std::fabs(ligint::meshArea(ligint::buildSurface(one, ses)) / exactSes - 1.0) < 0.05);
    }
    // a closed shell of atoms encloses a cavity that the solvent cannot reach: only the outer surface is built
    {
        std::vector<ligint::ProteinAtom> shell;
        const int n = 500;
        for (int i = 0; i < n; i++) {                             // Fibonacci sphere of radius 8 A
            const double y = 1 - 2.0 * (i + 0.5) / n, rr = std::sqrt(1 - y * y), phi = i * 2.399963229728653;
            ligint::ProteinAtom a;
            a.pos = {8 * rr * std::cos(phi), 8 * y, 8 * rr * std::sin(phi)};
            a.radius = 1.7f;
            shell.push_back(a);
        }
        ligint::SurfaceOptions o;
        o.spacing = 0.5;
        const ligint::Mesh m = ligint::buildSurface(shell, o);
        double rmin = 1e9;
        for (const auto& p : m.pos) rmin = std::min(rmin, radiusOf(p));
        CHECK(!m.empty() && rmin > 9.0);                          // the cavity surface would be at about 6 A
    }
}


static void testBuriednessPockets() {
    // A shell surrounding an empty cavity should have high directional buriedness.
    std::vector<ligint::PocketAtom> shell;
    for(int i=0;i<150;i++) {
        double y=1.0-2.0*(i+0.5)/150.0;
        double r=std::sqrt(1-y*y), a=i*2.399963229728653;
        shell.push_back({{5*r*std::cos(a),5*y,5*r*std::sin(a)},1.7,"A","SHELL"});
    }
    ligint::PocketOptions o;
    o.minClusterPoints=1;
    o.minBuriedness=0.8;
    o.maxGridPoints=250000;
    auto found=ligint::findPockets(shell,o);
    CHECK(found.selectedAtoms==shell.size());
    CHECK(found.sampledPoints>0);
    CHECK(!found.pockets.empty());
    // Disabling segmentation must reproduce the original connected components.
    ligint::PocketOptions raw=o;
    raw.segment=false;
    const auto unsegmented=ligint::findPockets(shell,raw);
    CHECK(!unsegmented.pockets.empty());
    std::size_t rawCount=0, splitCount=0;
    for(const auto& p:unsegmented.pockets)rawCount+=p.pointCount;
    for(const auto& p:found.pockets)splitCount+=p.pointCount;
    CHECK(splitCount<=rawCount);
    CHECK(found.pockets.size()>=unsegmented.pockets.size());
    ligint::PocketOptions invalid=o;
    invalid.seedSpacing=0;
    bool rejectedSeed=false;
    try {(void)ligint::findPockets(shell,invalid);} catch(const std::invalid_argument&){rejectedSeed=true;}
    CHECK(rejectedSeed);

    // The reader must work without any ligand-selection workflow.
    const auto temp = std::filesystem::temp_directory_path() / "ligint_pocket_input_test.pdb";
    {
        std::ofstream f(temp);
        f << atomLine("ATOM", 1, "CA", "ALA", 'A', 1, 0, 0, 0, "C") << "\n";
        f << atomLine("ATOM", 2, "CA", "GLY", 'B', 2, 4, 0, 0, "C") << "\n";
        f << atomLine("HETATM", 3, "O", "HOH", 'A', 3, 8, 0, 0, "O") << "\nEND\n";
    }
    const auto loaded=ligint::readPocketAtoms(temp.string());
    CHECK(loaded.size()==2);
    CHECK(loaded.size()<2 || (loaded[0].chain=="A" && loaded[1].chain=="B"));
    std::filesystem::remove(temp);
    ligint::PocketOptions filtered=o;
    filtered.chains={"B"};
    CHECK(ligint::findPockets(shell,filtered).selectedAtoms==0);
    filtered=o;
    filtered.spacing=0;
    bool rejected=false;
    try { (void)ligint::findPockets(shell,filtered); } catch(const std::invalid_argument&) {rejected=true;}
    CHECK(rejected);
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
    testKlifsPocketMol2();
    testMatrixList();
    testGpcrdbNumbers();
    testMatrixBuild();
    testPocketDefinition();
    testFingerprintWithPocket();
    testFingerprintOutputs();
    testScene3d();
    testBuriednessPockets();
    if (failures) {
        std::cerr << failures << " checks failed\n";
        return 1;
    }
    std::cout << "ligint_core_tests: all tests passed\n";
    return 0;
}
