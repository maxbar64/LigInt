// ligint - non-covalent ligand-protein interactions from a PDB file
//
// Public API of the library. Typical use:
//
//     ligint::Options opt;            // thresholds, ligand to analyze, SDF/MOL2 templates...
//     opt.ligName = "ATP";
//     ligint::Result r = ligint::analyzeFile("1ATP.pdb", opt);
//     for (const auto& lig : r.ligands)
//         for (const auto& h : lig.hits) { /* h.type, h.ligandAtom, h.proteinAtom, h.distance ... */ }
//
// Without ligName only the list of ligands present is returned (Result::listingOnly).
// The library never prints: warnings and information are returned in Result.
#pragma once

#include <array>
#include <climits>
#include <iosfwd>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef LIGINT_VERSION
#define LIGINT_VERSION "0.10.1"
#endif

namespace ligint {

using V3 = std::array<double, 3>;

// ---------------------------------------------------------------------------------------------
// Interaction types (values of HitRecord::type)
// ---------------------------------------------------------------------------------------------
namespace type {
inline constexpr const char* HydrogenBond = "Hydrogen bond";
inline constexpr const char* SaltBridge = "Salt bridge";
inline constexpr const char* ElectrostaticRepulsion = "Electrostatic repulsion";
inline constexpr const char* Hydrophobic = "Hydrophobic contact";
inline constexpr const char* BackboneApolar = "Apolar contact (backbone)";
inline constexpr const char* PiStackingParallel = "pi-stacking (parallel)";
inline constexpr const char* PiStackingT = "pi-stacking (T-shaped)";
inline constexpr const char* CationPi = "Cation-pi";
inline constexpr const char* HalogenBond = "Halogen bond";
inline constexpr const char* MetalCoordination = "Metal coordination";
inline constexpr const char* WaterBridge = "Water bridge";
inline constexpr const char* StericClash = "Steric clash";
inline constexpr const char* PossibleCovalent = "Possible covalent bond";
}  // namespace type

// ---------------------------------------------------------------------------------------------
// Thresholds
// ---------------------------------------------------------------------------------------------
struct Params {
    static constexpr double kMaxGridRadius = 6.0;   // maximum supported search radius (Angstrom)

    double hbMin = 2.2, hbMax = 3.5, hbHA = 2.7, hbDHA = 120.0;          // hydrogen bonds
    double hydMin = 3.3, hydMax = 4.0, saltMax = 4.0, repMax = 4.5;       // hydrophobic, salt bridges, repulsions
    double halMax = 3.5, halAngle = 140.0, metalMax = 3.5;                // halogens, metals (global cap)
    double piMax = 5.5, piOffset = 3.5, piAngle = 30.0, cpMax = 6.0, clash = 0.5;
    double wbMin = 75.0, wbMax = 140.0;                                   // ligand-O(water)-protein angle
    double hbLong = 0.0, hbLongTol = 40.0;                                // long H-bonds (0 = off)
    double saltCentroid = 0.0;                                            // salt bridges by centroids (0 = off)
    double hbAnte = 90.0;                                                 // antecedent-atom...partner angle

    // Sets a threshold by name (as on the command line, without "--"); false if the name is unknown.
    bool set(const std::string& name, double value);
    // Preset for comparison with PLIP: hbLong 4.1, saltCentroid 5.5, halMax 4.0.
    void applyPlipLike();
    // Limits to kMaxGridRadius the thresholds used by the grid search; returns the warnings.
    std::vector<std::string> clampToGrid();
};

// ---------------------------------------------------------------------------------------------
// Chemical templates (SDF / MOL2)
// ---------------------------------------------------------------------------------------------
struct TAtom { std::string el, type; V3 p{0, 0, 0}; int fc = 0; };
struct TBond { int a, b, order; };          // 0-based indices; order 4 = aromatic
struct Tmpl { std::string name; std::vector<TAtom> at; std::vector<TBond> bd; };

// Read all molecules in the file. Return false if the file cannot be opened; problems with single
// molecules (V3000, malformed file) are appended to `warnings` and the molecule is skipped.
bool loadSdf(const std::string& path, std::vector<Tmpl>& out, std::vector<std::string>& warnings);
bool loadMol2(const std::string& path, std::vector<Tmpl>& out, std::vector<std::string>& warnings);

// ---------------------------------------------------------------------------------------------
// Input
// ---------------------------------------------------------------------------------------------
// File format of the structure. Auto recognizes it from the content: PDBx/mmCIF files start with "data_" (or a
// comment), PDB files with a record name. Compressed files (gzip) are not supported: decompress them first.
enum class StructureFormat { Auto, Pdb, Mmcif };

struct Options {
    StructureFormat format = StructureFormat::Auto;
    Params params;
    std::string ligName;            // HETATM residue name (e.g. "ATP"); empty = only list the ligands
    std::string chain;              // optional: ligand chain
    int resi = INT_MIN;             // optional: residue number (INT_MIN = any)
    char icode = 0;                 // optional: insertion code (0 = any, ' ' = none, otherwise that code)
    char altloc = 0;                // alternate conformation to use (e.g. 'B'); 0 = the first conformation of each residue.
                                    // A residue without that conformation keeps its first one.
    bool useWater = true;           // consider waters (water bridges)
    bool metalDonors = false;       // O/N coordinated to a metal remain possible donors
    bool backboneHyd = false;       // apolar contacts with backbone CA/C
    bool mergeLinked = true;        // merge covalently linked HETATM residues into one ligand
    bool protein3d = false;         // also return the backbone and the atoms of the protein for the 3D graphics (see scene3d.hpp)
    std::vector<Tmpl> templates;    // chemical annotation of the ligands (matched by coordinates)
};

// ---------------------------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------------------------
struct AtomRecord {
    std::string name, res, chain, element;
    char icode = ' ';
    int resi = 0;
    V3 pos{0, 0, 0};
    bool hetero = false;
};
struct BondRecord { int a, b, order; };     // indices into the atom vector; order 0 = unknown, 4 = aromatic

struct HitRecord {
    std::string type;               // one of ligint::type::*
    std::string ligandAtom;         // label "RES nnnC:ATOM" (or of the ring/group)
    std::string proteinAtom;
    std::string residue;            // protein residue, e.g. "ASP 189A"
    std::string note;
    double distance = 0;            // Angstrom
    bool favorable = true;          // false for clashes and repulsions
    V3 ligandPos{0, 0, 0};          // ligand point (atom, ring centroid or charged-group centroid)
    V3 proteinPos{0, 0, 0};
    V3 waterPos{0, 0, 0};           // water bridges only
    bool hasWater = false;
    int residueCharge = 0;          // charge of the protein atom in salt bridges, repulsions and charge-assisted H-bonds (+1/-1/0)
    char hbResidueRole = 0;         // hydrogen bonds: 'D' = the protein residue donates, 'A' = it accepts, 0 = not an H-bond
};

struct LigandResult {
    std::string id;                 // ligand label, e.g. "BZC 201A"
    std::string resName, chain;     // first residue
    int resi = 0;
    char icode = ' ';               // insertion code of the first residue (' ' = none)
    std::vector<std::string> residues;      // residues merged into the composite ligand (if more than one)
    std::string annotation;                 // how the chemistry was obtained (template, PDB hydrogens, heuristics)
    std::vector<std::string> info;          // informational lines (composite, annotation, B-factor)
    std::vector<AtomRecord> atoms;          // heavy atoms of the ligand
    std::vector<BondRecord> bonds;
    std::size_t nRings = 0;
    std::vector<HitRecord> hits;            // sorted by type, atom, distance
    std::vector<AtomRecord> contactAtoms;   // all atoms of the protein residues that interact
    std::vector<BondRecord> contactBonds;
    std::vector<V3> waters;                 // waters involved in water bridges
    std::string metalSphere;                // text: coordination sphere of the metals touched
};

struct LigandId {                           // a HETATM residue present in the file
    std::string res, chain;
    int resi = 0;
    char icode = ' ';
    int atoms = 0;
};

// Protein description for the 3D graphics (filled when Options::protein3d is set; geometry in scene3d.hpp).
struct BackboneResidue {
    std::string res, chain;
    int resi = 0;
    char icode = ' ';
    V3 n{0, 0, 0}, ca{0, 0, 0}, c{0, 0, 0}, o{0, 0, 0};
    bool hasN = false, hasC = false, hasO = false;   // the C-alpha is always present
    bool linkedToPrevious = false;                   // peptide bond with the previous residue of the vector
    char ss = 'C';                                   // secondary structure: 'H' helix, 'E' strand, 'C' coil
    int chainIndex = 0;                              // 0, 1, 2... in order of appearance
    double bfactor = 0;                              // of the C-alpha
};
// Heavy atom of the protein, for the surfaces (the radius is the van der Waals radius).
struct ProteinAtom {
    V3 pos{0, 0, 0};
    float radius = 1.7f;
    char element[3] = {'C', 0, 0};
    int residue = -1;                                // index in Result::backbone, -1 if the residue has no backbone
    bool hetero = false;                             // HETATM record (ions, cofactors, modified residues)
};

struct Result {
    bool listingOnly = false;               // no ligand requested or found: only the listing
    std::map<std::string, int> hetCounts;   // name -> number of residues (only if listingOnly)
    std::vector<LigandId> hetResidues;      // same, residue by residue
    StructureFormat format = StructureFormat::Pdb;   // format that was read (never Auto)
    std::string notFound;                   // requested name that is absent (only if listingOnly)
    std::vector<char> altlocs;              // alternate-conformation identifiers found in the first model (sorted)
    std::vector<std::string> warnings;      // warnings (invalid coordinates, limited thresholds...)
    std::vector<std::string> info;          // global information (protein hydrogens, templates loaded)
    std::vector<LigandResult> ligands;
    std::vector<BackboneResidue> backbone;  // only with Options::protein3d: protein residues in file order, with secondary structure
    std::vector<ProteinAtom> proteinAtoms;  // only with Options::protein3d: heavy atoms of the protein (without the analyzed ligand)
};

// Analyze a structure file (PDB or PDBx/mmCIF). Throws std::runtime_error if the file cannot be opened
// or is compressed.
Result analyzeFile(const std::string& path, const Options& opt);
Result analyzeStream(std::istream& in, const Options& opt);

// ---------------------------------------------------------------------------------------------
// Text reports (identical to the command-line output)
// ---------------------------------------------------------------------------------------------
void printText(std::ostream& out, const Result& r);                       // readable table
void printCsv(std::ostream& out, std::ostream& info, const Result& r);    // CSV on `out`, information on `info`

}  // namespace ligint
