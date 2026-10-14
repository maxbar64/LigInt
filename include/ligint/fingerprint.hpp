// Interaction fingerprints: one row per protein residue, one column per interaction class.
//
// A fingerprint summarizes the interactions of one ligand (LigandResult) as a table that can be exported to
// CSV/JSON or drawn as an SVG figure (residues grouped by pocket region, with a glyph per interaction class).
//
// The pocket definition is optional. It is a plain text file that maps residues to a position label and a region
// (see docs/FINGERPRINT.md); without it the fingerprint lists only the residues in contact, sorted by number.
// The definition is data, not code: the same mechanism serves any protein family.
#pragma once

#include <array>
#include <iosfwd>
#include <string>
#include <vector>

#include "ligint/ligint.hpp"

namespace ligint {

// ---------------------------------------------------------------------------------------------
// Interaction classes (columns of the fingerprint)
// ---------------------------------------------------------------------------------------------
enum FpClass {
    FpHydrophobic = 0,
    FpFaceToFace,        // pi-stacking, parallel
    FpEdgeToFace,        // pi-stacking, T-shaped
    FpHBondDonor,        // the protein residue donates the hydrogen
    FpHBondAcceptor,     // the protein residue accepts the hydrogen
    FpIonicPositive,     // salt bridge with a positively charged residue
    FpIonicNegative,     // salt bridge with a negatively charged residue
    FpCationPi,
    FpHalogen,
    FpMetal,
    FpWaterBridge,
    FpClassCount
};

const char* fpClassName(int c);       // "hydrophobic", "face_to_face", "edge_to_face", "hbond_donor", ...
const char* fpClassTitle(int c);      // human-readable: "Hydrophobic", "Aromatic face-to-face", ...

// ---------------------------------------------------------------------------------------------
// Pocket definition
// ---------------------------------------------------------------------------------------------
struct PocketPosition {
    std::string chain;                // empty = any chain
    int resi = 0;
    char icode = ' ';
    std::string label;                // position label, e.g. "hinge.46" or "7"
    std::string region;               // region name, e.g. "hinge"
    std::string resName;              // optional residue name (3 letters), used to label empty positions
};

struct PocketDefinition {
    std::string name;
    std::vector<PocketPosition> positions;                       // in drawing order
    std::vector<std::pair<std::string, std::string>> colors;     // optional region colors ("@color REGION #rrggbb")
};

// Reads a pocket definition. Malformed lines are skipped and reported in `warnings`.
// parsePocketDefinition reads from a stream; loadPocketDefinition returns false if the file cannot be opened.
bool parsePocketDefinition(std::istream& in, PocketDefinition& out, std::vector<std::string>& warnings);
bool loadPocketDefinition(const std::string& path, PocketDefinition& out, std::vector<std::string>& warnings);

// ---------------------------------------------------------------------------------------------
// Fingerprint
// ---------------------------------------------------------------------------------------------
struct FingerprintEntry {
    std::string residue;              // label "ASP 189A"
    std::string resName, chain;
    int resi = 0;
    char icode = ' ';
    std::string position;             // from the pocket definition (empty without one)
    std::string region;               // idem
    bool inPocket = true;             // false: residue in contact that is not in the pocket definition
    std::array<int, FpClassCount> counts{};          // number of contacts per class
    std::array<double, FpClassCount> minDistance{};  // shortest distance per class (0 = none)
    bool any() const;
};

struct Fingerprint {
    std::string ligandId;
    std::string pocketName;
    std::vector<std::pair<std::string, std::string>> regionColors;   // from the pocket definition ("@color")
    std::vector<FingerprintEntry> entries;           // pocket order, or sorted by chain/number without a definition
    std::vector<std::string> warnings;               // e.g. residue name mismatches with the definition
};

// Builds the fingerprint of a ligand. Only favorable interactions count (clashes, repulsions and possible covalent
// bonds are ignored). With a pocket definition every position is listed, also those without interactions, and the
// residues in contact that are not in the definition are appended at the end (inPocket = false).
Fingerprint makeFingerprint(const LigandResult& lig, const PocketDefinition* pocket = nullptr);

// One-letter residue code ("ALA" -> 'A'; nucleotides DA/A -> 'A'; unknown -> '?').
char oneLetterCode(const std::string& resName);

// ---------------------------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------------------------
// CSV: ligand,residue,chain,resi,position,region,<one count column per class>. One header, one row per entry and ligand.
void writeFingerprintCsv(std::ostream& out, const std::vector<Fingerprint>& fps);
// JSON: an array with one object per ligand (classes, residues with counts and minimum distances).
void writeFingerprintJson(std::ostream& out, const std::vector<Fingerprint>& fps);

struct SvgOptions {
    int cellsPerRow = 20;             // residues per row of the figure
    bool showLegend = true;
    bool showPositionLabels = true;   // position label above each residue (needs a pocket definition)
    std::string title;                // default: the ligand id
};
// SVG figure of the fingerprint (self-contained, no external resources).
std::string fingerprintSvg(const Fingerprint& fp, const SvgOptions& opt = SvgOptions());

}  // namespace ligint
