// Interaction matrix: the fingerprints of several ligands, one row each, on the same columns.
//
// A column is a pocket position and an interaction class; a cell is 1 if the residue at that position makes at least one
// favorable interaction of that class with the ligand, 0 otherwise. With the 85 positions of a KLIFS pocket definition
// and the 7 classes of the KLIFS interaction fingerprint (IFP) this is a 595-bit fingerprint per ligand, comparable
// across structures because the columns are the pocket positions, not residue numbers.
//
// The ligands come from a list file (one line each: structure, ligand, optional label, pocket definition...), see
// docs/MATRIX.md. All the entries use a pocket definition, or none does: without one the columns are the residues in contact.
#pragma once

#include <climits>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

#include "ligint/fingerprint.hpp"
#include "ligint/ligint.hpp"

namespace ligint {

// ---------------------------------------------------------------------------------------------
// The list of ligands
// ---------------------------------------------------------------------------------------------
struct MatrixEntry {
    std::string structure;            // path of the PDB/mmCIF file
    std::string ligand;               // residue name of the ligand (required)
    std::string chain;                // optional: chain of the ligand
    int resi = INT_MIN;               // optional: residue number (INT_MIN = any)
    char icode = 0;                   // optional: insertion code (0 = any, ' ' = none)
    std::string label;                // name of the row; default: "<file name without extension> <ligand id>"
    std::string pocket;               // pocket definition of this structure; empty = the common one (MatrixOptions::pocket)
    char altloc = 0;                  // alternate conformation; 0 = the common one (Options::altloc)
    int line = 0;                     // line of the list file, for the messages
};

// Reads a list: a TSV or CSV file with a header line. Columns (in any order, other columns are ignored):
//   structure (required; also: file, path), ligand (resname), chain, resseq (resi), label (name), pocket, altloc
// Blank lines and lines starting with # are skipped. Relative paths are relative to the folder of the list.
// Throws std::runtime_error with the line number if the list is not valid.
std::vector<MatrixEntry> readMatrixList(const std::string& path);
std::vector<MatrixEntry> parseMatrixList(std::istream& in, const std::string& baseDir = std::string());

// Residue number of a ligand as the command line writes it: "501" (any insertion code), "501A" (code A), "501-" (no code).
bool parseLigandResidueNumber(const std::string& text, int& resi, char& icode);

// ---------------------------------------------------------------------------------------------
// The matrix
// ---------------------------------------------------------------------------------------------
struct MatrixOptions {
    Options analysis;                 // thresholds and the other options of the analysis; the ligand comes from each entry
    std::string pocket;               // pocket definition for the entries that do not have their own
    std::vector<int> classes;         // FpClass values of the columns; empty = defaultMatrixClasses()
    bool hideEmptyPositions = false;  // leave out the positions where no ligand has any interaction
};

// The 7 classes of the KLIFS IFP: hydrophobic, face-to-face, edge-to-face, H-bond donor and acceptor of the protein, ionic +/-.
std::vector<int> defaultMatrixClasses();
// "klifs" (the 7 above), "all" (every class) or a comma-separated list of class names (fpClassName). The result is in the
// order of the classes. Returns false with a message in `error` if a name is not known or the list is empty.
bool parseMatrixClasses(const std::string& text, std::vector<int>& out, std::string& error);

struct MatrixPosition {
    std::string label;                // pocket position label, or the residue ("ASP 189A") without a pocket definition
    std::string region;               // region of the pocket; empty without a pocket definition
};

struct MatrixRow {
    std::string label, structure, ligandId;
    // For each position, the residue as the pocket definition of this entry names it (or as it is in contact with the ligand, if it is).
    // "" = the pocket definition of this entry has no such position (with a KLIFS definition: a residue that is missing in the crystal).
    // Without a pocket definition: the residue in contact, "" for the others.
    std::vector<std::string> residues;
    std::vector<unsigned char> bits;         // positions.size() * classes.size() cells: index = position * classes.size() + class
    int bitCount() const;                    // number of cells set
};

struct InteractionMatrix {
    bool pocketMode = false;                 // the columns are pocket positions (false: residues in contact)
    std::vector<int> classes;
    std::vector<MatrixPosition> positions;
    std::vector<std::pair<std::string, std::string>> regionColors;    // "@color" of the pocket definitions
    std::vector<MatrixRow> rows;
    std::vector<std::string> warnings;
    bool bit(std::size_t row, std::size_t position, std::size_t cls) const {
        return rows[row].bits[position * classes.size() + cls] != 0;
    }
};

// Analyzes every entry and builds the matrix. Throws std::runtime_error naming the entry (line and label) if a structure cannot be
// read, the ligand is not found or ambiguous, or a pocket definition cannot be opened; or if the entries mix pocket definitions with none.
InteractionMatrix makeMatrix(const std::vector<MatrixEntry>& entries, const MatrixOptions& options = MatrixOptions());

// ---------------------------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------------------------
// CSV: label,structure,ligand, then one 0/1 column "<position>:<class>" for each position and class, row by row.
void writeMatrixCsv(std::ostream& out, const InteractionMatrix& m);
// JSON: {"classes": [...], "pocket": bool, "positions": [{"position", "region"}], "rows": [{"label", "structure", "ligand",
// "residues": [...], "bits": [...]}], "warnings": [...]}; bits are in the order of the CSV columns.
void writeMatrixJson(std::ostream& out, const InteractionMatrix& m);

struct MatrixSvgOptions {
    std::string title;
    bool showLegend = true;
};
// Heat map: a row for each ligand, a block of one stripe per class for each position, the regions on top. Self-contained.
std::string matrixSvg(const InteractionMatrix& m, const MatrixSvgOptions& opt = MatrixSvgOptions());

}  // namespace ligint
