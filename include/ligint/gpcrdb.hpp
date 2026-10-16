// Pocket definitions for GPCRs from the generic residue numbers of GPCRdb.
//
// A pocket definition (see fingerprint.hpp) names positions and regions. For a GPCR the positions are generic residue numbers,
// which are the same in every receptor of a class: "3x32" (GPCRdb scheme) or "3.32" (Ballesteros-Weinstein, class A), the helix and
// the position in it. The regions are the segments: TM1-TM7, H8 and the loops (ICL1-3, ECL1-3).
//
// The numbers of a structure come from GPCRdb itself: its tool "Generic residue numbering (PDB)" annotates an uploaded PDB file with
// the numbers in the B-factor column (GPCRdb scheme in the CA atoms, Ballesteros-Weinstein in the N atoms), or from a table
// (residue, generic number). The pocket of interest is a list of generic numbers that the caller gives: this is not a fixed pocket
// like the 85 positions of KLIFS, it depends on the receptor class and on the ligand site.
#pragma once

#include <string>
#include <vector>

#include "ligint/fingerprint.hpp"

namespace ligint {

struct GpcrdbOptions {
    // A PDB or mmCIF file annotated by GPCRdb (numbers in the B-factor column). Used if mappingFile is empty.
    std::string numberedStructure;
    // Alternatively, a table with a header line (TSV or CSV): residue number (resi, resseq...), GPCRdb number (gpcrdb, generic...)
    // and/or Ballesteros-Weinstein number (bw), optionally chain, segment and aa. See docs/GPCRDB.md.
    std::string mappingFile;
    std::string chain;                        // the chain of the receptor; empty = the only chain that has generic numbers
    std::vector<std::string> positions;       // the pocket: generic numbers, "3x32" (GPCRdb) or "3.32" (Ballesteros-Weinstein)
    std::string name;                         // name of the definition; empty = "GPCR pocket (GPCRdb numbering, chain X)"
};

struct GpcrdbPocket {
    PocketDefinition pocket;                  // the positions that this structure has, in the order of the helices
    std::vector<std::string> missing;         // requested positions that this structure does not have
    std::vector<std::string> warnings;
};

// Reads a generic residue number: "3x32", "3X32", "3.32", a bulge "2x551". `canonical` is the number as it is written here ("3x32"),
// `scheme` is 'x' (GPCRdb) or '.' (Ballesteros-Weinstein). False if the text is not a generic number.
bool parseGenericNumber(const std::string& text, std::string& canonical, char& scheme);

// The segment of a generic number: "TM1".."TM7", "H8", "ICL1".."ICL3", "ECL1".."ECL3" ("" if the number does not belong to one).
std::string genericSegment(const std::string& canonical);

// A list of generic numbers separated by commas, blanks or new lines ('#' starts a comment); if the text is the name of a file,
// the content of the file. Returns canonical numbers without repetitions. Throws std::runtime_error for a text that is not a number.
std::vector<std::string> parseGpcrdbPositions(const std::string& listOrFile);

// Builds the pocket definition of a structure. Throws std::runtime_error if the structure has no generic numbers, the chain is not
// clear or none of the positions is there.
GpcrdbPocket generateGpcrdbPocket(const GpcrdbOptions& options);

// Writes a pocket definition in the format that loadPocketDefinition reads, with a color for each region it uses.
void saveGpcrdbPocket(const PocketDefinition& pocket, const std::string& filename);

}  // namespace ligint
