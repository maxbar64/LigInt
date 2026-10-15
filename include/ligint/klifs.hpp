#pragma once
#include <stdexcept>
#include <string>
#include <vector>
#include "ligint/fingerprint.hpp"
namespace ligint {
struct KlifsOptions {
    std::string pdbId, chain, model;
    int structureId = 0;
    std::string cacheDirectory;
    std::string mappingFile; // offline KLIFS pocket: a CSV/TSV mapping, or the pocket.mol2 file that KLIFS provides for a structure
    // The 85-character sequence of the KLIFS pocket ('_' for a residue that is missing in the crystal): the field `pocket` of the service or of
    // the KLIFS overview. It is needed to place the residues of a pocket.mol2 that has fewer than 85 of them (a pocket with all of them needs nothing).
    std::string pocketSequence;
};
// A KLIFS structure that matches a request. The same PDB entry and chain can have several structures, one for each alternate model.
struct KlifsCandidate {
    int id = 0;
    std::string chain, model;      // model is empty when the answer of the service does not say it
    std::string pocket;            // the 85-character pocket sequence of this structure, if the answer has it
};
// Several KLIFS structures match: the candidates tell which ones, so that the caller can ask which to use (KlifsOptions::structureId).
class KlifsAmbiguousError : public std::runtime_error {
public:
    KlifsAmbiguousError(const std::string& what, std::vector<KlifsCandidate> c) : std::runtime_error(what), candidates(std::move(c)) {}
    std::vector<KlifsCandidate> candidates;
};
// Chooses the KLIFS structure ID in the answer of the service about a PDB entry (the JSON text of `structures_pdb_list`): the structure of
// options.chain (and of options.model if given). Returns options.structureId without reading anything if it is set. Throws
// KlifsAmbiguousError if several structures match and std::runtime_error if none does or the answer is not understood.
// In a build without network support it always throws std::runtime_error. Exposed for the tests, which cannot reach the service.
int selectKlifsStructure(const std::string& structuresJson, const KlifsOptions& options);

// Raises std::runtime_error on missing/ambiguous structures, bad responses or invalid mapping.
// Offline mappingFile expects KLIFS columns residue.klifs_id, residue.id, residue.klifs_region.
PocketDefinition generateKlifsPocket(const KlifsOptions& options);
void saveKlifsPocket(const PocketDefinition& pocket, const std::string& filename);
bool klifsOnlineAvailable();
}
