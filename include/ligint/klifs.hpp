#pragma once
#include <string>
#include <vector>
#include "ligint/fingerprint.hpp"
namespace ligint {
struct KlifsOptions {
    std::string pdbId, chain, model;
    int structureId = 0;
    std::string cacheDirectory;
    std::string mappingFile; // offline KLIFS pocket CSV/TSV
};
// Raises std::runtime_error on missing/ambiguous structures, bad responses or invalid mapping.
// Offline mappingFile expects KLIFS columns residue.klifs_id, residue.id, residue.klifs_region.
PocketDefinition generateKlifsPocket(const KlifsOptions& options);
void saveKlifsPocket(const PocketDefinition& pocket, const std::string& filename);
bool klifsOnlineAvailable();
}
