#pragma once
#include <string>

namespace ligint {
struct RcsbDownloadOptions {
    std::string pdbId;
    std::string format = "pdb"; // pdb or cif
    std::string outputFile;      // optional; defaults to ID.format in cacheDir
    std::string cacheDir = ".";
    bool refresh = false;
};
// Returns absolute or caller-relative output path. Throws std::runtime_error on failure.
std::string downloadRcsbStructure(const RcsbDownloadOptions& options);
bool rcsbOnlineAvailable();
}
