#include "ligint/rcsb.hpp"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#ifdef LIGINT_RCSB_ONLINE
#include <curl/curl.h>
#endif
namespace ligint {
namespace {
namespace fs = std::filesystem;
bool validId(const std::string& id) {
    return id.size() == 4 && std::all_of(id.begin(), id.end(), [](unsigned char c){return std::isalnum(c) != 0 && c < 128;});
}
std::string upper(std::string s) {for(char& c:s)c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));return s;}
bool validContents(const fs::path& path, const std::string& format) {
    std::error_code ec;
    if(!fs::is_regular_file(path,ec) || fs::file_size(path,ec) < 100) return false;
    std::ifstream f(path, std::ios::binary);
    if(!f) return false;
    std::string line;
    for(int n=0;n<60 && std::getline(f,line);++n){
        if(format=="cif" && line.rfind("data_",0)==0) return true;
        if(format=="pdb" && (line.rfind("HEADER",0)==0 || line.rfind("ATOM  ",0)==0 || line.rfind("HETATM",0)==0 || line.rfind("TITLE ",0)==0)) return true;
    }
    return false;
}
#ifdef LIGINT_RCSB_ONLINE
size_t writeBytes(char* data, size_t size, size_t count, void* handle) {
    return std::fwrite(data,size,count,static_cast<FILE*>(handle));
}
#endif
}
bool rcsbOnlineAvailable() {
#ifdef LIGINT_RCSB_ONLINE
    return true;
#else
    return false;
#endif
}
std::string downloadRcsbStructure(const RcsbDownloadOptions& o) {
    if(!validId(o.pdbId)) throw std::runtime_error("RCSB requires a four-character alphanumeric PDB ID");
    if(o.format!="pdb" && o.format!="cif") throw std::runtime_error("RCSB format must be pdb or cif");
    const std::string id = upper(o.pdbId);
    const fs::path dest = o.outputFile.empty() ? fs::path(o.cacheDir)/(id+"."+o.format) : fs::path(o.outputFile);
    if(!o.refresh && validContents(dest,o.format)) return dest.string();
#ifndef LIGINT_RCSB_ONLINE
    throw std::runtime_error("RCSB HTTPS unavailable: rebuild with libcurl (CURL::libcurl). Cached structures remain usable.");
#else
    if(dest.has_parent_path()) fs::create_directories(dest.parent_path());
    fs::path part = dest;
    part += ".download";
    // An already present, invalid cache entry is never considered authoritative.
    const std::string url="https://files.rcsb.org/download/"+id+"."+o.format;
    FILE* file=std::fopen(part.string().c_str(),"wb");
    if(!file) throw std::runtime_error("Cannot create temporary RCSB file: "+part.string());
    CURL* curl=curl_easy_init();
    if(!curl) {std::fclose(file);fs::remove(part);throw std::runtime_error("Cannot initialize libcurl");}
    curl_easy_setopt(curl,CURLOPT_URL,url.c_str());
    curl_easy_setopt(curl,CURLOPT_FOLLOWLOCATION,1L);
    curl_easy_setopt(curl,CURLOPT_MAXREDIRS,3L);
    curl_easy_setopt(curl,CURLOPT_PROTOCOLS_STR,"https");
    curl_easy_setopt(curl,CURLOPT_REDIR_PROTOCOLS_STR,"https");
    curl_easy_setopt(curl,CURLOPT_CONNECTTIMEOUT,10L);
    curl_easy_setopt(curl,CURLOPT_TIMEOUT,90L);
    curl_easy_setopt(curl,CURLOPT_USERAGENT,"LigInt/0.7 RCSB-structure-downloader");
    curl_easy_setopt(curl,CURLOPT_WRITEFUNCTION,writeBytes);
    curl_easy_setopt(curl,CURLOPT_WRITEDATA,file);
    curl_easy_setopt(curl,CURLOPT_FAILONERROR,1L);
    CURLcode result=curl_easy_perform(curl);
    long status=0;curl_easy_getinfo(curl,CURLINFO_RESPONSE_CODE,&status);
    std::string error=curl_easy_strerror(result);
    curl_easy_cleanup(curl);std::fclose(file);
    if(result!=CURLE_OK || status!=200 || !validContents(part,o.format)) {
        fs::remove(part);
        throw std::runtime_error("RCSB download failed for "+id+" (HTTP "+std::to_string(status)+", "+error+"), or returned invalid structure data");
    }
    std::error_code ec;
    // Windows rename() cannot replace an existing file; preserve original until download validates.
    if(fs::exists(dest)) fs::remove(dest,ec);
    fs::rename(part,dest,ec);
    if(ec) {fs::remove(part);throw std::runtime_error("Cannot save downloaded structure: "+ec.message());}
    return dest.string();
#endif
}
}
