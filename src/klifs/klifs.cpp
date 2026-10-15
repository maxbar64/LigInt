#include "ligint/klifs.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <set>
#include <vector>
#ifdef LIGINT_KLIFS_ONLINE
#include <curl/curl.h>
#include <json-c/json.h>
#endif
namespace ligint {
namespace {
std::string trim(std::string s) { auto p=s.find_first_not_of(" \t\r\n"); if(p==std::string::npos)return {}; auto q=s.find_last_not_of(" \t\r\n"); return s.substr(p,q-p+1); }
std::string upper(std::string s) { for(auto& c:s)c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));return s; }
std::string clean(std::string s) {for(auto& c:s) if(c==' '||c=='\t')c='_'; return s;}
std::vector<std::string> csvSplit(const std::string& line, char sep) {std::vector<std::string> v;std::string part;bool quoted=false;for(size_t i=0;i<line.size();++i){char c=line[i];if(c=='"'){if(quoted&&i+1<line.size()&&line[i+1]=='"'){part+='"';++i;}else quoted=!quoted;}else if(c==sep&&!quoted){v.push_back(trim(part));part.clear();}else part+=c;}v.push_back(trim(part));return v;}
std::string field(const std::map<std::string,std::string>& row, std::initializer_list<const char*> keys) { for(auto k:keys){auto i=row.find(k);if(i!=row.end() && !i->second.empty() && i->second!="_"&&i->second!="null"&&i->second!="None")return i->second;}return {}; }
std::string aa3(std::string aa) {aa=upper(aa);static const char* a[]={"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};std::string order="ACDEFGHIKLMNPQRSTVWY";if(aa.size()==1){auto n=order.find(aa[0]);return n==std::string::npos?"":a[n];}return aa.size()==3?aa:"";}
std::string color(const std::string& s) {static const std::map<std::string,std::string> colors={{"I","#e7dc71"},{"g.l","#8acb7b"},{"II","#f4cf65"},{"III","#e6ab63"},{"aC","#e47875"},{"b.l","#7bc96f"},{"IV","#eee679"},{"V","#f0e877"},{"GK","#efcb64"},{"hinge","#cc83d8"},{"linker","#80dbe1"},{"aD","#df7f82"},{"aE","#d77a75"},{"aF","#da7778"},{"VI","#f4e97d"},{"VII","#f3e579"},{"VIII","#f3e579"},{"c.l","#edb164"},{"DFG","#8ea5ec"},{"a.l","#8ea5ec"}};auto i=colors.find(s);return i==colors.end()?"#bdc6d0":i->second;}
PocketDefinition convert(const std::vector<std::map<std::string,std::string>>& rows,const KlifsOptions& o) {
 PocketDefinition pocket; pocket.name="Kinase pocket (KLIFS, "+upper(o.pdbId)+" chain "+o.chain+")";
 std::map<int,PocketPosition> ordered;std::set<std::string> colors;
 for(auto& r:rows){auto p=field(r,{"residue.klifs_id","klifs_id","klifs_position","position"});if(p.empty())throw std::runtime_error("KLIFS pocket row has no position");int n=0;try{size_t used=0;n=std::stoi(p,&used);if(used!=p.size()||n<1||n>85)throw std::invalid_argument("range");}catch(...){throw std::runtime_error("Invalid KLIFS pocket position: "+p);}
 if(ordered.count(n))throw std::runtime_error("Duplicate KLIFS position: "+p);
 auto rid=field(r,{"residue.id","residue_id","pdb_residue","resi"});if(rid.empty())continue;
 size_t used=0;int resi=0;try{resi=std::stoi(rid,&used);}catch(...){throw std::runtime_error("Invalid PDB residue: "+rid);}if(!used||used<rid.size()-1 || (used<rid.size()&&!std::isalpha(static_cast<unsigned char>(rid[used]))))throw std::runtime_error("Invalid PDB residue: "+rid);
 auto region=clean(field(r,{"residue.klifs_region","klifs_region","region"}));if(region.empty())throw std::runtime_error("KLIFS region missing at position "+p);
 PocketPosition pp;pp.chain=o.chain;pp.resi=resi;pp.icode=used<rid.size()?rid[used]:' ';pp.label=std::to_string(n);pp.region=region;pp.resName=aa3(field(r,{"residue.aa","residue.amino_acid","residue_name","resname"}));ordered[n]=pp;
 if(colors.insert(region).second)pocket.colors.push_back({region,color(region)});
 }
 if(ordered.empty())throw std::runtime_error("No mapped KLIFS residues");
 for(auto& p:ordered)pocket.positions.push_back(p.second);
 return pocket;
}
std::vector<std::map<std::string,std::string>> readCsv(const std::string& filename) {
 std::ifstream in(filename);if(!in)throw std::runtime_error("Cannot open KLIFS mapping: "+filename);std::string line;if(!std::getline(in,line))throw std::runtime_error("Empty KLIFS mapping");if(line.size()>2&&(unsigned char)line[0]==0xef)line.erase(0,3);char sep=line.find('\t')!=std::string::npos?'\t':',';auto headers=csvSplit(line,sep);std::vector<std::map<std::string,std::string>> rows;while(std::getline(in,line)){if(trim(line).empty()||line[0]=='#')continue;auto v=csvSplit(line,sep);std::map<std::string,std::string> r;for(size_t i=0;i<headers.size()&&i<v.size();++i)r[headers[i]]=v[i];rows.push_back(std::move(r));}return rows;
}
#ifdef LIGINT_KLIFS_ONLINE
size_t receive(char* ptr,size_t sz,size_t nm,void* user){auto& b=*static_cast<std::string*>(user);if(b.size()+sz*nm>8000000)return 0;b.append(ptr,sz*nm);return sz*nm;}
std::string fetch(const std::string& url){ CURL* c=curl_easy_init();if(!c)throw std::runtime_error("curl initialization failed");std::string result;curl_easy_setopt(c,CURLOPT_URL,url.c_str());curl_easy_setopt(c,CURLOPT_FOLLOWLOCATION,1L);curl_easy_setopt(c,CURLOPT_CONNECTTIMEOUT,10L);curl_easy_setopt(c,CURLOPT_TIMEOUT,25L);curl_easy_setopt(c,CURLOPT_USERAGENT,"LigInt/0.7 KLIFS-pocket");curl_easy_setopt(c,CURLOPT_WRITEFUNCTION,receive);curl_easy_setopt(c,CURLOPT_WRITEDATA,&result);auto rc=curl_easy_perform(c);long code=0;curl_easy_getinfo(c,CURLINFO_RESPONSE_CODE,&code);curl_easy_cleanup(c);if(rc!=CURLE_OK||code!=200)throw std::runtime_error("KLIFS request failed (HTTP "+std::to_string(code)+", curl "+std::to_string(rc)+"): "+url);return result;}
struct Json {json_object* p;explicit Json(const std::string& s):p(json_tokener_parse(s.c_str())){if(!p)throw std::runtime_error("Invalid KLIFS JSON response");}~Json(){json_object_put(p);} };
std::string val(json_object* o,const char* k){json_object* v=nullptr;if(!o||!json_object_object_get_ex(o,k,&v)||!v||json_object_get_type(v)==json_type_null)return {};return json_object_get_string(v);}
std::vector<json_object*> objects(json_object* root){std::vector<json_object*> result;if(!root)return result;if(json_object_get_type(root)==json_type_array){for(size_t i=0;i<json_object_array_length(root);++i)result.push_back(json_object_array_get_idx(root,i));}else if(json_object_get_type(root)==json_type_object)result.push_back(root);return result;}
std::vector<std::map<std::string,std::string>> online(const KlifsOptions& o){
 int id=o.structureId;
 if(!id){Json j(fetch("https://klifs.net/api_v2/structures_pdb_list?pdb-codes="+upper(o.pdbId)));std::vector<std::string> choices;int matches=0;
 for(auto* r:objects(j.p)){auto chain=val(r,"chain");if(chain.empty())chain=val(r,"structure.chain");auto model=val(r,"alternate_model");if(model.empty())model=val(r,"structure.alternate_model");auto sid=val(r,"structure_ID");if(sid.empty())sid=val(r,"structure.klifs_id");if(sid.empty())sid=val(r,"structure_IDs");choices.push_back("ID="+sid+" chain="+chain+" model="+model);if(upper(chain)==upper(o.chain)&&(o.model.empty()||upper(o.model)==upper(model))){++matches;try{id=std::stoi(sid);}catch(...){throw std::runtime_error("KLIFS response missing numeric structure ID");}}}
 if(matches!=1){std::ostringstream s;s<<"KLIFS matching structures: "<<matches<<" (must be exactly 1). Candidates: ";for(const auto& x:choices)s<<"["<<x<<"] ";throw std::runtime_error(s.str());}}
 Json data(fetch("https://klifs.net/api_v2/structure_get_pocket?structure_ID="+std::to_string(id)));
 std::vector<std::map<std::string,std::string>> rows;
 for(auto* item:objects(data.p)){std::map<std::string,std::string> r;if(json_object_get_type(item)!=json_type_object)continue;
 json_object_object_foreach(item,key,value){if(value&&json_object_get_type(value)!=json_type_null)r[key]=json_object_get_string(value);}
 rows.push_back(std::move(r));}
 if(rows.empty())throw std::runtime_error("No pocket residue rows returned by KLIFS for structure ID "+std::to_string(id));
 return rows;
}
#endif
}
bool klifsOnlineAvailable(){
#ifdef LIGINT_KLIFS_ONLINE
return true;
#else
return false;
#endif
}
PocketDefinition generateKlifsPocket(const KlifsOptions& o){if(o.chain.empty())throw std::runtime_error("KLIFS chain is required");if(o.mappingFile.empty()&&o.structureId<=0&&(o.pdbId.size()!=4||!std::all_of(o.pdbId.begin(),o.pdbId.end(),[](unsigned char c){return std::isalnum(c)!=0;})))throw std::runtime_error("KLIFS requires a four-character PDB ID");if(!o.mappingFile.empty())return convert(readCsv(o.mappingFile),o);
#ifdef LIGINT_KLIFS_ONLINE
return convert(online(o),o);
#else
throw std::runtime_error("Online KLIFS unavailable: rebuild with libcurl and json-c development packages, or provide --klifs-mapping FILE");
#endif
}
void saveKlifsPocket(const PocketDefinition& p,const std::string& filename){std::ofstream out(filename);if(!out)throw std::runtime_error("Cannot write KLIFS pocket: "+filename);out<<"# Generated by LigInt native KLIFS provider\n@name "<<p.name<<'\n';for(auto& c:p.colors)out<<"@color "<<c.first<<' '<<c.second<<'\n';for(auto& x:p.positions){out<<x.chain<<'\t'<<x.resi;if(x.icode!=' ')out<<x.icode;out<<'\t'<<x.label<<'\t'<<x.region;if(!x.resName.empty())out<<'\t'<<x.resName;out<<'\n';}if(!out)throw std::runtime_error("Error writing KLIFS pocket: "+filename);}
}
