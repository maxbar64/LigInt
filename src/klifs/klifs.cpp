#include "ligint/klifs.hpp"
#include "ligint/ligint.hpp"
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
std::string color(const std::string& s) {static const std::map<std::string,std::string> colors={{"I","#e7dc71"},{"g.l","#8acb7b"},{"II","#f4cf65"},{"III","#e6ab63"},{"aC","#e47875"},{"αC","#e47875"},{"b.l","#7bc96f"},{"IV","#eee679"},{"V","#f0e877"},{"GK","#efcb64"},{"hinge","#cc83d8"},{"linker","#80dbe1"},{"aD","#df7f82"},{"αD","#df7f82"},{"aE","#d77a75"},{"αE","#d77a75"},{"aF","#da7778"},{"αF","#da7778"},{"VI","#f4e97d"},{"VII","#f3e579"},{"VIII","#f3e579"},{"c.l","#edb164"},{"DFG","#8ea5ec"},{"xDFG","#8ea5ec"},{"a.l","#8ea5ec"}};auto i=colors.find(s);return i==colors.end()?"#bdc6d0":i->second;}
// The 85 positions of the KLIFS pocket and their regions (the KLIFS alignment, as the OpenCADD client defines it)
const char* regionOf(int position){
    static const struct {int last;const char* name;} regions[]={{3,"I"},{9,"g.l"},{13,"II"},{19,"III"},{30,"αC"},{37,"b.l"},{41,"IV"},{44,"V"},{45,"GK"},
        {48,"hinge"},{52,"linker"},{59,"αD"},{64,"αE"},{67,"VI"},{75,"c.l"},{78,"VII"},{79,"VIII"},{83,"xDFG"},{85,"a.l"}};
    for(const auto& r:regions)if(position<=r.last)return r.name;
    return "";
}
char oneLetter(const std::string& res3){
    static const std::map<std::string,char> m={{"ALA",'A'},{"CYS",'C'},{"ASP",'D'},{"GLU",'E'},{"PHE",'F'},{"GLY",'G'},{"HIS",'H'},{"ILE",'I'},{"LYS",'K'},{"LEU",'L'},
        {"MET",'M'},{"ASN",'N'},{"PRO",'P'},{"GLN",'Q'},{"ARG",'R'},{"SER",'S'},{"THR",'T'},{"VAL",'V'},{"TRP",'W'},{"TYR",'Y'}};
    auto i=m.find(res3);return i==m.end()?'\0':i->second;
}
bool isPocketMol2(const std::string& text){
    size_t i=0;if(text.size()>=3&&(unsigned char)text[0]==0xef)i=3;
    while(i<text.size()&&std::isspace(static_cast<unsigned char>(text[i])))++i;
    return text.compare(i,10,"@<TRIPOS>M")==0;
}
// The pocket.mol2 of a KLIFS structure: its residues are, in order, the positions 1..85 that the structure has. The sequence says which positions are
// missing; without it the file must have all 85 residues.
std::vector<std::map<std::string,std::string>> pocketFromMol2(const std::string& text,const std::string& sequenceIn){
    std::string sequence;
    for(char c:sequenceIn)if(!std::isspace(static_cast<unsigned char>(c)))sequence+=c;
    if(!sequence.empty()){
        if(sequence.size()!=85)throw std::runtime_error("The KLIFS pocket sequence must have 85 characters, it has "+std::to_string(sequence.size()));
        for(char& c:sequence){if(c=='-')c='_';if(!(std::isalpha(static_cast<unsigned char>(c))||c=='_'))throw std::runtime_error(std::string("Invalid character in the KLIFS pocket sequence: '")+c+"'");}
    }
    struct Res{std::string name,number;};
    std::vector<Res> residues;
    std::istringstream in(text);std::string line;bool inSub=false;
    while(std::getline(in,line)){
        if(!line.empty()&&line.back()=='\r')line.pop_back();
        if(line.compare(0,9,"@<TRIPOS>")==0){inSub=line.compare(0,21,"@<TRIPOS>SUBSTRUCTURE")==0;continue;}
        if(!inSub||trim(line).empty()||trim(line)[0]=='#')continue;
        std::istringstream ls(line);std::vector<std::string> t;std::string w;while(ls>>w)t.push_back(w);
        if(t.size()<7)throw std::runtime_error("Invalid substructure line in the KLIFS pocket MOL2: '"+trim(line).substr(0,80)+"'");
        Res r;r.name=upper(t[6]);
        std::string label=t[1],number;
        if(upper(label).compare(0,r.name.size(),r.name)==0)number=label.substr(r.name.size());
        else{size_t k=label.find_first_of("0123456789-");number=k==std::string::npos?"":label.substr(k);}
        size_t k=0;if(k<number.size()&&number[k]=='-')++k;size_t d0=k;while(k<number.size()&&std::isdigit(static_cast<unsigned char>(number[k])))++k;
        if(k==d0||k-d0>6||(k<number.size()&&!(k+1==number.size()&&std::isalpha(static_cast<unsigned char>(number[k])))))
            throw std::runtime_error("Invalid residue in the KLIFS pocket MOL2: '"+label+"'");
        r.number=number;residues.push_back(r);
    }
    if(residues.empty())throw std::runtime_error("The KLIFS pocket MOL2 has no residues (no SUBSTRUCTURE section)");
    if(residues.size()>85)throw std::runtime_error("The KLIFS pocket MOL2 has "+std::to_string(residues.size())+" residues: a KLIFS pocket has at most 85");
    std::vector<int> positions;
    if(sequence.empty()){
        if(residues.size()!=85)throw std::runtime_error("The KLIFS pocket has "+std::to_string(residues.size())+" residues, so some of the 85 positions are missing and the file does not say which: "
            "give the 85-character pocket sequence ('_' marks a missing residue; --klifs-sequence, or the field `pocket` of the KLIFS service)");
        for(int p=1;p<=85;++p)positions.push_back(p);
    }else{
        for(int p=1;p<=85;++p)if(sequence[p-1]!='_')positions.push_back(p);
        if(positions.size()!=residues.size())throw std::runtime_error("The KLIFS pocket sequence has "+std::to_string(positions.size())+" residues but the pocket MOL2 has "+std::to_string(residues.size()));
        for(size_t i=0;i<positions.size();++i){
            const char expected=oneLetter(residues[i].name);
            if(expected&&std::toupper(static_cast<unsigned char>(sequence[positions[i]-1]))!=expected)
                throw std::runtime_error("Position "+std::to_string(positions[i])+": the pocket sequence has "+std::string(1,sequence[positions[i]-1])+" but the residue in the MOL2 is "+residues[i].name+residues[i].number);
        }
    }
    std::vector<std::map<std::string,std::string>> rows;
    for(size_t i=0;i<residues.size();++i)
        rows.push_back({{"residue.klifs_id",std::to_string(positions[i])},{"residue.id",residues[i].number},{"residue.klifs_region",regionOf(positions[i])},{"residue.aa",residues[i].name}});
    return rows;
}
std::string slurp(const std::string& filename){
    std::ifstream in(filename,std::ios::binary);if(!in)throw std::runtime_error("Cannot open KLIFS mapping: "+filename);
    std::ostringstream ss;ss<<in.rdbuf();return ss.str();
}
PocketDefinition convert(const std::vector<std::map<std::string,std::string>>& rows,const KlifsOptions& o) {
 PocketDefinition pocket; pocket.name="Kinase pocket (KLIFS, "+upper(o.pdbId)+" chain "+o.chain+")";
 std::map<int,PocketPosition> ordered;std::set<std::string> colors;
 for(auto& r:rows){auto p=field(r,{"residue.klifs_id","klifs_id","klifs_position","position"});if(p.empty()){std::string names;for(auto& kv:r)names+=(names.empty()?"":", ")+kv.first;throw std::runtime_error("KLIFS pocket row has no position (its fields are: "+names+")");}int n=0;try{size_t used=0;n=std::stoi(p,&used);if(used!=p.size()||n<1||n>85)throw std::invalid_argument("range");}catch(...){throw std::runtime_error("Invalid KLIFS pocket position: "+p);}
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
std::string fetch(const std::string& url){ CURL* c=curl_easy_init();if(!c)throw std::runtime_error("curl initialization failed");std::string result;curl_easy_setopt(c,CURLOPT_URL,url.c_str());curl_easy_setopt(c,CURLOPT_FOLLOWLOCATION,1L);curl_easy_setopt(c,CURLOPT_MAXREDIRS,3L);
#if LIBCURL_VERSION_NUM >= 0x075500
    curl_easy_setopt(c,CURLOPT_PROTOCOLS_STR,"https");curl_easy_setopt(c,CURLOPT_REDIR_PROTOCOLS_STR,"https");
#else
    curl_easy_setopt(c,CURLOPT_PROTOCOLS,static_cast<long>(CURLPROTO_HTTPS));curl_easy_setopt(c,CURLOPT_REDIR_PROTOCOLS,static_cast<long>(CURLPROTO_HTTPS));
#endif
    curl_easy_setopt(c,CURLOPT_NOSIGNAL,1L);curl_easy_setopt(c,CURLOPT_CONNECTTIMEOUT,10L);curl_easy_setopt(c,CURLOPT_TIMEOUT,25L);curl_easy_setopt(c,CURLOPT_USERAGENT,"LigInt/" LIGINT_VERSION " KLIFS-pocket");curl_easy_setopt(c,CURLOPT_WRITEFUNCTION,receive);curl_easy_setopt(c,CURLOPT_WRITEDATA,&result);auto rc=curl_easy_perform(c);long code=0;curl_easy_getinfo(c,CURLINFO_RESPONSE_CODE,&code);curl_easy_cleanup(c);if(rc!=CURLE_OK||code!=200)throw std::runtime_error("KLIFS request failed (HTTP "+std::to_string(code)+", curl "+std::to_string(rc)+"): "+url);return result;}
// the first characters of an answer, for the messages (control characters replaced by blanks)
std::string snippet(const std::string& body){std::string t=body.substr(0,120);for(char& c:t)if(static_cast<unsigned char>(c)<32)c=' ';return t+(body.size()>120?"...":"");}
struct Json {json_object* p;
    explicit Json(const std::string& s):p(json_tokener_parse(s.c_str())){if(!p)throw std::runtime_error("Invalid KLIFS JSON response (it starts with: '"+snippet(s)+"')");}
    ~Json(){json_object_put(p);}
    Json(const Json&)=delete;Json& operator=(const Json&)=delete;};
std::string val(json_object* o,const char* k){json_object* v=nullptr;
    if(!o||!json_object_object_get_ex(o,k,&v)||!v||json_object_get_type(v)==json_type_null)return {};
    return json_object_get_string(v);}
// the value of the first of these fields that exists (the service and the client libraries do not name them all alike)
std::string firstVal(json_object* o,std::initializer_list<const char*> keys){for(auto k:keys){auto v=val(o,k);if(!v.empty())return v;}return {};}
std::vector<json_object*> objects(json_object* root){std::vector<json_object*> result;
    if(!root)return result;
    if(json_object_get_type(root)==json_type_array){for(size_t i=0;i<json_object_array_length(root);++i)result.push_back(json_object_array_get_idx(root,i));}
    else if(json_object_get_type(root)==json_type_object)result.push_back(root);
    return result;}
std::string fieldsOf(json_object* o){std::string names;if(!o||json_object_get_type(o)!=json_type_object)return names;
    json_object_object_foreach(o,key,value){(void)value;names+=(names.empty()?"":", ")+std::string(key);}
    return names;}
std::string describe(const KlifsCandidate& c){return "[ID="+std::to_string(c.id)+" chain="+c.chain+" model="+c.model+"]";}

int choose(const std::string& body,const KlifsOptions& o){
    if(o.structureId>0)return o.structureId;
    Json j(body);
    std::vector<KlifsCandidate> all,matching;
    std::string fields;
    for(auto* r:objects(j.p)){
        if(json_object_get_type(r)!=json_type_object)continue;
        if(fields.empty())fields=fieldsOf(r);
        KlifsCandidate c;
        c.chain=firstVal(r,{"chain","structure.chain"});
        c.model=firstVal(r,{"alt","alternate_model","structure.alternate_model","alt_model","model","structure.alt"});
        c.pocket=firstVal(r,{"pocket","structure.pocket"});
        const std::string sid=firstVal(r,{"structure_ID","structure.klifs_id","structure_id","klifs_id","structure_IDs"});
        try{size_t used=0;c.id=std::stoi(sid,&used);if(used!=sid.size()||c.id<=0)c.id=0;}catch(...){c.id=0;}
        if(!c.id)throw std::runtime_error("KLIFS response has no numeric structure ID (its fields are: "+fields+")");
        all.push_back(c);
        if(upper(c.chain)==upper(o.chain)&&(o.model.empty()||upper(o.model)==upper(c.model)))matching.push_back(c);
    }
    if(matching.size()==1)return matching[0].id;
    std::ostringstream s;
    if(matching.empty()){
        s<<"KLIFS has no structure for chain "<<o.chain<<(o.model.empty()?"":" and model "+o.model)<<" of "<<upper(o.pdbId)<<". Its structures: ";
        for(const auto& c:all)s<<describe(c)<<' ';
        if(all.empty())s<<"(none)";
        throw std::runtime_error(s.str());
    }
    bool anyModel=false;for(const auto& c:all)if(!c.model.empty())anyModel=true;
    s<<"KLIFS matching structures: "<<matching.size()<<" (must be exactly 1). Candidates: ";
    for(const auto& c:matching)s<<describe(c)<<' ';
    s<<"Choose one with --klifs-id N";
    if(anyModel)s<<" or with --klifs-model M";
    else s<<". The alternate model is not in the answer under a name that ligint knows, so --klifs-model cannot choose (the fields of the answer are: "<<fields<<")";
    throw KlifsAmbiguousError(s.str(),matching);
}
// the pocket sequence of the structure `id` in the answer about the PDB entry, or "" (an answer that cannot be read is not an error here)
std::string pocketOf(const std::string& body,int id){
    json_object* root=json_tokener_parse(body.c_str());if(!root)return {};
    std::string found;
    for(auto* r:objects(root)){
        if(json_object_get_type(r)!=json_type_object)continue;
        const std::string sid=firstVal(r,{"structure_ID","structure.klifs_id","structure_id","klifs_id"});
        if(sid==std::to_string(id)){found=firstVal(r,{"pocket","structure.pocket"});break;}
    }
    json_object_put(root);return found;
}
// the names of the fields of the structure `id` in the answer about the PDB entry (of the first structure if `id` is not there), for the messages
std::string fieldsOfStructure(const std::string& body,int id){
    json_object* root=json_tokener_parse(body.c_str());if(!root)return {};
    std::string names,first;
    for(auto* r:objects(root)){
        if(json_object_get_type(r)!=json_type_object)continue;
        if(first.empty())first=fieldsOf(r);
        if(firstVal(r,{"structure_ID","structure.klifs_id","structure_id","klifs_id"})==std::to_string(id)){names=fieldsOf(r);break;}
    }
    json_object_put(root);return names.empty()?first:names;
}
bool validPdbId(const std::string& id){return id.size()==4&&std::all_of(id.begin(),id.end(),[](unsigned char c){return std::isalnum(c)!=0;});}
std::vector<std::map<std::string,std::string>> online(const KlifsOptions& o){
    // the list of the structures of the PDB entry gives the structure ID and the pocket sequence
    std::string list;
    if(validPdbId(o.pdbId))list=fetch("https://klifs.net/api_v2/structures_pdb_list?pdb-codes="+upper(o.pdbId));
    const int id=o.structureId>0?o.structureId:choose(list,o);
    std::string sequence=o.pocketSequence;
    if(sequence.empty()&&!list.empty())sequence=pocketOf(list,id);
    const std::string body=fetch("https://klifs.net/api_v2/structure_get_pocket?structure_ID="+std::to_string(id));
    if(isPocketMol2(body)){                                                    // what the service returns: the pocket as a MOL2 file
        try{return pocketFromMol2(body,sequence);}
        catch(const std::runtime_error& e){
            if(sequence.empty()&&!list.empty())throw std::runtime_error(std::string(e.what())+" The answer of the service about the structures has these fields: "+fieldsOfStructure(list,id));
            throw;
        }
    }
    Json data(body);
    std::vector<std::map<std::string,std::string>> rows;
    for(auto* item:objects(data.p)){std::map<std::string,std::string> r;
        if(json_object_get_type(item)!=json_type_object)continue;
        json_object_object_foreach(item,key,value){if(value&&json_object_get_type(value)!=json_type_null)r[key]=json_object_get_string(value);}
        rows.push_back(std::move(r));}
    if(rows.empty())throw std::runtime_error("No pocket residue rows returned by KLIFS for structure ID "+std::to_string(id)+" (the answer starts with: '"+snippet(body)+"')");
    return rows;
}
#endif
}
#ifdef LIGINT_KLIFS_ONLINE
int selectKlifsStructure(const std::string& structuresJson,const KlifsOptions& options){return choose(structuresJson,options);}
#else
int selectKlifsStructure(const std::string&,const KlifsOptions& options){
    if(options.structureId>0)return options.structureId;
    throw std::runtime_error("Online KLIFS unavailable: rebuild with libcurl and json-c development packages, or provide --klifs-mapping FILE");
}
#endif
bool klifsOnlineAvailable(){
#ifdef LIGINT_KLIFS_ONLINE
return true;
#else
return false;
#endif
}
PocketDefinition generateKlifsPocket(const KlifsOptions& o){if(o.chain.empty())throw std::runtime_error("KLIFS chain is required");if(o.mappingFile.empty()&&o.structureId<=0&&(o.pdbId.size()!=4||!std::all_of(o.pdbId.begin(),o.pdbId.end(),[](unsigned char c){return std::isalnum(c)!=0;})))throw std::runtime_error("KLIFS requires a four-character PDB ID");if(!o.mappingFile.empty()){const std::string text=slurp(o.mappingFile);
    if(isPocketMol2(text))return convert(pocketFromMol2(text,o.pocketSequence),o);
    return convert(readCsv(o.mappingFile),o);}
#ifdef LIGINT_KLIFS_ONLINE
return convert(online(o),o);
#else
throw std::runtime_error("Online KLIFS unavailable: rebuild with libcurl and json-c development packages, or provide --klifs-mapping FILE");
#endif
}
void saveKlifsPocket(const PocketDefinition& p,const std::string& filename){std::ofstream out(filename);if(!out)throw std::runtime_error("Cannot write KLIFS pocket: "+filename);out<<"# Generated by LigInt native KLIFS provider\n@name "<<p.name<<'\n';for(auto& c:p.colors)out<<"@color "<<c.first<<' '<<c.second<<'\n';for(auto& x:p.positions){out<<x.chain<<'\t'<<x.resi;if(x.icode!=' ')out<<x.icode;out<<'\t'<<x.label<<'\t'<<x.region;if(!x.resName.empty())out<<'\t'<<x.resName;out<<'\n';}if(!out)throw std::runtime_error("Error writing KLIFS pocket: "+filename);}
}
