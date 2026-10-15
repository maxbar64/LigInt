// Reader of PDBx/mmCIF files: only the _atom_site category (the atoms of the first model) is used.
//
// The file is read as a stream of CIF tokens (comments, quoted strings, semicolon text fields, loops) so that
// a structure with millions of atoms never has to be held in memory as text; the rows of _atom_site are
// converted to RawAtom one at a time. The chain and the residue number are the author's ones (auth_asym_id and
// auth_seq_id), as in the legacy PDB format.
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "internal.hpp"

namespace ligint {
namespace detail {
namespace {

enum class Tok { End, Tag, Value, Loop, Data, Other };

bool istartsWith(const std::string& s, const char* prefix) {
    size_t i = 0;
    for (; prefix[i]; i++)
        if (i >= s.size() || std::tolower(static_cast<unsigned char>(s[i])) != prefix[i]) return false;
    return true;
}

std::string lowered(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

class CifTokenizer {
public:
    explicit CifTokenizer(std::istream& in) : in_(in) {}

    Tok next(std::string& tok) {
        if (hasPending_) {
            hasPending_ = false;
            tok = std::move(pendingTok_);
            return pendingKind_;
        }
        for (;;) {
            if (!haveLine_) {
                if (!std::getline(in_, line_)) return Tok::End;
                if (!line_.empty() && line_.back() == '\r') line_.pop_back();
                pos_ = 0;
                haveLine_ = true;
                if (!line_.empty() && line_[0] == ';') {          // text field: up to the next line that starts with ';'
                    tok.assign(line_, 1, std::string::npos);
                    std::string l;
                    bool closed = false;
                    while (std::getline(in_, l)) {
                        if (!l.empty() && l.back() == '\r') l.pop_back();
                        if (!l.empty() && l[0] == ';') {
                            line_.assign(l, 1, std::string::npos);  // what follows the closing ';' is read as usual
                            pos_ = 0;
                            closed = true;
                            break;
                        }
                        tok += '\n';
                        tok += l;
                    }
                    if (!closed) haveLine_ = false;                  // unterminated: the rest of the file was the text
                    return Tok::Value;
                }
            }
            while (pos_ < line_.size() && (line_[pos_] == ' ' || line_[pos_] == '\t')) pos_++;
            if (pos_ >= line_.size() || line_[pos_] == '#') {        // end of line or comment
                haveLine_ = false;
                continue;
            }
            const char c = line_[pos_];
            if (c == '\'' || c == '"') {                             // quoted: ends at the quote followed by a blank or the line end
                size_t i = pos_ + 1;
                for (; i < line_.size(); i++)
                    if (line_[i] == c && (i + 1 == line_.size() || line_[i + 1] == ' ' || line_[i + 1] == '\t')) break;
                if (i >= line_.size()) {                             // unterminated: take the rest of the line
                    tok.assign(line_, pos_ + 1, std::string::npos);
                    pos_ = line_.size();
                } else {
                    tok.assign(line_, pos_ + 1, i - pos_ - 1);
                    pos_ = i + 1;
                }
                return Tok::Value;
            }
            size_t e = pos_;
            while (e < line_.size() && line_[e] != ' ' && line_[e] != '\t') e++;
            tok.assign(line_, pos_, e - pos_);
            pos_ = e;
            if (tok[0] == '_') return Tok::Tag;
            if (tok.size() == 5 && istartsWith(tok, "loop_")) return Tok::Loop;
            if (istartsWith(tok, "data_")) return Tok::Data;
            if (istartsWith(tok, "save_") || istartsWith(tok, "global_") || istartsWith(tok, "stop_")) return Tok::Other;
            return Tok::Value;
        }
    }

    void unread(Tok kind, std::string tok) {
        pendingKind_ = kind;
        pendingTok_ = std::move(tok);
        hasPending_ = true;
    }

private:
    std::istream& in_;
    std::string line_;
    size_t pos_ = 0;
    bool haveLine_ = false;
    bool hasPending_ = false;
    Tok pendingKind_ = Tok::End;
    std::string pendingTok_;
};

// column indices of the _atom_site items that are used (-1 = absent)
struct Columns {
    int group = -1, type = -1, labelAtom = -1, authAtom = -1, labelAlt = -1, labelComp = -1, authComp = -1;
    int labelAsym = -1, authAsym = -1, labelSeq = -1, authSeq = -1, ins = -1, x = -1, y = -1, z = -1, bf = -1, model = -1;
};

bool isNull(const std::string& s) { return s.empty() || s == "?" || s == "."; }

Columns mapColumns(const std::vector<std::string>& tags) {
    Columns c;
    for (size_t i = 0; i < tags.size(); i++) {
        const std::string& t = tags[i];       // lower case: CIF tags are case-insensitive
        const int k = static_cast<int>(i);
        if (t == "_atom_site.group_pdb") c.group = k;
        else if (t == "_atom_site.type_symbol") c.type = k;
        else if (t == "_atom_site.label_atom_id") c.labelAtom = k;
        else if (t == "_atom_site.auth_atom_id") c.authAtom = k;
        else if (t == "_atom_site.label_alt_id") c.labelAlt = k;
        else if (t == "_atom_site.label_comp_id") c.labelComp = k;
        else if (t == "_atom_site.auth_comp_id") c.authComp = k;
        else if (t == "_atom_site.label_asym_id") c.labelAsym = k;
        else if (t == "_atom_site.auth_asym_id") c.authAsym = k;
        else if (t == "_atom_site.label_seq_id") c.labelSeq = k;
        else if (t == "_atom_site.auth_seq_id") c.authSeq = k;
        else if (t == "_atom_site.pdbx_pdb_ins_code") c.ins = k;
        else if (t == "_atom_site.cartn_x") c.x = k;
        else if (t == "_atom_site.cartn_y") c.y = k;
        else if (t == "_atom_site.cartn_z") c.z = k;
        else if (t == "_atom_site.b_iso_or_equiv") c.bf = k;
        else if (t == "_atom_site.pdbx_pdb_model_num") c.model = k;
    }
    return c;
}

bool hasRequired(const Columns& c) {
    return c.x >= 0 && c.y >= 0 && c.z >= 0 && (c.authAtom >= 0 || c.labelAtom >= 0) && (c.authComp >= 0 || c.labelComp >= 0) &&
           (c.authAsym >= 0 || c.labelAsym >= 0);
}

double number(const std::string& s) {
    if (isNull(s)) return std::numeric_limits<double>::quiet_NaN();
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return end == s.c_str() ? std::numeric_limits<double>::quiet_NaN() : v;
}

// one row of values -> RawAtom; false if the row has no usable atom name or residue name
bool makeAtom(const std::vector<std::string>& v, const Columns& c, RawAtom& r) {
    auto at = [&](int i) -> const std::string& {
        static const std::string empty;
        return i >= 0 && i < static_cast<int>(v.size()) ? v[i] : empty;
    };
    auto pick = [&](int preferred, int fallback) -> const std::string& {
        return !isNull(at(preferred)) ? at(preferred) : at(fallback);
    };
    const std::string& name = pick(c.authAtom, c.labelAtom);
    const std::string& res = pick(c.authComp, c.labelComp);
    const std::string& chain = pick(c.authAsym, c.labelAsym);
    if (isNull(name) || isNull(res)) return false;
    r.name = name;
    r.res = res;
    r.chain = isNull(chain) ? " " : chain;
    const std::string& seq = pick(c.authSeq, c.labelSeq);
    r.resi = isNull(seq) ? 0 : std::atoi(seq.c_str());
    const std::string& ins = at(c.ins);
    r.icode = isNull(ins) ? ' ' : ins[0];
    const std::string& alt = at(c.labelAlt);
    r.alt = isNull(alt) ? ' ' : alt[0];
    r.x = number(at(c.x));
    r.y = number(at(c.y));
    r.z = number(at(c.z));
    const double bf = number(at(c.bf));
    r.bf = std::isfinite(bf) ? bf : 0.0;
    r.het = at(c.group) == "HETATM";
    std::string el = at(c.type);
    if (isNull(el)) {                                   // no element: first letter of the atom name
        el.clear();
        for (char ch : name) if (std::isalpha(static_cast<unsigned char>(ch))) { el = std::string(1, ch); break; }
    }
    r.el = el;
    return true;
}

}  // namespace

void readMmcifAtoms(std::istream& in, const AtomSink& sink, CifStats& stats) {
    CifTokenizer rd(in);
    std::string tok;
    bool done = false;
    long modelFirst = LONG_MIN;
    std::vector<std::pair<std::string, std::string>> single;   // _atom_site items written as "tag value" pairs (a one-atom file)

    // converts the rows of an _atom_site category; returns false when the first model is over (nothing else is needed)
    auto emitRow = [&](const std::vector<std::string>& vals, const Columns& cols) -> bool {
        if (cols.model >= 0) {
            const long m = std::strtol(vals[cols.model].c_str(), nullptr, 10);
            if (modelFirst == LONG_MIN) modelFirst = m;
            else if (m != modelFirst) return false;     // models are consecutive: the first one is complete
        }
        RawAtom r;
        if (!makeAtom(vals, cols, r)) { stats.badRows++; return true; }
        sink(r);
        return true;
    };

    Tok k;
    while (!done && (k = rd.next(tok)) != Tok::End) {
        if (k == Tok::Data) {
            if (stats.sawAtomSite) break;               // only the first data block with atoms
            single.clear();
            continue;
        }
        if (k == Tok::Loop) {
            std::vector<std::string> tags;
            std::string t2;
            Tok k2;
            while ((k2 = rd.next(t2)) == Tok::Tag) tags.push_back(lowered(t2));
            const bool atomSite = !tags.empty() && tags[0].compare(0, 11, "_atom_site.") == 0;
            if (!atomSite || stats.sawAtomSite) {
                while (k2 == Tok::Value) k2 = rd.next(t2);   // skip the values of another category
                rd.unread(k2, std::move(t2));
                continue;
            }
            stats.sawAtomSite = true;
            const Columns cols = mapColumns(tags);
            if (!hasRequired(cols)) {
                stats.warnings.push_back("mmCIF: the _atom_site category lacks the coordinates, the atom, the residue or the chain: no atoms read");
                done = true;
                break;
            }
            std::vector<std::string> vals(tags.size());
            size_t filled = 0;
            while (k2 == Tok::Value) {
                vals[filled++].swap(t2);
                if (filled == vals.size()) {
                    filled = 0;
                    if (!emitRow(vals, cols)) { done = true; break; }
                }
                k2 = rd.next(t2);
            }
            if (!done) {
                if (filled != 0) stats.warnings.push_back("mmCIF: incomplete last row of _atom_site ignored");
                rd.unread(k2, std::move(t2));
                done = true;                            // the atoms are all here: the rest (e.g. anisotropic factors) is not needed
            }
            continue;
        }
        if (k == Tok::Tag) {                             // "tag value" pair
            std::string v;
            const Tok kv = rd.next(v);
            if (kv != Tok::Value) { rd.unread(kv, std::move(v)); continue; }
            const std::string t = lowered(tok);
            if (t.compare(0, 11, "_atom_site.") == 0 && !stats.sawAtomSite) single.push_back({t, v});
        }
    }
    if (!stats.sawAtomSite && !single.empty()) {        // an _atom_site category written as pairs: a single atom
        std::vector<std::string> tags, vals;
        for (auto& kv : single) { tags.push_back(kv.first); vals.push_back(kv.second); }
        const Columns cols = mapColumns(tags);
        stats.sawAtomSite = true;
        if (hasRequired(cols)) emitRow(vals, cols);
        else stats.warnings.push_back("mmCIF: the _atom_site category lacks the coordinates, the atom, the residue or the chain: no atoms read");
    }
    if (!stats.sawAtomSite) stats.warnings.push_back("mmCIF: no _atom_site category found: no atoms read");
}

}  // namespace detail
}  // namespace ligint
