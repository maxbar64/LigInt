#include "SceneBuilder.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <QSet>
#include <algorithm>
#include <cmath>
#include <map>

#include "Style.h"

namespace gui {
namespace {

using ligint::V3;

std::uint32_t pc(const QColor& c) { return (static_cast<std::uint32_t>(c.red()) << 16) | (static_cast<std::uint32_t>(c.green()) << 8) | static_cast<std::uint32_t>(c.blue()); }

QColor mix(const QColor& a, const QColor& b, double t) {
    t = std::clamp(t, 0.0, 1.0);
    return QColor::fromRgbF(a.redF() + (b.redF() - a.redF()) * t, a.greenF() + (b.greenF() - a.greenF()) * t, a.blueF() + (b.blueF() - a.blueF()) * t);
}

V3 sub3(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V3 add3(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
V3 mul3(const V3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot3(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V3 cross3(const V3& a, const V3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
double len3(const V3& a) { return std::sqrt(dot3(a, a)); }
V3 unit3(const V3& a) { const double l = len3(a); return l > 1e-12 ? mul3(a, 1.0 / l) : V3{0, 0, 1}; }

// ---- colors by residue -----------------------------------------------------------------------
class ResidueColors {
public:
    ResidueColors(const std::vector<ligint::BackboneResidue>& bb, ProteinColor scheme) : bb_(bb), scheme_(scheme) {
        int chains = 0;
        for (const auto& r : bb) chains = std::max(chains, r.chainIndex + 1);
        first_.assign(chains, static_cast<int>(bb.size()));
        last_.assign(chains, -1);
        for (int i = 0; i < static_cast<int>(bb.size()); i++) {
            first_[bb[i].chainIndex] = std::min(first_[bb[i].chainIndex], i);
            last_[bb[i].chainIndex] = std::max(last_[bb[i].chainIndex], i);
        }
        std::vector<double> b;
        for (const auto& r : bb) b.push_back(r.bfactor);
        if (!b.empty()) {
            std::sort(b.begin(), b.end());
            lo_ = b[static_cast<std::size_t>(0.05 * (b.size() - 1))];
            hi_ = b[static_cast<std::size_t>(0.95 * (b.size() - 1))];
            if (hi_ - lo_ < 1e-6) hi_ = lo_ + 1;
        }
    }

    QColor atResidue(int k) const {
        if (k < 0 || k >= static_cast<int>(bb_.size())) return QColor("#b9c4d4");
        const auto& r = bb_[k];
        switch (scheme_) {
            case ProteinColor::SecondaryStructure:
                return r.ss == 'H' ? QColor("#e0527f") : (r.ss == 'E' ? QColor("#f2b94b") : QColor("#b4c2d6"));
            case ProteinColor::Chain: {
                static const char* pal[] = {"#5c9ded", "#f0a04b", "#7bc67b", "#d96ba0", "#9a7fd1", "#4fc3c0", "#e6c84f", "#c47a5a"};
                return QColor(pal[r.chainIndex % 8]);
            }
            case ProteinColor::Rainbow: {
                const int a = first_[r.chainIndex], b = last_[r.chainIndex];
                const double t = b > a ? static_cast<double>(k - a) / (b - a) : 0.0;
                return QColor::fromHsvF((1.0 - t) * 0.667, 0.62, 0.95);
            }
            case ProteinColor::BFactor: return bfactorColor(r.bfactor);
            case ProteinColor::Uniform: break;
        }
        return QColor("#b9c4d4");
    }

    // position along the chain, as in Mesh::tag
    QColor atFraction(float u) const {
        if (bb_.empty()) return QColor("#b9c4d4");
        const int n = static_cast<int>(bb_.size());
        if (scheme_ != ProteinColor::Rainbow && scheme_ != ProteinColor::BFactor)
            return atResidue(std::clamp(static_cast<int>(u + 0.5f), 0, n - 1));
        const int k = std::clamp(static_cast<int>(std::floor(u)), 0, n - 1);
        const int k2 = std::min(n - 1, k + 1);
        if (bb_[k].chainIndex != bb_[k2].chainIndex) return atResidue(k);
        return mix(atResidue(k), atResidue(k2), u - k);
    }

private:
    QColor bfactorColor(double b) const {
        const double t = std::clamp((b - lo_) / (hi_ - lo_), 0.0, 1.0);
        const QColor low("#3b6fd8"), mid("#f2f2f2"), high("#d83b3b");
        return t < 0.5 ? mix(low, mid, t * 2) : mix(mid, high, (t - 0.5) * 2);
    }
    const std::vector<ligint::BackboneResidue>& bb_;
    ProteinColor scheme_;
    std::vector<int> first_, last_;
    double lo_ = 0, hi_ = 1;
};

// Kyte-Doolittle hydropathy
double hydropathy(const std::string& res) {
    static const std::map<std::string, double> kd = {
        {"ILE", 4.5}, {"VAL", 4.2}, {"LEU", 3.8}, {"PHE", 2.8}, {"CYS", 2.5}, {"MET", 1.9}, {"ALA", 1.8}, {"GLY", -0.4}, {"THR", -0.7},
        {"SER", -0.8}, {"TRP", -0.9}, {"TYR", -1.3}, {"PRO", -1.6}, {"HIS", -3.2}, {"GLU", -3.5}, {"GLN", -3.5}, {"ASP", -3.5},
        {"ASN", -3.5}, {"LYS", -3.9}, {"ARG", -4.5}, {"MSE", 1.9}};
    auto it = kd.find(res);
    return it == kd.end() ? 0.0 : it->second;
}

QColor hydrophobicityColor(double kdValue) {
    const double t = std::clamp((kdValue + 4.5) / 9.0, 0.0, 1.0);
    const QColor polar("#5f9be0"), mid("#f3f1ec"), apolar("#eba23a");
    return t < 0.5 ? mix(polar, mid, t * 2) : mix(mid, apolar, (t - 0.5) * 2);
}

QColor surfaceElementColor(const char* el) {
    const std::string e = el;
    if (e == "N") return QColor("#6f8fe0");
    if (e == "O") return QColor("#e8665f");
    if (e == "S") return QColor("#e9c64f");
    return QColor("#cfd5dc");
}

// ---- geometry ----------------------------------------------------------------------------------
void appendCylinder(ligint::Mesh& m, std::vector<std::uint32_t>& colors, const V3& a, const V3& b, double radius, std::uint32_t color, int sides) {
    const V3 d = unit3(sub3(b, a));
    const V3 axis = std::fabs(d[0]) < 0.6 ? V3{1, 0, 0} : V3{0, 1, 0};
    const V3 u = unit3(cross3(d, axis));
    const V3 v = cross3(d, u);                        // u x v = d (right-handed)
    const int base = static_cast<int>(m.pos.size());
    for (int ring = 0; ring < 2; ring++) {
        const V3& c = ring == 0 ? a : b;
        for (int j = 0; j < sides; j++) {
            const double th = 2.0 * M_PI * j / sides;
            const V3 n = add3(mul3(u, std::cos(th)), mul3(v, std::sin(th)));
            const V3 p = add3(c, mul3(n, radius));
            m.pos.push_back({(float)p[0], (float)p[1], (float)p[2]});
            m.nrm.push_back({(float)n[0], (float)n[1], (float)n[2]});
            m.tag.push_back(0.f);
            colors.push_back(color);
        }
    }
    for (int j = 0; j < sides; j++) {
        const int j2 = (j + 1) % sides;
        const int a0 = base + j, a1 = base + j2, b0 = base + sides + j, b1 = base + sides + j2;
        m.tri.push_back({a0, a1, b0});
        m.tri.push_back({a1, b1, b0});
    }
}

// a bond as two half cylinders, one per atom color
void appendBond(AtomsPiece& piece, const V3& a, const V3& b, double radius, const QColor& ca, const QColor& cb) {
    const V3 mid = mul3(add3(a, b), 0.5);
    if (pc(ca) == pc(cb)) {
        appendCylinder(piece.sticks, piece.stickColors, a, b, radius, pc(ca), 10);
    } else {
        appendCylinder(piece.sticks, piece.stickColors, a, mid, radius, pc(ca), 10);
        appendCylinder(piece.sticks, piece.stickColors, mid, b, radius, pc(cb), 10);
    }
}

double vdwRadius(const QString& el) {
    if (el == "H") return 1.1;
    if (el == "N") return 1.55;
    if (el == "O") return 1.52;
    if (el == "S" || el == "P") return 1.8;
    if (el == "F") return 1.47;
    if (el == "Cl") return 1.75;
    if (el == "Br") return 1.85;
    if (el == "I") return 1.98;
    if (isMetalElement(el)) return 1.4;
    return 1.7;
}

std::shared_ptr<const AtomsPiece> atomsPiece(const std::vector<ligint::AtomRecord>& atoms, const std::vector<ligint::BondRecord>& bonds, bool ligand,
                                             LigandStyle style, double stickRadius, AromaticDisplay aromaticDisplay = AromaticDisplay::Dashed) {
    auto piece = std::make_shared<AtomsPiece>();
    std::vector<QColor> col(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); i++) col[i] = elementColor(QString::fromStdString(atoms[i].element), ligand);
    for (const auto& b : bonds) {
        if (b.a < 0 || b.b < 0 || b.a >= (int)atoms.size() || b.b >= (int)atoms.size()) continue;
        const V3 a=atoms[b.a].pos, z=atoms[b.b].pos;
        const V3 axis=sub3(z,a);
        const double length=len3(axis);
        if(length<1e-5)continue;
        // Choose a stable orthogonal displacement even for almost vertical bonds.
        const V3 reference=std::abs(axis[2]/length)<0.85?V3{0,0,1}:V3{0,1,0};
        const V3 perpendicular={axis[1]*reference[2]-axis[2]*reference[1],
                                axis[2]*reference[0]-axis[0]*reference[2],
                                axis[0]*reference[1]-axis[1]*reference[0]};
        const double norm=len3(perpendicular);
        if(norm<1e-7)continue;
        const V3 normal=mul3(perpendicular,1.0/norm);
        auto stroke=[&](double shift,double radius,const V3& start,const V3& end){
            const V3 delta=mul3(normal,shift);
            appendBond(*piece,add3(start,delta),add3(end,delta),radius,col[b.a],col[b.b]);
        };
        if(b.order==2||b.order==3){
            const double d=0.17;
            stroke(-d,stickRadius*0.72,a,z);
            stroke(+d,stickRadius*0.72,a,z);
            if(b.order==3)stroke(0,stickRadius*0.72,a,z);
        }else if(b.order==4 && aromaticDisplay!=AromaticDisplay::Single){
            // Aromatic connections: continuous outer bond and an inset dashed
            // parallel segment. Fused systems use the same per-edge convention.
            stroke(0,stickRadius,a,z);
            const V3 p=add3(a,mul3(axis,0.22)),q=add3(a,mul3(axis,0.46));
            const V3 u=add3(a,mul3(axis,0.56)),v=add3(a,mul3(axis,0.80));
            if(aromaticDisplay==AromaticDisplay::KekuleLike){
                // This is an illustrative alternating appearance, not a formally assigned Kekule structure.
                stroke(0.16,stickRadius*0.48,p,v);
            }else{
                stroke(0.16,stickRadius*0.48,p,q);
                stroke(0.16,stickRadius*0.48,u,v);
            }
        }else appendBond(*piece,a,z,stickRadius,col[b.a],col[b.b]);
    }
    for (std::size_t i = 0; i < atoms.size(); i++) {
        const QString el = QString::fromStdString(atoms[i].element);
        double r = stickRadius;                                   // rounds the joints of the sticks
        if (isMetalElement(el)) r = 0.8;
        else if (style == LigandStyle::BallStick) r = 0.32 * elementRadius(el);
        else if (style == LigandStyle::Spheres) r = vdwRadius(el);
        piece->spheres.push_back({(float)atoms[i].pos[0], (float)atoms[i].pos[1], (float)atoms[i].pos[2], (float)r, pc(col[i])});
    }
    return piece;
}

}  // namespace

std::uint32_t packColor(const QColor& c) { return pc(c); }

QColor backgroundColor(Background b) {
    switch (b) {
        case Background::Light: return QColor("#e9edf2");
        case Background::Dark: return QColor("#14181f");
        case Background::White: break;
    }
    return QColor("#ffffff");
}

std::shared_ptr<const MeshPiece> makeCartoonPiece(const ligint::Result& r, const ViewStyle& style) {
    if (style.protein == ProteinStyle::Off || r.backbone.size() < 2) return nullptr;
    ligint::CartoonOptions o;
    if (style.protein == ProteinStyle::Tube) { o.tubeOnly = true; o.tubeRadius = 0.45; }
    if (style.protein == ProteinStyle::Trace) { o.tubeOnly = true; o.tubeRadius = 0.17; }
    auto piece = std::make_shared<MeshPiece>();
    piece->mesh = ligint::buildCartoon(r.backbone, o);
    ResidueColors rc(r.backbone, style.proteinColor);
    piece->colors.resize(piece->mesh.pos.size());
    for (std::size_t i = 0; i < piece->colors.size(); i++) piece->colors[i] = pc(rc.atFraction(piece->mesh.tag[i]));
    return piece;
}

std::shared_ptr<const MeshPiece> makeSurfacePiece(const ligint::Result& r, const ligint::LigandResult& lig, const ViewStyle& style) {
    if (style.surface == SurfaceMode::Off || r.proteinAtoms.empty()) return nullptr;
    std::vector<V3> ligandPoints;
    for (const auto& a : lig.atoms) ligandPoints.push_back(a.pos);
    const double selectRadius = 17.0, showRadius = 10.0;

    // atoms of the protein: not the ligand itself (composite ligands remain in the list), not ions or cofactors
    std::vector<ligint::ProteinAtom> sel;
    for (const auto& a : r.proteinAtoms) {
        if (a.hetero) continue;
        bool isLigand = false, near = style.surface == SurfaceMode::Protein;
        for (const auto& p : ligandPoints) {
            const double d2 = dot3(sub3(a.pos, p), sub3(a.pos, p));
            if (d2 < 0.0025) { isLigand = true; break; }
            if (d2 < selectRadius * selectRadius) near = true;
        }
        if (!isLigand && near) sel.push_back(a);
    }
    if (sel.empty()) return nullptr;

    ligint::SurfaceOptions so;
    so.kind = style.surfaceType == SurfaceType::SolventExcluded ? ligint::SurfaceKind::SolventExcluded : ligint::SurfaceKind::SolventAccessible;
    const std::size_t n = sel.size();
    so.spacing = style.surface == SurfaceMode::Pocket ? 0.5 : (n <= 4000 ? 0.6 : (n <= 12000 ? 0.8 : (n <= 30000 ? 1.0 : 1.3)));
    auto piece = std::make_shared<MeshPiece>();
    piece->mesh = ligint::buildSurface(sel, so);
    if (style.surface == SurfaceMode::Pocket && !ligandPoints.empty()) piece->mesh = ligint::clipMeshNear(piece->mesh, ligandPoints, showRadius);

    ResidueColors rc(r.backbone, style.proteinColor);
    piece->colors.resize(piece->mesh.pos.size());
    for (std::size_t i = 0; i < piece->colors.size(); i++) {
        const int ai = std::clamp(static_cast<int>(piece->mesh.tag[i] + 0.5f), 0, static_cast<int>(sel.size()) - 1);
        const ligint::ProteinAtom& atom = sel[ai];
        QColor c("#b3c3d6");
        switch (style.surfaceColor) {
            case SurfaceColor::Neutral: break;
            case SurfaceColor::Atoms: c = surfaceElementColor(atom.element); break;
            case SurfaceColor::Hydrophobicity:
                c = atom.residue >= 0 ? hydrophobicityColor(hydropathy(r.backbone[atom.residue].res)) : QColor("#d9d9d9");
                break;
            case SurfaceColor::Protein: c = rc.atResidue(atom.residue); break;
        }
        piece->colors[i] = pc(c);
    }
    return piece;
}

std::shared_ptr<const AtomsPiece> makeLigandPiece(const ligint::LigandResult& lig, const ViewStyle& style) {
    if (!style.showLigand) return nullptr;
    const double radius = style.ligand == LigandStyle::Sticks ? 0.15 : 0.11;
    return atomsPiece(lig.atoms, lig.bonds, true, style.ligand, radius, style.aromaticDisplay);
}

std::shared_ptr<const AtomsPiece> makeResiduePiece(const ligint::LigandResult& lig, const ViewStyle& style) {
    if (style.residues == ResidueStyle::Off) return nullptr;
    return atomsPiece(lig.contactAtoms, lig.contactBonds, false, LigandStyle::Sticks, 0.095);
}

std::shared_ptr<const AtomsPiece> makeWaterPiece(const ligint::LigandResult& lig, const ViewStyle& style) {
    if (!style.showWaters) return nullptr;
    auto piece = std::make_shared<AtomsPiece>();
    std::vector<V3> seen;
    auto add = [&](const V3& p) {
        for (const auto& q : seen) if (len3(sub3(p, q)) < 0.2) return;
        seen.push_back(p);
        piece->spheres.push_back({(float)p[0], (float)p[1], (float)p[2], 0.5f, pc(QColor("#74c0fc"))});
    };
    for (const auto& w : lig.waters) add(w);
    for (const auto& h : lig.hits) if (h.hasWater) add(h.waterPos);
    return piece;
}

}  // namespace gui
