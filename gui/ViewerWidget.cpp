#include "ViewerWidget.h"

#include <QFontMetricsF>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>
#include <map>

#include "Style.h"

namespace gui {
namespace {

using Mat3 = std::array<double, 9>;
using V3 = ligint::V3;

QPointF eventPos(const QMouseEvent* e) {
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    return e->position();
#else
    return e->localPos();
#endif
}

Mat3 mul(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) r[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
    return r;
}

Mat3 rotY(double t) { return {std::cos(t), 0, std::sin(t), 0, 1, 0, -std::sin(t), 0, std::cos(t)}; }
Mat3 rotX(double t) { return {1, 0, 0, 0, std::cos(t), -std::sin(t), 0, std::sin(t), std::cos(t)}; }

V3 rotateVec(const Mat3& m, const V3& v) {
    return {m[0] * v[0] + m[1] * v[1] + m[2] * v[2], m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
            m[6] * v[0] + m[7] * v[1] + m[8] * v[2]};
}

// Eigenvectors of a symmetric 3x3 matrix (Jacobi). Returns the rows sorted by decreasing eigenvalue.
Mat3 principalAxes(const std::vector<V3>& pts) {
    V3 c{0, 0, 0};
    for (auto& p : pts) for (int i = 0; i < 3; i++) c[i] += p[i] / std::max<size_t>(1, pts.size());
    double a[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (auto& p : pts)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) a[i][j] += (p[i] - c[i]) * (p[j] - c[j]);
    double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 30; sweep++) {
        double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        if (off < 1e-12) break;
        for (int p = 0; p < 2; p++)
            for (int q = p + 1; q < 3; q++) {
                if (std::fabs(a[p][q]) < 1e-14) continue;
                double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
                double t = (theta >= 0 ? 1 : -1) / (std::fabs(theta) + std::sqrt(theta * theta + 1));
                double cs = 1 / std::sqrt(t * t + 1), sn = t * cs;
                for (int k = 0; k < 3; k++) {
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = cs * akp - sn * akq; a[k][q] = sn * akp + cs * akq;
                }
                for (int k = 0; k < 3; k++) {
                    double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = cs * apk - sn * aqk; a[q][k] = sn * apk + cs * aqk;
                }
                for (int k = 0; k < 3; k++) {
                    double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = cs * vkp - sn * vkq; v[k][q] = sn * vkp + cs * vkq;
                }
            }
    }
    int order[3] = {0, 1, 2};
    std::sort(order, order + 3, [&](int x, int y) { return a[x][x] > a[y][y]; });
    V3 e1{v[0][order[0]], v[1][order[0]], v[2][order[0]]};
    V3 e2{v[0][order[1]], v[1][order[1]], v[2][order[1]]};
    V3 e3{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]};   // right-handed
    return {e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2]};
}

}  // namespace

ViewerWidget::ViewerWidget(QWidget* parent) : QWidget(parent) {
    setMinimumSize(320, 240);
    setMouseTracking(false);
    setFocusPolicy(Qt::StrongFocus);
    setToolTip(tr("Drag: rotate · right/middle button: pan · wheel: zoom · double click: reset"));
}

void ViewerWidget::setLigand(const ligint::LigandResult& r, const std::vector<ligint::AtomRecord>& proteinAtoms,
                             const std::vector<ligint::BondRecord>& proteinBonds) {
    data_ = r;
    proteinAtoms_ = proteinAtoms;
    proteinBonds_ = proteinBonds;
    hasData_ = !data_.atoms.empty();
    highlight_ = -1;
    orient();
    update();
}

void ViewerWidget::clear() {
    data_ = ligint::LigandResult();
    proteinAtoms_.clear();
    proteinBonds_.clear();
    hasData_ = false;
    highlight_ = -1;
    update();
}

void ViewerWidget::orient(bool includeProtein) {
    std::vector<V3> pts;
    if (showLigand_) for (auto& a : data_.atoms) pts.push_back(a.pos);
    if (includeProtein && showProtein_) for (auto& a : proteinAtoms_) pts.push_back(a.pos);
    for (auto& h : data_.hits) { pts.push_back(h.ligandPos); pts.push_back(h.proteinPos); }
    if (pts.empty()) return;
    V3 lo = pts[0], hi = pts[0];
    for (auto& p : pts) for (int i = 0; i < 3; i++) { lo[i] = std::min(lo[i], p[i]); hi[i] = std::max(hi[i], p[i]); }
    for (int i = 0; i < 3; i++) center_[i] = 0.5 * (lo[i] + hi[i]);
    radius_ = 0;
    for (auto& p : pts) {
        double d = std::sqrt(std::pow(p[0] - center_[0], 2) + std::pow(p[1] - center_[1], 2) + std::pow(p[2] - center_[2], 2));
        radius_ = std::max(radius_, d);
    }
    const double hitRadius = std::max(radius_, 4.0);
    for (auto& a : data_.contactAtoms) {
        double d = std::sqrt(std::pow(a.pos[0] - center_[0], 2) + std::pow(a.pos[1] - center_[1], 2) + std::pow(a.pos[2] - center_[2], 2));
        radius_ = std::max(radius_, d);
    }
    radius_ = std::min(radius_, 1.35 * hitRadius);       // the farthest residues may fall outside the frame
    radius_ = std::max(radius_, 4.0);
    // initial orientation: the ligand (with its contacts) is shown with the maximum extent on screen
    rot_ = principalAxes(pts);
    zoom_ = 1.0;
    pan_ = QPointF(0, 0);
}

void ViewerWidget::resetView() {
    if (hasData_) orient(showProtein_);
    update();
}

void ViewerWidget::focusLigand() {
    if (hasData_) orient(false);
    update();
}

void ViewerWidget::fitProtein() {
    showProtein_ = true;
    if (hasData_) orient(true);
    update();
}

void ViewerWidget::setProteinStyle(int style) {
    proteinStyle_ = std::clamp(style, 0, 2);
    update();
}

void ViewerWidget::setVisibleTypes(const QSet<QString>& types) {
    visible_ = types;
    allVisible_ = false;
    update();
}
void ViewerWidget::setHighlight(int i) { highlight_ = i; update(); }
void ViewerWidget::setShowLabels(bool on) { showLabels_ = on; update(); }
void ViewerWidget::setShowResidues(bool on) { showResidues_ = on; update(); }
void ViewerWidget::setShowLigand(bool on) { showLigand_ = on; update(); }
void ViewerWidget::setShowProtein(bool on) { showProtein_ = on; if (hasData_) orient(showProtein_); update(); }
void ViewerWidget::setShowDistances(bool on) { showDistances_ = on; update(); }
void ViewerWidget::setShowLegend(bool on) { showLegend_ = on; update(); }

bool ViewerWidget::typeVisible(const std::string& type) const {
    return allVisible_ || visible_.contains(QString::fromStdString(type));
}

double ViewerWidget::fitScale(const QRectF& rect) const {
    return 0.47 * std::min(rect.width(), rect.height()) / radius_ * zoom_;
}

ViewerWidget::Proj ViewerWidget::project(const V3& v, const QRectF& rect) const {
    V3 q{v[0] - center_[0], v[1] - center_[1], v[2] - center_[2]};
    V3 r = rotateVec(rot_, q);
    double s = fitScale(rect);
    return {QPointF(rect.center().x() + pan_.x() + s * r[0], rect.center().y() + pan_.y() - s * r[1]), r[2]};
}

void ViewerWidget::render(QPainter& p, const QRectF& rect) const {
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setRenderHint(QPainter::TextAntialiasing, true);
    p.fillRect(rect, QColor("#fbfbfd"));
    if (!hasData_) {
        p.setPen(QColor("#868e96"));
        p.drawText(rect, Qt::AlignCenter, QObject::tr("Open a structure file, choose a ligand and press \u00abAnalyze\u00bb"));
        return;
    }
    const double s = fitScale(rect);

    // ---- depth, for color attenuation
    double zmin = 1e9, zmax = -1e9;
    auto track = [&](const V3& v) { double z = project(v, rect).z; zmin = std::min(zmin, z); zmax = std::max(zmax, z); };
    if (showLigand_) for (auto& a : data_.atoms) track(a.pos);
    if (showProtein_) for (auto& a : proteinAtoms_) track(a.pos);
    for (auto& a : data_.contactAtoms) track(a.pos);
    if (zmax - zmin < 1e-6) zmax = zmin + 1;
    auto depthAlpha = [&](double z) { return 0.55 + 0.45 * (z - zmin) / (zmax - zmin); };

    struct Item { double z; int kind; int a, b; bool lig; };   // kind 0 = bond, 1 = atom, 2/3 = protein bond/CA
    std::vector<Item> items;
    auto zOf = [&](const ligint::AtomRecord& a) { return project(a.pos, rect).z; };
    if (showLigand_) for (size_t i = 0; i < data_.bonds.size(); i++) {
        const auto& b = data_.bonds[i];
        items.push_back({(zOf(data_.atoms[b.a]) + zOf(data_.atoms[b.b])) / 2, 0, (int)i, 0, true});
    }
    if (showLigand_) for (size_t i = 0; i < data_.atoms.size(); i++) items.push_back({zOf(data_.atoms[i]), 1, (int)i, 0, true});
    if (showProtein_) for (size_t i = 0; i < proteinBonds_.size(); i++) {
        const auto& b = proteinBonds_[i];
        items.push_back({(zOf(proteinAtoms_[b.a]) + zOf(proteinAtoms_[b.b])) / 2, 2, (int)i, 0, false});
    }
    if (showProtein_ && proteinStyle_ == 2)
        for (size_t i = 0; i < proteinAtoms_.size(); i++) items.push_back({zOf(proteinAtoms_[i]), 3, (int)i, 0, false});
    if (showResidues_) {
        for (size_t i = 0; i < data_.contactBonds.size(); i++) {
            const auto& b = data_.contactBonds[i];
            items.push_back({(zOf(data_.contactAtoms[b.a]) + zOf(data_.contactAtoms[b.b])) / 2, 0, (int)i, 0, false});
        }
        for (size_t i = 0; i < data_.contactAtoms.size(); i++) items.push_back({zOf(data_.contactAtoms[i]), 1, (int)i, 0, false});
    }
    std::sort(items.begin(), items.end(), [](const Item& x, const Item& y) { return x.z < y.z; });

    auto elColor = [&](const ligint::AtomRecord& a, bool lig) { return elementColor(QString::fromStdString(a.element), lig); };
    for (const Item& it : items) {
        if (it.kind == 2) {
            const auto& bd = proteinBonds_[it.a];
            const Proj p1 = project(proteinAtoms_[bd.a].pos, rect), p2 = project(proteinAtoms_[bd.b].pos, rect);
            const double alpha = proteinStyle_ == 1 ? 0.72 : 0.60;
            const double width = proteinStyle_ == 1 ? std::clamp(0.18 * s, 3.0, 7.0) : std::clamp(0.07 * s, 0.9, 2.0);
            p.setOpacity(alpha * depthAlpha(it.z));
            p.setPen(QPen(QColor("#1971c2"), width, Qt::SolidLine, Qt::RoundCap));
            p.drawLine(p1.pt, p2.pt);
        } else if (it.kind == 3) {
            const Proj pr = project(proteinAtoms_[it.a].pos, rect);
            const double rad = std::clamp(0.16 * s, 1.5, 4.0);
            p.setOpacity(0.78 * depthAlpha(it.z));
            p.setPen(QPen(QColor("#1864ab"), 0.8));
            p.setBrush(QColor("#74c0fc"));
            p.drawEllipse(pr.pt, rad, rad);
        } else if (it.kind == 0) {
            const auto& A = it.lig ? data_.atoms : data_.contactAtoms;
            const auto& B = it.lig ? data_.bonds : data_.contactBonds;
            const auto& bd = B[it.a];
            Proj p1 = project(A[bd.a].pos, rect), p2 = project(A[bd.b].pos, rect);
            QColor c1 = elColor(A[bd.a], it.lig), c2 = elColor(A[bd.b], it.lig);
            double w = it.lig ? std::clamp(0.17 * s, 2.0, 6.0) : std::clamp(0.10 * s, 1.2, 3.0);
            p.setOpacity(depthAlpha(it.z));
            QPointF mid = (p1.pt + p2.pt) / 2;
            QPointF d = p2.pt - p1.pt;
            double len = std::hypot(d.x(), d.y());
            QPointF nrm = len > 1e-6 ? QPointF(-d.y() / len, d.x() / len) : QPointF(0, 0);
            std::vector<double> offs;
            if (it.lig && bd.order == 2) offs = {-0.55 * w, 0.55 * w};
            else if (it.lig && bd.order == 3) offs = {-1.1 * w, 0, 1.1 * w};
            else offs = {0};
            double lw = offs.size() > 1 ? w * 0.6 : w;
            for (double o : offs) {
                QPointF sh = nrm * o;
                p.setPen(QPen(c1, lw, Qt::SolidLine, Qt::RoundCap));
                p.drawLine(p1.pt + sh, mid + sh);
                p.setPen(QPen(c2, lw, Qt::SolidLine, Qt::RoundCap));
                p.drawLine(mid + sh, p2.pt + sh);
            }
        } else {
            const auto& a = (it.lig ? data_.atoms : data_.contactAtoms)[it.a];
            Proj pr = project(a.pos, rect);
            QString el = QString::fromStdString(a.element);
            double rad = isMetalElement(el) ? 0.9 * s : (it.lig ? 0.30 : 0.20) * elementRadius(el) * s;
            rad = std::clamp(rad, 2.0, isMetalElement(el) ? 22.0 : 14.0);
            p.setOpacity(depthAlpha(it.z));
            QColor col = elColor(a, it.lig);
            p.setPen(QPen(col.darker(150), 1.0));
            p.setBrush(col);
            p.drawEllipse(pr.pt, rad, rad);
        }
    }
    p.setOpacity(1.0);

    // ---- interactions
    QFont small = p.font();
    small.setPointSizeF(std::max(7.5, small.pointSizeF() * 0.9));
    std::map<std::string, QPointF> labelPos;
    std::vector<std::pair<std::string, int>> typesUsed;
    for (size_t i = 0; i < data_.hits.size(); i++) {
        const auto& h = data_.hits[i];
        if (!typeVisible(h.type)) continue;
        const QString typeName = QString::fromStdString(h.type);
        const TypeStyle& st = styleForType(typeName);
        QPointF a = project(h.ligandPos, rect).pt, b = project(h.proteinPos, rect).pt;
        const bool hl = static_cast<int>(i) == highlight_;
        auto strokeLine = [&](QPointF from, QPointF to) {
            if (hl) {
                p.setPen(QPen(QColor(255, 255, 255, 230), st.width + 5, Qt::SolidLine, Qt::RoundCap));
                p.drawLine(from, to);
            }
            QPen pen(st.color, hl ? st.width + 1.5 : st.width, st.pen, Qt::RoundCap);
            p.setPen(pen);
            p.drawLine(from, to);
        };
        if (h.hasWater) {
            QPointF w = project(h.waterPos, rect).pt;
            strokeLine(a, w);
            strokeLine(w, b);
            const double wr = std::clamp(0.45 * s, 3.0, 9.0);
            p.setPen(QPen(QColor("#1971c2"), 1.0));
            p.setBrush(QColor("#a5d8ff"));
            p.drawEllipse(w, wr, wr);
        } else {
            strokeLine(a, b);
        }
        if (h.type.rfind("pi-stacking", 0) == 0 || h.type == ligint::type::CationPi) {
            p.setPen(QPen(st.color, 1.6));
            p.setBrush(Qt::NoBrush);
            p.drawEllipse(a, 6, 6);
            p.drawEllipse(b, 6, 6);
        }
        if (!labelPos.count(h.residue)) labelPos[h.residue] = b;
        if (showDistances_ && (hl || (h.type != ligint::type::Hydrophobic && h.type != ligint::type::BackboneApolar))) {
            QPointF mid = (a + b) / 2;
            QString t = QString::number(h.distance, 'f', 2);
            p.setFont(small);
            QFontMetricsF fm(small);
            QRectF r(mid.x() - fm.horizontalAdvance(t) / 2 - 2, mid.y() - fm.height() / 2, fm.horizontalAdvance(t) + 4, fm.height());
            p.setPen(Qt::NoPen);
            p.setBrush(QColor(255, 255, 255, 200));
            p.drawRoundedRect(r, 2, 2);
            p.setPen(st.color.darker(130));
            p.drawText(r, Qt::AlignCenter, t);
        }
        bool seen = false;
        for (auto& u : typesUsed) if (u.first == h.type) { u.second++; seen = true; }
        if (!seen) typesUsed.push_back({h.type, 1});
    }

    // ---- residue labels
    if (showLabels_) {
        p.setFont(small);
        QFontMetricsF fm(small);
        for (auto& kv : labelPos) {
            QString t = QString::fromStdString(kv.first);
            QPointF at = kv.second + QPointF(7, -7);
            QRectF r(at.x(), at.y() - fm.ascent(), fm.horizontalAdvance(t) + 6, fm.height() + 2);
            p.setPen(Qt::NoPen);
            p.setBrush(QColor(255, 255, 255, 215));
            p.drawRoundedRect(r, 3, 3);
            p.setPen(QColor("#212529"));
            p.drawText(r, Qt::AlignCenter, t);
        }
    }

    // ---- legend and title
    p.setFont(small);
    QFontMetricsF fm(small);
    if (showLegend_ && !typesUsed.empty()) {
        double w = 0;
        for (auto& u : typesUsed) w = std::max(w, fm.horizontalAdvance(QString::fromStdString(u.first) + QString(" (%1)").arg(u.second)));
        const double lineH = fm.height() + 3;
        QRectF box(rect.left() + 8, rect.top() + 8, w + 52, lineH * typesUsed.size() + 8);
        p.setPen(QPen(QColor("#dee2e6"), 1));
        p.setBrush(QColor(255, 255, 255, 225));
        p.drawRoundedRect(box, 4, 4);
        double y = box.top() + 4;
        for (auto& u : typesUsed) {
            const QString legendType = QString::fromStdString(u.first);
            const TypeStyle& st = styleForType(legendType);
            p.setPen(QPen(st.color, st.width, st.pen, Qt::RoundCap));
            p.drawLine(QPointF(box.left() + 6, y + lineH / 2), QPointF(box.left() + 34, y + lineH / 2));
            p.setPen(QColor("#343a40"));
            p.drawText(QPointF(box.left() + 42, y + fm.ascent() + 1), QString::fromStdString(u.first) + QString(" (%1)").arg(u.second));
            y += lineH;
        }
    }
    p.setPen(QColor("#495057"));
    p.drawText(QPointF(rect.left() + 10, rect.bottom() - 10), QString::fromStdString(data_.id));
}

void ViewerWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    render(p, QRectF(rect()));
}

QImage ViewerWidget::renderImage(const QSize& size) {
    QImage img(size, QImage::Format_ARGB32);
    img.fill(Qt::white);
    QPainter p(&img);
    render(p, QRectF(QPointF(0, 0), QSizeF(size)));
    p.end();
    return img;
}

void ViewerWidget::mousePressEvent(QMouseEvent* e) { last_ = eventPos(e).toPoint(); }

void ViewerWidget::mouseMoveEvent(QMouseEvent* e) {
    const QPoint cur = eventPos(e).toPoint();
    QPoint d = cur - last_;
    last_ = cur;
    if (!hasData_) return;
    if (e->buttons() & Qt::LeftButton) {
        rot_ = mul(rotX(d.y() * 0.01), mul(rotY(d.x() * 0.01), rot_));
    } else if (e->buttons() & (Qt::RightButton | Qt::MiddleButton)) {
        pan_ += QPointF(d);
    } else {
        return;
    }
    update();
}

void ViewerWidget::wheelEvent(QWheelEvent* e) {
    if (!hasData_) return;
    zoom_ = std::clamp(zoom_ * std::pow(1.0015, e->angleDelta().y()), 0.2, 25.0);
    update();
}

void ViewerWidget::mouseDoubleClickEvent(QMouseEvent*) { resetView(); }

}  // namespace gui
