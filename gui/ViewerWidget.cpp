#include "ViewerWidget.h"
#include "PocketSmoothMesh.h"
#include <map>

#include <QFontMetricsF>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QPolygonF>
#include <QResizeEvent>
#include <QWheelEvent>
#include <QtConcurrent/QtConcurrent>
#include <algorithm>
#include <cmath>
#include <cstring>
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

double dist3(const V3& a, const V3& b) { return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])); }

bool darkBackground(Background b) { return b == Background::Dark; }

}  // namespace

// ------------------------------------------------------------------------------------------------
ViewerWidget::ViewerWidget(QWidget* parent) : QWidget(parent) {
    setMinimumSize(320, 240);
    setMouseTracking(false);
    setFocusPolicy(Qt::StrongFocus);
    setToolTip(tr("Drag: rotate · right/middle button: pan · wheel: zoom · double click: reset"));
    hqTimer_.setSingleShot(true);
    hqTimer_.setInterval(160);
    connect(&hqTimer_, &QTimer::timeout, this, &ViewerWidget::startHighQuality);

    hqWatcher_ = new QFutureWatcher<Frame>(this);
    connect(hqWatcher_, &QFutureWatcher<Frame>::finished, this, [this]() {
        if (hqWatcher_->future().resultCount() == 0) return;
        const Frame f = hqWatcher_->result();
        if (f.key == stateKey_ && f.size == size()) {
            frame_ = f;
            update();
        }
        if (hqPending_) { hqPending_ = false; startHighQuality(); }
    });

    surfaceWatcher_ = new QFutureWatcher<std::shared_ptr<const MeshPiece>>(this);
    connect(surfaceWatcher_, &QFutureWatcher<std::shared_ptr<const MeshPiece>>::finished, this, [this]() {
        // resetting the watcher with an empty future also emits finished(): there is nothing to apply then
        if (surfaceWatcher_->future().resultCount() == 0) return;
        pieces_.surface = surfaceWatcher_->result();
        surfaceBusy_ = false;
        emit busyChanged(false);
        invalidate(false);
    });
}

ViewerWidget::~ViewerWidget() {
    hqWatcher_->waitForFinished();
    surfaceWatcher_->waitForFinished();
}

// ---- data ----------------------------------------------------------------------------------------
void ViewerWidget::setLigand(const ligint::LigandResult& r) {
    data_ = r;
    protein_.reset();
    hasData_ = !data_.atoms.empty();
    highlight_ = -1;
    rebuildAll();
    orient();
    invalidate(false);
}

void ViewerWidget::setStructure(const ligint::Result& result, int ligandIndex) {
    const bool withLigand = ligandIndex >= 0 && ligandIndex < static_cast<int>(result.ligands.size());
    data_ = withLigand ? result.ligands[ligandIndex] : ligint::LigandResult();
    protein_ = std::make_shared<ligint::Result>();
    protein_->backbone = result.backbone;
    protein_->proteinAtoms = result.proteinAtoms;
    hasData_ = !data_.atoms.empty() || !protein_->backbone.empty();
    highlight_ = -1;
    rebuildAll();
    if (withLigand) orient();
    else { rot_ = Mat3{1,0,0,0,1,0,0,0,1}; fitProtein(); }
    invalidate(false);
}

void ViewerWidget::setPocketAppearance(bool visible,double opacity,double specular){
    pocketVisible_=visible;
    pocketOpacity_=std::clamp(opacity,0.05,1.0);
    pocketSpecular_=std::clamp(specular,0.0,1.0);
    invalidate(false);
}

void ViewerWidget::setPocketOverlay(const std::vector<ligint::PocketPoint>& points,
                                    double spacing,int smoothing) {
    pocketPoints_=points;
    pocketSpacing_=spacing;
    pocketSmoothing_=smoothing;
    pieces_.pocket.reset();
    if(!points.empty()){
        const auto triangles=makeSmoothPocketMesh(points,spacing,smoothing);
        auto piece=std::make_shared<MeshPiece>();
        auto& mesh=piece->mesh;
        // Shared position indices allow averaged vertex normals and smooth per-pixel lighting.
        using Key=std::array<long long,3>;
        std::map<Key,int> vertices;
        for(const auto& tri:triangles){
            ligint::Tri indices{};
            for(int j=0;j<3;j++){
                const auto& v=tri.vertices[j];
                Key key={std::llround(v[0]*100000),std::llround(v[1]*100000),
                         std::llround(v[2]*100000)};
                auto it=vertices.find(key);
                if(it==vertices.end()){
                    const int id=static_cast<int>(mesh.pos.size());
                    vertices.emplace(key,id);
                    mesh.pos.push_back({static_cast<float>(v[0]),static_cast<float>(v[1]),
                                        static_cast<float>(v[2])});
                    mesh.nrm.push_back({0,0,0});
                    indices[j]=id;
                }else indices[j]=it->second;
            }
            const auto& a=mesh.pos[indices[0]];
            const auto& b=mesh.pos[indices[1]];
            const auto& d=mesh.pos[indices[2]];
            const float ux=b[0]-a[0],uy=b[1]-a[1],uz=b[2]-a[2];
            const float vx=d[0]-a[0],vy=d[1]-a[1],vz=d[2]-a[2];
            const ligint::V3f normal={uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx};
            for(int id:indices)for(int k=0;k<3;k++)mesh.nrm[id][k]+=normal[k];
            mesh.tri.push_back(indices);
        }
        for(auto& n:mesh.nrm){
            const float norm=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            if(norm>1e-8f)for(auto& x:n)x/=norm;
            else n={0,0,1};
        }
        piece->colors.assign(mesh.pos.size(),0x40c4df);
        pieces_.pocket=std::move(piece);
    }
    invalidate(false);
}

void ViewerWidget::clear() {
    data_ = ligint::LigandResult();
    protein_.reset();
    pieces_ = Pieces();
    pocketPoints_.clear();
    hasData_ = false;
    highlight_ = -1;
    surfaceWatcher_->setFuture(QFuture<std::shared_ptr<const MeshPiece>>());
    if (surfaceBusy_) { surfaceBusy_ = false; emit busyChanged(false); }
    invalidate(false);
}

void ViewerWidget::rebuildAll() {
    rebuildCartoon();
    rebuildAtoms();
    startSurface();
}

void ViewerWidget::rebuildCartoon() {
    pieces_.cartoon = protein_ ? makeCartoonPiece(*protein_, style_) : nullptr;
}

void ViewerWidget::rebuildAtoms() {
    pieces_.ligand = makeLigandPiece(data_, style_);
    pieces_.residues = makeResiduePiece(data_, style_);
    pieces_.water = makeWaterPiece(data_, style_);
}

void ViewerWidget::startSurface() {
    if (style_.surface == SurfaceMode::Off || !protein_ || !hasData_) {
        surfaceWatcher_->setFuture(QFuture<std::shared_ptr<const MeshPiece>>());
        pieces_.surface.reset();
        if (surfaceBusy_) { surfaceBusy_ = false; emit busyChanged(false); }
        return;
    }
    if (!surfaceBusy_) { surfaceBusy_ = true; emit busyChanged(true); }
    const std::shared_ptr<ligint::Result> prot = protein_;
    const ligint::LigandResult lig = data_;
    const ViewStyle st = style_;
    surfaceWatcher_->setFuture(QtConcurrent::run([prot, lig, st]() { return makeSurfacePiece(*prot, lig, st); }));
}

void ViewerWidget::setStyle(const ViewStyle& n) {
    const ViewStyle o = style_;
    style_ = n;
    if (n.protein != o.protein || n.proteinColor != o.proteinColor) rebuildCartoon();
    if (n.ligand != o.ligand || n.showLigand != o.showLigand || n.residues != o.residues || n.showWaters != o.showWaters) rebuildAtoms();
    if (n.surface != o.surface || n.surfaceType != o.surfaceType || n.surfaceColor != o.surfaceColor ||
        (n.surfaceColor == SurfaceColor::Protein && n.proteinColor != o.proteinColor))
        startSurface();
    invalidate(false);
    emit styleChanged();
}

void ViewerWidget::setShowLabels(bool on) { ViewStyle s = style_; s.showLabels = on; setStyle(s); }
void ViewerWidget::setShowResidues(bool on) { ViewStyle s = style_; s.residues = on ? ResidueStyle::Sticks : ResidueStyle::Off; setStyle(s); }
void ViewerWidget::setShowDistances(bool on) { ViewStyle s = style_; s.showDistances = on; setStyle(s); }
void ViewerWidget::setShowLegend(bool on) { ViewStyle s = style_; s.showLegend = on; setStyle(s); }

void ViewerWidget::setVisibleTypes(const QSet<QString>& types) {
    visible_ = types;
    allVisible_ = false;
    invalidate(false);
}
void ViewerWidget::setHighlight(int i) { highlight_ = i; invalidate(false); }

bool ViewerWidget::typeVisible(const std::string& type) const {
    return allVisible_ || visible_.contains(QString::fromStdString(type));
}

std::size_t ViewerWidget::triangleCount() const {
    std::size_t n = 0;
    if (pieces_.cartoon) n += pieces_.cartoon->mesh.tri.size();
    if (pieces_.surface) n += pieces_.surface->mesh.tri.size();
    if (pieces_.pocket) n += pieces_.pocket->mesh.tri.size();
    return n;
}

// ---- camera --------------------------------------------------------------------------------------
void ViewerWidget::setFrame(const V3& center, double radius) {
    center_ = center;
    radius_ = std::max(radius, 4.0);
    zoom_ = 1.0;
    panFrac_ = QPointF(0, 0);
}

void ViewerWidget::orient() {
    std::vector<V3> pts;
    for (auto& a : data_.atoms) pts.push_back(a.pos);
    for (auto& h : data_.hits) { pts.push_back(h.ligandPos); pts.push_back(h.proteinPos); }
    if (pts.empty()) return;
    focusLigand();
    rot_ = principalAxes(pts);      // the ligand (with its contacts) is shown with the maximum extent on screen
}

void ViewerWidget::focusLigand() {
    std::vector<V3> pts;
    for (auto& a : data_.atoms) pts.push_back(a.pos);
    for (auto& h : data_.hits) { pts.push_back(h.ligandPos); pts.push_back(h.proteinPos); }
    if (pts.empty()) return;
    V3 lo = pts[0], hi = pts[0];
    for (auto& p : pts) for (int i = 0; i < 3; i++) { lo[i] = std::min(lo[i], p[i]); hi[i] = std::max(hi[i], p[i]); }
    V3 c{0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
    double radius = 0;
    for (auto& p : pts) radius = std::max(radius, dist3(p, c));
    const double hitRadius = std::max(radius, 4.0);
    for (auto& a : data_.contactAtoms) radius = std::max(radius, dist3(a.pos, c));
    radius = std::min(radius, 1.35 * hitRadius);       // the farthest residues may fall outside the frame
    setFrame(c, radius);
    invalidate(false);
}

void ViewerWidget::fitProtein() {
    if (!hasProtein()) { focusLigand(); return; }
    V3 c{0, 0, 0};
    for (auto& r : protein_->backbone) for (int i = 0; i < 3; i++) c[i] += r.ca[i] / protein_->backbone.size();
    double radius = 0;
    for (auto& r : protein_->backbone) radius = std::max(radius, dist3(r.ca, c));
    setFrame(c, radius + 4.0);
    invalidate(false);
}

void ViewerWidget::resetView() {
    if (!data_.atoms.empty()) orient();
    else if (hasProtein()) { rot_ = Mat3{1,0,0,0,1,0,0,0,1}; fitProtein(); }
    invalidate(false);
}

raster::Camera ViewerWidget::makeCamera(const QSizeF& size, double ratio, bool withClip) const {
    raster::Camera c;
    c.rot = rot_;
    c.center = center_;
    const double m = std::min(size.width(), size.height());
    c.width = std::max(1, static_cast<int>(std::lround(size.width() * ratio)));
    c.height = std::max(1, static_cast<int>(std::lround(size.height() * ratio)));
    c.scale = 0.47 * m / radius_ * zoom_ * ratio;
    c.panX = panFrac_.x() * m * ratio;
    c.panY = panFrac_.y() * m * ratio;
    c.distance = style_.perspective ? 5.0 * radius_ : 0.0;
    c.clipZ = 1e9;
    if (withClip && style_.clipMargin > 0 && (pieces_.cartoon || pieces_.surface)) {
        // the protein in front of the ligand is cut away, so that the ligand is not hidden by the protein around it
        double zLig = -1e9;
        auto z = [&](const V3& p) {
            const double dx = p[0] - center_[0], dy = p[1] - center_[1], dz = p[2] - center_[2];
            return rot_[6] * dx + rot_[7] * dy + rot_[8] * dz;
        };
        for (auto& a : data_.atoms) zLig = std::max(zLig, z(a.pos));
        for (auto& h : data_.hits) { zLig = std::max(zLig, z(h.ligandPos)); zLig = std::max(zLig, z(h.proteinPos)); }
        if (zLig > -1e8) c.clipZ = zLig + style_.clipMargin;
    }
    return c;
}

// ---- rendering -----------------------------------------------------------------------------------
ViewerWidget::Frame ViewerWidget::renderScene(const Pieces& pieces, const ViewStyle& style, const raster::Camera& full, const raster::Camera& clipped,
                                              const QSize& pixelSize, int supersample, double fogNear, double fogFar) {
    auto ren = std::make_shared<raster::Renderer>(pixelSize.width(), pixelSize.height(), supersample);
    raster::Lighting light;
    light.background = packColor(backgroundColor(style.background));
    light.fogStrength = style.depthCue ? (darkBackground(style.background) ? 0.65f : 0.5f) : 0.0f;
    light.fogNear = fogNear;
    light.fogFar = fogFar;
    ren->setLighting(light);

    ren->setCamera(clipped);                       // the protein layers are cut in front of the ligand...
    if (pieces.cartoon) {
        raster::MeshView mv;
        mv.mesh = &pieces.cartoon->mesh;
        mv.colors = &pieces.cartoon->colors;
        ren->drawMesh(mv);
    }
    if (pieces.surface) {
        raster::MeshView mv;
        mv.mesh = &pieces.surface->mesh;
        mv.colors = &pieces.surface->colors;
        mv.alpha = static_cast<float>(std::clamp(style.surfaceOpacity, 0.05, 1.0));
        mv.cullBack = mv.alpha < 0.999f;
        mv.specular = 0.18f;
        ren->drawMesh(mv);
    }
    // The pocket shares the protein's world coordinates, camera and z-buffer.
    // Render with the full camera to preserve the selected cavity.
    ren->setCamera(full);
    if(pieces.pocket && style.pocketVisible){
        raster::MeshView pv;
        pv.mesh=&pieces.pocket->mesh;
        pv.colors=&pieces.pocket->colors;
        pv.alpha=static_cast<float>(style.pocketOpacity);
        // Display only the front shell: avoids double blending of entry/exit walls.
        pv.cullBack=pv.alpha<0.999f;
        pv.specular=static_cast<float>(style.pocketSpecular);
        ren->drawMesh(pv);
    }
    // ...ligand and residues remain in the same depth-buffered composition.

    for (const auto& piece : {pieces.residues, pieces.ligand, pieces.water}) {
        if (!piece) continue;
        if (!piece->sticks.empty()) {
            raster::MeshView mv;
            mv.mesh = &piece->sticks;
            mv.colors = &piece->stickColors;
            mv.specular = 0.3f;
            ren->drawMesh(mv);
        }
        ren->drawSpheres(piece->spheres);
    }
    ren->render();

    Frame f;
    const std::vector<std::uint32_t> px = ren->resolve();
    f.image = QImage(pixelSize, QImage::Format_ARGB32);
    for (int y = 0; y < pixelSize.height(); y++)
        std::memcpy(f.image.scanLine(y), px.data() + static_cast<std::size_t>(y) * pixelSize.width(), static_cast<std::size_t>(pixelSize.width()) * 4);
    f.renderer = ren;
    return f;
}

ViewerWidget::Frame ViewerWidget::renderCurrent(const QSize& logicalSize, double ratio, int supersample) const {
    const raster::Camera full = makeCamera(logicalSize, ratio, false);
    const raster::Camera clipped = makeCamera(logicalSize, ratio, true);
    const double base = style_.perspective ? 5.0 * radius_ : 1000.0;
    ViewStyle renderStyle=style_;
    renderStyle.pocketVisible=pocketVisible_;
    renderStyle.pocketOpacity=pocketOpacity_;
    renderStyle.pocketSpecular=pocketSpecular_;
    Frame f = renderScene(pieces_, renderStyle, full, clipped, QSize(full.width, full.height), supersample, base - radius_, base + 2.5 * radius_);
    f.size = logicalSize;
    f.ratio = ratio;
    f.key = stateKey_;
    return f;
}

void ViewerWidget::invalidate(bool interactive) {
    stateKey_++;
    interacting_ = interactive;
    if (interactive) hqTimer_.stop();
    update();
}

void ViewerWidget::startHighQuality() {
    if (!hasData_ || size().isEmpty()) return;
    if (frame_.key == stateKey_ && frame_.highQuality && frame_.size == size()) return;
    if (hqWatcher_->isRunning()) { hqPending_ = true; return; }
    const Pieces pieces = pieces_;
    ViewStyle style = style_;
    style.pocketVisible=pocketVisible_;
    style.pocketOpacity=pocketOpacity_;
    style.pocketSpecular=pocketSpecular_;
    const double ratio = devicePixelRatioF();
    const QSize sz = size();
    const raster::Camera full = makeCamera(sz, ratio, false);
    const raster::Camera clipped = makeCamera(sz, ratio, true);
    const double base = style_.perspective ? 5.0 * radius_ : 1000.0;
    const double fogNear = base - radius_, fogFar = base + 2.5 * radius_;
    const quint64 key = stateKey_;
    const int ss = ratio >= 1.9 ? 1 : 2;
    hqWatcher_->setFuture(QtConcurrent::run([=]() {
        Frame f = renderScene(pieces, style, full, clipped, QSize(full.width, full.height), ss, fogNear, fogFar);
        f.size = sz;
        f.ratio = ratio;
        f.key = key;
        f.highQuality = true;
        return f;
    }));
}

QImage ViewerWidget::renderImage(const QSize& size) {
    if (size.isEmpty()) return QImage();
    Frame f = renderCurrent(size, 1.0, 3);
    QImage img = f.image.convertToFormat(QImage::Format_ARGB32);
    if (!hasData_) {
        img.fill(backgroundColor(style_.background));
        return img;
    }
    QPainter p(&img);
    drawOverlay(p, QRectF(QPointF(0, 0), QSizeF(size)), f);
    p.end();
    return img;
}

// ---- painting ------------------------------------------------------------------------------------
void ViewerWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    const QColor bg = backgroundColor(style_.background);
    if (!hasData_) {
        p.fillRect(rect(), bg);
        p.setPen(darkBackground(style_.background) ? QColor("#adb5bd") : QColor("#868e96"));
        p.drawText(rect(), Qt::AlignCenter, tr("Open a structure file, choose a ligand and press \u00abAnalyze\u00bb"));
        return;
    }
    const double ratio = devicePixelRatioF();
    if (frame_.key != stateKey_ || frame_.size != size()) {
        const double lowFactor = (interacting_ && triangleCount() > 500000) ? 0.5 : 1.0;
        frame_ = renderCurrent(size(), ratio * lowFactor, 1);
        if (!interacting_) hqTimer_.start();
    }
    p.setRenderHint(QPainter::SmoothPixmapTransform, true);
    p.drawImage(rect(), frame_.image);
    p.setRenderHint(QPainter::Antialiasing, true);
    drawOverlay(p, QRectF(rect()), frame_);
    if (surfaceBusy_) {
        p.setPen(darkBackground(style_.background) ? QColor("#ced4da") : QColor("#495057"));
        p.drawText(QPointF(rect().right() - 150, rect().top() + 18), tr("Computing the surface…"));
    }
}

void ViewerWidget::resizeEvent(QResizeEvent* e) {
    QWidget::resizeEvent(e);
    hqTimer_.start();
}

void ViewerWidget::drawOverlay(QPainter& p, const QRectF& rect, const Frame& frame) const {
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setRenderHint(QPainter::TextAntialiasing, true);
    const bool dark = darkBackground(style_.background);
    const QColor fg = dark ? QColor("#e9ecef") : QColor("#212529");
    const QColor boxFill = dark ? QColor(30, 36, 46, 215) : QColor(255, 255, 255, 215);
    const QColor textColor = dark ? QColor("#dee2e6") : QColor("#343a40");

    const raster::Camera cam = makeCamera(rect.size(), 1.0, false);
    const double s = cam.scale;     // pixels per Angstrom at the centre
    auto project = [&](const V3& v) { return cam.project(v); };
    auto toPt = [&](const raster::Camera::Point& q) { return QPointF(rect.left() + q.x, rect.top() + q.y); };
    // hidden by the protein? (compares with the depth of the image)
    auto visible = [&](const raster::Camera::Point& q) {
        if (!frame.renderer || !q.ok) return q.ok;
        const float d = frame.renderer->depthAt(q.x * frame.ratio, q.y * frame.ratio);
        return q.w <= d + 1.3;
    };

    QFont small = p.font();
    small.setPointSizeF(std::max(7.5, small.pointSizeF() * 0.9));
    p.setFont(small);
    QFontMetricsF fm(small);

    std::map<std::string, QPointF> labelPos;
    std::vector<std::pair<std::string, int>> typesUsed;

    // dashed line, only the parts that are not behind the protein
    auto strokeVisible = [&](const V3& a3, const V3& b3, const QPen& pen) {
        const auto pa = project(a3), pb = project(b3);
        const double len = std::hypot(pb.x - pa.x, pb.y - pa.y);
        const int n = std::clamp(static_cast<int>(len / 3.0), 2, 240);
        QPolygonF run;
        p.setPen(pen);
        auto flush = [&]() {
            if (run.size() >= 2) p.drawPolyline(run);
            run.clear();
        };
        for (int i = 0; i <= n; i++) {
            const double t = static_cast<double>(i) / n;
            const V3 w{a3[0] + (b3[0] - a3[0]) * t, a3[1] + (b3[1] - a3[1]) * t, a3[2] + (b3[2] - a3[2]) * t};
            const auto q = project(w);
            if (visible(q)) run << toPt(q); else flush();
        }
        flush();
    };

    for (size_t i = 0; i < data_.hits.size(); i++) {
        const auto& h = data_.hits[i];
        if (!typeVisible(h.type)) continue;
        const QString typeName = QString::fromStdString(h.type);
        const TypeStyle& st = styleForType(typeName);
        const bool hl = static_cast<int>(i) == highlight_;
        const double widthScale = std::clamp(s / 22.0, 0.8, 1.6);
        auto drawSeg = [&](const V3& a, const V3& b) {
            if (hl) strokeVisible(a, b, QPen(dark ? QColor(255, 255, 255, 120) : QColor(255, 255, 255, 230), (st.width + 5) * widthScale, Qt::SolidLine, Qt::RoundCap));
            strokeVisible(a, b, QPen(st.color, (hl ? st.width + 1.5 : st.width) * widthScale, st.pen, Qt::RoundCap));
        };
        if (h.hasWater) {
            drawSeg(h.ligandPos, h.waterPos);
            drawSeg(h.waterPos, h.proteinPos);
        } else {
            drawSeg(h.ligandPos, h.proteinPos);
        }
        const auto qa = project(h.ligandPos), qb = project(h.proteinPos);
        if (h.type.rfind("pi-stacking", 0) == 0 || h.type == ligint::type::CationPi) {
            p.setPen(QPen(st.color, 1.6));
            p.setBrush(Qt::NoBrush);
            if (visible(qa)) p.drawEllipse(toPt(qa), 6, 6);
            if (visible(qb)) p.drawEllipse(toPt(qb), 6, 6);
        }
        if (!labelPos.count(h.residue) && visible(qb)) labelPos[h.residue] = toPt(qb);
        const auto mid = project({(h.ligandPos[0] + h.proteinPos[0]) / 2, (h.ligandPos[1] + h.proteinPos[1]) / 2, (h.ligandPos[2] + h.proteinPos[2]) / 2});
        if (style_.showDistances && visible(mid) && (hl || (h.type != ligint::type::Hydrophobic && h.type != ligint::type::BackboneApolar))) {
            const QPointF at = toPt(mid);
            const QString t = QString::number(h.distance, 'f', 2);
            QRectF r(at.x() - fm.horizontalAdvance(t) / 2 - 2, at.y() - fm.height() / 2, fm.horizontalAdvance(t) + 4, fm.height());
            p.setPen(Qt::NoPen);
            p.setBrush(boxFill);
            p.drawRoundedRect(r, 2, 2);
            p.setPen(dark ? st.color.lighter(130) : st.color.darker(130));
            p.drawText(r, Qt::AlignCenter, t);
        }
        bool seen = false;
        for (auto& u : typesUsed) if (u.first == h.type) { u.second++; seen = true; }
        if (!seen) typesUsed.push_back({h.type, 1});
    }

    if (style_.showLabels) {
        for (auto& kv : labelPos) {
            const QString t = QString::fromStdString(kv.first);
            const QPointF at = kv.second + QPointF(7, -7);
            QRectF r(at.x(), at.y() - fm.ascent(), fm.horizontalAdvance(t) + 6, fm.height() + 2);
            p.setPen(Qt::NoPen);
            p.setBrush(boxFill);
            p.drawRoundedRect(r, 3, 3);
            p.setPen(fg);
            p.drawText(r, Qt::AlignCenter, t);
        }
    }

    if (style_.showLegend && !typesUsed.empty()) {
        double w = 0;
        for (auto& u : typesUsed) w = std::max(w, fm.horizontalAdvance(QString::fromStdString(u.first) + QString(" (%1)").arg(u.second)));
        const double lineH = fm.height() + 3;
        QRectF box(rect.left() + 8, rect.top() + 8, w + 52, lineH * typesUsed.size() + 8);
        p.setPen(QPen(dark ? QColor("#495057") : QColor("#dee2e6"), 1));
        p.setBrush(boxFill);
        p.drawRoundedRect(box, 4, 4);
        double y = box.top() + 4;
        for (auto& u : typesUsed) {
            const QString legendType = QString::fromStdString(u.first);
            const TypeStyle& st = styleForType(legendType);
            p.setPen(QPen(st.color, st.width, st.pen, Qt::RoundCap));
            p.drawLine(QPointF(box.left() + 6, y + lineH / 2), QPointF(box.left() + 34, y + lineH / 2));
            p.setPen(textColor);
            p.drawText(QPointF(box.left() + 42, y + fm.ascent() + 1), QString::fromStdString(u.first) + QString(" (%1)").arg(u.second));
            y += lineH;
        }
    }
    p.setPen(dark ? QColor("#adb5bd") : QColor("#495057"));
    p.drawText(QPointF(rect.left() + 10, rect.bottom() - 10), QString::fromStdString(data_.id));
}

// ---- mouse ---------------------------------------------------------------------------------------
void ViewerWidget::mousePressEvent(QMouseEvent* e) { last_ = eventPos(e).toPoint(); }

void ViewerWidget::mouseMoveEvent(QMouseEvent* e) {
    const QPoint cur = eventPos(e).toPoint();
    const QPoint d = cur - last_;
    last_ = cur;
    if (!hasData_) return;
    if (e->buttons() & Qt::LeftButton) {
        rot_ = mul(rotX(d.y() * 0.01), mul(rotY(d.x() * 0.01), rot_));
    } else if (e->buttons() & (Qt::RightButton | Qt::MiddleButton)) {
        const double m = std::max(1, std::min(width(), height()));
        panFrac_ += QPointF(d.x() / m, d.y() / m);
    } else {
        return;
    }
    invalidate(true);
}

void ViewerWidget::mouseReleaseEvent(QMouseEvent*) {
    if (interacting_) {
        interacting_ = false;
        stateKey_++;          // a fresh frame at the normal size, then the high quality one
        update();
    }
}

void ViewerWidget::wheelEvent(QWheelEvent* e) {
    if (!hasData_) return;
    zoom_ = std::clamp(zoom_ * std::pow(1.0015, e->angleDelta().y()), 0.2, 25.0);
    invalidate(false);
}

}  // namespace gui
