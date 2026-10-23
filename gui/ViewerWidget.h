// 3D view of a ligand in its protein: cartoon, molecular surface, ligand, residues in contact and interactions.
// Rendered by a software rasterizer (no OpenGL, no WebEngine): smooth at low quality while the mouse moves, refined in the
// background when it stops. The interactions, labels and legend are drawn over the image and hidden where the protein hides them.
#pragma once

#include <QFutureWatcher>
#include <QImage>
#include <QColor>
#include <QPointF>
#include <QSet>
#include <QString>
#include <QTimer>
#include <QWidget>
#include <array>
#include <memory>
#include <vector>

#include "Raster.h"
#include "SceneBuilder.h"
#include "ViewStyle.h"
#include "ligint/ligint.hpp"
#include "ligint/pocket.hpp"

namespace gui {

class ViewerWidget : public QWidget {
    Q_OBJECT
public:
    explicit ViewerWidget(QWidget* parent = nullptr);
    ~ViewerWidget() override;

    // The ligand alone (no protein layers)...
    void setLigand(const ligint::LigandResult& r);
    // ...or the ligand with the protein of the analysis (Result::backbone and Result::proteinAtoms, see Options::protein3d).
    void setStructure(const ligint::Result& result, int ligandIndex = 0);
    void setImportPreview(const ligint::Result& result, int selectedIndex = -1);
    void setPreviewDonorHydrogens(bool enabled);
    bool previewDonorHydrogens() const { return previewDonorHydrogens_; }
    void clear();
    // Overlay the raw pocket volume in the same depth-buffered scene as protein and ligand.
    void setPocketOverlay(const std::vector<ligint::PocketPoint>& points, double spacing, int smoothing=2,
                          const std::vector<std::vector<ligint::PocketPoint>>& groups = {});
    void setPocketDistinctColors(bool enabled);
    QColor pocketDisplayColor(const std::vector<ligint::PocketPoint>& source) const;
    void setPocketAppearance(bool visible, double opacity, double specular);
    bool hasData() const { return hasData_; }
    bool hasProtein() const { return protein_ && !protein_->backbone.empty(); }

    void setVisibleTypes(const QSet<QString>& types);
    void setHighlight(int hitIndex);                 // index in LigandResult::hits, -1 = none
    void setShowLabels(bool on);
    void setShowResidues(bool on);
    void setShowDistances(bool on);
    void setShowLegend(bool on);

    const ViewStyle& style() const { return style_; }
    void setStyle(const ViewStyle& style);           // rebuilds only what changed

    // Draws the scene on an image at full quality (same view as on the screen): "Save image" and the tests.
    QImage renderImage(const QSize& size);

    bool surfaceBusy() const { return surfaceBusy_; }
    std::size_t surfaceTriangles() const { return pieces_.surface ? pieces_.surface->mesh.tri.size() : 0; }   // 0 while it is computed
    std::size_t cartoonTriangles() const { return pieces_.cartoon ? pieces_.cartoon->mesh.tri.size() : 0; }

    QSize sizeHint() const override { return {720, 540}; }

signals:
    void styleChanged();
    void busyChanged(bool busy);
    void importLigandPicked(int index);
    void pocketSurfacePicked(int visibleGroupIndex, bool additive);

public slots:
    void resetView();        // initial orientation, framing the ligand
    void focusLigand();      // frames the ligand and its pocket, keeps the orientation
    void fitProtein();       // frames the whole protein, keeps the orientation

protected:
    void paintEvent(QPaintEvent* e) override;
    void resizeEvent(QResizeEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent*) override { resetView(); }

private:
    using Mat3 = std::array<double, 9>;

    struct Pieces {
        std::shared_ptr<const MeshPiece> cartoon, surface, pocket;
        std::shared_ptr<const AtomsPiece> ligand, residues, water, previewLigands;
    };
    struct Frame {                                   // result of a render
        QImage image;
        std::shared_ptr<raster::Renderer> renderer;  // keeps the depth, used to hide the interactions behind the protein
        QSize size;
        double ratio = 1;
        quint64 key = 0;
        bool highQuality = false;
    };

    raster::Camera makeCamera(const QSizeF& size, double pixelRatio, bool withClip) const;
    static Frame renderScene(const Pieces& pieces, const ViewStyle& style, const raster::Camera& full, const raster::Camera& clipped,
                             const QSize& pixelSize, int supersample, double fogNear, double fogFar);
    Frame renderCurrent(const QSize& logicalSize, double ratio, int supersample) const;
    void drawOverlay(QPainter& p, const QRectF& rect, const Frame& frame) const;
    bool typeVisible(const std::string& type) const;
    void invalidate(bool interactive);
    void startHighQuality();
    void rebuildAll();
    void rebuildCartoon();
    void startSurface();
    void rebuildAtoms();
    void orient();
    void setFrame(const ligint::V3& center, double radius);
    std::size_t triangleCount() const;

    ligint::LigandResult data_;
    std::shared_ptr<ligint::Result> protein_;
    std::vector<ligint::LigandResult> importedLigands_;
    QPoint pressPoint_;
    bool previewDonorHydrogens_ = false;        // backbone and atoms only
    bool hasData_ = false;
    ViewStyle style_;
    Pieces pieces_;
    std::vector<ligint::PocketPoint> pocketPoints_;
    double pocketSpacing_=1.0;
    int pocketSmoothing_=2;
    bool pocketVisible_=true;
    double pocketOpacity_=0.38, pocketSpecular_=0.12;
    bool pocketDistinctColors_=false;
    std::vector<std::vector<ligint::PocketPoint>> pocketGroups_;

    // camera
    Mat3 rot_{1, 0, 0, 0, 1, 0, 0, 0, 1};
    ligint::V3 center_{0, 0, 0};
    double radius_ = 5;                              // radius of what is framed
    double zoom_ = 1.0;
    QPointF panFrac_{0, 0};                          // fraction of the smaller side of the view
    QPoint last_;
    bool interacting_ = false;

    QSet<QString> visible_;
    bool allVisible_ = true;
    int highlight_ = -1;

    // frames
    Frame frame_;
    quint64 stateKey_ = 1;                           // changes with the camera, the style or the data
    QTimer hqTimer_;
    QFutureWatcher<Frame>* hqWatcher_ = nullptr;
    bool hqPending_ = false;
    QFutureWatcher<std::shared_ptr<const MeshPiece>>* surfaceWatcher_ = nullptr;
    quint64 surfaceGeneration_ = 0;
    bool surfaceBusy_ = false;
};

}  // namespace gui
