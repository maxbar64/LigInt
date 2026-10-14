// 3D graphic view (orthographic projection, software rendering with QPainter).
// Shows the ligand, the residues in contact, the ions, the waters and the interactions found.
#pragma once

#include <QImage>
#include <QPointF>
#include <QSet>
#include <QString>
#include <QWidget>
#include <array>
#include <vector>

#include "ligint/ligint.hpp"

namespace gui {

class ViewerWidget : public QWidget {
    Q_OBJECT
public:
    explicit ViewerWidget(QWidget* parent = nullptr);

    void setLigand(const ligint::LigandResult& r);   // copies the data and orients the view
    void clear();
    bool hasData() const { return hasData_; }

    void setVisibleTypes(const QSet<QString>& types);
    void setHighlight(int hitIndex);                 // index in LigandResult::hits, -1 = none
    void setShowLabels(bool on);
    void setShowResidues(bool on);
    void setShowDistances(bool on);
    void setShowLegend(bool on);

    // Draws the scene on an image (same rendering as the screen): used for "Save image" and in the tests.
    QImage renderImage(const QSize& size);

    QSize sizeHint() const override { return {640, 480}; }

public slots:
    void resetView();

protected:
    void paintEvent(QPaintEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;

private:
    using Mat3 = std::array<double, 9>;   // rotation matrix, row-major
    struct Proj { QPointF pt; double z; };

    void render(QPainter& p, const QRectF& rect) const;
    Proj project(const ligint::V3& v, const QRectF& rect) const;
    double fitScale(const QRectF& rect) const;
    void orient();
    bool typeVisible(const std::string& type) const;

    ligint::LigandResult data_;
    bool hasData_ = false;
    ligint::V3 center_{0, 0, 0};
    double radius_ = 5;

    Mat3 rot_{1, 0, 0, 0, 1, 0, 0, 0, 1};
    double zoom_ = 1.0;
    QPointF pan_{0, 0};
    QPoint last_;

    QSet<QString> visible_;
    bool allVisible_ = true;
    int highlight_ = -1;
    bool showLabels_ = true, showResidues_ = true, showDistances_ = true, showLegend_ = true;
};

}  // namespace gui
