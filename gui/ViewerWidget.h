// Vista grafica 3D (proiezione ortografica, rendering software con QPainter).
// Mostra il ligando, i residui a contatto, gli ioni, le acque e le interazioni trovate.
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

    void setLigand(const ligint::LigandResult& r);   // copia i dati e orienta la vista
    void clear();
    bool hasData() const { return hasData_; }

    void setVisibleTypes(const QSet<QString>& types);
    void setHighlight(int hitIndex);                 // indice in LigandResult::hits, -1 = nessuno
    void setShowLabels(bool on);
    void setShowResidues(bool on);
    void setShowDistances(bool on);
    void setShowLegend(bool on);

    // Disegna la scena su un'immagine (stesso rendering dello schermo): usata per "Salva immagine" e nei test.
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
    using Mat3 = std::array<double, 9>;   // matrice di rotazione per righe
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
