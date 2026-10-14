// Colors and styles shared between the graphic view, the tables and the legend.
#pragma once

#include <QColor>
#include <QString>
#include <QVector>
#include <Qt>

namespace gui {

struct TypeStyle {
    QString type;        // string returned by the library (one of ligint::type::*)
    QColor color;
    Qt::PenStyle pen;
    double width;        // thickness in pixels
};

// All known interaction types, in display order.
const QVector<TypeStyle>& typeStyles();
const TypeStyle& styleForType(const QString& type);

QColor elementColor(const QString& element, bool ligand);
double elementRadius(const QString& element);          // relative radius (carbon = 1)
bool isMetalElement(const QString& element);

}  // namespace gui
