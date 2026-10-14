// Colori e stili condivisi tra vista grafica, tabelle e legenda.
#pragma once

#include <QColor>
#include <QString>
#include <QVector>
#include <Qt>

namespace gui {

struct TypeStyle {
    QString type;        // stringa restituita dalla libreria (es. "Legame a idrogeno")
    QColor color;
    Qt::PenStyle pen;
    double width;        // spessore in pixel
};

// Tutti i tipi di interazione noti, nell'ordine di visualizzazione.
const QVector<TypeStyle>& typeStyles();
const TypeStyle& styleForType(const QString& type);

QColor elementColor(const QString& element, bool ligand);
double elementRadius(const QString& element);          // raggio relativo (carbonio = 1)
bool isMetalElement(const QString& element);

}  // namespace gui
