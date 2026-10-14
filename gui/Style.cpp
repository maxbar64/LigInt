#include "Style.h"

#include "ligint/ligint.hpp"

#include <QSet>

namespace gui {

const QVector<TypeStyle>& typeStyles() {
    static const QVector<TypeStyle> styles = {
        {QString(ligint::type::HydrogenBond),             QColor("#1c7ed6"), Qt::DashLine,       2.2},
        {QString(ligint::type::SaltBridge),                 QColor("#e8590c"), Qt::DashDotLine,    2.4},
        {QString(ligint::type::ElectrostaticRepulsion),  QColor("#c92a2a"), Qt::DashDotDotLine, 1.6},
        {QString(ligint::type::Hydrophobic),           QColor("#74b816"), Qt::DotLine,        1.2},
        {QString(ligint::type::BackboneApolar),    QColor("#a9b79c"), Qt::DotLine,        1.0},
        {QString(ligint::type::PiStackingParallel),     QColor("#2f9e44"), Qt::DashLine,       2.6},
        {QString(ligint::type::PiStackingT),           QColor("#099268"), Qt::DashLine,       2.6},
        {QString(ligint::type::CationPi),                     QColor("#f08c00"), Qt::DashLine,       2.4},
        {QString(ligint::type::HalogenBond),             QColor("#0ca678"), Qt::DashLine,       2.2},
        {QString(ligint::type::MetalCoordination),       QColor("#7048e8"), Qt::SolidLine,      2.4},
        {QString(ligint::type::WaterBridge),                 QColor("#4dabf7"), Qt::DashLine,       1.8},
        {QString(ligint::type::StericClash),                 QColor("#fa5252"), Qt::SolidLine,      2.6},
        {QString(ligint::type::PossibleCovalent),     QColor("#212529"), Qt::SolidLine,      3.0},
    };
    return styles;
}

const TypeStyle& styleForType(const QString& type) {
    for (const auto& s : typeStyles())
        if (s.type == type) return s;
    static const TypeStyle fallback{"?", QColor("#868e96"), Qt::SolidLine, 1.5};
    return fallback;
}

QColor fpClassColor(int c) {
    static const char* col[] = {"#868e96", "#2f9e44", "#94d82d", "#e03131", "#1c7ed6", "#1971c2", "#c92a2a",
                                "#f08c00", "#0ca678", "#7048e8", "#4dabf7"};
    return QColor(c >= 0 && c < static_cast<int>(sizeof col / sizeof col[0]) ? col[c] : "#868e96");
}

bool isMetalElement(const QString& e) {
    static const QSet<QString> m = {"Zn", "Mg", "Ca", "Fe", "Mn", "Cu", "Ni", "Co", "Cd", "Hg", "Na", "K"};
    return m.contains(e);
}

QColor elementColor(const QString& e, bool ligand) {
    if (e == "C") return ligand ? QColor("#d6336c") : QColor("#40a060");
    if (e == "N") return QColor("#3b5bdb");
    if (e == "O") return QColor("#e03131");
    if (e == "S") return QColor("#f59f00");
    if (e == "P") return QColor("#e8590c");
    if (e == "F") return QColor("#82c91e");
    if (e == "Cl") return QColor("#37b24d");
    if (e == "Br") return QColor("#a52a2a");
    if (e == "I") return QColor("#6f2dbd");
    if (e == "H") return QColor("#e9ecef");
    if (e == "Zn") return QColor("#868e96");
    if (e == "Mg") return QColor("#51cf66");
    if (e == "Ca") return QColor("#4dabf7");
    if (e == "Fe") return QColor("#d9480f");
    return QColor("#adb5bd");
}

double elementRadius(const QString& e) {
    if (e == "H") return 0.6;
    if (isMetalElement(e)) return 1.7;
    if (e == "S" || e == "P" || e == "Cl") return 1.25;
    if (e == "Br" || e == "I") return 1.4;
    return 1.0;
}

}  // namespace gui
