#include "HitsModel.h"

#include <QBrush>

#include "Style.h"

namespace gui {

HitsModel::HitsModel(QObject* parent) : QAbstractTableModel(parent) {}

void HitsModel::setHits(std::vector<ligint::HitRecord> hits) {
    beginResetModel();
    hits_ = std::move(hits);
    endResetModel();
}

int HitsModel::rowCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : static_cast<int>(hits_.size()); }
int HitsModel::columnCount(const QModelIndex& parent) const { return parent.isValid() ? 0 : ColumnCount; }

QVariant HitsModel::data(const QModelIndex& index, int role) const {
    if (!index.isValid() || index.row() < 0 || index.row() >= static_cast<int>(hits_.size())) return {};
    const ligint::HitRecord& h = hits_[index.row()];
    if (role == HitIndexRole) return index.row();
    if (role == SortRole) {
        if (index.column() == Distance) return h.distance;
        if (index.column() == Favorable) return h.favorable ? 1 : 0;
    }
    if (role == Qt::ForegroundRole && index.column() == Type) return QBrush(styleForType(QString::fromStdString(h.type)).color.darker(115));
    if (role == Qt::TextAlignmentRole && index.column() == Distance) return int(Qt::AlignRight | Qt::AlignVCenter);
    if (role != Qt::DisplayRole && role != SortRole) return {};
    switch (index.column()) {
        case Type: return QString::fromStdString(h.type);
        case Favorable: return h.favorable ? QStringLiteral("+") : QStringLiteral("-");
        case LigandAtom: return QString::fromStdString(h.ligandAtom);
        case ProteinAtom: return QString::fromStdString(h.proteinAtom);
        case Distance: return role == SortRole ? QVariant(h.distance) : QVariant(QString::number(h.distance, 'f', 2));
        case Note: return QString::fromStdString(h.note);
    }
    return {};
}

QVariant HitsModel::headerData(int section, Qt::Orientation orientation, int role) const {
    if (orientation != Qt::Horizontal || role != Qt::DisplayRole) return {};
    static const char* names[] = {"Type", "Fav", "Ligand atom", "Protein atom", "d (Å)", "Notes"};
    return section >= 0 && section < ColumnCount ? QString::fromUtf8(names[section]) : QVariant();
}

HitsFilterProxy::HitsFilterProxy(QObject* parent) : QSortFilterProxyModel(parent) {
    setSortRole(HitsModel::SortRole);
    for (const auto& s : typeStyles()) visible_.insert(s.type);
}

void HitsFilterProxy::setVisibleTypes(const QSet<QString>& types) {
    visible_ = types;
    invalidateFilter();
}

bool HitsFilterProxy::filterAcceptsRow(int row, const QModelIndex& parent) const {
    const QModelIndex idx = sourceModel()->index(row, HitsModel::Type, parent);
    return visible_.contains(sourceModel()->data(idx, Qt::DisplayRole).toString());
}

}  // namespace gui
