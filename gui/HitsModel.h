// Modello della tabella delle interazioni e filtro per tipo.
#pragma once

#include <QAbstractTableModel>
#include <QSet>
#include <QSortFilterProxyModel>
#include <QString>
#include <vector>

#include "ligint/ligint.hpp"

namespace gui {

class HitsModel : public QAbstractTableModel {
    Q_OBJECT
public:
    enum Column { Type = 0, Favorable, LigandAtom, ProteinAtom, Distance, Note, ColumnCount };
    enum Role { SortRole = Qt::UserRole + 1, HitIndexRole };

    explicit HitsModel(QObject* parent = nullptr);
    void setHits(std::vector<ligint::HitRecord> hits);
    const std::vector<ligint::HitRecord>& hits() const { return hits_; }

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    int columnCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role) const override;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const override;

private:
    std::vector<ligint::HitRecord> hits_;
};

class HitsFilterProxy : public QSortFilterProxyModel {
    Q_OBJECT
public:
    explicit HitsFilterProxy(QObject* parent = nullptr);
    void setVisibleTypes(const QSet<QString>& types);   // insieme vuoto = nessun tipo visibile
    const QSet<QString>& visibleTypes() const { return visible_; }

protected:
    bool filterAcceptsRow(int row, const QModelIndex& parent) const override;

private:
    QSet<QString> visible_;
};

}  // namespace gui
