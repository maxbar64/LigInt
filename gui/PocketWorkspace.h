#pragma once
#include <QWidget>
#include <QFutureWatcher>
#include <QPoint>
#include "ligint/pocket.hpp"
#include "ligint/pocket_groups.hpp"
class QComboBox;class QDoubleSpinBox;class QSpinBox;class QListWidget;class QLabel;class QPushButton;
namespace gui {
struct PocketGuiOutcome {
    ligint::PocketSearchResult result;
    std::vector<ligint::PocketGroup> groups;
    std::vector<ligint::PocketAtom> atoms;
    QString error;
};
class PocketWorkspace : public QWidget {
    Q_OBJECT
public:
    explicit PocketWorkspace(QWidget* parent=nullptr);
    ~PocketWorkspace() override;
    void setStructurePath(const QString& path);
private:
    QWidget* canvas_=nullptr;
    QComboBox *renderMode_=nullptr,*ligand_=nullptr;
    QDoubleSpinBox *seedSpacing_=nullptr,*seedProminence_=nullptr,*mergeDistance_=nullptr,*mergeDrop_=nullptr;
    QSpinBox *maxSeeds_=nullptr,*noiseNeighbors_=nullptr,*noiseComponent_=nullptr;
    std::vector<ligint::V3> ligandCoords_;
    QComboBox *morph_=nullptr,*display_=nullptr,*chain_=nullptr,*palette_=nullptr;
    QDoubleSpinBox *spacing_=nullptr,*buried_=nullptr,*ray_=nullptr,*surface_=nullptr,*probe_=nullptr;
    QSpinBox *radius_=nullptr,*minPts_=nullptr,*opacity_=nullptr,*pointSize_=nullptr;
    QListWidget* entries_=nullptr;
    QLabel* status_=nullptr;
    QPushButton* run_=nullptr;
    QString path_;
    ligint::PocketSearchResult result_;
    std::vector<ligint::PocketAtom> atoms_;
    std::vector<ligint::PocketGroup> groups_;
    QFutureWatcher<PocketGuiOutcome> watcher_;
    std::vector<std::size_t> visibleMembers() const;
    friend class PocketCanvas;
};
} // namespace gui
