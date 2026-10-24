#pragma once
#include <QWidget>
#include <QFutureWatcher>
#include <QPoint>
#include <QStringList>
#include <QList>
#include "ligint/pocket.hpp"
#include "ligint/pocket_groups.hpp"
#include "ligint/pocket_alpha.hpp"
class QCheckBox;class QComboBox;class QDoubleSpinBox;class QSpinBox;class QListWidget;class QLabel;class QPushButton;
namespace gui {
struct PocketGuiOutcome {
    ligint::PocketSearchResult result;
    std::vector<ligint::PocketGroup> groups;
    std::vector<ligint::PocketAtom> atoms;
    std::vector<ligint::AlphaPocket> alphaPockets;
    ligint::AlphaMergeAudit alphaAudit;
    bool alphaMode = false;
    QString error;
};
class PocketWorkspace : public QWidget {
    Q_OBJECT
public:
    explicit PocketWorkspace(QWidget* parent=nullptr);
    ~PocketWorkspace() override;
    void setStructurePath(const QString& path);
    QStringList pocketChoices() const;
    std::vector<ligint::PocketPoint> pocketColorSample(int row) const;
    int selectedPocketRow() const;
    int rowForVisiblePocket(int visibleIndex) const;
    void selectPocketRow(int row);
    QList<int> selectedPocketRows() const;
    void setSelectedPocketRows(const QList<int>& rows);
    int pocketOpacityPercent() const;
    void setPocketOpacityPercent(int percent);
    bool pocketOverlayVisible() const;
    void setPocketOverlayVisible(bool visible);
    int pocketMaterialIndex() const;
    void setPocketMaterialIndex(int index);
signals:
    void pocketChoicesChanged();
    void pocketSelectionChanged(int row);
    void pocketSelectedRowsChanged(const QList<int>& rows);
    void pocketOpacityPercentChanged(int percent);
    void pocketMaterialChanged(int index);
    void pocketVisibilityChanged(bool visible);
    void pocketAppearanceChanged(bool visible,double opacity,double specular);
    void pocketOverlayChanged(const std::vector<ligint::PocketPoint>& points, double spacing, int smoothing, const std::vector<std::vector<ligint::PocketPoint>>& groups);
private:
    QWidget* canvas_=nullptr;
    QCheckBox* overlayVisible_=nullptr;
    QComboBox* overlayMaterial_=nullptr;
    QSpinBox *overlayOpacity_=nullptr,*overlaySpecular_=nullptr;
    QComboBox *renderMode_=nullptr,*ligand_=nullptr,*alphaRender_=nullptr;
    QDoubleSpinBox *seedSpacing_=nullptr,*seedProminence_=nullptr,*mergeDistance_=nullptr,*mergeDrop_=nullptr;
    QSpinBox *maxSeeds_=nullptr,*noiseNeighbors_=nullptr,*noiseComponent_=nullptr,*smoothSteps_=nullptr,*lightStrength_=nullptr,*shine_=nullptr;
    std::vector<ligint::V3> ligandCoords_;
    QComboBox *algorithm_=nullptr,*morph_=nullptr,*display_=nullptr,*chain_=nullptr,*palette_=nullptr;
    QDoubleSpinBox *spacing_=nullptr,*buried_=nullptr,*ray_=nullptr,*surface_=nullptr,*probe_=nullptr;
    QDoubleSpinBox *alphaMinRadius_=nullptr,*alphaMaxRadius_=nullptr,*alphaNeighborCutoff_=nullptr,
                   *alphaClusterDistance_=nullptr,*alphaMaxClusterRadius_=nullptr,
                   *alphaMacroLink_=nullptr,*alphaMacroDiameter_=nullptr,*alphaMacroClearance_=nullptr;
    QSpinBox *alphaMinContacts_=nullptr,*alphaMinSpheres_=nullptr,*alphaOpacity_=nullptr;
    QCheckBox *alphaMerge_=nullptr;
    QComboBox *alphaRank_=nullptr;
    QSpinBox *radius_=nullptr,*minPts_=nullptr,*opacity_=nullptr,*pointSize_=nullptr;
    QListWidget* entries_=nullptr;
    QLabel* status_=nullptr;
    QPushButton *run_=nullptr,*exportAlpha_=nullptr,*copyAlphaCommand_=nullptr;
    QString path_;
    ligint::PocketSearchResult result_;
    std::vector<ligint::PocketAtom> atoms_;
    std::vector<ligint::PocketGroup> groups_;
    std::vector<ligint::AlphaPocket> alphaPockets_;
    ligint::AlphaMergeAudit alphaAudit_;
    QFutureWatcher<PocketGuiOutcome> watcher_;
    std::vector<std::size_t> visibleMembers() const;
    void publishSelection();
    friend class PocketCanvas;
};
} // namespace gui
