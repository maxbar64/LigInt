// Finestra principale di ligint-gui.
#pragma once

#include <QFutureWatcher>
#include <QMainWindow>
#include <QString>
#include <QStringList>
#include <vector>

#include "ligint/ligint.hpp"

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
class QTableView;
class QTableWidget;
class QTabWidget;

namespace gui {

class HitsFilterProxy;
class HitsModel;
class ViewerWidget;

struct AnalysisOutcome {
    ligint::Result result;
    QString error;
};

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

    // Con interactive = false gli errori non aprono finestre di dialogo (usato dai test).
    void setInteractive(bool on) { interactive_ = on; }

    bool openStructure(const QString& path);      // legge il PDB ed elenca i ligandi
    int ligandCount() const;
    void selectLigand(int index);
    void analyze(bool wait = false);              // wait = true: attende il risultato (test)
    bool addTemplate(const QString& path);        // SDF o MOL2

    const ligint::Result& result() const { return result_; }
    ViewerWidget* viewer() const { return viewer_; }
    HitsModel* hitsModel() const { return model_; }
    HitsFilterProxy* hitsProxy() const { return proxy_; }
    QString logText() const;

signals:
    void analysisFinished(bool ok);

protected:
    void dragEnterEvent(QDragEnterEvent* e) override;
    void dropEvent(QDropEvent* e) override;

private slots:
    void onOpen();
    void onExportCsv();
    void onSaveImage();
    void onAbout();
    void onAddTemplate();
    void onPresetToggled(bool on);
    void onTypeFilterChanged();
    void onRowChanged();
    void onFinished();

private:
    void buildUi();
    void refreshLigandCombo();
    ligint::Options buildOptions() const;
    void applyOutcome(const AnalysisOutcome& o);
    void fillResidueTable();
    void fillTypeList();
    void log(const QString& line);
    void showError(const QString& title, const QString& text);

    bool interactive_ = true;
    QString path_;
    std::vector<ligint::LigandId> allLigands_, shownLigands_;
    std::vector<ligint::Tmpl> templates_;
    QStringList templateNames_;
    ligint::Result result_;

    // pannello di controllo
    QLabel* fileLabel_ = nullptr;
    QComboBox* ligandCombo_ = nullptr;
    QCheckBox* showSmall_ = nullptr;
    QCheckBox *water_ = nullptr, *merge_ = nullptr, *plipLike_ = nullptr, *backbone_ = nullptr, *metalDonors_ = nullptr;
    QDoubleSpinBox *hbMax_ = nullptr, *hbLong_ = nullptr, *saltMax_ = nullptr, *saltCentroid_ = nullptr,
                   *hydMax_ = nullptr, *piAngle_ = nullptr;
    QLabel* templateLabel_ = nullptr;
    QPushButton* runButton_ = nullptr;
    QListWidget* typeList_ = nullptr;

    // risultati
    ViewerWidget* viewer_ = nullptr;
    QTabWidget* tabs_ = nullptr;
    QTableView* table_ = nullptr;
    QTableWidget* residueTable_ = nullptr;
    QPlainTextEdit* log_ = nullptr;
    HitsModel* model_ = nullptr;
    HitsFilterProxy* proxy_ = nullptr;
    QFutureWatcher<AnalysisOutcome>* watcher_ = nullptr;
    QAction *exportAct_ = nullptr, *imageAct_ = nullptr;
};

}  // namespace gui
