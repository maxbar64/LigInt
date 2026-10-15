// Main window of ligint-gui.
#pragma once

#include <QFutureWatcher>
#include <QMainWindow>
#include <QString>
#include <QStringList>
#include <vector>

#include "ligint/fingerprint.hpp"
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
class QScrollArea;

namespace gui {

class HitsFilterProxy;
class HitsModel;
class SvgView;
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

    // With interactive = false errors do not open dialog boxes (used by the tests).
    void setInteractive(bool on) { interactive_ = on; }

    bool openStructure(const QString& path);      // reads the PDB and lists the ligands
    int ligandCount() const;
    void selectLigand(int index);
    void analyze(bool wait = false);              // wait = true: waits for the result (tests)
    bool addTemplate(const QString& path);        // SDF or MOL2
    bool loadPocket(const QString& path);         // pocket definition for the fingerprint
    void clearPocket();
    bool exportFingerprint(const QString& path);  // format from the suffix: .svg, .csv or .json

    const ligint::Result& result() const { return result_; }
    ViewerWidget* viewer() const { return viewer_; }
    HitsModel* hitsModel() const { return model_; }
    HitsFilterProxy* hitsProxy() const { return proxy_; }
    const ligint::Fingerprint& fingerprint() const { return fp_; }
    QTableWidget* fingerprintTable() const { return fpTable_; }
    QString logText() const;

    // The `ligint` command line equivalent to the GUI settings (POSIX shell quoting; program name `ligint`).
    // commandLine(): the current settings; lastCommandLine(): the settings of the last analysis;
    // fingerprintCommandLine(): the last analysis with --fingerprint and the pocket definition currently loaded.
    QString commandLine() const;
    QString lastCommandLine() const;
    QString fingerprintCommandLine(const QString& format = QStringLiteral("svg")) const;

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
    void onLoadPocket();
    void onExportFingerprint();
    void onCopyCommand();
    void onCopyFingerprintCommand();
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
    void updateFingerprint();
    QStringList currentArgs() const;
    void log(const QString& line);
    void showError(const QString& title, const QString& text);

    bool interactive_ = true;
    QString path_;
    std::vector<ligint::LigandId> allLigands_, shownLigands_;
    std::vector<ligint::Tmpl> templates_;
    QStringList templateNames_;
    ligint::Result result_;

    // control panel
    QLabel* fileLabel_ = nullptr;
    QComboBox* ligandCombo_ = nullptr;
    QCheckBox* showSmall_ = nullptr;
    QCheckBox *water_ = nullptr, *merge_ = nullptr, *plipLike_ = nullptr, *backbone_ = nullptr, *metalDonors_ = nullptr;
    QDoubleSpinBox *hbMax_ = nullptr, *hbLong_ = nullptr, *saltMax_ = nullptr, *saltCentroid_ = nullptr,
                   *hydMax_ = nullptr, *piAngle_ = nullptr;
    QLabel* templateLabel_ = nullptr;
    QLabel* pocketLabel_ = nullptr;
    QPushButton* runButton_ = nullptr;
    QListWidget* typeList_ = nullptr;

    // results
    ViewerWidget* viewer_ = nullptr;
    QTabWidget* tabs_ = nullptr;
    QTableView* table_ = nullptr;
    QTableWidget* residueTable_ = nullptr;
    QPlainTextEdit* log_ = nullptr;
    HitsModel* model_ = nullptr;
    HitsFilterProxy* proxy_ = nullptr;
    QFutureWatcher<AnalysisOutcome>* watcher_ = nullptr;
    QAction *exportAct_ = nullptr, *imageAct_ = nullptr, *fpAct_ = nullptr;
    ligint::PocketDefinition pocket_;
    QString pocketPath_;
    QStringList templatePaths_;
    QStringList pendingArgs_, lastArgs_;
    QAction *cmdAct_ = nullptr, *fpCmdAct_ = nullptr;
    bool havePocket_ = false;
    ligint::Fingerprint fp_;
    QTableWidget* fpTable_ = nullptr;
    SvgView* fpFigure_ = nullptr;
    QLabel* fpNote_ = nullptr;
};

}  // namespace gui
