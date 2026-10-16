// Main window of ligint-gui.
#pragma once

#include <QFutureWatcher>
#include <functional>
#include <QMainWindow>
#include <QString>
#include <QStringList>
#include <vector>

#include "ligint/fingerprint.hpp"
#include "ligint/matrix.hpp"
#include <QPointer>
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

class HelpWindow;
class HitsFilterProxy;
class HitsModel;
class SvgView;
class ViewerWidget;

struct AnalysisOutcome {
    ligint::Result result;
    QString error;
};

// Result of a task that runs in the background and ends by opening a file: a download or the generation of a pocket definition.
struct TaskOutcome {
    bool ok = false;
    QString error;
    QString path;          // the file that was downloaded or written
    int positions = 0;     // pocket positions (KLIFS)
    QStringList candidateLabels;   // KLIFS: several structures match; the user can choose one
    QList<int> candidateIds;
};

// Result of building the interaction matrix in the background.
struct MatrixOutcome {
    bool ok = false;
    QString error;
    ligint::InteractionMatrix matrix;
};

struct KlifsRequest {      // the last KLIFS generation asked, to repeat it with a structure that the user chooses
    QString pdb, chain, model, mapping, output;
};

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

    // With interactive = false errors do not open dialog boxes (used by the tests).
    void setInteractive(bool on) { interactive_ = on; }

    bool openStructure(const QString& path);      // reads the structure (PDB or mmCIF) and lists the ligands
    // Downloads a structure from the RCSB (format "pdb" or "cif") into cacheDir (default: the application data folder) and opens it.
    // A structure that is already in the cache is not downloaded again. Runs in a background thread; wait = true waits for it (tests).
    // False if the request is refused (invalid ID, another task running); the result of the transfer is reported in the log.
    bool downloadRcsb(const QString& id, const QString& format, const QString& cacheDir = QString(), bool wait = false);
    // Builds a pocket definition from a KLIFS mapping file (mappingFile) or, when it is empty, from the KLIFS service; saves it
    // in `output` and loads it. Same rules for the thread and the return value.
    // structureId > 0 names the KLIFS structure to use (several can match the same PDB entry and chain: one for each alternate model).
    bool generateKlifsPocket(const QString& pdb, const QString& chain, const QString& model, const QString& mappingFile,
                             const QString& output, bool wait = false, int structureId = 0);
    bool taskRunning() const;                     // a download or a KLIFS generation is in progress

    // Interaction matrix: the fingerprints of the ligands of a list file on the same columns (docs/MATRIX.md).
    // loadMatrixList reads the list (and shows the Matrix tab); runMatrix builds the matrix in a background thread with the
    // analysis settings of the panel and the pocket definition that is loaded (wait = true: waits for it, for the tests);
    // exportMatrix writes .svg, .csv or .json. False if the request is refused; the details are in the log.
    bool loadMatrixList(const QString& path);
    bool runMatrix(bool wait = false);
    bool exportMatrix(const QString& path);
    bool matrixRunning() const;
    const ligint::InteractionMatrix& matrix() const { return matrix_; }
    QTableWidget* matrixTable() const { return matrixTable_; }
    // The `ligint --matrix` command that reproduces the last matrix (format: csv, json or svg).
    QString matrixCommandLine(const QString& format = QStringLiteral("csv")) const;
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
    // The manual (Help menu, F1): `target` is a page of the manual, e.g. "gui.html#the-log"; empty = the page of the current tab.
    void showManual(const QString& target = QString());
    HelpWindow* helpWindow() const;
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
    void onDownloadRcsb();
    void onExportCsv();
    void onSaveImage();
    void onAbout();
    void onAddTemplate();
    void onLoadPocket();
    void onGenerateKlifsPocket();
    void onExportFingerprint();
    void onHelp();
    void onCopyCommand();
    void onCopyFingerprintCommand();
    void onOpenMatrixList();
    void onExportMatrix();
    void onCopyMatrixCommand();
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
    QStringList optionArgs() const;               // the options that differ from the defaults (without structure and ligand)
    void applyMatrix(const MatrixOutcome& o);
    void fillMatrixTable();
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
    QComboBox* altlocCombo_ = nullptr;
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
    QFutureWatcher<TaskOutcome>* taskWatcher_ = nullptr;
    QFutureWatcher<MatrixOutcome>* matrixWatcher_ = nullptr;
    QString matrixListPath_;
    std::vector<ligint::MatrixEntry> matrixEntries_;
    ligint::InteractionMatrix matrix_;
    bool haveMatrix_ = false;
    QStringList lastMatrixArgs_;
    QWidget* matrixTab_ = nullptr;
    QLabel *matrixInfo_ = nullptr, *matrixNote_ = nullptr;
    QPushButton *matrixOpenButton_ = nullptr, *matrixRunButton_ = nullptr;
    QCheckBox *matrixAllClasses_ = nullptr, *matrixHide_ = nullptr;
    QTableWidget* matrixTable_ = nullptr;
    SvgView* matrixFigure_ = nullptr;
    QAction *matrixExportAct_ = nullptr, *matrixCmdAct_ = nullptr;
    QString taskKind_;                            // "rcsb" or "klifs" while a task runs
    QString taskLabel_;                           // PDB ID of the running task, for the messages
    KlifsRequest lastKlifs_;
    void applyTask(const TaskOutcome& o);
    void startTask(const QString& kind, const QString& label, const QString& message, std::function<TaskOutcome()> job, bool wait);
    QAction *exportAct_ = nullptr, *imageAct_ = nullptr, *fpAct_ = nullptr;
    QPointer<HelpWindow> help_;
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
