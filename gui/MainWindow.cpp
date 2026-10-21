#include "MainWindow.h"
#include "ligint/gpcrdb.hpp"
#include "ligint/klifs.hpp"
#include "ligint/rcsb.hpp"
#include <QStandardPaths>
#include <QInputDialog>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QLineEdit>

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QDragEnterEvent>
#include <QDesktopServices>
#include <QDropEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QListWidget>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QPlainTextEdit>
#include <QProcess>
#include <QProgressDialog>
#include <QTimer>
#include <QPushButton>
#include <QScrollArea>
#include <QSettings>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QTableView>
#include <QTableWidget>
#include <QUrl>
#include <QVBoxLayout>
#include <QtConcurrent/QtConcurrent>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "CommandLine.h"
#include "HelpWindow.h"
#include "HitsModel.h"
#include "Style.h"
#include "SvgView.h"
#include "ViewerOptionsWidget.h"
#include "ViewerWidget.h"
#include "PocketWorkspace.h"

namespace gui {
namespace {

AnalysisOutcome runAnalysis(const QString& path, const ligint::Options& opt) {
    AnalysisOutcome out;
    try {
        out.result = ligint::analyzeFile(path.toStdString(), opt);
    } catch (const std::exception& e) {
        out.error = QString::fromUtf8(e.what());
    }
    return out;
}

QDoubleSpinBox* spin(double lo, double hi, double value, double step = 0.1, int decimals = 1) {
    auto* s = new QDoubleSpinBox;
    s->setRange(lo, hi);
    s->setSingleStep(step);
    s->setDecimals(decimals);
    s->setValue(value);
    s->setSuffix(QStringLiteral(" Å"));
    return s;
}

}  // namespace

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle(tr("ligint — ligand-protein interactions"));
    setAcceptDrops(true);
    model_ = new HitsModel(this);
    proxy_ = new HitsFilterProxy(this);
    proxy_->setSourceModel(model_);
    taskWatcher_ = new QFutureWatcher<TaskOutcome>(this);
    connect(taskWatcher_, &QFutureWatcher<TaskOutcome>::finished, this, [this] { applyTask(taskWatcher_->result()); });
    matrixWatcher_ = new QFutureWatcher<MatrixOutcome>(this);
    connect(matrixWatcher_, &QFutureWatcher<MatrixOutcome>::finished, this, [this] { applyMatrix(matrixWatcher_->result()); });
    watcher_ = new QFutureWatcher<AnalysisOutcome>(this);
    connect(watcher_, &QFutureWatcher<AnalysisOutcome>::finished, this, &MainWindow::onFinished);
    buildUi();
    resize(1280, 820);
    statusBar()->showMessage(tr("Open a structure file, PDB or mmCIF (File → Open, or drag it here)"));
}

MainWindow::~MainWindow() {
    if (watcher_ && watcher_->isRunning()) watcher_->waitForFinished();
    if (taskWatcher_ && taskWatcher_->isRunning()) taskWatcher_->waitForFinished();
    if (matrixWatcher_ && matrixWatcher_->isRunning()) matrixWatcher_->waitForFinished();
}

void MainWindow::buildUi() {
    // ---- menu
    auto* file = menuBar()->addMenu(tr("&File"));
    auto makeAction = [this](QMenu* menu, const QString& text, auto slot, QKeySequence key = QKeySequence()) {
        auto* a = new QAction(text, this);
        if (!key.isEmpty()) a->setShortcut(key);
        connect(a, &QAction::triggered, this, slot);
        menu->addAction(a);
        return a;
    };
    makeAction(file, tr("&Open structure…"), &MainWindow::onOpen, QKeySequence::Open);
    makeAction(file, tr("Download from &RCSB…"), &MainWindow::onDownloadRcsb);
    exportAct_ = makeAction(file, tr("Export interactions (&CSV)…"), &MainWindow::onExportCsv);
    imageAct_ = makeAction(file, tr("&Save image…"), &MainWindow::onSaveImage);
    fpAct_ = makeAction(file, tr("Export &fingerprint…"), &MainWindow::onExportFingerprint);
    file->addSeparator();
    cmdAct_ = makeAction(file, tr("Copy &command line"), &MainWindow::onCopyCommand, QKeySequence(Qt::CTRL | Qt::SHIFT | Qt::Key_C));
    fpCmdAct_ = makeAction(file, tr("Copy fingerprint command line"), &MainWindow::onCopyFingerprintCommand);
    cmdAct_->setToolTip(tr("The ligint command that reproduces the last analysis outside the GUI"));
    cmdAct_->setEnabled(false);
    fpCmdAct_->setEnabled(false);
    file->addSeparator();
    makeAction(file, tr("Open ligand &list for the matrix…"), &MainWindow::onOpenMatrixList);
    matrixExportAct_ = makeAction(file, tr("Export &matrix…"), &MainWindow::onExportMatrix);
    matrixCmdAct_ = makeAction(file, tr("Copy matrix command line"), &MainWindow::onCopyMatrixCommand);
    matrixExportAct_->setEnabled(false);
    matrixCmdAct_->setEnabled(false);
    file->addSeparator();
    exportAct_->setEnabled(false);
    imageAct_->setEnabled(false);
    fpAct_->setEnabled(false);
    file->addSeparator();
    makeAction(file, tr("&Quit"), &QWidget::close, QKeySequence::Quit);
    auto* help = menuBar()->addMenu(tr("&Help"));
    makeAction(help, tr("&User manual"), &MainWindow::onHelp, QKeySequence::HelpContents);
    {
        QAction* cli = help->addAction(tr("&Command-line reference"));
        connect(cli, &QAction::triggered, this, [this] { showManual(QStringLiteral("cli.html")); });
        QAction* web = help->addAction(tr("Manual in &web browser"));
        web->setToolTip(tr("Open the HTML files of the manual installed with the program in your web browser"));
        web->setEnabled(!HelpWindow::installedManualDir().isEmpty());
        connect(web, &QAction::triggered, this, [] {
            const QString dir = HelpWindow::installedManualDir();
            if (!dir.isEmpty()) QDesktopServices::openUrl(QUrl::fromLocalFile(dir + QStringLiteral("/index.html")));
        });
        help->addSeparator();
    }
    makeAction(help, tr("&About ligint…"), &MainWindow::onAbout);

    // ---- control panel (left)
    auto* left = new QWidget;
    auto* lv = new QVBoxLayout(left);

    auto* gStruct = new QGroupBox(tr("Structure"));
    auto* sv = new QVBoxLayout(gStruct);
    fileLabel_ = new QLabel(tr("(no file)"));
    fileLabel_->setWordWrap(true);
    auto* openBtn = new QPushButton(tr("Open structure…"));
    connect(openBtn, &QPushButton::clicked, this, &MainWindow::onOpen);
    ligandCombo_ = new QComboBox;
    ligandCombo_->setToolTip(tr("Ligands (HETATM residues) present in the file"));
    showSmall_ = new QCheckBox(tr("Also show ions and small molecules"));
    connect(showSmall_, &QCheckBox::toggled, this, [this] { refreshLigandCombo(); });
    sv->addWidget(openBtn);
    auto* rcsbBtn = new QPushButton(tr("Download from RCSB…"));
    connect(rcsbBtn, &QPushButton::clicked, this, &MainWindow::onDownloadRcsb);
    sv->addWidget(rcsbBtn);
    sv->addWidget(fileLabel_);
    sv->addWidget(new QLabel(tr("Ligand:")));
    sv->addWidget(ligandCombo_);
    sv->addWidget(showSmall_);
    lv->addWidget(gStruct);

    auto* gOpt = new QGroupBox(tr("Options"));
    auto* ov = new QVBoxLayout(gOpt);
    water_ = new QCheckBox(tr("Consider waters (water bridges)"));
    water_->setChecked(true);
    merge_ = new QCheckBox(tr("Merge linked residues (glycans, dipeptides)"));
    merge_->setChecked(true);
    plipLike_ = new QCheckBox(tr("\u00abPLIP-like\u00bb preset"));
    plipLike_->setToolTip(tr("H-bonds up to 4.1 \u00c5 if the geometry is possible, salt bridges by centroids at 5.5 \u00c5, halogens at 4.0 \u00c5"));
    backbone_ = new QCheckBox(tr("Apolar contacts with the backbone"));
    metalDonors_ = new QCheckBox(tr("O/N coordinated to metals may donate H"));
    for (auto* c : {water_, merge_, plipLike_, backbone_, metalDonors_}) ov->addWidget(c);
    {
        auto* altRow = new QHBoxLayout;
        altRow->addWidget(new QLabel(tr("Alternate conformation:")));
        altlocCombo_ = new QComboBox;
        altlocCombo_->addItem(tr("first (default)"), QString());
        altlocCombo_->setToolTip(tr("Conformation to use for the residues that have alternate locations (altloc). "
                                    "Residues without it keep their first conformation."));
        altlocCombo_->setEnabled(false);
        altRow->addWidget(altlocCombo_, 1);
        ov->addLayout(altRow);
    }
    connect(plipLike_, &QCheckBox::toggled, this, &MainWindow::onPresetToggled);
    lv->addWidget(gOpt);

    auto* gThr = new QGroupBox(tr("Thresholds"));
    auto* form = new QFormLayout(gThr);
    hbMax_ = spin(2.5, 6.0, 3.5);
    hbLong_ = spin(0.0, 6.0, 0.0);
    hbLong_->setSpecialValueText(tr("off"));
    saltMax_ = spin(2.5, 6.0, 4.0);
    saltCentroid_ = spin(0.0, 8.0, 0.0);
    saltCentroid_->setSpecialValueText(tr("off"));
    hydMax_ = spin(3.3, 6.0, 4.0);
    piMax_ = spin(0.0, 12.0, 5.5);
    piOffset_ = spin(0.0, 12.0, 3.5);
    piMax_->setSuffix(QStringLiteral(" Å"));
    piOffset_->setSuffix(QStringLiteral(" Å"));
    piAngle_ = spin(10.0, 60.0, 30.0, 1.0, 0);
    piAngle_->setSuffix(QStringLiteral("°"));
    form->addRow(tr("H-bond, D···A max"), hbMax_);
    form->addRow(tr("Long H-bond, up to"), hbLong_);
    form->addRow(tr("Salt bridge (atoms)"), saltMax_);
    form->addRow(tr("Salt bridge (centroids)"), saltCentroid_);
    form->addRow(tr("Hydrophobic max"), hydMax_);
    form->addRow(tr("π-stacking, max distance"), piMax_);
    form->addRow(tr("π-stacking, max offset"), piOffset_);
    form->addRow(tr("π-stacking, max angle"), piAngle_);
    lv->addWidget(gThr);

    auto* gTmpl = new QGroupBox(tr("Chemical templates (SDF / MOL2)"));
    auto* tv = new QVBoxLayout(gTmpl);
    templateLabel_ = new QLabel(tr("none: chemistry inferred from the geometry"));
    templateLabel_->setWordWrap(true);
    auto* tbtn = new QPushButton(tr("Add template…"));
    connect(tbtn, &QPushButton::clicked, this, &MainWindow::onAddTemplate);
    tv->addWidget(tbtn);
    tv->addWidget(templateLabel_);
    lv->addWidget(gTmpl);

    auto* gPocket = new QGroupBox(tr("Pocket definition (fingerprint)"));
    auto* pkv = new QVBoxLayout(gPocket);
    pocketLabel_ = new QLabel(tr("none: residues in contact, sorted by number"));
    pocketLabel_->setWordWrap(true);
    auto* pkRow = new QHBoxLayout;
    auto* pkLoad = new QPushButton(tr("Load definition…"));
    auto* pkClear = new QPushButton(tr("Clear"));
    auto* pkKlifs = new QPushButton(tr("Generate from KLIFS…"));
    pkKlifs->setToolTip(tr("Kinases: the 85 positions of the KLIFS pocket"));
    connect(pkKlifs, &QPushButton::clicked, this, &MainWindow::onGenerateKlifsPocket);
    auto* pkGpcr = new QPushButton(tr("Generate from GPCRdb…"));
    pkGpcr->setObjectName("pocket_gpcrdb");
    pkGpcr->setToolTip(tr("GPCRs: the positions that you choose, by generic residue number (3x32...), read from a structure annotated by GPCRdb or from a table"));
    connect(pkGpcr, &QPushButton::clicked, this, &MainWindow::onGenerateGpcrdbPocket);
    connect(pkLoad, &QPushButton::clicked, this, &MainWindow::onLoadPocket);
    connect(pkClear, &QPushButton::clicked, this, [this] { clearPocket(); });
    pkRow->addWidget(pkLoad);
    pkRow->addWidget(pkClear);
    pkv->addLayout(pkRow);
    pkv->addWidget(pkKlifs);
    pkv->addWidget(pkGpcr);
    auto* pkShared = new QPushButton(tr("Build shared GPCRdb panel…"));
    pkShared->setObjectName("gpcrdb_shared_panel");
    pkShared->setToolTip(tr("Build matched GPCRdb positions across multiple receptor structures; missing positions are audited and never treated as zero contacts."));
    connect(pkShared, &QPushButton::clicked, this, [this] {
        QSettings settings("ligint", "ligint-gui");
        const QStringList files = QFileDialog::getOpenFileNames(
            this, tr("Select GPCRdb mapping tables (two or more)"),
            settings.value("lastGpcrMappingDir", path_.isEmpty() ? QString() : QFileInfo(path_).absolutePath()).toString(),
            tr("GPCRdb mapping TSV (*.tsv);;All files (*)"));
        if (files.isEmpty()) return;
        if (files.size() < 2) {
            showError(tr("Shared GPCRdb fingerprints"), tr("Choose at least two receptor mapping tables.")); return;
        }
        settings.setValue("lastGpcrMappingDir", QFileInfo(files.front()).absolutePath());
        bool accepted = false;
        const QStringList modes = {tr("Intersection (only positions mapped in all structures; recommended)"),
                                   tr("Curated GPCRdb position list (reject incomplete coverage)"),
                                   tr("Union (reject incomplete coverage)")};
        const QString mode = QInputDialog::getItem(this, tr("Common GPCRdb positions"),
             tr("How should the common fingerprint positions be selected?"),
             modes, 0, false, &accepted);
        if (!accepted) return;
        QString curated;
        if (mode == modes[1]) {
            curated = QFileDialog::getOpenFileName(this, tr("Curated generic-position list"),
                        QFileInfo(files.front()).absolutePath(), tr("Text files (*.txt);;All files (*)"));
            if (curated.isEmpty()) return;
        }
        const QString output = QFileDialog::getExistingDirectory(this, tr("Output directory for shared GPCRdb pocket definitions"),
                              QFileInfo(files.front()).absolutePath(), QFileDialog::ShowDirsOnly);
        if (output.isEmpty()) return;
        const QString script = QStringLiteral(LIGINT_SOURCE_ROOT) + "/tools/gpcrdb_shared_fingerprints.py";
        if (!QFileInfo::exists(script)) {
            showError(tr("Shared GPCRdb fingerprints"), tr("Cannot locate GPCRdb panel tool: ") + script); return;
        }
        QString python = QStandardPaths::findExecutable("python3");
        if (python.isEmpty()) python = QStandardPaths::findExecutable("python");
        if (python.isEmpty()) {
            showError(tr("Shared GPCRdb fingerprints"), tr("Python 3 is required.")); return;
        }
        QStringList args = {script, "--output-dir", output};
        if (!curated.isEmpty()) args << "--positions" << curated;
        else args << "--mode" << (mode == modes[0] ? "intersection" : "union");
        args << files;
        QProcess process;
        process.start(python, args);
        if (!process.waitForStarted(3000)) {
            showError(tr("Shared GPCRdb fingerprints"), tr("Failed to start Python: ") + process.errorString()); return;
        }
        QProgressDialog progress(tr("Checking GPCRdb position coverage across structures…"), tr("Cancel"), 0, 0, this);
        progress.setWindowModality(Qt::WindowModal);
        progress.setMinimumDuration(0);
        bool done = false;
        connect(&process, QOverload<int,QProcess::ExitStatus>::of(&QProcess::finished),
                &progress, [&](int, QProcess::ExitStatus) { done = true; progress.accept(); });
        connect(&progress, &QProgressDialog::canceled, &process, &QProcess::kill);
        QTimer guard;
        guard.setSingleShot(true);
        connect(&guard, &QTimer::timeout, &process, &QProcess::kill);
        guard.start(30000);
        if (process.state() != QProcess::NotRunning) progress.exec();
        if (process.state() != QProcess::NotRunning) process.waitForFinished(2000);
        const QString reportPath = output + "/coverage_audit.json";
        QFile reportFile(reportPath);
        QString summary;
        if (reportFile.open(QIODevice::ReadOnly)) {
            const auto doc = QJsonDocument::fromJson(reportFile.readAll());
            const auto object = doc.object();
            const auto positions = object.value("positions").toArray();
            const auto structures = object.value("structures").toArray();
            summary = tr("Common generic positions: %1\\nStructures: %2\\n\\n")
                          .arg(positions.size()).arg(structures.size());
            for (const auto& record : structures) {
                const auto value = record.toObject();
                summary += tr("%1: %2/%3 mapped; %4 missing\\n")
                               .arg(value.value("mapping").toString())
                               .arg(value.value("covered").toInt())
                               .arg(value.value("total").toInt())
                               .arg(value.value("missing").toArray().size());
            }
        }
        const bool success = process.exitStatus() == QProcess::NormalExit && process.exitCode() == 0;
        if (!success) {
            showError(tr("Shared GPCRdb panel not accepted"),
                      QString::fromUtf8(process.readAllStandardError()) +
                      tr("\\n\\nCoverage report: ") + reportPath +
                      (summary.isEmpty() ? QString() : "\\n\\n" + summary));
            return;
        }
        log(tr("Shared GPCRdb pocket definitions: ") + output);
        log(tr("Coverage audit: ") + reportPath);
        QMessageBox::information(this, tr("Shared GPCRdb panel created"),
                                 summary + tr("\\nDefinitions saved in: ") + output +
                                 tr("\\n\\nUse the per-structure *.pocket.tsv files in the Matrix input list."));
    });
    pkv->addWidget(pkShared);
    pkv->addWidget(pocketLabel_);
    lv->addWidget(gPocket);

    runButton_ = new QPushButton(tr("Analyze"));
    runButton_->setEnabled(false);
    runButton_->setDefault(true);
    connect(runButton_, &QPushButton::clicked, this, [this] { analyze(false); });
    lv->addWidget(runButton_);

    auto* gTypes = new QGroupBox(tr("Interaction types"));
    auto* yv = new QVBoxLayout(gTypes);
    typeList_ = new QListWidget;
    typeList_->setMinimumHeight(120);
    connect(typeList_, &QListWidget::itemChanged, this, [this] { onTypeFilterChanged(); });
    yv->addWidget(typeList_);
    lv->addWidget(gTypes, 1);

    // ---- view and tables (right)
    viewer_ = new ViewerWidget;
    connect(viewer_, &ViewerWidget::importLigandPicked, this, [this](int previewIndex) {
        if (!result_.ligands.empty()) return; // only import preview has pickable HETATM residues
        if (previewIndex < 0 || previewIndex >= static_cast<int>(importPreview_.ligands.size())) return;
        const auto& picked = importPreview_.ligands[previewIndex];
        for (int i=0;i<static_cast<int>(shownLigands_.size());++i) {
            const auto& id=shownLigands_[i];
            if(id.chain==picked.chain && id.resi==picked.resi &&
               id.icode==picked.icode && id.res==picked.resName) {
                ligandCombo_->setCurrentIndex(i);
                return;
            }
        }
        if (!showSmall_->isChecked()) {
            showSmall_->setChecked(true);
            for (int i=0;i<static_cast<int>(shownLigands_.size());++i) {
                const auto& id=shownLigands_[i];
                if(id.chain==picked.chain && id.resi==picked.resi &&
                   id.icode==picked.icode && id.res==picked.resName) {
                    ligandCombo_->setCurrentIndex(i);return;
                }
            }
        }
    });
    connect(ligandCombo_, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, [this](int index) {
        if (!result_.ligands.empty() || importPreview_.ligands.empty()) return;
        if(index<0 || index>=static_cast<int>(shownLigands_.size())) return;
        const auto& selected=shownLigands_[index];
        for (int n=0;n<static_cast<int>(importPreview_.ligands.size());++n) {
            const auto& lig=importPreview_.ligands[n];
            if(lig.chain==selected.chain && lig.resi==selected.resi &&
               lig.icode==selected.icode && lig.resName==selected.res) {
                viewer_->setImportPreview(importPreview_,n);
                return;
            }
        }
    });
    graphicsOptions_ = new ViewerOptionsWidget(viewer_);
    auto* graphicsScroll = new QScrollArea;
    graphicsScroll->setObjectName("graphicsOptionsScroll");
    graphicsScroll->setWidgetResizable(true);
    graphicsScroll->setFrameShape(QFrame::NoFrame);
    graphicsScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    graphicsScroll->setMinimumWidth(225);
    graphicsScroll->setMaximumWidth(280);
    graphicsScroll->setWidget(graphicsOptions_);
    auto* viewerArea = new QWidget;
    auto* viewerLayout = new QHBoxLayout(viewerArea);
    viewerLayout->setContentsMargins(0, 0, 0, 0);
    viewerLayout->setSpacing(4);
    viewerLayout->addWidget(viewer_, 1);
    viewerLayout->addWidget(graphicsScroll);
    auto* modes3d = new QTabWidget;
    auto* pocketWorkspace = new PocketWorkspace;
    pocketWorkspace->setObjectName("pocketWorkspace");
    connect(pocketWorkspace,&PocketWorkspace::pocketAppearanceChanged,viewer_,
        &ViewerWidget::setPocketAppearance);
    connect(pocketWorkspace,&PocketWorkspace::pocketOverlayChanged,viewer_,
        [this](const std::vector<ligint::PocketPoint>& points,double spacing,int smoothing){
            viewer_->setPocketOverlay(points,spacing,smoothing);
        });
    modes3d->addTab(viewerArea, tr("Interactions 3D"));
    modes3d->addTab(pocketWorkspace, tr("Pocket Finder 3D"));
    tabs_ = new QTabWidget;
    table_ = new QTableView;
    table_->setModel(proxy_);
    table_->setSortingEnabled(true);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->setAlternatingRowColors(true);
    table_->verticalHeader()->setVisible(false);
    table_->horizontalHeader()->setStretchLastSection(true);
    connect(table_->selectionModel(), &QItemSelectionModel::currentRowChanged, this, [this] { onRowChanged(); });
    residueTable_ = new QTableWidget(0, 3);
    residueTable_->setHorizontalHeaderLabels({tr("Residue"), tr("Contacts"), tr("Types")});
    residueTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    residueTable_->horizontalHeader()->setStretchLastSection(true);
    residueTable_->verticalHeader()->setVisible(false);
    log_ = new QPlainTextEdit;
    log_->setReadOnly(true);
    tabs_->addTab(table_, tr("Interactions"));
    tabs_->addTab(residueTable_, tr("Residues"));

    // fingerprint tab: figure on top (if Qt Svg is available) and table below
    auto* fpTab = new QSplitter(Qt::Vertical);
    fpFigure_ = new SvgView;
    if (SvgView::available()) {
        auto* scroll = new QScrollArea;
        scroll->setWidget(fpFigure_);
        scroll->setWidgetResizable(false);
        scroll->setBackgroundRole(QPalette::Base);
        fpTab->addWidget(scroll);
    } else {
        fpNote_ = new QLabel(tr("The figure needs Qt Svg: this build shows the fingerprint as a table only."));
        fpNote_->setAlignment(Qt::AlignCenter);
        fpTab->addWidget(fpNote_);
        fpFigure_->hide();
    }
    fpTable_ = new QTableWidget(0, 3 + ligint::FpClassCount);
    {
        QStringList heads = {tr("Residue"), tr("Position"), tr("Region")};
        for (int c = 0; c < ligint::FpClassCount; c++) heads << QString::fromUtf8(ligint::fpClassTitle(c));
        fpTable_->setHorizontalHeaderLabels(heads);
    }
    fpTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    fpTable_->verticalHeader()->setVisible(false);
    fpTable_->setAlternatingRowColors(true);
    fpTab->addWidget(fpTable_);
    fpTab->setStretchFactor(0, 3);
    fpTab->setStretchFactor(1, 2);
    tabs_->addTab(fpTab, tr("Fingerprint"));

    // matrix tab: controls, the heat map (if Qt Svg is available) and a table with a row for each ligand
    matrixTab_ = new QWidget;
    {
        auto* box = new QVBoxLayout(matrixTab_);
        auto* bar = new QHBoxLayout;
        matrixOpenButton_ = new QPushButton(tr("Open list…"));
        matrixOpenButton_->setObjectName("matrix_open");
        matrixOpenButton_->setToolTip(tr("A TSV or CSV file with one line for each ligand: structure, ligand, label, pocket… (see the manual)"));
        matrixRunButton_ = new QPushButton(tr("Compute matrix"));
        matrixRunButton_->setObjectName("matrix_run");
        matrixRunButton_->setEnabled(false);
        matrixAllClasses_ = new QCheckBox(tr("All interaction classes"));
        matrixAllClasses_->setObjectName("matrix_all");
        matrixAllClasses_->setToolTip(tr("Unchecked: the 7 classes of the KLIFS interaction fingerprint. Checked: also cation-π, halogen, metal and water bridge."));
        matrixHide_ = new QCheckBox(tr("Hide positions without interactions"));
        matrixHide_->setObjectName("matrix_hide");
        matrixInfo_ = new QLabel(tr("No list: File → Open ligand list for the matrix…"));
        bar->addWidget(matrixOpenButton_);
        bar->addWidget(matrixRunButton_);
        bar->addWidget(matrixAllClasses_);
        bar->addWidget(matrixHide_);
        bar->addWidget(matrixInfo_, 1);
        box->addLayout(bar);
        auto* split = new QSplitter(Qt::Vertical);
        matrixFigure_ = new SvgView;
        if (SvgView::available()) {
            auto* scroll = new QScrollArea;
            scroll->setWidget(matrixFigure_);
            scroll->setWidgetResizable(false);
            scroll->setBackgroundRole(QPalette::Base);
            split->addWidget(scroll);
        } else {
            matrixNote_ = new QLabel(tr("The heat map needs Qt Svg: this build shows the matrix as a table only."));
            matrixNote_->setAlignment(Qt::AlignCenter);
            split->addWidget(matrixNote_);
            matrixFigure_->hide();
        }
        matrixTable_ = new QTableWidget(0, 5);
        matrixTable_->setHorizontalHeaderLabels({tr("Label"), tr("Ligand"), tr("Structure"), tr("Interactions"), tr("Positions")});
        matrixTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
        matrixTable_->verticalHeader()->setVisible(false);
        matrixTable_->setAlternatingRowColors(true);
        matrixTable_->horizontalHeader()->setStretchLastSection(true);
        split->addWidget(matrixTable_);
        split->setStretchFactor(0, 3);
        split->setStretchFactor(1, 2);
        box->addWidget(split, 1);
        connect(matrixOpenButton_, &QPushButton::clicked, this, &MainWindow::onOpenMatrixList);
        connect(matrixRunButton_, &QPushButton::clicked, this, [this] { runMatrix(); });
    }
    tabs_->addTab(matrixTab_, tr("Matrix"));
    tabs_->addTab(log_, tr("Log"));

    auto* right = new QSplitter(Qt::Vertical);
    right->addWidget(modes3d);
    right->addWidget(tabs_);
    right->setStretchFactor(0, 3);
    right->setStretchFactor(1, 2);

    // the control panel scrolls when the window is too small for all its groups
    auto* leftScroll = new QScrollArea;
    leftScroll->setWidget(left);
    leftScroll->setWidgetResizable(true);
    leftScroll->setFrameShape(QFrame::NoFrame);
    leftScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    auto* split = new QSplitter(Qt::Horizontal);
    split->addWidget(leftScroll);
    split->addWidget(right);
    split->setStretchFactor(1, 1);
    leftScroll->setMaximumWidth(380);
    setCentralWidget(split);

    // object names: used by the tests to drive the controls
    water_->setObjectName("opt_water");
    merge_->setObjectName("opt_merge");
    plipLike_->setObjectName("opt_plip");
    backbone_->setObjectName("opt_backbone");
    metalDonors_->setObjectName("opt_metalDonors");
    altlocCombo_->setObjectName("opt_altloc");
    hbMax_->setObjectName("thr_hbMax");
    hbLong_->setObjectName("thr_hbLong");
    saltMax_->setObjectName("thr_saltMax");
    saltCentroid_->setObjectName("thr_saltCentroid");
    hydMax_->setObjectName("thr_hydMax");
    piMax_->setObjectName("thr_piMax");
    piOffset_->setObjectName("thr_piOffset");
    piAngle_->setObjectName("thr_piAngle");
}

// ------------------------------------------------------------------------------------------------
// opening the file and choosing the ligand
// ------------------------------------------------------------------------------------------------
void MainWindow::log(const QString& line) { log_->appendPlainText(line); }

QString MainWindow::logText() const { return log_->toPlainText(); }

void MainWindow::showError(const QString& title, const QString& text) {
    log(title + ": " + text);
    statusBar()->showMessage(text, 8000);
    if (interactive_) QMessageBox::warning(this, title, text);
}

bool MainWindow::openStructure(const QString& path) {
    ligint::Result listing;
    try {
        ligint::Options previewOptions;
        previewOptions.protein3d = true;
        previewOptions.previewLigands3d = true;
        listing = ligint::analyzeFile(path.toStdString(), previewOptions);
    } catch (const std::exception& e) {
        showError(tr("Cannot open the file"), QString::fromUtf8(e.what()));
        return false;
    }
    path_ = path;
    if(auto* pockets=findChild<PocketWorkspace*>("pocketWorkspace"))pockets->setStructurePath(path);
    allLigands_ = listing.hetResidues;
    {   // alternate conformations present in the file; the previous choice is kept if it still exists
        const QString previous = altlocCombo_->currentData().toString();
        altlocCombo_->clear();
        altlocCombo_->addItem(tr("first (default)"), QString());
        for (char c : listing.altlocs) altlocCombo_->addItem(QString(QChar(c)), QString(QChar(c)));
        const int keep = altlocCombo_->findData(previous);
        altlocCombo_->setCurrentIndex(keep >= 0 ? keep : 0);
        altlocCombo_->setEnabled(!listing.altlocs.empty());
    }
    result_ = ligint::Result();
    importPreview_ = listing;
    model_->setHits({});
    if (!listing.backbone.empty() || !listing.ligands.empty()) viewer_->setImportPreview(importPreview_, -1);
    else viewer_->clear();
    residueTable_->setRowCount(0);
    typeList_->clear();
    fp_ = ligint::Fingerprint();
    fpTable_->setRowCount(0);
    fpAct_->setEnabled(false);
    lastArgs_.clear();
    cmdAct_->setEnabled(false);
    fpCmdAct_->setEnabled(false);
    fileLabel_->setText(QFileInfo(path).fileName());
    setWindowTitle(tr("ligint — %1").arg(QFileInfo(path).fileName()));
    for (const auto& w : listing.warnings) log(QString::fromStdString(w));
    log(tr("Opened %1: %2 HETATM residues").arg(path).arg(allLigands_.size()));
    refreshLigandCombo();
    // Keep the full-structure framing initially; combo selection can focus a ligand later.
    if (viewer_->hasData()) viewer_->setImportPreview(importPreview_, -1);
    runButton_->setEnabled(ligandCombo_->count() > 0);
    exportAct_->setEnabled(false);
    imageAct_->setEnabled(viewer_->hasData());
    if (ligandCombo_->count() == 0) statusBar()->showMessage(tr("No HETATM ligand in this file"));
    else statusBar()->showMessage(tr("%1 ligands: choose one and press \u00abAnalyze\u00bb").arg(ligandCombo_->count()));
    return true;
}

void MainWindow::refreshLigandCombo() {
    shownLigands_.clear();
    for (const auto& l : allLigands_)
        if (showSmall_->isChecked() || l.atoms >= 6) shownLigands_.push_back(l);
    if (shownLigands_.empty()) shownLigands_ = allLigands_;       // only ions or small molecules
    std::stable_sort(shownLigands_.begin(), shownLigands_.end(), [](const ligint::LigandId& a, const ligint::LigandId& b) { return a.atoms > b.atoms; });
    ligandCombo_->clear();
    for (const auto& l : shownLigands_)
        ligandCombo_->addItem(tr("%1 %2 %3%4  (%5 atoms)")
                                  .arg(QString::fromStdString(l.res))
                                  .arg(QString::fromStdString(l.chain))
                                  .arg(l.resi)
                                  .arg(l.icode == ' ' ? QString() : QString(QChar(l.icode)))
                                  .arg(l.atoms));
    runButton_->setEnabled(ligandCombo_->count() > 0 && !path_.isEmpty());
}

int MainWindow::ligandCount() const { return ligandCombo_->count(); }
void MainWindow::selectLigand(int index) { ligandCombo_->setCurrentIndex(index); }

bool MainWindow::addTemplate(const QString& path) {
    std::vector<ligint::Tmpl> t;
    std::vector<std::string> warnings;
    const bool mol2 = path.endsWith(".mol2", Qt::CaseInsensitive);
    bool ok = mol2 ? ligint::loadMol2(path.toStdString(), t, warnings) : ligint::loadSdf(path.toStdString(), t, warnings);
    for (const auto& w : warnings) log(QString::fromStdString(w));
    if (!ok || t.empty()) {
        showError(tr("Invalid template"), tr("No readable molecule in %1").arg(path));
        return false;
    }
    templates_.insert(templates_.end(), t.begin(), t.end());
    templatePaths_ << QFileInfo(path).absoluteFilePath();
    templateNames_ << QFileInfo(path).fileName() + tr(" (%1 molecules)").arg(t.size());
    templateLabel_->setText(templateNames_.join("\n"));
    log(tr("Template loaded: %1").arg(path));
    return true;
}

void MainWindow::onAddTemplate() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getOpenFileName(this, tr("SDF or MOL2 template"), s.value("lastDir").toString(),
                                             tr("Molecules (*.sdf *.mol2);;All files (*)"));
    if (!f.isEmpty()) addTemplate(f);
}

void MainWindow::onPresetToggled(bool on) {
    if (on) {
        hbLong_->setValue(4.1);
        saltCentroid_->setValue(5.5);
    } else {
        hbLong_->setValue(0.0);
        saltCentroid_->setValue(0.0);
    }
}

ligint::Options MainWindow::buildOptions() const {
    ligint::Options opt;
    const int idx = ligandCombo_->currentIndex();
    if (idx >= 0 && idx < static_cast<int>(shownLigands_.size())) {
        const auto& l = shownLigands_[idx];
        opt.ligName = l.res;
        opt.chain = l.chain;
        opt.resi = l.resi;
        opt.icode = l.icode;
    }
    opt.useWater = water_->isChecked();
    opt.mergeLinked = merge_->isChecked();
    opt.backboneHyd = backbone_->isChecked();
    opt.metalDonors = metalDonors_->isChecked();
    {
        const QString alt = altlocCombo_->currentData().toString();
        opt.altloc = alt.isEmpty() ? 0 : alt[0].toLatin1();
    }
    opt.templates = templates_;
    auto& p = opt.params;
    p.hbMax = hbMax_->value();
    p.hbLong = hbLong_->value();
    p.saltMax = saltMax_->value();
    p.saltCentroid = saltCentroid_->value();
    p.hydMax = hydMax_->value();
    p.piMax = piMax_->value();
    p.piOffset = piOffset_->value();
    p.piAngle = piAngle_->value();
    if (plipLike_->isChecked()) p.halMax = 4.0;
    return opt;
}

// ------------------------------------------------------------------------------------------------
// analysis
// ------------------------------------------------------------------------------------------------
void MainWindow::analyze(bool wait) {
    if (path_.isEmpty() || ligandCombo_->currentIndex() < 0) return;
    if (watcher_->isRunning()) return;
    ligint::Options opt = buildOptions();
    opt.protein3d = true;                   // backbone and atoms for the cartoon and the surfaces (not part of the command line)
    pendingArgs_ = currentArgs();
    const QString path = path_;
    runButton_->setEnabled(false);
    statusBar()->showMessage(tr("Analysis in progress…"));
    if (wait) {
        applyOutcome(runAnalysis(path, opt));
        return;
    }
    watcher_->setFuture(QtConcurrent::run([path, opt]() { return runAnalysis(path, opt); }));
}

void MainWindow::onFinished() { applyOutcome(watcher_->result()); }

void MainWindow::applyOutcome(const AnalysisOutcome& o) {
    runButton_->setEnabled(true);
    if (!o.error.isEmpty()) {
        showError(tr("Analysis failed"), o.error);
        emit analysisFinished(false);
        return;
    }
    result_ = o.result;
    for (const auto& w : result_.warnings) log(QString::fromStdString(w));
    for (const auto& s : result_.info) log(QString::fromStdString(s));
    if (result_.ligands.empty()) {
        model_->setHits({});
        viewer_->clear();
        statusBar()->showMessage(tr("No ligand analyzed"));
        emit analysisFinished(false);
        return;
    }
    const ligint::LigandResult& lig = result_.ligands.front();
    lastArgs_ = pendingArgs_;
    for (const auto& s : lig.info) log(QString::fromStdString(s));
    log(tr("Command line: ") + lastCommandLine());
    if (!lig.metalSphere.empty()) log(tr("Coordination sphere of the metals:\n") + QString::fromStdString(lig.metalSphere));
    model_->setHits(lig.hits);
    viewer_->setStructure(result_, 0);
    graphicsOptions_->updateAvailability();
    fillTypeList();
    fillResidueTable();
    updateFingerprint();
    fpAct_->setEnabled(true);
    cmdAct_->setEnabled(true);
    fpCmdAct_->setEnabled(true);
    table_->resizeColumnsToContents();
    table_->sortByColumn(HitsModel::Type, Qt::AscendingOrder);
    exportAct_->setEnabled(true);
    imageAct_->setEnabled(true);
    statusBar()->showMessage(tr("%1: %2 interactions, %3 residues").arg(QString::fromStdString(lig.id)).arg(lig.hits.size()).arg(residueTable_->rowCount()));
    emit analysisFinished(true);
}

void MainWindow::fillTypeList() {
    typeList_->blockSignals(true);
    typeList_->clear();
    std::map<QString, int> counts;
    for (const auto& h : result_.ligands.front().hits) counts[QString::fromStdString(h.type)]++;
    QSet<QString> all;
    for (const auto& st : typeStyles()) {
        auto it = counts.find(st.type);
        if (it == counts.end()) continue;
        auto* item = new QListWidgetItem(QStringLiteral("%1 (%2)").arg(st.type).arg(it->second));
        item->setData(Qt::UserRole, st.type);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Checked);
        item->setForeground(st.color.darker(120));
        typeList_->addItem(item);
        all.insert(st.type);
    }
    typeList_->blockSignals(false);
    proxy_->setVisibleTypes(all);
    viewer_->setVisibleTypes(all);
}

void MainWindow::onTypeFilterChanged() {
    QSet<QString> on;
    for (int i = 0; i < typeList_->count(); i++)
        if (typeList_->item(i)->checkState() == Qt::Checked) on.insert(typeList_->item(i)->data(Qt::UserRole).toString());
    proxy_->setVisibleTypes(on);
    viewer_->setVisibleTypes(on);
}

void MainWindow::fillResidueTable() {
    std::map<std::string, std::pair<int, std::set<std::string>>> by;
    for (const auto& h : result_.ligands.front().hits) { by[h.residue].first++; by[h.residue].second.insert(h.type); }
    std::vector<std::pair<std::string, std::pair<int, std::set<std::string>>>> rows(by.begin(), by.end());
    std::stable_sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) { return a.second.first > b.second.first; });
    residueTable_->setRowCount(static_cast<int>(rows.size()));
    for (int i = 0; i < static_cast<int>(rows.size()); i++) {
        QString types;
        for (const auto& t : rows[i].second.second) types += (types.isEmpty() ? "" : ", ") + QString::fromStdString(t);
        residueTable_->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(rows[i].first)));
        residueTable_->setItem(i, 1, new QTableWidgetItem(QString::number(rows[i].second.first)));
        residueTable_->setItem(i, 2, new QTableWidgetItem(types));
    }
    residueTable_->resizeColumnsToContents();
}

// ------------------------------------------------------------------------------------------------
// command line equivalent to the GUI settings
// ------------------------------------------------------------------------------------------------
QStringList MainWindow::currentArgs() const {
    QStringList a;
    a << QFileInfo(path_).absoluteFilePath();
    const int idx = ligandCombo_->currentIndex();
    if (idx >= 0 && idx < static_cast<int>(shownLigands_.size())) {
        const auto& l = shownLigands_[idx];
        a << QString::fromStdString(l.res) << QString::fromStdString(l.chain);
        // the insertion code is written only when it is needed to tell this residue from another one with the same number
        bool ambiguous = false;
        for (const auto& o : allLigands_)
            if (o.res == l.res && o.chain == l.chain && o.resi == l.resi && o.icode != l.icode) ambiguous = true;
        QString number = QString::number(l.resi);
        if (ambiguous) number += l.icode == ' ' ? QStringLiteral("-") : QString(QChar(l.icode));
        a << number;
    }
    a << optionArgs();
    return a;
}

QStringList MainWindow::optionArgs() const {
    QStringList a;
    // only the settings that differ from the program defaults
    const ligint::Params d;
    auto thr = [&](const char* name, double value, double def) {
        if (std::fabs(value - def) > 1e-9) a << QStringLiteral("--%1").arg(name) << QString::number(value, 'g', 6);
    };
    if (!water_->isChecked()) a << "--no-water";
    if (!merge_->isChecked()) a << "--no-merge";
    if (backbone_->isChecked()) a << "--backbone-hyd";
    if (metalDonors_->isChecked()) a << "--metal-donors";
    if (!altlocCombo_->currentData().toString().isEmpty()) a << "--altloc" << altlocCombo_->currentData().toString();
    thr("hb-max", hbMax_->value(), d.hbMax);
    thr("hb-long", hbLong_->value(), d.hbLong);
    thr("salt-max", saltMax_->value(), d.saltMax);
    thr("salt-centroid", saltCentroid_->value(), d.saltCentroid);
    thr("hyd-max", hydMax_->value(), d.hydMax);
    thr("pi-max", piMax_->value(), d.piMax);
    thr("pi-offset", piOffset_->value(), d.piOffset);
    thr("pi-angle", piAngle_->value(), d.piAngle);
    if (plipLike_->isChecked()) thr("hal-max", 4.0, d.halMax);     // the preset also lowers the halogen distance
    for (const QString& t : templatePaths_) a << (t.endsWith(".mol2", Qt::CaseInsensitive) ? "--mol2" : "--sdf") << t;
    return a;
}

QString MainWindow::commandLine() const { return joinCommand("ligint", currentArgs()); }

QString MainWindow::lastCommandLine() const { return lastArgs_.isEmpty() ? QString() : joinCommand("ligint", lastArgs_); }

QString MainWindow::fingerprintCommandLine(const QString& format) const {
    if (lastArgs_.isEmpty()) return QString();
    QStringList a = lastArgs_;
    a << "--fingerprint" << format;
    if (havePocket_ && !pocketPath_.isEmpty()) a << "--pocket" << pocketPath_;
    return joinCommand("ligint", a);
}

void MainWindow::onCopyCommand() {
    const QString c = lastCommandLine();
    if (c.isEmpty()) return;
    QApplication::clipboard()->setText(c);
    statusBar()->showMessage(tr("Command line copied"), 4000);
}

void MainWindow::onCopyFingerprintCommand() {
    const QString c = fingerprintCommandLine();
    if (c.isEmpty()) return;
    QApplication::clipboard()->setText(c);
    statusBar()->showMessage(tr("Fingerprint command line copied"), 4000);
}

// ------------------------------------------------------------------------------------------------
// interaction matrix
// ------------------------------------------------------------------------------------------------
bool MainWindow::matrixRunning() const { return matrixWatcher_ && matrixWatcher_->isRunning(); }

bool MainWindow::loadMatrixList(const QString& path) {
    if (matrixRunning()) return false;
    std::vector<ligint::MatrixEntry> entries;
    try {
        entries = ligint::readMatrixList(path.toStdString());
    } catch (const std::exception& e) {
        showError(tr("Matrix list"), QString::fromUtf8(e.what()));
        return false;
    }
    matrixEntries_ = std::move(entries);
    matrixListPath_ = path;
    haveMatrix_ = false;
    matrix_ = ligint::InteractionMatrix();
    lastMatrixArgs_.clear();
    matrixExportAct_->setEnabled(false);
    matrixCmdAct_->setEnabled(false);
    matrixRunButton_->setEnabled(true);
    matrixInfo_->setText(tr("%1: %2 ligands").arg(QFileInfo(path).fileName()).arg(matrixEntries_.size()));
    if (SvgView::available()) matrixFigure_->load(QByteArray());
    fillMatrixTable();
    log(tr("Matrix list %1: %2 ligands").arg(path).arg(matrixEntries_.size()));
    tabs_->setCurrentWidget(matrixTab_);
    return true;
}

void MainWindow::fillMatrixTable() {
    matrixTable_->setRowCount(0);
    auto cell = [&](int row, int col, const QString& text) { matrixTable_->setItem(row, col, new QTableWidgetItem(text)); };
    if (haveMatrix_) {
        matrixTable_->setRowCount(static_cast<int>(matrix_.rows.size()));
        const std::size_t nc = matrix_.classes.size();
        for (std::size_t r = 0; r < matrix_.rows.size(); ++r) {
            const auto& row = matrix_.rows[r];
            int positions = 0;
            for (std::size_t p = 0; p < matrix_.positions.size(); ++p) {
                bool any = false;
                for (std::size_t c = 0; c < nc; ++c) any = any || row.bits[p * nc + c] != 0;
                positions += any;
            }
            const int i = static_cast<int>(r);
            cell(i, 0, QString::fromStdString(row.label));
            cell(i, 1, QString::fromStdString(row.ligandId));
            cell(i, 2, QFileInfo(QString::fromStdString(row.structure)).fileName());
            cell(i, 3, QString::number(row.bitCount()));
            cell(i, 4, QString::number(positions));
        }
    } else {
        matrixTable_->setRowCount(static_cast<int>(matrixEntries_.size()));
        for (std::size_t r = 0; r < matrixEntries_.size(); ++r) {
            const auto& e = matrixEntries_[r];
            const int i = static_cast<int>(r);
            cell(i, 0, QString::fromStdString(e.label));
            cell(i, 1, QString::fromStdString(e.ligand));
            cell(i, 2, QFileInfo(QString::fromStdString(e.structure)).fileName());
            cell(i, 3, QStringLiteral("—"));
            cell(i, 4, QStringLiteral("—"));
        }
    }
    matrixTable_->resizeColumnsToContents();
}

bool MainWindow::runMatrix(bool wait) {
    if (matrixEntries_.empty()) {
        showError(tr("Matrix"), tr("Open a ligand list first (File → Open ligand list for the matrix…)."));
        return false;
    }
    if (matrixRunning()) return false;
    ligint::MatrixOptions options;
    options.analysis = buildOptions();                              // the ligand of each row comes from the list
    if (havePocket_ && !pocketPath_.isEmpty()) options.pocket = pocketPath_.toStdString();
    options.hideEmptyPositions = matrixHide_->isChecked();
    if (matrixAllClasses_->isChecked()) {
        std::string error;
        ligint::parseMatrixClasses("all", options.classes, error);
    }
    // the command line of this run
    QStringList args;
    args << "--matrix" << QFileInfo(matrixListPath_).absoluteFilePath();
    if (havePocket_ && !pocketPath_.isEmpty()) args << "--pocket" << pocketPath_;
    if (matrixAllClasses_->isChecked()) args << "--matrix-classes" << "all";
    if (matrixHide_->isChecked()) args << "--matrix-hide-empty";
    args << optionArgs();
    lastMatrixArgs_ = args;
    const std::vector<ligint::MatrixEntry> entries = matrixEntries_;
    auto job = [entries, options]() {
        MatrixOutcome o;
        try {
            o.matrix = ligint::makeMatrix(entries, options);
            o.ok = true;
        } catch (const std::exception& e) {
            o.error = QString::fromUtf8(e.what());
        }
        return o;
    };
    statusBar()->showMessage(tr("Computing the matrix…"));
    matrixRunButton_->setEnabled(false);
    if (wait) {
        applyMatrix(job());
        return true;
    }
    QApplication::setOverrideCursor(Qt::BusyCursor);
    matrixWatcher_->setFuture(QtConcurrent::run(job));
    return true;
}

void MainWindow::applyMatrix(const MatrixOutcome& o) {
    while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
    statusBar()->clearMessage();
    matrixRunButton_->setEnabled(!matrixEntries_.empty());
    if (!o.ok) {
        // no matrix is better than one that no command reproduces: go back to the list of the entries
        lastMatrixArgs_.clear();
        haveMatrix_ = false;
        matrix_ = ligint::InteractionMatrix();
        matrixExportAct_->setEnabled(false);
        matrixCmdAct_->setEnabled(false);
        if (SvgView::available()) matrixFigure_->load(QByteArray());
        fillMatrixTable();
        showError(tr("Matrix failed"), o.error);
        return;
    }
    matrix_ = o.matrix;
    haveMatrix_ = true;
    for (const auto& w : matrix_.warnings) log(tr("Warning: ") + QString::fromStdString(w));
    int bits = 0;
    for (const auto& row : matrix_.rows) bits += row.bitCount();
    log(tr("Matrix: %1 ligands, %2 %3, %4 classes, %5 interactions")
            .arg(matrix_.rows.size())
            .arg(matrix_.positions.size())
            .arg(matrix_.pocketMode ? tr("pocket positions") : tr("residues"))
            .arg(matrix_.classes.size())
            .arg(bits));
    fillMatrixTable();
    if (SvgView::available()) {
        ligint::MatrixSvgOptions so;
        so.title = QFileInfo(matrixListPath_).completeBaseName().toStdString();
        if (!matrixFigure_->load(QByteArray::fromStdString(ligint::matrixSvg(matrix_, so)))) log(tr("The matrix figure could not be drawn"));
    }
    matrixExportAct_->setEnabled(true);
    matrixCmdAct_->setEnabled(true);
    log(tr("Command line: ") + matrixCommandLine());
    matrixInfo_->setText(tr("%1: %2 ligands, %3 interactions").arg(QFileInfo(matrixListPath_).fileName()).arg(matrix_.rows.size()).arg(bits));
}

bool MainWindow::exportMatrix(const QString& path) {
    if (!haveMatrix_) return false;
    const QString suffix = QFileInfo(path).suffix().toLower();
    std::ofstream file(path.toStdString(), std::ios::binary);
    if (!file) return false;
    if (suffix == "svg") {
        ligint::MatrixSvgOptions so;
        so.title = QFileInfo(matrixListPath_).completeBaseName().toStdString();
        file << ligint::matrixSvg(matrix_, so);
    } else if (suffix == "csv") ligint::writeMatrixCsv(file, matrix_);
    else if (suffix == "json") ligint::writeMatrixJson(file, matrix_);
    else return false;
    file.flush();
    return static_cast<bool>(file);
}

QString MainWindow::matrixCommandLine(const QString& format) const {
    if (lastMatrixArgs_.isEmpty()) return QString();
    QStringList a = lastMatrixArgs_;
    if (format != QLatin1String("csv")) a << "--matrix-format" << format;
    // the figure of the window has the name of the list as its title
    if (format == QLatin1String("svg")) a << "--matrix-title" << QFileInfo(matrixListPath_).completeBaseName();
    return joinCommand("ligint", a);
}

void MainWindow::onOpenMatrixList() {
    QSettings s("ligint", "ligint-gui");
    const QString f = QFileDialog::getOpenFileName(this, tr("Open ligand list for the matrix"), s.value("lastDir").toString(),
                                                   tr("Ligand lists (*.tsv *.csv *.txt);;All files (*)"));
    if (f.isEmpty()) return;
    s.setValue("lastDir", QFileInfo(f).absolutePath());
    loadMatrixList(f);
}

void MainWindow::onExportMatrix() {
    if (!haveMatrix_) return;
    QSettings s("ligint", "ligint-gui");
    QString selected;
    QString f = QFileDialog::getSaveFileName(this, tr("Export matrix"), s.value("lastDir").toString() + "/matrix.svg",
                                             tr("SVG heat map (*.svg);;CSV table (*.csv);;JSON (*.json)"), &selected);
    if (f.isEmpty()) return;
    if (QFileInfo(f).suffix().isEmpty()) f += selected.contains("CSV") ? ".csv" : selected.contains("JSON") ? ".json" : ".svg";
    if (!exportMatrix(f)) showError(tr("Export failed"), f);
    else statusBar()->showMessage(tr("Exported %1").arg(f), 5000);
}

void MainWindow::onCopyMatrixCommand() {
    const QString c = matrixCommandLine();
    if (c.isEmpty()) return;
    QApplication::clipboard()->setText(c);
    statusBar()->showMessage(tr("Matrix command line copied"), 4000);
}

// ------------------------------------------------------------------------------------------------
// fingerprint and pocket definition
// ------------------------------------------------------------------------------------------------
bool MainWindow::loadPocket(const QString& path) {
    ligint::PocketDefinition def;
    std::vector<std::string> warnings;
    if (!ligint::loadPocketDefinition(path.toStdString(), def, warnings)) {
        showError(tr("Cannot open the pocket definition"), path);
        return false;
    }
    for (const auto& w : warnings) log(tr("Pocket definition: ") + QString::fromStdString(w));
    if (def.positions.empty()) {
        showError(tr("Invalid pocket definition"), tr("No position found in %1").arg(path));
        return false;
    }
    pocket_ = def;
    pocketPath_ = QFileInfo(path).absoluteFilePath();
    havePocket_ = true;
    const QString name = def.name.empty() ? QFileInfo(path).fileName() : QString::fromStdString(def.name);
    pocketLabel_->setText(tr("%1 (%2 positions)").arg(name).arg(def.positions.size()));
    log(tr("Pocket definition loaded: %1").arg(path));
    updateFingerprint();
    return true;
}

void MainWindow::clearPocket() {
    pocket_ = ligint::PocketDefinition();
    pocketPath_.clear();
    havePocket_ = false;
    pocketLabel_->setText(tr("none: residues in contact, sorted by number"));
    updateFingerprint();
}

// ------------------------------------------------------------------------------------------------
// downloads from the RCSB and pocket definitions from KLIFS (background tasks)
// ------------------------------------------------------------------------------------------------
bool MainWindow::taskRunning() const { return taskWatcher_ && taskWatcher_->isRunning(); }

void MainWindow::startTask(const QString& kind, const QString& label, const QString& message, std::function<TaskOutcome()> job, bool wait) {
    taskKind_ = kind;
    taskLabel_ = label;
    statusBar()->showMessage(message);
    if (wait) {
        applyTask(job());
        return;
    }
    QApplication::setOverrideCursor(Qt::BusyCursor);
    taskWatcher_->setFuture(QtConcurrent::run(job));
}

void MainWindow::applyTask(const TaskOutcome& o) {
    while (QApplication::overrideCursor()) QApplication::restoreOverrideCursor();
    statusBar()->clearMessage();
    if (!o.ok) {
        if (taskKind_ == QLatin1String("klifs") && !o.candidateIds.isEmpty() && interactive_) {
            log(tr("KLIFS: ") + o.error);
            bool ok = false;
            const QString pick = QInputDialog::getItem(this, tr("KLIFS structure"),
                                                       tr("Several KLIFS structures match this chain; choose the one to use:"), o.candidateLabels, 0, false, &ok);
            const int index = o.candidateLabels.indexOf(pick);
            if (ok && index >= 0) {
                const KlifsRequest r = lastKlifs_;
                generateKlifsPocket(r.pdb, r.chain, r.model, r.mapping, r.output, false, o.candidateIds[index]);
            }
            return;
        }
        showError(taskKind_ == QLatin1String("rcsb")     ? tr("RCSB download failed")
                  : taskKind_ == QLatin1String("gpcrdb") ? tr("GPCRdb generation failed")
                                                         : tr("KLIFS generation failed"),
                  o.error);
        return;
    }
    if (taskKind_ == QLatin1String("rcsb")) {
        log(tr("RCSB: structure %1 saved to %2").arg(taskLabel_, o.path));
        if (openStructure(o.path)) statusBar()->showMessage(tr("Loaded %1 from RCSB; select a ligand to analyze.").arg(taskLabel_), 10000);
    } else if (taskKind_ == QLatin1String("gpcrdb")) {
        for (const QString& w : o.warnings) log(tr("Warning: ") + w);
        if (loadPocket(o.path)) {
            log(tr("GPCRdb: generated %1 positions").arg(o.positions));
            log(tr("Command line: ") + lastGpcrdbCommand_);
        }
    } else {
        if (loadPocket(o.path)) log(tr("KLIFS: generated %1 mapped residues").arg(o.positions));
    }
}

bool MainWindow::generateGpcrdbPocket(const QString& numberedStructure, const QString& mappingFile, const QString& chain, const QString& positions,
                                      const QString& output, bool wait) {
    if (output.isEmpty()) {
        showError(tr("GPCRdb"), tr("No file to save the pocket definition in."));
        return false;
    }
    ligint::GpcrdbOptions options;
    try {
        options.positions = ligint::parseGpcrdbPositions(positions.toStdString());          // reads the list (or the file) here: a mistake is reported at once
    } catch (const std::exception& e) {
        showError(tr("GPCRdb"), QString::fromUtf8(e.what()));
        return false;
    }
    if (taskRunning()) return false;
    QString numbered = numberedStructure;
    if (numbered.isEmpty() && mappingFile.isEmpty()) numbered = path_;                        // the open structure
    if (numbered.isEmpty() && mappingFile.isEmpty()) {
        showError(tr("GPCRdb"), tr("Open a structure annotated by GPCRdb, choose one, or give a table of generic numbers."));
        return false;
    }
    options.mappingFile = mappingFile.toStdString();
    options.numberedStructure = numbered.toStdString();
    options.chain = chain.trimmed().toStdString();
    // the command that makes the same pocket
    QStringList list, args;
    for (const auto& p : options.positions) list << QString::fromStdString(p);
    if (!mappingFile.isEmpty()) args << "--gpcrdb-mapping" << mappingFile;
    else args << "--gpcrdb-pdb" << numbered;
    if (!options.chain.empty()) args << "--gpcrdb-chain" << chain.trimmed();
    args << "--gpcrdb-positions" << list.join(",") << "--gpcrdb-output" << output;
    lastGpcrdbCommand_ = joinCommand("ligint", args);
    const std::string out = output.toStdString();
    startTask(QStringLiteral("gpcrdb"), QFileInfo(mappingFile.isEmpty() ? numbered : mappingFile).fileName(), tr("Reading the generic residue numbers…"), [options, out]() {
        TaskOutcome o;
        try {
            const ligint::GpcrdbPocket g = ligint::generateGpcrdbPocket(options);
            ligint::saveGpcrdbPocket(g.pocket, out);
            o.path = QString::fromStdString(out);
            o.positions = static_cast<int>(g.pocket.positions.size());
            for (const auto& w : g.warnings) o.warnings << QString::fromStdString(w);
            o.ok = true;
        } catch (const std::exception& e) {
            o.error = QString::fromUtf8(e.what());
        }
        return o;
    }, wait);
    return true;
}

void MainWindow::onGenerateGpcrdbPocket() {
    QSettings s("ligint", "ligint-gui");
    QStringList sources = {tr("Automatically fetch GPCRdb numbering and align receptor chain (recommended)"), tr("The open structure (annotated by GPCRdb)"), tr("Another structure annotated by GPCRdb…"), tr("A table of generic numbers (TSV or CSV)…")};
    if (path_.isEmpty()) { sources.removeFirst(); sources.removeFirst(); } // requires a structure for automatic numbering
    bool ok = false;
    const QString source = QInputDialog::getItem(this, tr("GPCR pocket"), tr("How should numbering be obtained?"), sources, 0, false, &ok);
    if (!ok) return;
    QString numbered, mapping;
    const QString folder = path_.isEmpty() ? s.value("lastDir").toString() : QFileInfo(path_).absolutePath();
    const bool automatic = source.startsWith(tr("Automatically fetch GPCRdb numbering"));
    const int index = sources.indexOf(source) + (path_.isEmpty() ? 1 : -1); // -1 auto, 0 open, 1 another, 2 table
    if (index == 1) {
        numbered = QFileDialog::getOpenFileName(this, tr("Structure annotated by GPCRdb"), folder, tr("Structures (*.pdb *.ent *.cif *.mmcif);;All files (*)"));
        if (numbered.isEmpty()) return;
    } else if (index == 2) {
        mapping = QFileDialog::getOpenFileName(this, tr("Table of generic residue numbers"), folder, tr("Tables (*.tsv *.csv *.txt);;All files (*)"));
        if (mapping.isEmpty()) return;
    }
    QString automaticChain, reviewedPositions;
    if(automatic){
        if(path_.isEmpty()){showError(tr("GPCRdb"),tr("Open a receptor-ligand PDB first."));return;}
        const QString defaultEntry=QStringLiteral("auto");
        const QString entry=QInputDialog::getText(this,tr("Identify receptor"),
            tr("GPCRdb receptor entry (auto = identify from PDB ID; or enter adrb2_human):"),QLineEdit::Normal,defaultEntry,&ok).trimmed();
        if(!ok||entry.isEmpty())return;
        const QString autoChain=QInputDialog::getText(this,tr("Receptor chain"),
            tr("PDB receptor chain (auto = find best receptor chain by sequence):"),
            QLineEdit::Normal,QStringLiteral("auto"),&ok).trimmed();
        if(!ok||autoChain.isEmpty())return;
        const QString out=QFileInfo(path_).absolutePath()+"/"+
                          QFileInfo(path_).completeBaseName()+"-gpcrdb-auto.tsv";
        const QString audit=out+".audit.json";
        const QString script=QStringLiteral(LIGINT_SOURCE_ROOT)+"/tools/gpcrdb_auto_mapping.py";
        if(!QFileInfo::exists(script)){
            showError(tr("GPCRdb"),tr("The GPCRdb alignment script was not found: ")+script);return;
        }
        QProcess process;
        QProgressDialog progress(tr("Downloading GPCRdb reference and validating alignment…"),
                                 tr("Cancel"),0,0,this);
        progress.setWindowModality(Qt::WindowModal);
        progress.setMinimumDuration(0);
        bool finished=false;
        connect(&process,QOverload<int,QProcess::ExitStatus>::of(&QProcess::finished),
                &progress,[&](int,QProcess::ExitStatus){finished=true;progress.accept();});
        connect(&progress,&QProgressDialog::canceled,&process,&QProcess::kill);
        QTimer timeout;
        timeout.setSingleShot(true);
        connect(&timeout,&QTimer::timeout,&process,&QProcess::kill);
        QString python=QStandardPaths::findExecutable("python3");
        if(python.isEmpty())python=QStandardPaths::findExecutable("python");
        if(python.isEmpty()){
            showError(tr("GPCRdb"),tr("Python 3 was not found. Install Python 3 to enable automatic GPCRdb alignment."));return;
        }
        process.start(python,{"-B",script,"--pdb",path_,"--chain",autoChain,
                                 "--entry",entry,"--pdb-id",QFileInfo(path_).completeBaseName().left(4),
                                 "--output",out,"--audit",audit});
        if(!process.waitForStarted(3000)){
            showError(tr("GPCRdb"),tr("Python 3 could not be started: ")+process.errorString());return;
        }
        timeout.start(45000);
        if(!finished)progress.exec();
        if(!finished||process.exitStatus()!=QProcess::NormalExit||process.exitCode()!=0){
            showError(tr("GPCRdb automatic mapping failed"),
                QString::fromUtf8(process.readAllStandardError())+
                tr("\nSee alignment audit: ")+audit);
            return;
        }
        QFile auditFile(audit);
        if(!auditFile.open(QIODevice::ReadOnly)){
            showError(tr("GPCRdb"),tr("Cannot read the generated alignment audit: ")+audit);return;
        }
        const auto auditDoc=QJsonDocument::fromJson(auditFile.readAll());
        if(!auditDoc.isObject()){
            showError(tr("GPCRdb"),tr("Malformed alignment audit."));return;
        }
        const auto auditObj=auditDoc.object();
        automaticChain=auditObj.value("chain").toString();
        if(automaticChain.isEmpty()){
            showError(tr("GPCRdb"),tr("The alignment did not identify a receptor chain."));return;
        }
        QFile mappingFile(out);
        if(!mappingFile.open(QIODevice::ReadOnly|QIODevice::Text)){
            showError(tr("GPCRdb"),tr("Cannot read generated mapping: ")+out);return;
        }
        QDialog review(this);
        review.setWindowTitle(tr("Review GPCRdb numbering — %1 (%2)").arg(entry,automaticChain));
        review.resize(770,570);
        auto* layout=new QVBoxLayout(&review);
        auto* summary=new QLabel(tr("Chain %1 | identity %2% | coverage %3% | GPCRdb positions %4")
            .arg(automaticChain)
            .arg(auditObj.value("identity").toDouble()*100.0,0,'f',1)
            .arg(auditObj.value("coverage").toDouble()*100.0,0,'f',1)
            .arg(auditObj.value("numbered_positions").toInt()),&review);
        layout->addWidget(summary);
        auto* table=new QTableWidget(&review);
        table->setColumnCount(6);
        table->setHorizontalHeaderLabels({tr("Chain"),tr("PDB residue"),tr("AA"),
                                        tr("GPCRdb"),tr("Segment"),tr("Ligand contact")});
        std::set<std::pair<QString,int>> contacted;
        if(!result_.ligands.empty())
            for(const auto& atom:result_.ligands.front().contactAtoms)
                contacted.insert({QString::fromStdString(atom.chain),atom.resi});
        table->setEditTriggers(QAbstractItemView::NoEditTriggers);
        table->setSelectionBehavior(QAbstractItemView::SelectRows);
        table->horizontalHeader()->setStretchLastSection(true);
        while(!mappingFile.atEnd()){
            const QString line=QString::fromUtf8(mappingFile.readLine()).trimmed();
            if(line.startsWith("chain")||line.isEmpty())continue;
            const auto fields=line.split('\t');
            if(fields.size()<5)continue;
            const int row=table->rowCount();
            table->insertRow(row);
            for(int col=0;col<5;++col)
                table->setItem(row,col,new QTableWidgetItem(fields[col]));
            const bool isContact=contacted.count({fields[0],fields[1].toInt()})>0;
            auto* choice=table->item(row,3);
            choice->setFlags(choice->flags()|Qt::ItemIsUserCheckable);
            choice->setCheckState(isContact?Qt::Checked:Qt::Unchecked);
            table->setItem(row,5,new QTableWidgetItem(isContact?tr("Yes"):tr("No")));
        }
        layout->addWidget(new QLabel(tr("Check GPCRdb positions to include in the fingerprint. Contact positions are preselected."),&review));
        layout->addWidget(table);
        auto* buttons=new QDialogButtonBox(QDialogButtonBox::Ok|QDialogButtonBox::Cancel,&review);
        connect(buttons,&QDialogButtonBox::accepted,&review,&QDialog::accept);
        connect(buttons,&QDialogButtonBox::rejected,&review,&QDialog::reject);
        layout->addWidget(buttons);
        if(review.exec()!=QDialog::Accepted)return;
        QStringList selected;
        for(int row=0;row<table->rowCount();++row)
            if(table->item(row,3)->checkState()==Qt::Checked)
                selected<<table->item(row,3)->text();
        if(selected.isEmpty()){
            showError(tr("GPCRdb"),tr("Select at least one generic position."));return;
        }
        reviewedPositions=selected.join(" ");
        s.setValue("gpcrdbEntry",auditObj.value("entry").toString(entry));
        mapping=out;
        log(tr("GPCRdb automatic mapping created: ")+out);
        log(tr("Alignment audit: ")+audit);
    }
    const QString chain = automatic ? automaticChain : QInputDialog::getText(this, tr("GPCR pocket"),
        tr("Receptor chain (e.g. R for 3SN6; do not select a G-protein chain):"),
        QLineEdit::Normal, QString(), &ok);
    if (!automatic && !ok) return;
    QString positions=reviewedPositions;
    if(!reviewedPositions.isEmpty()){
        const QString from=mapping;
        const QString output=QFileDialog::getSaveFileName(this,tr("Save GPCR pocket"),
            QFileInfo(from).absolutePath()+"/"+QFileInfo(path_).completeBaseName()+"-gpcr-pocket.tsv",
            tr("TSV files (*.tsv)"));
        if(!output.isEmpty())generateGpcrdbPocket(numbered,mapping,chain,positions,output);
        return;
    }
    const QStringList methods={tr("Suggest GPCRdb positions from analyzed ligand contacts (recommended)"),
                               tr("Enter a custom list of generic positions")};
    const QString method=QInputDialog::getItem(this,tr("Define GPCR fingerprint pocket"),
        tr("How should LigInt choose the pocket positions?"),methods,0,false,&ok);
    if(!ok)return;
    if(method==methods.front()){
        if(result_.ligands.empty()){
            showError(tr("GPCRdb"),tr("Analyze the selected ligand first; its contacts are needed for the suggested pocket."));
            return;
        }
        if(chain.trimmed().isEmpty()){
            showError(tr("GPCRdb"),tr("Choose the receptor chain to map ligand contacts."));
            return;
        }
        const QString numberedSource=numbered.isEmpty()?path_:numbered;
        std::vector<std::pair<std::string,int>> contacts;
        std::set<std::pair<std::string,int>> unique;
        for(const auto& atom:result_.ligands.front().contactAtoms){
            const auto pair=std::make_pair(atom.chain,atom.resi);
            if(unique.insert(pair).second)contacts.push_back(pair);
        }
        try{
            const auto mappingResult=ligint::mapGpcrdbContacts(
                numberedSource.toStdString(),mapping.toStdString(),
                chain.trimmed().toStdString(),contacts);
            QStringList proposed,details;
            for(const auto& p:mappingResult.positions)proposed<<QString::fromStdString(p);
            for(const auto& item:mappingResult.details)details<<QString::fromStdString(item);
            positions=proposed.join(" ");
            const QString audit=tr("%1 contact positions were mapped. %2 could not be mapped.\n\n%3")
                .arg(proposed.size()).arg(mappingResult.missing.size()).arg(details.join("\n"));
            log(tr("GPCRdb contact mapping:\n")+audit);
            if(!mappingResult.missing.empty()){
                QStringList missing;
                for(const auto& item:mappingResult.missing)missing<<QString::fromStdString(item);
                log(tr("Unmapped contact residues: ")+missing.join(", "));
            }
            if(QMessageBox::question(this,tr("Confirm contact-derived GPCR pocket"),
                audit+tr("\n\nUse these positions for the fingerprint?"),
                QMessageBox::Yes|QMessageBox::No,QMessageBox::Yes)!=QMessageBox::Yes)return;
        }catch(const std::exception& e){
            showError(tr("GPCRdb mapping failed"),QString::fromUtf8(e.what()));return;
        }
    }else{
        positions = QInputDialog::getMultiLineText(
            this,tr("GPCR pocket"),
            tr("Enter GPCRdb positions such as 3x32 3x33 6x48, or a file path."),
            s.value("gpcrdbPositions").toString(),&ok);
        if(!ok||positions.trimmed().isEmpty())return;
    }
    const QString from = !mapping.isEmpty() ? mapping : (!numbered.isEmpty() ? numbered : path_);
    const QString output = QFileDialog::getSaveFileName(this, tr("Save GPCR pocket"), QFileInfo(from).absolutePath() + "/" + QFileInfo(from).completeBaseName() + "-gpcr-pocket.tsv",
                                                        tr("TSV files (*.tsv)"));
    if (output.isEmpty()) return;
    s.setValue("gpcrdbPositions", positions);
    generateGpcrdbPocket(numbered, mapping, chain, positions, output);
}

bool MainWindow::downloadRcsb(const QString& idIn, const QString& format, const QString& cacheDir, bool wait) {
    const QString id = idIn.trimmed().toUpper();
    if (id.size() != 4 || !std::all_of(id.begin(), id.end(), [](QChar c) { return c.isLetterOrNumber() && c.unicode() < 128; })) {
        showError(tr("RCSB"), tr("Enter a valid four-character PDB ID."));
        return false;
    }
    if (format != QLatin1String("pdb") && format != QLatin1String("cif")) {
        showError(tr("RCSB"), tr("The format must be pdb or cif."));
        return false;
    }
    if (taskRunning()) return false;
    ligint::RcsbDownloadOptions options;
    options.pdbId = id.toStdString();
    options.format = format.toStdString();
    options.cacheDir = (cacheDir.isEmpty() ? QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) + "/rcsb" : cacheDir).toStdString();
    startTask(QStringLiteral("rcsb"), id, tr("Downloading %1 from RCSB…").arg(id), [options]() {
        TaskOutcome o;
        try {
            o.path = QString::fromStdString(ligint::downloadRcsbStructure(options));
            o.ok = true;
        } catch (const std::exception& e) {
            o.error = QString::fromUtf8(e.what());
        }
        return o;
    }, wait);
    return true;
}

bool MainWindow::generateKlifsPocket(const QString& pdb, const QString& chain, const QString& model, const QString& mappingFile,
                                     const QString& output, bool wait, int structureId) {
    if (output.isEmpty()) {
        showError(tr("KLIFS"), tr("No file to save the pocket definition in."));
        return false;
    }
    if (taskRunning()) return false;
    ligint::KlifsOptions options;
    options.mappingFile = mappingFile.toStdString();
    options.pdbId = pdb.trimmed().toStdString();
    options.chain = chain.trimmed().toStdString();
    options.model = model.trimmed().toStdString();
    options.structureId = structureId > 0 ? structureId : 0;
    lastKlifs_ = {pdb, chain, model, mappingFile, output};
    const std::string out = output.toStdString();
    startTask(QStringLiteral("klifs"), pdb, mappingFile.isEmpty() ? tr("Asking KLIFS for the pocket of %1…").arg(pdb) : tr("Reading the KLIFS mapping…"),
              [options, out]() {
                  TaskOutcome o;
                  try {
                      const ligint::PocketDefinition pocket = ligint::generateKlifsPocket(options);
                      ligint::saveKlifsPocket(pocket, out);
                      o.path = QString::fromStdString(out);
                      o.positions = static_cast<int>(pocket.positions.size());
                      o.ok = true;
                  } catch (const ligint::KlifsAmbiguousError& e) {                       // several structures match: let the user choose
                      o.error = QString::fromUtf8(e.what());
                      for (const auto& c : e.candidates) {
                          o.candidateIds << c.id;
                          o.candidateLabels << QStringLiteral("ID %1 (chain %2, model %3)").arg(c.id).arg(QString::fromStdString(c.chain),
                              c.model.empty() ? QObject::tr("unknown") : QString::fromStdString(c.model));
                      }
                  } catch (const std::exception& e) {
                      o.error = QString::fromUtf8(e.what());
                  }
                  return o;
              }, wait);
    return true;
}

void MainWindow::onGenerateKlifsPocket() {
    if (path_.isEmpty()) { showError(tr("KLIFS"), tr("Open a PDB/mmCIF structure first.")); return; }
    const QString suggestion = QFileInfo(path_).completeBaseName().left(4).toUpper();
    bool ok = false;
    const QString pdb = QInputDialog::getText(this, tr("KLIFS PDB"), tr("Four-character PDB ID:"), QLineEdit::Normal, suggestion, &ok);
    if (!ok) return;
    const QString chain = QInputDialog::getText(this, tr("KLIFS chain"), tr("PDB author chain:"), QLineEdit::Normal, "A", &ok);
    if (!ok) return;
    QStringList modes = {tr("Local KLIFS mapping CSV/TSV")};
    if (ligint::klifsOnlineAvailable()) modes.prepend(tr("Online (KLIFS service)"));      // the online choice only exists if this build can do it
    const QString mode = modes.size() == 1 ? modes[0] : QInputDialog::getItem(this, tr("KLIFS source"), tr("Mapping source:"), modes, 0, false, &ok);
    if (!ok) return;
    const bool online = ligint::klifsOnlineAvailable() && mode == modes.first();
    QString mapping, model;
    if (online) {
        model = QInputDialog::getText(this, tr("KLIFS alternate model"), tr("KLIFS model (blank = auto if unique):"), QLineEdit::Normal, QString(), &ok);
        if (!ok) return;
    } else {
        mapping = QFileDialog::getOpenFileName(this, tr("KLIFS mapping export"), QFileInfo(path_).absolutePath(), tr("KLIFS mapping or pocket (*.tsv *.csv *.mol2);;All files (*)"));
        if (mapping.isEmpty()) return;
    }
    const QString output = QFileDialog::getSaveFileName(this, tr("Save KLIFS pocket"), QFileInfo(path_).absolutePath() + "/" + pdb + "-klifs-pocket.tsv",
                                                        tr("TSV files (*.tsv)"));
    if (output.isEmpty()) return;
    generateKlifsPocket(pdb, chain, model, mapping, output);
}

void MainWindow::onLoadPocket() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getOpenFileName(this, tr("Pocket definition"), s.value("lastDir").toString(),
                                             tr("Pocket definitions (*.tsv *.csv *.txt);;All files (*)"));
    if (!f.isEmpty()) loadPocket(f);
}

void MainWindow::updateFingerprint() {
    if (result_.ligands.empty()) return;
    fp_ = ligint::makeFingerprint(result_.ligands.front(), havePocket_ ? &pocket_ : nullptr);
    for (const auto& w : fp_.warnings) log(tr("Fingerprint: ") + QString::fromStdString(w));

    fpTable_->setRowCount(static_cast<int>(fp_.entries.size()));
    for (int i = 0; i < static_cast<int>(fp_.entries.size()); i++) {
        const ligint::FingerprintEntry& e = fp_.entries[i];
        fpTable_->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(e.residue)));
        fpTable_->setItem(i, 1, new QTableWidgetItem(QString::fromStdString(e.position)));
        fpTable_->setItem(i, 2, new QTableWidgetItem(QString::fromStdString(e.region)));
        for (int c = 0; c < ligint::FpClassCount; c++) {
            auto* item = new QTableWidgetItem(e.counts[c] > 0 ? QString::number(e.counts[c]) : QString());
            item->setTextAlignment(Qt::AlignCenter);
            if (e.counts[c] > 0) {
                QColor col = fpClassColor(c);
                col.setAlpha(70);
                item->setBackground(col);
            }
            fpTable_->setItem(i, 3 + c, item);
        }
    }
    fpTable_->resizeColumnsToContents();
    if (SvgView::available()) {
        const std::string svg = ligint::fingerprintSvg(fp_);
        if (!fpFigure_->load(QByteArray::fromStdString(svg))) log(tr("Fingerprint figure could not be drawn"));
    }
}

bool MainWindow::exportFingerprint(const QString& path) {
    if (result_.ligands.empty()) return false;
    const QString suffix = QFileInfo(path).suffix().toLower();
    std::ofstream file(path.toStdString());
    if (!file) return false;
    if (suffix == "svg") file << ligint::fingerprintSvg(fp_);
    else if (suffix == "csv") ligint::writeFingerprintCsv(file, {fp_});
    else if (suffix == "json") ligint::writeFingerprintJson(file, {fp_});
    else return false;
    file.flush();
    return static_cast<bool>(file);
}

void MainWindow::onExportFingerprint() {
    if (result_.ligands.empty()) return;
    QSettings s("ligint", "ligint-gui");
    QString selected;
    QString f = QFileDialog::getSaveFileName(this, tr("Export fingerprint"), s.value("lastDir").toString() + "/fingerprint.svg",
                                             tr("SVG figure (*.svg);;CSV table (*.csv);;JSON (*.json)"), &selected);
    if (f.isEmpty()) return;
    if (QFileInfo(f).suffix().isEmpty()) f += selected.contains("CSV") ? ".csv" : selected.contains("JSON") ? ".json" : ".svg";
    if (!exportFingerprint(f)) showError(tr("Export failed"), f);
    else statusBar()->showMessage(tr("Exported %1").arg(f), 5000);
}

void MainWindow::onRowChanged() {
    const QModelIndex cur = table_->selectionModel()->currentIndex();
    if (!cur.isValid()) { viewer_->setHighlight(-1); return; }
    viewer_->setHighlight(proxy_->mapToSource(cur).data(HitsModel::HitIndexRole).toInt());
}

// ------------------------------------------------------------------------------------------------
// files and dialogs
// ------------------------------------------------------------------------------------------------
void MainWindow::onDownloadRcsb() {
    bool ok = false;
    const QString label = ligint::rcsbOnlineAvailable() ? tr("PDB ID (four characters):")
                                                         : tr("PDB ID (four characters). This build has no network support: only structures already downloaded can be opened.");
    const QString id = QInputDialog::getText(this, tr("Download from RCSB"), label, QLineEdit::Normal, "3EYG", &ok);
    if (!ok) return;
    const QStringList formats = {"PDB (.pdb)", "mmCIF (.cif)"};
    const QString format = QInputDialog::getItem(this, tr("RCSB format"), tr("Download format:"), formats, 0, false, &ok);
    if (!ok) return;
    downloadRcsb(id, format == formats[0] ? QStringLiteral("pdb") : QStringLiteral("cif"));
}

void MainWindow::onOpen() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getOpenFileName(this, tr("Open structure"), s.value("lastDir").toString(),
                                             tr("Structures (*.pdb *.ent *.cif *.mmcif);;All files (*)"));
    if (f.isEmpty()) return;
    s.setValue("lastDir", QFileInfo(f).absolutePath());
    openStructure(f);
}

void MainWindow::onExportCsv() {
    if (result_.ligands.empty()) return;
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getSaveFileName(this, tr("Export interactions"), s.value("lastDir").toString() + "/interactions.csv", tr("CSV (*.csv)"));
    if (f.isEmpty()) return;
    std::ostringstream out, info;
    ligint::printCsv(out, info, result_);
    std::ofstream file(f.toStdString());
    file << out.str();
    statusBar()->showMessage(tr("Exported %1").arg(f), 5000);
}

void MainWindow::onSaveImage() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getSaveFileName(this, tr("Save image"), s.value("lastDir").toString() + "/ligint.png", tr("PNG image (*.png)"));
    if (f.isEmpty()) return;
    if (!viewer_->renderImage(QSize(1600, 1200)).save(f)) showError(tr("Saving failed"), f);
    else statusBar()->showMessage(tr("Saved %1").arg(f), 5000);
}

void MainWindow::showManual(const QString& target) {
    if (!help_) help_ = new HelpWindow(this);
    help_->setWindowFlag(Qt::Window, true);
    QString page = target;
    if (page.isEmpty()) page = HelpWindow::topicForTab(tabs_->tabText(tabs_->currentIndex()));
    help_->showPage(page);
    help_->show();
    help_->raise();
    help_->activateWindow();
}

void MainWindow::onHelp() { showManual(); }

HelpWindow* MainWindow::helpWindow() const { return help_.data(); }

void MainWindow::onAbout() {
    QMessageBox::about(this, tr("About ligint"),
                       tr("<b>ligint %1</b><br>Non-covalent ligand-protein interactions from PDB and mmCIF files.<br><br>"
                          "The interactions are computed from the geometry; the chemistry of the ligand is inferred "
                          "(or read from SDF/MOL2). The comparison with PLIP is described in the documentation.")
                           .arg(LIGINT_VERSION));
}

void MainWindow::dragEnterEvent(QDragEnterEvent* e) {
    if (e->mimeData()->hasUrls()) e->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent* e) {
    const auto urls = e->mimeData()->urls();
    if (!urls.isEmpty()) openStructure(urls.first().toLocalFile());
}

}  // namespace gui
