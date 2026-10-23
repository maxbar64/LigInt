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
#include <QPixmap>
#include <QIcon>
#include <QListWidget>
#include <QMenu>
#include <QActionGroup>
#include <QSpinBox>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QPlainTextEdit>
#include <QProcess>
#include <QProgressDialog>
#include <QTimer>
#include <QPushButton>
#include <QScrollArea>
#include <QSlider>
#include <QSettings>
#include <QSignalBlocker>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QTableView>
#include <QTableWidget>
#include <QColor>
#include <QGuiApplication>
#include <QScreen>
#include <QToolButton>
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
    poseWatcher_ = new QFutureWatcher<std::vector<ligand_pose_batch::Row>>(this);
    connect(poseWatcher_, &QFutureWatcher<std::vector<ligand_pose_batch::Row>>::finished, this, [this] {
        poseRows_ = poseWatcher_->result();
        const auto poseGroups = ligand_pose_batch::exactFingerprintGroups(poseRows_);
        poseTable_->setSortingEnabled(false);
        poseTable_->setRowCount(static_cast<int>(poseRows_.size()));
        for (int i = 0; i < static_cast<int>(poseRows_.size()); ++i) {
            const auto& r = poseRows_[i];
            const QStringList fields = {QString::number(r.index), QString::fromStdString(r.name),
                r.error.empty() ? tr("OK") : tr("Error"), QString::number(r.favorable),
                QString::number(r.hbonds), QString::number(r.hydrophobic), QString::number(r.pi),
                QString::number(r.water), QString::number(r.clashes), QString::number(r.residues),
                poseGroups[i] ? tr("G%1").arg(poseGroups[i]) : tr("Error")};
            for (int j = 0; j < fields.size(); ++j) {
                auto* item = new QTableWidgetItem(fields[j]);
                if (j != 1 && j != 2 && j != 10) item->setData(Qt::EditRole, j == 0 ? QVariant(static_cast<int>(r.index)) : QVariant(fields[j].toInt()));
                poseTable_->setItem(i, j, item);
            }
            if (!r.error.empty()) poseTable_->item(i, 2)->setToolTip(QString::fromStdString(r.error));
        }
        poseTable_->setSortingEnabled(true);
        std::map<std::string, int> frequency;
        int successful = 0;
        for (const auto& r : poseRows_) {
            if (!r.error.empty()) continue;
            ++successful;
            for (const auto& residue : r.contactResidues) ++frequency[residue];
        }
        poseResidueTable_->setSortingEnabled(false);
        poseResidueTable_->setColumnCount(3 + static_cast<int>(poseRows_.size()));
        QStringList headers = {tr("Residue"), tr("Poses"), tr("Frequency (%)")};
        for (const auto& pose : poseRows_) headers << tr("Pose %1").arg(pose.index);
        poseResidueTable_->setHorizontalHeaderLabels(headers);
        poseResidueTable_->setRowCount(static_cast<int>(frequency.size()));
        int residueRow = 0;
        for (const auto& entry : frequency) {
            auto* name = new QTableWidgetItem(QString::fromStdString(entry.first));
            auto* count = new QTableWidgetItem;
            count->setData(Qt::EditRole, entry.second);
            auto* percent = new QTableWidgetItem;
            percent->setData(Qt::EditRole, successful ? 100.0 * entry.second / successful : 0.0);
            poseResidueTable_->setItem(residueRow, 0, name);
            poseResidueTable_->setItem(residueRow, 1, count);
            poseResidueTable_->setItem(residueRow, 2, percent);
            for (int pi = 0; pi < static_cast<int>(poseRows_.size()); ++pi) {
                const auto& pose = poseRows_[pi];
                QStringList codes, descriptions;
                const auto typesIt = pose.residueTypes.find(entry.first);
                if (pose.error.empty() && typesIt != pose.residueTypes.end()) {
                    for (const auto& type : typesIt->second) {
                        descriptions << QString::fromStdString(type);
                        if (type == ligint::type::HydrogenBond) codes << QStringLiteral("HB");
                        else if (type == ligint::type::Hydrophobic || type == ligint::type::BackboneApolar) codes << QStringLiteral("HY");
                        else if (type == ligint::type::PiStackingParallel || type == ligint::type::PiStackingT || type == ligint::type::CationPi) codes << QStringLiteral("PI");
                        else if (type == ligint::type::SaltBridge) codes << QStringLiteral("SB");
                        else if (type == ligint::type::WaterBridge) codes << QStringLiteral("WB");
                        else if (type == ligint::type::HalogenBond) codes << QStringLiteral("XB");
                        else if (type == ligint::type::MetalCoordination) codes << QStringLiteral("ME");
                        else codes << QStringLiteral("OT");
                    }
                    codes.removeDuplicates();
                }
                auto* cell = new QTableWidgetItem(pose.error.empty() ? codes.join(QStringLiteral("+")) : QStringLiteral("—"));
                cell->setToolTip(descriptions.join(QStringLiteral("\n")));
                cell->setTextAlignment(Qt::AlignCenter);
                // Heatmap: color reflects interaction class, not strength or binding affinity.
                // With multiple classes, use the priority below and retain every class in the text/tooltip.
                if (!pose.error.empty()) {
                    cell->setBackground(QColor(220, 220, 220));
                } else if (codes.contains(QStringLiteral("HB"))) {
                    cell->setBackground(QColor(175, 212, 245));
                } else if (codes.contains(QStringLiteral("SB"))) {
                    cell->setBackground(QColor(224, 182, 241));
                } else if (codes.contains(QStringLiteral("PI"))) {
                    cell->setBackground(QColor(248, 216, 154));
                } else if (codes.contains(QStringLiteral("HY"))) {
                    cell->setBackground(QColor(189, 225, 176));
                } else if (codes.contains(QStringLiteral("WB"))) {
                    cell->setBackground(QColor(175, 230, 232));
                } else if (codes.contains(QStringLiteral("XB"))) {
                    cell->setBackground(QColor(246, 192, 180));
                } else if (codes.contains(QStringLiteral("ME"))) {
                    cell->setBackground(QColor(223, 204, 161));
                } else if (!codes.isEmpty()) {
                    cell->setBackground(QColor(220, 220, 220));
                }
                if (!codes.isEmpty()) cell->setForeground(QColor(25, 35, 45));
                poseResidueTable_->setItem(residueRow, 3 + pi, cell);
            }
            ++residueRow;
        }
        poseResidueTable_->setSortingEnabled(true);
        poseResidueTable_->sortItems(1, Qt::DescendingOrder);
        poseRunButton_->setEnabled(true);
        poseExportButton_->setEnabled(!poseRows_.empty());
        poseMatrixExportButton_->setEnabled(!poseRows_.empty());
        statusBar()->showMessage(tr("Pose comparison: %1 poses analyzed").arg(poseRows_.size()));
    });
    watcher_ = new QFutureWatcher<AnalysisOutcome>(this);
    connect(watcher_, &QFutureWatcher<AnalysisOutcome>::finished, this, &MainWindow::onFinished);
    buildUi();
    // Keep the initial window within the usable desktop area (including display scaling).
    const QRect desktop = QGuiApplication::primaryScreen() ? QGuiApplication::primaryScreen()->availableGeometry() : QRect(0, 0, 1280, 800);
    QSettings layoutSettings("ligint", "ligint-gui");
    const QByteArray savedGeometry = layoutSettings.value("layout/geometry").toByteArray();
    if (!savedGeometry.isEmpty()) restoreGeometry(savedGeometry);
    else resize(qMin(1280, qMax(640, desktop.width() - 40)),
                qMin(820, qMax(480, desktop.height() - 60)));
    if (mainSplitter_ && !layoutSettings.value("layout/mainSplitter").toByteArray().isEmpty())
        mainSplitter_->restoreState(layoutSettings.value("layout/mainSplitter").toByteArray());
    if (resultsSplitter_ && !layoutSettings.value("layout/resultsSplitter").toByteArray().isEmpty())
        resultsSplitter_->restoreState(layoutSettings.value("layout/resultsSplitter").toByteArray());
    if (graphicsToggle_) graphicsToggle_->setChecked(layoutSettings.value("layout/graphicsVisible", false).toBool());
    // Restore controls after they have been constructed. Every graphics widget
    // has a stable objectName and updates the renderer through its normal signal.
    if (graphicsOptions_) {
        const auto controls = graphicsOptions_->findChildren<QWidget*>();
        for (QWidget* w : controls) {
            const QString name = w->objectName();
            if (!name.startsWith("graphics_")) continue;
            const QString key = "graphics/" + name;
            if (!layoutSettings.contains(key)) continue;
            if (auto* combo = qobject_cast<QComboBox*>(w))
                combo->setCurrentIndex(layoutSettings.value(key).toInt());
            else if (auto* check = qobject_cast<QCheckBox*>(w))
                check->setChecked(layoutSettings.value(key).toBool());
            else if (auto* slider = qobject_cast<QSlider*>(w))
                slider->setValue(layoutSettings.value(key).toInt());
        }
    }
    if (auto* visible = findChild<QCheckBox*>("quickPocketVisible"))
        visible->setChecked(layoutSettings.value("graphics/pocketVisible", true).toBool());
    if (auto* colors = findChild<QCheckBox*>("quickPocketColors"))
        colors->setChecked(layoutSettings.value("graphics/pocketColors", false).toBool());
    if (auto* opacity = findChild<QSpinBox*>("quickPocketOpacity"))
        opacity->setValue(layoutSettings.value("graphics/pocketOpacity", 38).toInt());
    if (auto* material = findChild<QComboBox*>("quickPocketMaterial"))
        material->setCurrentIndex(layoutSettings.value("graphics/pocketMaterial", 0).toInt());
    if (auto* legend = findChild<QToolButton*>("pocketLegendToggle"))
        legend->setChecked(layoutSettings.value("graphics/pocketLegendVisible", false).toBool());
    if (controlPages_) controlPages_->setCurrentIndex(0);
    statusBar()->showMessage(tr("Open a structure file, PDB or mmCIF (File → Open, or drag it here)"));
}

MainWindow::~MainWindow() {
    // A pending CCD lookup is a child QProcess. Stop and reap it before QObject
    // tears down its children; otherwise Qt warns (and may crash on teardown).
    ++ccdGeneration_;
    if (ccdProcess_) {
        ccdProcess_->disconnect(this);
        if (ccdProcess_->state() != QProcess::NotRunning) {
            ccdProcess_->kill();
            ccdProcess_->waitForFinished(3000);
        }
        ccdProcess_ = nullptr;
    }
    if (watcher_ && watcher_->isRunning()) watcher_->waitForFinished();
    if (taskWatcher_ && taskWatcher_->isRunning()) taskWatcher_->waitForFinished();
    if (matrixWatcher_ && matrixWatcher_->isRunning()) matrixWatcher_->waitForFinished();
    if (poseWatcher_ && poseWatcher_->isRunning()) poseWatcher_->waitForFinished();
    QSettings layoutSettings("ligint", "ligint-gui");
    layoutSettings.setValue("layout/geometry", saveGeometry());
    if (mainSplitter_) layoutSettings.setValue("layout/mainSplitter", mainSplitter_->saveState());
    if (resultsSplitter_) layoutSettings.setValue("layout/resultsSplitter", resultsSplitter_->saveState());
    if (graphicsToggle_) layoutSettings.setValue("layout/graphicsVisible", graphicsToggle_->isChecked());
    if (graphicsOptions_) {
        const auto controls = graphicsOptions_->findChildren<QWidget*>();
        for (QWidget* w : controls) {
            const QString name = w->objectName();
            if (!name.startsWith("graphics_")) continue;
            const QString key = "graphics/" + name;
            if (auto* combo = qobject_cast<QComboBox*>(w))
                layoutSettings.setValue(key, combo->currentIndex());
            else if (auto* check = qobject_cast<QCheckBox*>(w))
                layoutSettings.setValue(key, check->isChecked());
            else if (auto* slider = qobject_cast<QSlider*>(w))
                layoutSettings.setValue(key, slider->value());
        }
    }
    if (auto* visible = findChild<QCheckBox*>("quickPocketVisible"))
        layoutSettings.setValue("graphics/pocketVisible", visible->isChecked());
    if (auto* colors = findChild<QCheckBox*>("quickPocketColors"))
        layoutSettings.setValue("graphics/pocketColors", colors->isChecked());
    if (auto* opacity = findChild<QSpinBox*>("quickPocketOpacity"))
        layoutSettings.setValue("graphics/pocketOpacity", opacity->value());
    if (auto* material = findChild<QComboBox*>("quickPocketMaterial"))
        layoutSettings.setValue("graphics/pocketMaterial", material->currentIndex());
    if (auto* legend = findChild<QToolButton*>("pocketLegendToggle"))
        layoutSettings.setValue("graphics/pocketLegendVisible", legend->isChecked());

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
    ligandCombo_->setToolTip(tr("Choose a crystallographic ligand from the receptor or an external SDF/MOL2 ligand pose."));
    showSmall_ = new QCheckBox(tr("Show ions / small molecules"));
    connect(showSmall_, &QCheckBox::toggled, this, [this] { refreshLigandCombo(); });
    sv->addWidget(openBtn);
    auto* rcsbBtn = new QPushButton(tr("Get from RCSB…"));
    connect(rcsbBtn, &QPushButton::clicked, this, &MainWindow::onDownloadRcsb);
    sv->addWidget(rcsbBtn);
    sv->addWidget(fileLabel_);
    sv->addWidget(new QLabel(tr("Ligand to analyze:")));
    sv->addWidget(ligandCombo_);
    sv->addWidget(showSmall_);
    auto* previewHydrogens=new QCheckBox(tr("Preview donor H (CCD)"));
    previewHydrogens->setObjectName("preview_donor_hydrogens");
    previewHydrogens->setToolTip(tr("Draw N-H/O-H only at sites explicitly defined by the selected CCD component; coordinates are estimated and do not establish experimental protonation."));
    connect(previewHydrogens,&QCheckBox::toggled,this,[this](bool checked){
        viewer_->setPreviewDonorHydrogens(checked);
    });
    sv->addWidget(previewHydrogens);

    auto* ccdButton=new QPushButton(tr("Enrich from CCD…"));
    ccdButton->setObjectName("preview_ccd_enrich");
    ccdButton->setToolTip(tr("Fetch cached or official RCSB Chemical Component Dictionary bonds. Only exact, unambiguous atom-name matches are accepted."));
    connect(ccdButton,&QPushButton::clicked,this,&MainWindow::startCcdPreview);
    sv->addWidget(ccdButton);

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

    auto* gTmpl = new QGroupBox(tr("Ligand poses (SDF/MOL2)"));
    auto* tv = new QVBoxLayout(gTmpl);
    templateLabel_ = new QLabel(tr("none: ligand chemistry inferred from geometry"));
    templateLabel_->setWordWrap(true);
    auto* tbtn = new QPushButton(tr("Add ligand structure…"));
    connect(tbtn, &QPushButton::clicked, this, &MainWindow::onAddTemplate);
    tv->addWidget(tbtn);
    tv->addWidget(templateLabel_);
    excludePdbLigand_ = new QLineEdit;
    excludePdbLigand_->setPlaceholderText(tr("e.g. BZC (optional)"));
    excludePdbLigand_->setToolTip(tr("When analyzing an external pose, omit all HETATM atoms with this residue name from receptor interactions. Other cofactors and waters remain."));
    auto* excludeForm = new QFormLayout;
    excludeForm->addRow(tr("Exclude PDB:"), excludePdbLigand_);
    tv->addLayout(excludeForm);
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
    graphicsToggle_ = new QToolButton;
    auto* graphicsToggle = graphicsToggle_;
    graphicsToggle->setText(tr("3D options"));
    graphicsToggle->setCheckable(true);
    graphicsToggle->setChecked(false);
    graphicsToggle->setToolTip(tr("Show or hide the 3D appearance controls to enlarge the viewer"));
    graphicsToggle->setToolButtonStyle(Qt::ToolButtonTextOnly);
    auto* viewerColumn = new QVBoxLayout;
    viewerColumn->setContentsMargins(0, 0, 0, 0);
    viewerColumn->setSpacing(2);
    viewerColumn->addWidget(graphicsToggle, 0, Qt::AlignRight);
    viewerColumn->addWidget(viewer_, 1);
    viewerLayout->addLayout(viewerColumn, 1);
    graphicsScroll->setVisible(false);
    viewerLayout->addWidget(graphicsScroll);
    connect(graphicsToggle, &QToolButton::toggled, graphicsScroll, &QWidget::setVisible);
    auto* modes3d = new QTabWidget;
    auto* pocketWorkspace = new PocketWorkspace;
    pocketWorkspace->setObjectName("pocketWorkspace");
    connect(pocketWorkspace,&PocketWorkspace::pocketAppearanceChanged,viewer_,
        &ViewerWidget::setPocketAppearance);
    connect(pocketWorkspace,&PocketWorkspace::pocketOverlayChanged,viewer_,
        [this](const std::vector<ligint::PocketPoint>& points,double spacing,int smoothing,
               const std::vector<std::vector<ligint::PocketPoint>>& groups){
            viewer_->setPocketOverlay(points,spacing,smoothing,groups);
        });
    connect(viewer_, &ViewerWidget::pocketSurfacePicked, pocketWorkspace,
            [this, pocketWorkspace](int visibleIndex, Qt::KeyboardModifiers modifiers, bool isolate) {
                const int target = pocketWorkspace->rowForVisiblePocket(visibleIndex);
                if (target < 0) return;
                viewer_->highlightPocket(visibleIndex);
                QList<int> selected = pocketWorkspace->selectedPocketRows();
                if (isolate) selected = {target};
                else if (modifiers & Qt::ControlModifier) selected.removeAll(target);
                else if ((modifiers & Qt::ShiftModifier) && !selected.contains(target))
                    selected.append(target);
                // Normal click highlights, but does not change pocket visibility.
                if (isolate || (modifiers & (Qt::ControlModifier | Qt::ShiftModifier)))
                    pocketWorkspace->setSelectedPocketRows(selected);
                if (auto* legend = findChild<QListWidget*>("pocketLegend")) {
                    legend->setCurrentRow(target);
                    legend->scrollToItem(legend->item(target));
                }
                statusBar()->showMessage(tr("Pocket %1 — %2")
                    .arg(target + 1)
                    .arg(isolate ? tr("isolated") : tr("selected")), 3500);
            });
    // Compact shared controls: change pocket, material and visibility without tab switching.
    auto* pocketQuickBar = new QWidget;
    pocketQuickBar->setObjectName("pocketQuickBar");
    auto* pocketQuickLayout = new QHBoxLayout(pocketQuickBar);
    pocketQuickLayout->setContentsMargins(2, 0, 2, 0);
    pocketQuickLayout->setSpacing(5);
    auto* quickPocket = new QToolButton;
    quickPocket->setObjectName("quickPocketSelector");
    quickPocket->setText(tr("Select pockets"));
    quickPocket->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    quickPocket->setPopupMode(QToolButton::InstantPopup);
    auto* pocketMenu = new QMenu(quickPocket);
    quickPocket->setMenu(pocketMenu);
    auto* quickMaterial = new QComboBox;
    quickMaterial->setObjectName("quickPocketMaterial");
    quickMaterial->addItems({"Satin translucent", "Matte", "Glossy", "Custom"});
    auto* quickVisible = new QCheckBox(tr("Show pocket"));
    quickVisible->setObjectName("quickPocketVisible");
    quickVisible->setChecked(pocketWorkspace->pocketOverlayVisible());
    auto* quickOpacity = new QSpinBox;
    quickOpacity->setObjectName("quickPocketOpacity");
    quickOpacity->setRange(5, 100);
    quickOpacity->setSuffix("%");
    quickOpacity->setValue(pocketWorkspace->pocketOpacityPercent());
    quickOpacity->setToolTip(tr("3D pocket opacity"));
    auto* distinctColors = new QCheckBox(tr("Different colors"));
    distinctColors->setObjectName("quickPocketColors");
    pocketQuickLayout->addWidget(new QLabel(tr("Pockets")));
    pocketQuickLayout->addWidget(quickPocket);
    auto* showAllPockets = new QToolButton;
    showAllPockets->setObjectName("quickPocketShowAll");
    showAllPockets->setText(tr("Show all"));
    showAllPockets->setToolTip(tr("Display every pocket"));
    auto* hideAllPockets = new QToolButton;
    hideAllPockets->setObjectName("quickPocketHideAll");
    hideAllPockets->setText(tr("Hide all"));
    hideAllPockets->setToolTip(tr("Hide all pockets; restore them with Show all or the legend"));
    pocketQuickLayout->addWidget(showAllPockets);
    pocketQuickLayout->addWidget(hideAllPockets);
    connect(showAllPockets, &QToolButton::clicked, pocketWorkspace, [pocketWorkspace] {
        QList<int> rows;
        const int count = pocketWorkspace->pocketChoices().size();
        for (int i = 0; i < count; ++i) rows.append(i);
        pocketWorkspace->setSelectedPocketRows(rows);
    });
    connect(hideAllPockets, &QToolButton::clicked, pocketWorkspace, [pocketWorkspace] {
        pocketWorkspace->setSelectedPocketRows({});
    });
    pocketQuickLayout->addWidget(new QLabel(tr("Material")));
    pocketQuickLayout->addWidget(quickMaterial);
    pocketQuickLayout->addWidget(new QLabel(tr("Opacity")));
    pocketQuickLayout->addWidget(quickOpacity);
    pocketQuickLayout->addWidget(distinctColors);
    pocketQuickLayout->addWidget(quickVisible);
    viewerColumn->insertWidget(1, pocketQuickBar);
    // Collapsible interactive legend: checked items mirror the Pocket Finder selection.
    auto* legendToggle = new QToolButton;
    legendToggle->setObjectName("pocketLegendToggle");
    legendToggle->setText(tr("Pocket legend"));
    legendToggle->setCheckable(true);
    legendToggle->setArrowType(Qt::RightArrow);
    legendToggle->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    auto* pocketLegend = new QListWidget;
    pocketLegend->setObjectName("pocketLegend");
    pocketLegend->setMaximumHeight(175);
    pocketLegend->setVisible(false);
    pocketLegend->setSelectionMode(QAbstractItemView::NoSelection);
    pocketLegend->setToolTip(tr("Check pocket candidates to display; the color swatch matches the 3D surface."));
    viewerColumn->insertWidget(2, legendToggle);
    viewerColumn->insertWidget(3, pocketLegend);
    connect(legendToggle, &QToolButton::toggled, this, [legendToggle, pocketLegend](bool on) {
        pocketLegend->setVisible(on);
        legendToggle->setArrowType(on ? Qt::DownArrow : Qt::RightArrow);
    });
    auto updateLegend = [pocketWorkspace, pocketLegend, viewer = viewer_]() {
        const QSignalBlocker blocker(pocketLegend);
        pocketLegend->clear();
        const QStringList choices = pocketWorkspace->pocketChoices();
        const QList<int> selected = pocketWorkspace->selectedPocketRows();
        for (int i = 0; i < choices.size(); ++i) {
            auto* item = new QListWidgetItem(choices[i], pocketLegend);
            item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
            item->setCheckState(selected.contains(i) ? Qt::Checked : Qt::Unchecked);
            const QColor color = viewer->pocketDisplayColor(pocketWorkspace->pocketColorSample(i));
            QPixmap swatch(18, 14);
            swatch.fill(color);
            item->setIcon(QIcon(swatch));
            item->setToolTip(choices[i]);
        }
        pocketLegend->setEnabled(!choices.isEmpty());
    };
    connect(pocketLegend, &QListWidget::itemChanged, pocketWorkspace,
            [pocketLegend, pocketWorkspace](QListWidgetItem* item) {
                const int index = pocketLegend->row(item);
                QList<int> selected = pocketWorkspace->selectedPocketRows();
                if (item->checkState() == Qt::Checked && !selected.contains(index))
                    selected.append(index);
                else if (item->checkState() != Qt::Checked)
                    selected.removeAll(index);
                pocketWorkspace->setSelectedPocketRows(selected);
            });
    connect(pocketWorkspace, &PocketWorkspace::pocketChoicesChanged,
            this, updateLegend);
    connect(pocketWorkspace, &PocketWorkspace::pocketSelectedRowsChanged,
            this, [pocketLegend](const QList<int>& selected) {
                const QSignalBlocker blocker(pocketLegend);
                for (int i = 0; i < pocketLegend->count(); ++i)
                    pocketLegend->item(i)->setCheckState(selected.contains(i)
                        ? Qt::Checked : Qt::Unchecked);
            });
    connect(viewer_, &ViewerWidget::styleChanged, this, updateLegend);
    connect(distinctColors, &QCheckBox::toggled, this, updateLegend);

    auto updateQuickPockets = [pocketWorkspace, quickPocket, pocketMenu]() {
        pocketMenu->clear();
        const QStringList choices = pocketWorkspace->pocketChoices();
        const QList<int> selected = pocketWorkspace->selectedPocketRows();
        for (int i = 0; i < choices.size(); ++i) {
            QAction* action = pocketMenu->addAction(choices[i]);
            action->setCheckable(true);
            action->setChecked(selected.contains(i));
            connect(action, &QAction::toggled, pocketWorkspace,
                    [pocketWorkspace, i](bool checked) {
                        QList<int> rows = pocketWorkspace->selectedPocketRows();
                        if (checked && !rows.contains(i)) rows.append(i);
                        else if (!checked) rows.removeAll(i);
                        pocketWorkspace->setSelectedPocketRows(rows);
                    });
        }
        quickPocket->setEnabled(!choices.isEmpty());
        quickPocket->setText(choices.isEmpty()
            ? QObject::tr("No pockets")
            : QObject::tr("Pockets (%1)").arg(selected.size()));
    };
    connect(pocketWorkspace, &PocketWorkspace::pocketChoicesChanged,
            this, updateQuickPockets);
    connect(pocketWorkspace, &PocketWorkspace::pocketSelectedRowsChanged,
            this, [pocketMenu, quickPocket](const QList<int>& selected) {
                const auto actions = pocketMenu->actions();
                for (int i = 0; i < actions.size(); ++i) {
                    const QSignalBlocker blocker(actions[i]);
                    actions[i]->setChecked(selected.contains(i));
                }
                quickPocket->setText(actions.isEmpty()
                    ? QObject::tr("No pockets")
                    : QObject::tr("Pockets (%1)").arg(selected.size()));
            });
    connect(quickOpacity, QOverload<int>::of(&QSpinBox::valueChanged),
            pocketWorkspace, &PocketWorkspace::setPocketOpacityPercent);
    connect(pocketWorkspace, &PocketWorkspace::pocketOpacityPercentChanged,
            this, [quickOpacity](int value) {
                const QSignalBlocker blocker(quickOpacity);
                quickOpacity->setValue(value);
            });
    connect(distinctColors, &QCheckBox::toggled,
            viewer_, &ViewerWidget::setPocketDistinctColors);
    connect(quickVisible, &QCheckBox::toggled,
            pocketWorkspace, &PocketWorkspace::setPocketOverlayVisible);
    connect(pocketWorkspace, &PocketWorkspace::pocketVisibilityChanged,
            this, [quickVisible](bool visible) {
                const QSignalBlocker blocker(quickVisible);
                quickVisible->setChecked(visible);
            });
    connect(quickMaterial, QOverload<int>::of(&QComboBox::currentIndexChanged),
            pocketWorkspace, &PocketWorkspace::setPocketMaterialIndex);
    connect(pocketWorkspace, &PocketWorkspace::pocketMaterialChanged,
            this, [quickMaterial](int index) {
                const QSignalBlocker blocker(quickMaterial);
                quickMaterial->setCurrentIndex(index);
            });
    {
        const QSignalBlocker blocker(quickMaterial);
        quickMaterial->setCurrentIndex(pocketWorkspace->pocketMaterialIndex());
    }
    updateQuickPockets();
    updateLegend();
    modes3d->addTab(viewerArea, tr("Interactions 3D"));
    modes3d->addTab(pocketWorkspace, tr("Pocket Finder 3D"));
    tabs_ = new QTabWidget;
    tabs_->setObjectName(QStringLiteral("resultsTabs"));
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
    auto* posePanel = new QWidget;
    auto* poseLayout = new QVBoxLayout(posePanel);
    auto* poseTools = new QHBoxLayout;
    poseRunButton_ = new QPushButton(tr("Analyze all imported poses"));
    poseExportButton_ = new QPushButton(tr("Export CSV…"));
    poseExportButton_->setEnabled(false);
    poseMatrixExportButton_ = new QPushButton(tr("Export residue matrix…"));
    poseMatrixExportButton_->setEnabled(false);
    poseTools->addWidget(poseRunButton_);
    poseTools->addWidget(poseExportButton_);
    poseTools->addWidget(poseMatrixExportButton_);
    poseLayout->addLayout(poseTools);
    poseTable_ = new QTableWidget(0, 11);
    poseTable_->setHorizontalHeaderLabels({tr("Index"), tr("Pose"), tr("Status"), tr("Favorable"), tr("H-bonds"), tr("Hydrophobic"), tr("Pi"), tr("Water"), tr("Clashes"), tr("Residues"), tr("IFP group")});
    poseTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    poseTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
    poseTable_->setSelectionMode(QAbstractItemView::SingleSelection);
    poseTable_->setSortingEnabled(true);
    poseTable_->horizontalHeader()->setStretchLastSection(true);
    auto* poseSplit = new QSplitter(Qt::Vertical);
    poseSplit->addWidget(poseTable_);
    poseResidueTable_ = new QTableWidget(0, 3);
    poseResidueTable_->setHorizontalHeaderLabels({tr("Residue"), tr("Poses"), tr("Frequency (%)")});
    poseResidueTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    poseResidueTable_->setSortingEnabled(true);
    poseResidueTable_->horizontalHeader()->setStretchLastSection(true);
    poseSplit->addWidget(poseResidueTable_);
    poseSplit->setStretchFactor(0, 3);
    poseSplit->setStretchFactor(1, 2);
    poseLayout->addWidget(poseSplit);
    auto* heatmapLegend = new QLabel(tr("Interaction abbreviations: HB = Hydrogen bond; HY = Hydrophobic contact; PI = Pi-stacking / cation-pi; SB = Salt bridge; WB = Water bridge; XB = Halogen bond; ME = Metal coordination; OT = Other favorable contact. Multiple codes indicate multiple types. Empty = no favorable contact; dash = failed pose."));
    heatmapLegend->setWordWrap(true);
    poseLayout->addWidget(heatmapLegend);
    auto* groupNote = new QLabel(tr("IFP group: poses with identical residue + favorable-interaction-type fingerprints share a group (G1, G2, ...). Groups do not imply identical geometry or binding affinity."));
    groupNote->setWordWrap(true);
    poseLayout->addWidget(groupNote);
    connect(poseRunButton_, &QPushButton::clicked, this, &MainWindow::runPoseComparison);
    connect(poseTable_, &QTableWidget::cellDoubleClicked, this, [this](int row, int) { selectComparedPose(row); });
    connect(poseExportButton_, &QPushButton::clicked, this, [this] {
        const QString file = QFileDialog::getSaveFileName(this, tr("Export pose comparison"), QStringLiteral("pose_comparison.csv"), tr("CSV (*.csv)"));
        if (file.isEmpty()) return;
        std::ofstream out(file.toStdString(), std::ios::binary | std::ios::trunc);
        if (!out) { showError(tr("Export failed"), tr("Cannot open output file")); return; }
        ligand_pose_batch::write(out, poseRows_, "csv");
        if (!out) showError(tr("Export failed"), tr("Could not write comparison table"));
    });
    connect(poseMatrixExportButton_, &QPushButton::clicked, this, [this] {
        const QString file = QFileDialog::getSaveFileName(this, tr("Export conserved residue matrix"), QStringLiteral("pose_residue_matrix.csv"), tr("CSV (*.csv);;JSON (*.json)"));
        if (file.isEmpty()) return;
        const std::string format = file.endsWith(QStringLiteral(".json"), Qt::CaseInsensitive) ? "json" : "csv";
        std::ofstream out(file.toStdString(), std::ios::binary | std::ios::trunc);
        if (!out) { showError(tr("Export failed"), tr("Cannot open output file")); return; }
        ligand_pose_batch::writeResidueMatrix(out, poseRows_, format);
        if (!out) showError(tr("Export failed"), tr("Could not write residue matrix"));
    });
    tabs_->addTab(posePanel, tr("Pose Comparison"));
    tabs_->addTab(log_, tr("Log"));

    auto* right = new QSplitter(Qt::Vertical);
    resultsSplitter_ = right;
    modes3d->setMinimumHeight(0);
    tabs_->setMinimumHeight(0);
    modes3d->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Ignored);
    tabs_->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Ignored);
    right->addWidget(modes3d);
    right->addWidget(tabs_);
    right->setStretchFactor(0, 3);
    right->setStretchFactor(1, 2);

    // Compact navigation: keep primary workflow controls visible without scrolling
    // through thresholds, pocket tools and interaction filters.
    for (QWidget* widget : std::initializer_list<QWidget*>{gStruct, gOpt, gThr, gTmpl, gPocket, runButton_, gTypes})
        lv->removeWidget(widget);
    auto* inputPage = new QWidget;
    auto* inputLayout = new QVBoxLayout(inputPage);
    inputLayout->addWidget(gStruct);
    inputLayout->addWidget(gTmpl);
    inputLayout->addWidget(runButton_);
    inputLayout->addStretch(1);
    auto* settingsPage = new QWidget;
    auto* settingsLayout = new QVBoxLayout(settingsPage);
    settingsLayout->addWidget(gOpt);
    settingsLayout->addWidget(gThr);
    settingsLayout->addStretch(1);
    auto* advancedPage = new QWidget;
    auto* advancedLayout = new QVBoxLayout(advancedPage);
    advancedLayout->addWidget(gPocket);
    advancedLayout->addWidget(gTypes);
    advancedLayout->addStretch(1);
    delete left; // All control widgets have been reparented to their new page layouts.
    auto* controlPages = new QTabWidget;
    controlPages_ = controlPages;
    controlPages->setObjectName("controlPages");
    auto makeScrollPage = [](QWidget* content) {
        auto* scroll = new QScrollArea;
        scroll->setWidgetResizable(true);
        scroll->setFrameShape(QFrame::NoFrame);
        scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        scroll->setWidget(content);
        return scroll;
    };
    controlPages->addTab(makeScrollPage(inputPage), tr("Input"));
    controlPages->addTab(makeScrollPage(settingsPage), tr("Settings"));
    controlPages->addTab(makeScrollPage(advancedPage), tr("Advanced"));
    // Each tab already owns its vertical scroll area: avoid another nested scroll
    // widget, whose size hint previously inflated the window's minimum height.
    controlPages->setMinimumSize(0, 0);
    controlPages->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Ignored);
    auto* split = new QSplitter(Qt::Horizontal);
    mainSplitter_ = split;
    split->addWidget(controlPages);
    split->addWidget(right);
    split->setStretchFactor(1, 1);
    controlPages->setMinimumWidth(215);
    controlPages->setMaximumWidth(340);
    split->setChildrenCollapsible(false);
    right->setChildrenCollapsible(false);
    split->setMinimumHeight(0);
    right->setMinimumHeight(0);
    split->setSizes({265, 1015});
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


void MainWindow::applyCcdPreview(const QJsonObject& results) {
        int updated=0,skipped=0;
        QStringList failedComponents;
        for(auto& ligand:importPreview_.ligands){
            const QString component=QString::fromStdString(ligand.resName);
            const auto entry=results.value(QString::fromStdString(ligand.resName).toUpper()).toObject();
            if(!entry.contains("bonds")){
                ++skipped;failedComponents<<component;
                log(tr("CCD %1: %2").arg(QString::fromStdString(ligand.resName),
                    entry.value("error").toString(tr("No CCD bond data returned"))));
                continue;
            }
            const auto bonds=entry.value("bonds").toArray();
            std::map<QString,int> byName;
            bool collision=false;
            for(int i=0;i<static_cast<int>(ligand.atoms.size());++i){
                const QString name=QString::fromStdString(ligand.atoms[i].name).trimmed();
                if(byName.count(name)){collision=true;break;}
                byName[name]=i;
            }
            if(collision){++skipped;failedComponents<<component;log(tr("CCD %1: duplicate atom names in PDB ligand").arg(QString::fromStdString(ligand.resName)));continue;}
            // Accept CCD chemistry only for a complete, uniquely matched heavy-atom
            // inventory; partial components can assign unsafe connectivity/protonation.
            const auto elements=entry.value("elements").toObject();
            int ccdHeavy=0;
            bool compatible=!elements.isEmpty();
            for(auto it=elements.begin();it!=elements.end();++it){
                const QString element=it.value().toString().toUpper();
                if(element=="H"||element=="D")continue;
                ++ccdHeavy;
                if(!byName.count(it.key()) ||
                   QString::fromStdString(ligand.atoms[byName[it.key()]].element).toUpper()!=element){
                    compatible=false;break;
                }
            }
            // A crystallographic ligand may have unresolved atoms. Preserve
            // CCD bonds whose endpoints are both present and element-validated;
            // allow CCD-derived H counts only for fully matching components.
            bool validObserved=true;
            for(const auto& atom:ligand.atoms){
                const QString name=QString::fromStdString(atom.name).trimmed();
                if(!elements.contains(name) ||
                   elements.value(name).toString().toUpper()!=QString::fromStdString(atom.element).toUpper()){
                    validObserved=false;break;
                }
            }
            if(!validObserved || !compatible){
                ++skipped;failedComponents<<component;
                log(tr("CCD %1: atom identity mismatch; geometric fallback retained.")
                    .arg(QString::fromStdString(ligand.resName)));
                continue;
            }
            const bool complete=ccdHeavy==static_cast<int>(ligand.atoms.size());
            std::vector<ligint::BondRecord> newBonds;
            std::set<std::pair<int,int>> seen;
            bool invalid=false;
            for(const auto& bondValue:bonds){
                const auto fields=bondValue.toArray();
                if(fields.size()!=3){invalid=true;break;}
                const QString a=fields[0].toString(),b=fields[1].toString();
                // CCD hydrogen atoms are frequently absent from deposited structures.
                if(!byName.count(a)||!byName.count(b))continue;
                const int ia=byName[a],ib=byName[b],order=fields[2].toInt();
                const auto key=std::minmax(ia,ib);
                if(ia==ib||order<1||order>4||!seen.insert(key).second){invalid=true;break;}
                newBonds.push_back({ia,ib,order});
            }
            if(invalid||newBonds.empty()){
                ++skipped;failedComponents<<component;log(tr("CCD %1: invalid or empty matched bond table").arg(QString::fromStdString(ligand.resName)));
                continue;
            }
            ligand.bonds=std::move(newBonds);
            ligand.ccdDonorHydrogens.clear();
            const auto donorSites=complete ? entry.value("donor_hydrogens").toObject() : QJsonObject();
            for(auto it=donorSites.begin();it!=donorSites.end();++it){
                const QString name=it.key();
                const int n=it.value().toInt();
                if(byName.count(name) && n>=1 && n<=3){
                    const int idx=byName[name];
                    const auto& element=ligand.atoms[idx].element;
                    if(element=="N"||element=="O")ligand.ccdDonorHydrogens.push_back({idx,n});
                }
            }
            ligand.annotation=complete ? "CCD atom-name matched connectivity" :
                "Partial CCD connectivity (no donor H assignment)";
            if(!complete)log(tr("CCD %1: partial topology applied; estimated hydrogens disabled.")
                .arg(QString::fromStdString(ligand.resName)));
            ++updated;
        }
        if(updated){
            // Keep the selected ligand when rebuilding the import preview.
            const auto index=ligandCombo_->currentIndex();
            int selected=-1;
            if(index>=0&&index<static_cast<int>(shownLigands_.size())){
                const auto& id=shownLigands_[index];
                for(int i=0;i<static_cast<int>(importPreview_.ligands.size());++i){
                    const auto& lig=importPreview_.ligands[i];
                    if(lig.resName==id.res&&lig.chain==id.chain&&lig.resi==id.resi&&lig.icode==id.icode){
                        selected=i;break;
                    }
                }
            }
            viewer_->setImportPreview(importPreview_,selected);
        }
        log(tr("CCD ligand preview: %1 updated, %2 unchanged (details above).").arg(updated).arg(skipped));
        statusBar()->showMessage(tr("CCD: %1 enriched, %2 unchanged; see Log").arg(updated).arg(skipped),12000);
        if(skipped>0){
            failedComponents.removeDuplicates();
            const QString sample=failedComponents.mid(0,6).join(", ") + (failedComponents.size()>6 ? ", …" : "");
            const QString warning=tr("CCD preview could not enrich %1 ligand(s) (%2). Geometric fallback is retained. See the Log for details.").arg(skipped).arg(sample);
            if(interactive_) QMessageBox::warning(this,tr("CCD ligand preview warning"),warning);
        }

}

void MainWindow::startCcdPreview() {
    if(importPreview_.ligands.empty())return;
    ++ccdGeneration_;
    if(ccdProcess_){
        ccdProcess_->disconnect(this);
        if(ccdProcess_->state()!=QProcess::NotRunning){
            ccdProcess_->kill();
            ccdProcess_->waitForFinished(3000);
        }
        ccdProcess_->deleteLater();
        ccdProcess_=nullptr;
    }
    QStringList ids;
    for(const auto& lig:importPreview_.ligands){
        const QString id=QString::fromStdString(lig.resName).toUpper();
        if(!ids.contains(id)) ids<<id;
    }
    const QString script=QStringLiteral(LIGINT_SOURCE_ROOT)+"/tools/ccd_preview.py";
    if(!QFileInfo::exists(script)){ const QString msg=tr("CCD: helper not found: ")+script;log(msg);if(interactive_)QMessageBox::warning(this,tr("CCD ligand preview warning"),msg);return; }
    QString python=QStandardPaths::findExecutable("python3");
    if(python.isEmpty())python=QStandardPaths::findExecutable("python");
    if(python.isEmpty()){const QString msg=tr("CCD: Python 3 not found");log(msg);if(interactive_)QMessageBox::warning(this,tr("CCD ligand preview warning"),msg);return;}
    QProcess* process=new QProcess(this);
    ccdProcess_=process;
    const quint64 generation=ccdGeneration_;
    connect(process,QOverload<int,QProcess::ExitStatus>::of(&QProcess::finished),this,
        [this,process,generation](int code,QProcess::ExitStatus status) {
            const QByteArray output=process->readAllStandardOutput();
            const QString error=QString::fromUtf8(process->readAllStandardError());
            if(ccdProcess_==process) ccdProcess_=nullptr;
            process->deleteLater();
            if(generation!=ccdGeneration_)return; // another PDB was opened
            // Analyze may finish while CCD networking is still in progress.
            // Never replace an analyzed ligand scene with import-preview data.
            if(!result_.ligands.empty() || !importPreview_.listingOnly)return;
            if(status!=QProcess::NormalExit || code!=0){
                const QString msg=tr("CCD automatic preview failed: ")+error;
                log(msg);
                if(interactive_)QMessageBox::warning(this,tr("CCD ligand preview warning"),msg);
                return;
            }
            QJsonParseError parseError{};
            const auto json=QJsonDocument::fromJson(output,&parseError);
            if(!json.isObject()){
                const QString msg=tr("CCD: invalid helper output: ")+parseError.errorString();
                log(msg);
                if(interactive_)QMessageBox::warning(this,tr("CCD ligand preview warning"),msg);
                return;
            }
            applyCcdPreview(json.object());
        });
    connect(process,&QProcess::errorOccurred,this,
        [this,process,generation](QProcess::ProcessError){
            if(generation==ccdGeneration_)log(tr("CCD: ")+process->errorString());
        });
    QPointer<QProcess> safeProcess(process);
    QTimer::singleShot(90000,this,[safeProcess] { if(safeProcess && safeProcess->state()!=QProcess::NotRunning) safeProcess->kill(); });
    QStringList args{script}; args<<ids;
    process->start(python,args);
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
    ++ccdGeneration_;
    if(ccdProcess_){
        ccdProcess_->disconnect(this);
        if(ccdProcess_->state()!=QProcess::NotRunning){
            ccdProcess_->kill();
            ccdProcess_->waitForFinished(3000);
        }
        ccdProcess_->deleteLater();
        ccdProcess_=nullptr;
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
    graphicsOptions_->updateAvailability();
    runButton_->setEnabled(ligandCombo_->count() > 0);
    exportAct_->setEnabled(false);
    imageAct_->setEnabled(viewer_->hasData());
    startCcdPreview();
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
        ligandCombo_->addItem(tr("Crystal: %1 %2 %3%4  (%5 atoms)")
                                  .arg(QString::fromStdString(l.res))
                                  .arg(QString::fromStdString(l.chain))
                                  .arg(l.resi)
                                  .arg(l.icode == ' ' ? QString() : QString(QChar(l.icode)))
                                  .arg(l.atoms));
    for (int i = 0; i < static_cast<int>(shownLigands_.size()); ++i)
        ligandCombo_->setItemData(i, tr("Crystallographic ligand (HETATM) from the receptor structure."), Qt::ToolTipRole);
    for (size_t ti = 0; ti < templates_.size(); ++ti) {
        ligint::LigandId ext;
        ext.res = "EXT"; ext.chain = "__EXTERNAL__"; ext.resi = static_cast<int>(ti);
        for (const auto& atom : templates_[ti].at) if (atom.el != "H" && atom.el != "D") ++ext.atoms;
        shownLigands_.push_back(ext);
        ligandCombo_->addItem(tr("External pose %1: %2 (%3 heavy atoms)")
                                  .arg(ti + 1).arg(QString::fromStdString(templates_[ti].name)).arg(ext.atoms));
        ligandCombo_->setItemData(ligandCombo_->count() - 1,
                                  tr("Imported SDF/MOL2 molecule %1. Coordinates must share the receptor frame.").arg(ti + 1),
                                  Qt::ToolTipRole);
    }
    runButton_->setEnabled(ligandCombo_->count() > 0 && !path_.isEmpty());
}

int MainWindow::ligandCount() const { return ligandCombo_->count(); }
void MainWindow::selectLigand(int index) { ligandCombo_->setCurrentIndex(index); }

void MainWindow::runPoseComparison() {
    if (path_.isEmpty() || templates_.empty() || poseWatcher_->isRunning()) {
        if (path_.isEmpty() || templates_.empty()) showError(tr("Pose comparison"), tr("Load a receptor and one or more SDF/MOL2 ligand structures first."));
        return;
    }
    poseRunButton_->setEnabled(false);
    poseExportButton_->setEnabled(false);
    poseMatrixExportButton_->setEnabled(false);
    poseRows_.clear();
    poseTable_->setRowCount(0);
    poseResidueTable_->setRowCount(0);
    ligint::Options opt = buildOptions();
    opt.protein3d = false;
    opt.excludePdbLigand = excludePdbLigand_->text().trimmed().toStdString();
    const std::string receptor = path_.toStdString();
    statusBar()->showMessage(tr("Analyzing imported poses…"));
    poseWatcher_->setFuture(QtConcurrent::run([receptor, opt]() {
        std::vector<ligand_pose_batch::Row> rows;
        rows.reserve(opt.templates.size());
        for (size_t i = 0; i < opt.templates.size(); ++i) rows.push_back(ligand_pose_batch::analyze(receptor, opt, i));
        return rows;
    }));
}

void MainWindow::selectComparedPose(int row) {
    if (row < 0 || row >= poseTable_->rowCount()) return;
    const auto* indexItem = poseTable_->item(row, 0);
    if (!indexItem) return;
    const int index = indexItem->data(Qt::EditRole).toInt() - 1;
    if (index < 0 || index >= static_cast<int>(poseRows_.size()) || !poseRows_[index].error.empty()) return;
    for (int i = 0; i < static_cast<int>(shownLigands_.size()); ++i) {
        if (shownLigands_[i].chain == "__EXTERNAL__" && shownLigands_[i].resi == index) {
            ligandCombo_->setCurrentIndex(i);
            analyze(false);
            return;
        }
    }
}

bool MainWindow::addTemplate(const QString& path) {
    std::vector<ligint::Tmpl> t;
    std::vector<std::string> warnings;
    const bool mol2 = path.endsWith(".mol2", Qt::CaseInsensitive);
    bool ok = mol2 ? ligint::loadMol2(path.toStdString(), t, warnings) : ligint::loadSdf(path.toStdString(), t, warnings);
    for (const auto& w : warnings) log(QString::fromStdString(w));
    if (!ok || t.empty()) {
        showError(tr("Invalid ligand structure"), tr("No readable molecule in %1").arg(path));
        return false;
    }
    templates_.insert(templates_.end(), t.begin(), t.end());
    templatePaths_ << QFileInfo(path).absoluteFilePath();
    templateNames_ << QFileInfo(path).fileName() + tr(" (%1 molecules)").arg(t.size());
    templateLabel_->setText(templateNames_.join("\n"));
    log(tr("Ligand structure loaded: %1").arg(path));
    refreshLigandCombo();
    return true;
}

void MainWindow::onAddTemplate() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getOpenFileName(this, tr("Ligand structure (SDF or MOL2)"), s.value("lastDir").toString(),
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
        if (l.chain == "__EXTERNAL__") opt.externalLigandIndex = l.resi;
        else {
            opt.ligName = l.res;
            opt.chain = l.chain;
            opt.resi = l.resi;
            opt.icode = l.icode;
        }
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
    if (opt.externalLigandIndex >= 0) opt.excludePdbLigand = excludePdbLigand_->text().trimmed().toStdString();
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
        if (l.chain != "__EXTERNAL__") {
        a << QString::fromStdString(l.res) << QString::fromStdString(l.chain);
        // the insertion code is written only when it is needed to tell this residue from another one with the same number
        bool ambiguous = false;
        for (const auto& o : allLigands_)
            if (o.res == l.res && o.chain == l.chain && o.resi == l.resi && o.icode != l.icode) ambiguous = true;
        QString number = QString::number(l.resi);
        if (ambiguous) number += l.icode == ' ' ? QStringLiteral("-") : QString(QChar(l.icode));
        a << number;
        }
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
    const int idx = ligandCombo_->currentIndex();
    if (idx >= 0 && idx < static_cast<int>(shownLigands_.size()) && shownLigands_[idx].chain == "__EXTERNAL__") {
        a << "--external-index" << QString::number(shownLigands_[idx].resi + 1);
        if (!excludePdbLigand_->text().trimmed().isEmpty())
            a << "--exclude-pdb-ligand" << excludePdbLigand_->text().trimmed();
    }
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
