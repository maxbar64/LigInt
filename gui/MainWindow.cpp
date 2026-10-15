#include "MainWindow.h"

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QDragEnterEvent>
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
#include <QPushButton>
#include <QScrollArea>
#include <QSettings>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QTableView>
#include <QTableWidget>
#include <QToolBar>
#include <QUrl>
#include <QVBoxLayout>
#include <QtConcurrent/QtConcurrent>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>

#include "CommandLine.h"
#include "HitsModel.h"
#include "Style.h"
#include "SvgView.h"
#include "ViewerWidget.h"

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
    watcher_ = new QFutureWatcher<AnalysisOutcome>(this);
    connect(watcher_, &QFutureWatcher<AnalysisOutcome>::finished, this, &MainWindow::onFinished);
    buildUi();
    resize(1280, 820);
    statusBar()->showMessage(tr("Open a structure file, PDB or mmCIF (File → Open, or drag it here)"));
}

MainWindow::~MainWindow() {
    if (watcher_ && watcher_->isRunning()) watcher_->waitForFinished();
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
    exportAct_->setEnabled(false);
    imageAct_->setEnabled(false);
    fpAct_->setEnabled(false);
    file->addSeparator();
    makeAction(file, tr("&Quit"), &QWidget::close, QKeySequence::Quit);
    auto* help = menuBar()->addMenu(tr("&Help"));
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
    piAngle_ = spin(10.0, 60.0, 30.0, 1.0, 0);
    piAngle_->setSuffix(QStringLiteral("°"));
    form->addRow(tr("H-bond, D···A max"), hbMax_);
    form->addRow(tr("Long H-bond, up to"), hbLong_);
    form->addRow(tr("Salt bridge (atoms)"), saltMax_);
    form->addRow(tr("Salt bridge (centroids)"), saltCentroid_);
    form->addRow(tr("Hydrophobic max"), hydMax_);
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
    connect(pkLoad, &QPushButton::clicked, this, &MainWindow::onLoadPocket);
    connect(pkClear, &QPushButton::clicked, this, [this] { clearPocket(); });
    pkRow->addWidget(pkLoad);
    pkRow->addWidget(pkClear);
    pkv->addLayout(pkRow);
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
    tabs_->addTab(log_, tr("Log"));

    auto* right = new QSplitter(Qt::Vertical);
    right->addWidget(viewer_);
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
    piAngle_->setObjectName("thr_piAngle");

    // ---- view toolbar
    auto* tb = addToolBar(tr("View"));
    tb->setMovable(false);
    auto addToggle = [&](const QString& text, bool on, auto setter) {
        QAction* a = tb->addAction(text);
        a->setCheckable(true);
        a->setChecked(on);
        connect(a, &QAction::toggled, viewer_, setter);
        return a;
    };
    addToggle(tr("Labels"), true, &ViewerWidget::setShowLabels);
    addToggle(tr("Residues"), true, &ViewerWidget::setShowResidues);
    addToggle(tr("Distances"), true, &ViewerWidget::setShowDistances);
    addToggle(tr("Legend"), true, &ViewerWidget::setShowLegend);
    tb->addSeparator();
    QAction* reset = tb->addAction(tr("Reset view"));
    connect(reset, &QAction::triggered, viewer_, &ViewerWidget::resetView);
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
        listing = ligint::analyzeFile(path.toStdString(), ligint::Options());
    } catch (const std::exception& e) {
        showError(tr("Cannot open the file"), QString::fromUtf8(e.what()));
        return false;
    }
    path_ = path;
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
    model_->setHits({});
    viewer_->clear();
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
    runButton_->setEnabled(ligandCombo_->count() > 0);
    exportAct_->setEnabled(false);
    imageAct_->setEnabled(false);
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
    const ligint::Options opt = buildOptions();
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
    viewer_->setLigand(lig);
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
