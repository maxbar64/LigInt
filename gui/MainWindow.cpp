#include "MainWindow.h"

#include <QAction>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHeaderView>
#include <QLabel>
#include <QListWidget>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QPlainTextEdit>
#include <QPushButton>
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
#include <fstream>
#include <map>
#include <sstream>

#include "HitsModel.h"
#include "Style.h"
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
    setWindowTitle(tr("ligint — interazioni ligando-proteina"));
    setAcceptDrops(true);
    model_ = new HitsModel(this);
    proxy_ = new HitsFilterProxy(this);
    proxy_->setSourceModel(model_);
    watcher_ = new QFutureWatcher<AnalysisOutcome>(this);
    connect(watcher_, &QFutureWatcher<AnalysisOutcome>::finished, this, &MainWindow::onFinished);
    buildUi();
    resize(1280, 820);
    statusBar()->showMessage(tr("Apri un file PDB (File → Apri, oppure trascinalo qui)"));
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
    makeAction(file, tr("&Apri PDB…"), &MainWindow::onOpen, QKeySequence::Open);
    exportAct_ = makeAction(file, tr("Esporta interazioni (&CSV)…"), &MainWindow::onExportCsv);
    imageAct_ = makeAction(file, tr("&Salva immagine…"), &MainWindow::onSaveImage);
    exportAct_->setEnabled(false);
    imageAct_->setEnabled(false);
    file->addSeparator();
    makeAction(file, tr("&Esci"), &QWidget::close, QKeySequence::Quit);
    auto* help = menuBar()->addMenu(tr("&Aiuto"));
    makeAction(help, tr("&Informazioni su ligint…"), &MainWindow::onAbout);

    // ---- pannello di controllo (sinistra)
    auto* left = new QWidget;
    auto* lv = new QVBoxLayout(left);

    auto* gStruct = new QGroupBox(tr("Struttura"));
    auto* sv = new QVBoxLayout(gStruct);
    fileLabel_ = new QLabel(tr("(nessun file)"));
    fileLabel_->setWordWrap(true);
    auto* openBtn = new QPushButton(tr("Apri PDB…"));
    connect(openBtn, &QPushButton::clicked, this, &MainWindow::onOpen);
    ligandCombo_ = new QComboBox;
    ligandCombo_->setToolTip(tr("Ligandi (residui HETATM) presenti nel file"));
    showSmall_ = new QCheckBox(tr("Mostra anche ioni e molecole piccole"));
    connect(showSmall_, &QCheckBox::toggled, this, [this] { refreshLigandCombo(); });
    sv->addWidget(openBtn);
    sv->addWidget(fileLabel_);
    sv->addWidget(new QLabel(tr("Ligando:")));
    sv->addWidget(ligandCombo_);
    sv->addWidget(showSmall_);
    lv->addWidget(gStruct);

    auto* gOpt = new QGroupBox(tr("Opzioni"));
    auto* ov = new QVBoxLayout(gOpt);
    water_ = new QCheckBox(tr("Considera le acque (ponti d'acqua)"));
    water_->setChecked(true);
    merge_ = new QCheckBox(tr("Unisci residui legati (glicani, dipeptidi)"));
    merge_->setChecked(true);
    plipLike_ = new QCheckBox(tr("Preset «alla PLIP»"));
    plipLike_->setToolTip(tr("H-bond fino a 4,1 Å se la geometria è possibile, ponti salini per centroidi a 5,5 Å, alogeni a 4,0 Å"));
    backbone_ = new QCheckBox(tr("Contatti apolari con la catena principale"));
    metalDonors_ = new QCheckBox(tr("O/N coordinati a metalli possono donare H"));
    for (auto* c : {water_, merge_, plipLike_, backbone_, metalDonors_}) ov->addWidget(c);
    connect(plipLike_, &QCheckBox::toggled, this, &MainWindow::onPresetToggled);
    lv->addWidget(gOpt);

    auto* gThr = new QGroupBox(tr("Soglie"));
    auto* form = new QFormLayout(gThr);
    hbMax_ = spin(2.5, 6.0, 3.5);
    hbLong_ = spin(0.0, 6.0, 0.0);
    hbLong_->setSpecialValueText(tr("spento"));
    saltMax_ = spin(2.5, 6.0, 4.0);
    saltCentroid_ = spin(0.0, 8.0, 0.0);
    saltCentroid_->setSpecialValueText(tr("spento"));
    hydMax_ = spin(3.3, 6.0, 4.0);
    piAngle_ = spin(10.0, 60.0, 30.0, 1.0, 0);
    piAngle_->setSuffix(QStringLiteral("°"));
    form->addRow(tr("H-bond, D···A max"), hbMax_);
    form->addRow(tr("H-bond lungo, fino a"), hbLong_);
    form->addRow(tr("Ponte salino (atomi)"), saltMax_);
    form->addRow(tr("Ponte salino (centroidi)"), saltCentroid_);
    form->addRow(tr("Idrofobico max"), hydMax_);
    form->addRow(tr("π-stacking, angolo max"), piAngle_);
    lv->addWidget(gThr);

    auto* gTmpl = new QGroupBox(tr("Template chimici (SDF / MOL2)"));
    auto* tv = new QVBoxLayout(gTmpl);
    templateLabel_ = new QLabel(tr("nessuno: chimica dedotta dalla geometria"));
    templateLabel_->setWordWrap(true);
    auto* tbtn = new QPushButton(tr("Aggiungi template…"));
    connect(tbtn, &QPushButton::clicked, this, &MainWindow::onAddTemplate);
    tv->addWidget(tbtn);
    tv->addWidget(templateLabel_);
    lv->addWidget(gTmpl);

    runButton_ = new QPushButton(tr("Analizza"));
    runButton_->setEnabled(false);
    runButton_->setDefault(true);
    connect(runButton_, &QPushButton::clicked, this, [this] { analyze(false); });
    lv->addWidget(runButton_);

    auto* gTypes = new QGroupBox(tr("Tipi di interazione"));
    auto* yv = new QVBoxLayout(gTypes);
    typeList_ = new QListWidget;
    typeList_->setMinimumHeight(120);
    connect(typeList_, &QListWidget::itemChanged, this, [this] { onTypeFilterChanged(); });
    yv->addWidget(typeList_);
    lv->addWidget(gTypes, 1);

    // ---- vista e tabelle (destra)
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
    residueTable_->setHorizontalHeaderLabels({tr("Residuo"), tr("Contatti"), tr("Tipi")});
    residueTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    residueTable_->horizontalHeader()->setStretchLastSection(true);
    residueTable_->verticalHeader()->setVisible(false);
    log_ = new QPlainTextEdit;
    log_->setReadOnly(true);
    tabs_->addTab(table_, tr("Interazioni"));
    tabs_->addTab(residueTable_, tr("Residui"));
    tabs_->addTab(log_, tr("Registro"));

    auto* right = new QSplitter(Qt::Vertical);
    right->addWidget(viewer_);
    right->addWidget(tabs_);
    right->setStretchFactor(0, 3);
    right->setStretchFactor(1, 2);

    auto* split = new QSplitter(Qt::Horizontal);
    split->addWidget(left);
    split->addWidget(right);
    split->setStretchFactor(1, 1);
    left->setMaximumWidth(360);
    setCentralWidget(split);

    // ---- barra degli strumenti della vista
    auto* tb = addToolBar(tr("Vista"));
    tb->setMovable(false);
    auto addToggle = [&](const QString& text, bool on, auto setter) {
        QAction* a = tb->addAction(text);
        a->setCheckable(true);
        a->setChecked(on);
        connect(a, &QAction::toggled, viewer_, setter);
        return a;
    };
    addToggle(tr("Etichette"), true, &ViewerWidget::setShowLabels);
    addToggle(tr("Residui"), true, &ViewerWidget::setShowResidues);
    addToggle(tr("Distanze"), true, &ViewerWidget::setShowDistances);
    addToggle(tr("Legenda"), true, &ViewerWidget::setShowLegend);
    tb->addSeparator();
    QAction* reset = tb->addAction(tr("Ripristina vista"));
    connect(reset, &QAction::triggered, viewer_, &ViewerWidget::resetView);
}

// ------------------------------------------------------------------------------------------------
// apertura del file e scelta del ligando
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
        showError(tr("Impossibile aprire il file"), QString::fromUtf8(e.what()));
        return false;
    }
    path_ = path;
    allLigands_ = listing.hetResidues;
    result_ = ligint::Result();
    model_->setHits({});
    viewer_->clear();
    residueTable_->setRowCount(0);
    typeList_->clear();
    fileLabel_->setText(QFileInfo(path).fileName());
    setWindowTitle(tr("ligint — %1").arg(QFileInfo(path).fileName()));
    for (const auto& w : listing.warnings) log(QString::fromStdString(w));
    log(tr("Aperto %1: %2 residui HETATM").arg(path).arg(allLigands_.size()));
    refreshLigandCombo();
    runButton_->setEnabled(ligandCombo_->count() > 0);
    exportAct_->setEnabled(false);
    imageAct_->setEnabled(false);
    if (ligandCombo_->count() == 0) statusBar()->showMessage(tr("Nessun ligando HETATM in questo file"));
    else statusBar()->showMessage(tr("%1 ligandi: scegli e premi «Analizza»").arg(ligandCombo_->count()));
    return true;
}

void MainWindow::refreshLigandCombo() {
    shownLigands_.clear();
    for (const auto& l : allLigands_)
        if (showSmall_->isChecked() || l.atoms >= 6) shownLigands_.push_back(l);
    if (shownLigands_.empty()) shownLigands_ = allLigands_;       // solo ioni o molecole piccole
    std::stable_sort(shownLigands_.begin(), shownLigands_.end(), [](const ligint::LigandId& a, const ligint::LigandId& b) { return a.atoms > b.atoms; });
    ligandCombo_->clear();
    for (const auto& l : shownLigands_)
        ligandCombo_->addItem(tr("%1 %2 %3%4  (%5 atomi)")
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
        showError(tr("Template non valido"), tr("Nessuna molecola leggibile in %1").arg(path));
        return false;
    }
    templates_.insert(templates_.end(), t.begin(), t.end());
    templateNames_ << QFileInfo(path).fileName() + tr(" (%1 molecole)").arg(t.size());
    templateLabel_->setText(templateNames_.join("\n"));
    log(tr("Template caricato: %1").arg(path));
    return true;
}

void MainWindow::onAddTemplate() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getOpenFileName(this, tr("Template SDF o MOL2"), s.value("lastDir").toString(),
                                             tr("Molecole (*.sdf *.mol2);;Tutti i file (*)"));
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
    }
    opt.useWater = water_->isChecked();
    opt.mergeLinked = merge_->isChecked();
    opt.backboneHyd = backbone_->isChecked();
    opt.metalDonors = metalDonors_->isChecked();
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
// analisi
// ------------------------------------------------------------------------------------------------
void MainWindow::analyze(bool wait) {
    if (path_.isEmpty() || ligandCombo_->currentIndex() < 0) return;
    if (watcher_->isRunning()) return;
    const ligint::Options opt = buildOptions();
    const QString path = path_;
    runButton_->setEnabled(false);
    statusBar()->showMessage(tr("Analisi in corso…"));
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
        showError(tr("Analisi fallita"), o.error);
        emit analysisFinished(false);
        return;
    }
    result_ = o.result;
    for (const auto& w : result_.warnings) log(QString::fromStdString(w));
    for (const auto& s : result_.info) log(QString::fromStdString(s));
    if (result_.ligands.empty()) {
        model_->setHits({});
        viewer_->clear();
        statusBar()->showMessage(tr("Nessun ligando analizzato"));
        emit analysisFinished(false);
        return;
    }
    const ligint::LigandResult& lig = result_.ligands.front();
    for (const auto& s : lig.info) log(QString::fromStdString(s));
    if (!lig.metalSphere.empty()) log(tr("Sfera di coordinazione dei metalli:\n") + QString::fromStdString(lig.metalSphere));
    model_->setHits(lig.hits);
    viewer_->setLigand(lig);
    fillTypeList();
    fillResidueTable();
    table_->resizeColumnsToContents();
    table_->sortByColumn(HitsModel::Type, Qt::AscendingOrder);
    exportAct_->setEnabled(true);
    imageAct_->setEnabled(true);
    statusBar()->showMessage(tr("%1: %2 interazioni, %3 residui").arg(QString::fromStdString(lig.id)).arg(lig.hits.size()).arg(residueTable_->rowCount()));
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

void MainWindow::onRowChanged() {
    const QModelIndex cur = table_->selectionModel()->currentIndex();
    if (!cur.isValid()) { viewer_->setHighlight(-1); return; }
    viewer_->setHighlight(proxy_->mapToSource(cur).data(HitsModel::HitIndexRole).toInt());
}

// ------------------------------------------------------------------------------------------------
// file e dialoghi
// ------------------------------------------------------------------------------------------------
void MainWindow::onOpen() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getOpenFileName(this, tr("Apri struttura"), s.value("lastDir").toString(),
                                             tr("Strutture PDB (*.pdb *.ent);;Tutti i file (*)"));
    if (f.isEmpty()) return;
    s.setValue("lastDir", QFileInfo(f).absolutePath());
    openStructure(f);
}

void MainWindow::onExportCsv() {
    if (result_.ligands.empty()) return;
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getSaveFileName(this, tr("Esporta interazioni"), s.value("lastDir").toString() + "/interazioni.csv", tr("CSV (*.csv)"));
    if (f.isEmpty()) return;
    std::ostringstream out, info;
    ligint::printCsv(out, info, result_);
    std::ofstream file(f.toStdString());
    file << out.str();
    statusBar()->showMessage(tr("Esportato %1").arg(f), 5000);
}

void MainWindow::onSaveImage() {
    QSettings s("ligint", "ligint-gui");
    QString f = QFileDialog::getSaveFileName(this, tr("Salva immagine"), s.value("lastDir").toString() + "/ligint.png", tr("Immagine PNG (*.png)"));
    if (f.isEmpty()) return;
    if (!viewer_->renderImage(QSize(1600, 1200)).save(f)) showError(tr("Salvataggio fallito"), f);
    else statusBar()->showMessage(tr("Salvata %1").arg(f), 5000);
}

void MainWindow::onAbout() {
    QMessageBox::about(this, tr("Informazioni su ligint"),
                       tr("<b>ligint %1</b><br>Interazioni non covalenti ligando-proteina da file PDB.<br><br>"
                          "Le interazioni sono calcolate dalla geometria; la chimica del ligando è dedotta "
                          "(o letta da SDF/MOL2). Il confronto con PLIP è descritto nella documentazione.")
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
