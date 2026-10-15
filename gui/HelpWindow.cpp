#include "HelpWindow.h"

#include <QApplication>
#include <QDesktopServices>
#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QKeySequence>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QRegularExpression>
#include <QSplitter>
#include <QStatusBar>
#include <QTabWidget>
#include <QTextBrowser>
#include <QTextCursor>
#include <QTextDocument>
#include <QToolBar>
#include <QTreeWidget>
#include <QVBoxLayout>
#include <algorithm>

// Qt resources of a static library must be initialized from outside any namespace
static void initManualResources() { Q_INIT_RESOURCE(manual); }

namespace gui {

namespace {

QByteArray resourceBytes(const QString& file) {
    QFile f(QStringLiteral(":/manual/") + file);
    return f.open(QIODevice::ReadOnly) ? f.readAll() : QByteArray();
}

QString htmlToText(const QByteArray& html) {
    QTextDocument doc;
    doc.setHtml(QString::fromUtf8(html));
    return doc.toPlainText();
}

}  // namespace

HelpWindow::HelpWindow(QWidget* parent) : QMainWindow(parent) {
    initManualResources();
    loadToc();
    setWindowTitle(title_);
    resize(1050, 720);

    // ---- toolbar
    QToolBar* tb = addToolBar(tr("Navigation"));
    tb->setMovable(false);
    backAct_ = tb->addAction(tr("Back"));
    backAct_->setShortcut(QKeySequence::Back);
    forwardAct_ = tb->addAction(tr("Forward"));
    forwardAct_->setShortcut(QKeySequence::Forward);
    QAction* home = tb->addAction(tr("Home"));
    home->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Home));
    tb->addSeparator();
    QAction* zin = tb->addAction(tr("Zoom in"));
    zin->setShortcut(QKeySequence::ZoomIn);
    QAction* zout = tb->addAction(tr("Zoom out"));
    zout->setShortcut(QKeySequence::ZoomOut);
    tb->addSeparator();
    browserAct_ = tb->addAction(tr("Open in browser"));
    browserAct_->setToolTip(tr("Open the HTML files of the manual installed with the program in your web browser"));
    browserAct_->setEnabled(!installedManualDir().isEmpty());

    // ---- left: contents and search
    auto* side = new QTabWidget;
    contents_ = new QTreeWidget;
    contents_->setHeaderHidden(true);
    side->addTab(contents_, tr("Contents"));
    auto* searchPage = new QWidget;
    auto* sl = new QVBoxLayout(searchPage);
    searchEdit_ = new QLineEdit;
    searchEdit_->setPlaceholderText(tr("Search the manual..."));
    searchEdit_->setClearButtonEnabled(true);
    results_ = new QListWidget;
    results_->setWordWrap(true);
    sl->addWidget(searchEdit_);
    sl->addWidget(results_);
    side->addTab(searchPage, tr("Search"));

    // ---- right: page and find bar
    browser_ = new QTextBrowser;
    browser_->setOpenLinks(false);
    browser_->setOpenExternalLinks(false);
    browser_->setSearchPaths({QStringLiteral(":/manual")});
    findEdit_ = new QLineEdit;
    findEdit_->setPlaceholderText(tr("Find in this page (Ctrl+F)"));
    findEdit_->setClearButtonEnabled(true);
    auto* right = new QWidget;
    auto* rl = new QVBoxLayout(right);
    rl->setContentsMargins(0, 0, 0, 0);
    rl->addWidget(browser_, 1);
    rl->addWidget(findEdit_);

    auto* split = new QSplitter(Qt::Horizontal);
    split->addWidget(side);
    split->addWidget(right);
    split->setStretchFactor(0, 0);
    split->setStretchFactor(1, 1);
    split->setSizes({280, 770});
    setCentralWidget(split);

    // ---- contents tree
    for (const Page& p : pages_) {
        auto* top = new QTreeWidgetItem(contents_, QStringList(p.title));
        top->setData(0, Qt::UserRole, p.file);
        for (const Heading& h : p.headings) {
            if (h.level > 2) continue;
            auto* it = new QTreeWidgetItem(top, QStringList(h.text));
            it->setData(0, Qt::UserRole, p.file + QLatin1Char('#') + h.id);
        }
    }

    // ---- connections
    connect(backAct_, &QAction::triggered, browser_, &QTextBrowser::backward);
    connect(forwardAct_, &QAction::triggered, browser_, &QTextBrowser::forward);
    connect(home, &QAction::triggered, this, [this] { showPage(QStringLiteral("index.html")); });
    connect(zin, &QAction::triggered, this, [this] { browser_->zoomIn(1); });
    connect(zout, &QAction::triggered, this, [this] { browser_->zoomOut(1); });
    connect(browserAct_, &QAction::triggered, this, &HelpWindow::openInBrowser);
    connect(browser_, &QTextBrowser::backwardAvailable, backAct_, &QAction::setEnabled);
    connect(browser_, &QTextBrowser::forwardAvailable, forwardAct_, &QAction::setEnabled);
    connect(browser_, &QTextBrowser::anchorClicked, this, &HelpWindow::onAnchorClicked);
    connect(browser_, &QTextBrowser::sourceChanged, this, &HelpWindow::onSourceChanged);
    connect(contents_, &QTreeWidget::itemActivated, this, &HelpWindow::onContentsActivated);
    connect(contents_, &QTreeWidget::itemClicked, this, &HelpWindow::onContentsActivated);
    connect(searchEdit_, &QLineEdit::returnPressed, this, [this] { search(searchEdit_->text()); });
    connect(results_, &QListWidget::itemActivated, this, &HelpWindow::onResultActivated);
    connect(results_, &QListWidget::itemClicked, this, &HelpWindow::onResultActivated);
    connect(findEdit_, &QLineEdit::returnPressed, this, [this] { findInPage(false); });
    connect(findEdit_, &QLineEdit::textEdited, this, [this] {
        browser_->moveCursor(QTextCursor::Start);
        findInPage(false);
    });
    auto* find = new QAction(this);
    find->setShortcut(QKeySequence::Find);
    addAction(find);
    connect(find, &QAction::triggered, this, [this] { findEdit_->setFocus(); findEdit_->selectAll(); });
    auto* findPrev = new QAction(this);
    findPrev->setShortcut(QKeySequence::FindPrevious);
    addAction(findPrev);
    connect(findPrev, &QAction::triggered, this, [this] { findInPage(true); });
    backAct_->setEnabled(false);
    forwardAct_->setEnabled(false);

    showPage(QStringLiteral("index.html"));
}

void HelpWindow::loadToc() {
    const QJsonObject root = QJsonDocument::fromJson(resourceBytes(QStringLiteral("toc.json"))).object();
    title_ = root.value(QStringLiteral("title")).toString(QStringLiteral("ligint manual"));
    for (const QJsonValue& pv : root.value(QStringLiteral("pages")).toArray()) {
        const QJsonObject po = pv.toObject();
        Page p;
        p.file = po.value(QStringLiteral("file")).toString();
        p.title = po.value(QStringLiteral("title")).toString();
        for (const QJsonValue& hv : po.value(QStringLiteral("headings")).toArray()) {
            const QJsonObject ho = hv.toObject();
            p.headings.append({ho.value(QStringLiteral("level")).toInt(), ho.value(QStringLiteral("text")).toString(),
                               ho.value(QStringLiteral("id")).toString()});
        }
        if (!p.file.isEmpty()) pages_.append(p);
    }
}

QStringList HelpWindow::pageFiles() const {
    QStringList files;
    for (const Page& p : pages_) files << p.file;
    return files;
}

bool HelpWindow::topicExists(const QString& target) {
    initManualResources();
    const QString file = target.section(QLatin1Char('#'), 0, 0);
    const QString frag = target.section(QLatin1Char('#'), 1);
    const QByteArray bytes = resourceBytes(file);
    if (bytes.isEmpty()) return false;
    return frag.isEmpty() || bytes.contains(("id=\"" + frag + "\"").toUtf8());
}

QString HelpWindow::topicForTab(const QString& tab) {
    if (tab == QLatin1String("Interactions")) return QStringLiteral("gui.html#the-interactions-table");
    if (tab == QLatin1String("Residues")) return QStringLiteral("gui.html#the-residues-table");
    if (tab == QLatin1String("Fingerprint")) return QStringLiteral("gui.html#the-fingerprint-tab");
    if (tab == QLatin1String("Log")) return QStringLiteral("gui.html#the-log");
    return QStringLiteral("index.html");
}

QString HelpWindow::installedManualDir() {
    const QDir app(QApplication::applicationDirPath());
    for (const char* rel : {"../share/doc/ligint/manual", "share/doc/ligint/manual"}) {
        const QString dir = QDir::cleanPath(app.filePath(QString::fromLatin1(rel)));
        if (QFile::exists(dir + QStringLiteral("/index.html"))) return dir;
    }
    return QString();
}

bool HelpWindow::showPage(const QString& target) {
    const QString file = target.section(QLatin1Char('#'), 0, 0);
    const QString frag = target.section(QLatin1Char('#'), 1);
    if (resourceBytes(file).isEmpty()) return false;
    QUrl url(QStringLiteral("qrc:/manual/") + file);
    if (!frag.isEmpty()) url.setFragment(frag);
    browser_->setSource(url);
    return true;
}

QString HelpWindow::currentPage() const {
    return browser_->source().path().section(QLatin1Char('/'), -1);
}

void HelpWindow::onAnchorClicked(const QUrl& url) {
    if (url.scheme() == QLatin1String("http") || url.scheme() == QLatin1String("https") || url.scheme() == QLatin1String("mailto")) {
        QDesktopServices::openUrl(url);
        return;
    }
    const QUrl abs = browser_->source().resolved(url);
    showPage(abs.path().section(QLatin1Char('/'), -1) + (abs.hasFragment() ? QStringLiteral("#") + abs.fragment() : QString()));
}

void HelpWindow::onSourceChanged(const QUrl& url) {
    const QString file = url.path().section(QLatin1Char('/'), -1);
    contents_->blockSignals(true);
    for (int i = 0; i < contents_->topLevelItemCount(); i++) {
        QTreeWidgetItem* it = contents_->topLevelItem(i);
        if (it->data(0, Qt::UserRole).toString() == file) {
            contents_->setCurrentItem(it);
            setWindowTitle(it->text(0) + QStringLiteral(" - ") + title_);
        }
    }
    contents_->blockSignals(false);
}

void HelpWindow::onContentsActivated(QTreeWidgetItem* item) {
    if (item) showPage(item->data(0, Qt::UserRole).toString());
}

void HelpWindow::buildSearchIndex() {
    if (indexBuilt_) return;
    for (Page& p : pages_) p.text = htmlToText(resourceBytes(p.file));
    indexBuilt_ = true;
}

QStringList HelpWindow::search(const QString& query) {
    buildSearchIndex();
    results_->clear();
    const QStringList words = query.toLower().split(QRegularExpression(QStringLiteral("\\s+")), Qt::SkipEmptyParts);
    QStringList files;
    if (words.isEmpty()) return files;
    struct Hit { int score, order; const Page* page; };
    QList<Hit> hits;
    int order = 0;
    for (const Page& p : pages_) {
        const QString lower = p.text.toLower(), title = p.title.toLower();
        int score = 0;
        bool all = true;
        for (const QString& w : words) {
            const int n = lower.count(w);
            if (n == 0 && !title.contains(w)) { all = false; break; }
            score += n + (title.contains(w) ? 20 : 0);
        }
        if (all) hits.append({score, order, &p});
        order++;
    }
    std::stable_sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.score > b.score; });
    for (const Hit& h : hits) {
        const QString lower = h.page->text.toLower();
        int at = lower.indexOf(words.first());
        QString snippet;
        if (at >= 0) {
            const int from = std::max(0, at - 60);
            snippet = h.page->text.mid(from, 150).simplified();
            snippet = (from > 0 ? QStringLiteral("...") : QString()) + snippet + QStringLiteral("...");
        }
        auto* item = new QListWidgetItem(h.page->title + (snippet.isEmpty() ? QString() : QStringLiteral("\n") + snippet), results_);
        item->setData(Qt::UserRole, h.page->file);
        item->setData(Qt::UserRole + 1, words.first());
        files << h.page->file;
    }
    if (files.isEmpty()) results_->addItem(tr("No page contains all these words."));
    return files;
}

void HelpWindow::onResultActivated() {
    QListWidgetItem* it = results_->currentItem();
    if (!it || it->data(Qt::UserRole).toString().isEmpty()) return;
    showPage(it->data(Qt::UserRole).toString());
    browser_->moveCursor(QTextCursor::Start);
    browser_->find(it->data(Qt::UserRole + 1).toString());
}

void HelpWindow::findInPage(bool backward) {
    const QString text = findEdit_->text();
    if (text.isEmpty()) return;
    QTextDocument::FindFlags flags;
    if (backward) flags |= QTextDocument::FindBackward;
    if (!browser_->find(text, flags)) {                                  // wrap around
        browser_->moveCursor(backward ? QTextCursor::End : QTextCursor::Start);
        browser_->find(text, flags);
    }
}

void HelpWindow::openInBrowser() {
    const QString dir = installedManualDir();
    if (dir.isEmpty()) return;
    const QString page = currentPage();
    QUrl url = QUrl::fromLocalFile(dir + QLatin1Char('/') + (page.isEmpty() ? QStringLiteral("index.html") : page));
    if (browser_->source().hasFragment()) url.setFragment(browser_->source().fragment());
    QDesktopServices::openUrl(url);
}

}  // namespace gui
