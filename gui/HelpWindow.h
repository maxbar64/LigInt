// Window with the user manual: the HTML pages are embedded in the program (docs/manual/html, built by tools/build_manual.py).
#pragma once

#include <QList>
#include <QMainWindow>
#include <QString>
#include <QStringList>
#include <QUrl>

class QAction;
class QLineEdit;
class QListWidget;
class QTextBrowser;
class QTreeWidget;
class QTreeWidgetItem;

namespace gui {

class HelpWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit HelpWindow(QWidget* parent = nullptr);

    // Shows a page of the manual: "gui.html" or "gui.html#the-log". False if the page does not exist.
    bool showPage(const QString& target);
    QString currentPage() const;                    // file name of the page shown, for example "gui.html"
    int pageCount() const { return pages_.size(); }
    QStringList pageFiles() const;

    // Looks for the words of `query` in all the pages: returns the files of the pages that contain all of them, best first.
    // The results are also shown in the Search tab.
    QStringList search(const QString& query);

    QTextBrowser* browser() const { return browser_; }
    QTreeWidget* contents() const { return contents_; }
    QListWidget* results() const { return results_; }

    // The place of the manual that describes the results tab with this title ("Interactions", "Residues", "Fingerprint", "Log").
    static QString topicForTab(const QString& tabTitle);
    // True if the manual has this page ("gui.html") and, when given, this section ("gui.html#the-log").
    static bool topicExists(const QString& target);
    // Folder with the HTML files of the manual installed next to the program (share/doc/ligint/manual), or an empty string.
    static QString installedManualDir();

private slots:
    void onAnchorClicked(const QUrl& url);
    void onSourceChanged(const QUrl& url);
    void onContentsActivated(QTreeWidgetItem* item);
    void onResultActivated();
    void findInPage(bool backward = false);
    void openInBrowser();

private:
    struct Heading { int level; QString text, id; };
    struct Page { QString file, title; QList<Heading> headings; QString text; };    // text: plain text, built on the first search

    void loadToc();
    void buildSearchIndex();

    QList<Page> pages_;
    QString title_;
    QTextBrowser* browser_ = nullptr;
    QTreeWidget* contents_ = nullptr;
    QLineEdit* searchEdit_ = nullptr;
    QListWidget* results_ = nullptr;
    QLineEdit* findEdit_ = nullptr;
    QAction* backAct_ = nullptr;
    QAction* forwardAct_ = nullptr;
    QAction* browserAct_ = nullptr;
    bool indexBuilt_ = false;
};

}  // namespace gui
