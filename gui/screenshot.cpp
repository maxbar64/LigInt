// Documentation tool: opens a PDB, analyzes the largest ligand and saves a screenshot.
//   ligint-screenshot file.pdb out.png [--viewer] [--pocket FILE] [--tab NAME]   (use QT_QPA_PLATFORM=offscreen without a display)
//   NAME is the title of a result tab (Interactions, Residues, Fingerprint, Log)
#include <QApplication>
#include <QImage>
#include <QString>
#include <QTabWidget>
#include <QTextBrowser>
#include "HelpWindow.h"
#include <iostream>

#include "MainWindow.h"
#include "ViewerWidget.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    QString pdb, out, pocket, tab, manualPage, manualSearch;
    bool viewerOnly = false;
    for (int i = 1; i < argc; i++) {
        QString a = argv[i];
        if (a == "--viewer") viewerOnly = true;
        else if (a == "--pocket" && i + 1 < argc) pocket = argv[++i];
        else if (a == "--tab" && i + 1 < argc) tab = argv[++i];
        else if (a == "--manual" && i + 1 < argc) manualPage = argv[++i];
        else if (a == "--search" && i + 1 < argc) manualSearch = argv[++i];
        else if (pdb.isEmpty()) pdb = a;
        else if (out.isEmpty()) out = a;
    }
    if (!manualPage.isEmpty() || !manualSearch.isEmpty()) {      // the manual window: ligint-screenshot x out.png --manual gui.html [--search words]
        if (out.isEmpty() && !pdb.isEmpty()) out = pdb;
        gui::HelpWindow help;
        help.resize(1100, 760);
        help.show();
        QApplication::processEvents();
        if (!manualSearch.isEmpty()) {
            const QStringList found = help.search(manualSearch);
            std::cout << "search '" << manualSearch.toStdString() << "': " << found.size() << " pages\n";
        }
        if (!manualPage.isEmpty() && !help.showPage(manualPage)) { std::cerr << "no such page\n"; return 1; }
        QApplication::processEvents();
        if (!help.grab().toImage().save(out)) { std::cerr << "saving failed\n"; return 1; }
        std::cout << "saved " << out.toStdString() << "\n";
        return 0;
    }
    if (pdb.isEmpty() || out.isEmpty()) { std::cerr << "usage: ligint-screenshot file.pdb out.png [--viewer]\n"; return 2; }
    gui::MainWindow w;
    w.setInteractive(false);
    w.resize(1280, 800);
    w.show();
    if (!w.openStructure(pdb) || w.ligandCount() == 0) { std::cerr << "no ligand\n"; return 1; }
    w.selectLigand(0);
    if (!pocket.isEmpty() && !w.loadPocket(pocket)) { std::cerr << "cannot load the pocket definition\n"; return 1; }
    w.analyze(true);
    if (!tab.isEmpty()) {
        QTabWidget* tabs = w.findChild<QTabWidget*>();
        for (int i = 0; tabs && i < tabs->count(); i++)
            if (tabs->tabText(i) == tab) tabs->setCurrentIndex(i);
    }
    if (w.result().ligands.empty()) { std::cerr << "analysis failed\n"; return 1; }
    QImage img = viewerOnly ? w.viewer()->renderImage(QSize(1100, 760)) : w.grab().toImage();
    if (!img.save(out)) { std::cerr << "saving failed\n"; return 1; }
    std::cout << "saved " << out.toStdString() << "\n";
    return 0;
}
