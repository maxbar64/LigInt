// Documentation tool: opens a PDB, analyzes the largest ligand and saves a screenshot.
//   ligint-screenshot file.pdb out.png [--viewer] [--plip-like]     (use QT_QPA_PLATFORM=offscreen without a display)
#include <QApplication>
#include <QImage>
#include <QString>
#include <iostream>

#include "MainWindow.h"
#include "ViewerWidget.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    QString pdb, out;
    bool viewerOnly = false;
    for (int i = 1; i < argc; i++) {
        QString a = argv[i];
        if (a == "--viewer") viewerOnly = true;
        else if (pdb.isEmpty()) pdb = a;
        else if (out.isEmpty()) out = a;
    }
    if (pdb.isEmpty() || out.isEmpty()) { std::cerr << "usage: ligint-screenshot file.pdb out.png [--viewer]\n"; return 2; }
    gui::MainWindow w;
    w.setInteractive(false);
    w.resize(1280, 800);
    w.show();
    if (!w.openStructure(pdb) || w.ligandCount() == 0) { std::cerr << "no ligand\n"; return 1; }
    w.selectLigand(0);
    w.analyze(true);
    if (w.result().ligands.empty()) { std::cerr << "analysis failed\n"; return 1; }
    QImage img = viewerOnly ? w.viewer()->renderImage(QSize(1100, 760)) : w.grab().toImage();
    if (!img.save(out)) { std::cerr << "saving failed\n"; return 1; }
    std::cout << "saved " << out.toStdString() << "\n";
    return 0;
}
