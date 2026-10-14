// Strumento per la documentazione: apre un PDB, analizza il ligando piu' grande e salva uno screenshot.
//   ligint-screenshot file.pdb out.png [--viewer] [--plip-like]     (usare QT_QPA_PLATFORM=offscreen senza schermo)
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
    if (pdb.isEmpty() || out.isEmpty()) { std::cerr << "uso: ligint-screenshot file.pdb out.png [--viewer]\n"; return 2; }
    gui::MainWindow w;
    w.setInteractive(false);
    w.resize(1280, 800);
    w.show();
    if (!w.openStructure(pdb) || w.ligandCount() == 0) { std::cerr << "nessun ligando\n"; return 1; }
    w.selectLigand(0);
    w.analyze(true);
    if (w.result().ligands.empty()) { std::cerr << "analisi fallita\n"; return 1; }
    QImage img = viewerOnly ? w.viewer()->renderImage(QSize(1100, 760)) : w.grab().toImage();
    if (!img.save(out)) { std::cerr << "salvataggio fallito\n"; return 1; }
    std::cout << "salvato " << out.toStdString() << "\n";
    return 0;
}
