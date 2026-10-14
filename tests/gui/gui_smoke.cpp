// GUI smoke test: opens a PDB, analyzes, checks tables, filters and rendering ("offscreen" platform).
#include <QApplication>
#include <QFile>
#include <QImage>
#include <QSet>
#include <QTableWidget>
#include <QTemporaryDir>
#include <iostream>

#include "HitsModel.h"
#include "MainWindow.h"
#include "ViewerWidget.h"

static int failures = 0;
#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            failures++;                                                                \
        }                                                                              \
    } while (0)

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    if (argc < 2) { std::cerr << "usage: ligint_gui_smoke file.pdb\n"; return 2; }

    gui::MainWindow w;
    w.setInteractive(false);                          // no dialog boxes
    w.resize(1200, 800);
    w.show();

    CHECK(!w.openStructure("/no/such.pdb"));       // error handled without crashing
    CHECK(w.openStructure(argv[1]));
    CHECK(w.ligandCount() > 0);
    CHECK(!w.viewer()->hasData());                    // nothing to show before the analysis

    w.selectLigand(0);
    w.analyze(true);
    CHECK(!w.result().ligands.empty());
    CHECK(w.hitsModel()->rowCount() > 0);
    CHECK(w.hitsProxy()->rowCount() == w.hitsModel()->rowCount());
    CHECK(w.viewer()->hasData());
    CHECK(w.logText().contains("Ligand annotation"));

    // the type filter acts on the table
    w.hitsProxy()->setVisibleTypes(QSet<QString>());
    CHECK(w.hitsProxy()->rowCount() == 0);
    w.hitsProxy()->setVisibleTypes(QSet<QString>{QString(ligint::type::HydrogenBond)});
    CHECK(w.hitsProxy()->rowCount() > 0 && w.hitsProxy()->rowCount() <= w.hitsModel()->rowCount());

    // rendering: the image is not uniform and is saved as an artifact
    w.viewer()->setVisibleTypes(QSet<QString>{QString(ligint::type::HydrogenBond), QString(ligint::type::MetalCoordination), QString(ligint::type::Hydrophobic)});
    QImage img = w.viewer()->renderImage(QSize(900, 700));
    CHECK(img.width() == 900 && img.height() == 700);
    QSet<QRgb> colors;
    for (int y = 0; y < img.height(); y += 3)
        for (int x = 0; x < img.width(); x += 3) colors.insert(img.pixel(x, y));
    CHECK(colors.size() > 20);
    CHECK(img.save("ligint_gui_smoke.png"));
    QImage window = w.grab().toImage();
    CHECK(!window.isNull() && window.width() > 600);
    window.save("ligint_gui_window.png");

    // fingerprint: table, pocket definition and exports
    CHECK(!w.fingerprint().entries.empty());
    CHECK(w.fingerprintTable()->rowCount() == static_cast<int>(w.fingerprint().entries.size()));
    const size_t contacts = w.fingerprint().entries.size();
    QTemporaryDir dir;
    CHECK(dir.isValid());
    const QString pocketPath = dir.filePath("pocket.tsv");
    {
        QFile f(pocketPath);
        CHECK(f.open(QIODevice::WriteOnly));
        f.write("@name Smoke pocket\nA\t91\t1\tr1\tPHE\nA\t94\t2\tr1\tHIS\nA\t199\t3\tr2\tTHR\n");
    }
    CHECK(!w.loadPocket("/no/such/pocket.tsv"));
    CHECK(w.loadPocket(pocketPath));
    CHECK(w.fingerprint().pocketName == "Smoke pocket");
    // 3 positions + the contact residues that are not in the definition
    CHECK(w.fingerprint().entries.size() == 3 + (contacts - 2));
    CHECK(w.fingerprint().entries[1].residue.find("HIS 94") == 0 && !w.fingerprint().entries[1].any());
    CHECK(w.fingerprintTable()->rowCount() == static_cast<int>(w.fingerprint().entries.size()));
    for (const char* name : {"fp.svg", "fp.csv", "fp.json"}) {
        const QString path = dir.filePath(name);
        CHECK(w.exportFingerprint(path));
        QFile f(path);
        CHECK(f.exists() && f.size() > 100);
    }
    CHECK(!w.exportFingerprint(dir.filePath("fp.txt")));          // unknown format
    w.clearPocket();
    CHECK(w.fingerprint().entries.size() == contacts);

    if (failures) { std::cerr << failures << " checks failed\n"; return 1; }
    std::cout << "ligint_gui_smoke: ok (" << w.hitsModel()->rowCount() << " interactions, image saved)\n";
    return 0;
}
