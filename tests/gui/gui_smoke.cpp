// Test di fumo della GUI: apre un PDB, analizza, controlla tabelle, filtri e rendering (piattaforma "offscreen").
#include <QApplication>
#include <QImage>
#include <QSet>
#include <iostream>

#include "HitsModel.h"
#include "MainWindow.h"
#include "ViewerWidget.h"

static int failures = 0;
#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK fallito: " #cond "\n"; \
            failures++;                                                                \
        }                                                                              \
    } while (0)

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    if (argc < 2) { std::cerr << "uso: ligint_gui_smoke file.pdb\n"; return 2; }

    gui::MainWindow w;
    w.setInteractive(false);                          // nessuna finestra di dialogo
    w.resize(1200, 800);
    w.show();

    CHECK(!w.openStructure("/non/esiste.pdb"));       // errore gestito senza crash
    CHECK(w.openStructure(argv[1]));
    CHECK(w.ligandCount() > 0);
    CHECK(!w.viewer()->hasData());                    // nulla da mostrare prima dell'analisi

    w.selectLigand(0);
    w.analyze(true);
    CHECK(!w.result().ligands.empty());
    CHECK(w.hitsModel()->rowCount() > 0);
    CHECK(w.hitsProxy()->rowCount() == w.hitsModel()->rowCount());
    CHECK(w.viewer()->hasData());
    CHECK(w.logText().contains("Annotazione ligando"));

    // il filtro per tipo agisce sulla tabella
    w.hitsProxy()->setVisibleTypes(QSet<QString>());
    CHECK(w.hitsProxy()->rowCount() == 0);
    w.hitsProxy()->setVisibleTypes(QSet<QString>{"Legame a idrogeno"});
    CHECK(w.hitsProxy()->rowCount() > 0 && w.hitsProxy()->rowCount() <= w.hitsModel()->rowCount());

    // rendering: l'immagine non e' uniforme e viene salvata come artefatto
    w.viewer()->setVisibleTypes(QSet<QString>{"Legame a idrogeno", "Coordinazione metallica", "Contatto idrofobico"});
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

    if (failures) { std::cerr << failures << " controlli falliti\n"; return 1; }
    std::cout << "ligint_gui_smoke: ok (" << w.hitsModel()->rowCount() << " interazioni, immagine salvata)\n";
    return 0;
}
