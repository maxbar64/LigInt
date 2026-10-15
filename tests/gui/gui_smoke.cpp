// GUI smoke test: opens a PDB, analyzes, checks tables, filters and rendering ("offscreen" platform).
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QClipboard>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileInfo>
#include <QAction>
#include <QProcess>
#include <QImage>
#include <QSet>
#include <QTableWidget>
#include <QTemporaryDir>
#include <iostream>

#include "CommandLine.h"
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


#ifndef Q_OS_WIN
// runs a command line through /bin/sh (so that the quoting is really exercised); returns stdout
static QString runShell(const QString& cmd, bool* ok) {
    QProcess p;
    p.start("/bin/sh", {"-c", cmd});
    const bool done = p.waitForFinished(60000);
    if (ok) *ok = done && p.exitStatus() == QProcess::NormalExit && p.exitCode() == 0;
    return QString::fromUtf8(p.readAllStandardOutput());
}

// the command printed by the GUI, with `ligint` replaced by the real binary and --csv appended
static QString withBinary(const QString& cmd, const QString& binary, const QString& extra) {
    return gui::shellQuote(binary) + cmd.mid(QString("ligint").size()) + extra;
}

static QStringList csvKeys(const QString& csv) {
    QStringList keys;
    const QStringList lines = csv.split('\n', Qt::SkipEmptyParts);
    for (int i = 1; i < lines.size(); i++) {
        const QStringList f = lines[i].split(',');
        if (f.size() >= 7) keys << (f[1] + "|" + f[3] + "|" + f[4] + "|" + f[5] + "|" + f[6]);
    }
    keys.sort();
    return keys;
}

static QStringList hitKeys(const ligint::LigandResult& lig) {
    QStringList keys;
    for (const auto& h : lig.hits)
        keys << (QString::fromStdString(h.type) + "|" + QString::fromStdString(h.ligandAtom) + "|" + QString::fromStdString(h.proteinAtom) + "|" +
                 QString::number(h.distance, 'f', 2) + "|" + QString::fromStdString(h.note));
    keys.sort();
    return keys;
}

// runs the last command of the window with the real CLI and compares it with what the window shows
static void checkEquivalence(gui::MainWindow& w, const QString& binary, const char* what, bool allowEmpty = false) {
    bool ok = false;
    const QString out = runShell(withBinary(w.lastCommandLine(), binary, " --csv"), &ok);
    const QStringList a = csvKeys(out), b = hitKeys(w.result().ligands.front());
    if (!ok || (a.isEmpty() && !allowEmpty) || a != b) {
        std::cerr << "command line not equivalent (" << what << "): " << w.lastCommandLine().toStdString() << "\n  cli " << a.size()
                  << " rows, gui " << b.size() << " rows, ok=" << ok << "\n";
        failures++;
    }
}
#endif

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


#ifndef Q_OS_WIN
    // ---- command line equivalent to the GUI settings (needs the real CLI binary as 2nd argument)
    if (argc > 2) {
        const QString binary = argv[2];
        const QString dataDir = QFileInfo(argv[1]).absolutePath();
        const QString absPdb = QFileInfo(argv[1]).absoluteFilePath();
        auto box = [&](const char* name) { return w.findChild<QCheckBox*>(name); };
        auto spin = [&](const char* name) { return w.findChild<QDoubleSpinBox*>(name); };
        CHECK(box("opt_water") && box("opt_plip") && spin("thr_hbMax"));
        auto resetUi = [&] {
            box("opt_plip")->setChecked(false);
            box("opt_water")->setChecked(true);
            box("opt_merge")->setChecked(true);
            box("opt_backbone")->setChecked(false);
            box("opt_metalDonors")->setChecked(false);
            spin("thr_hbMax")->setValue(3.5); spin("thr_hbLong")->setValue(0); spin("thr_saltMax")->setValue(4.0);
            spin("thr_saltCentroid")->setValue(0); spin("thr_hydMax")->setValue(4.0); spin("thr_piAngle")->setValue(30);
        };
        CHECK(w.openStructure(argv[1]));
        w.selectLigand(0);
        resetUi();
        CHECK(w.lastCommandLine().isEmpty());                                   // nothing before an analysis
        const QString base = "ligint " + gui::shellQuote(absPdb) + " AZM A 262";
        CHECK(w.commandLine() == base);                                         // defaults: no options at all
        w.analyze(true);
        CHECK(w.lastCommandLine() == base);
        CHECK(w.logText().contains("Command line: " + base));
        checkEquivalence(w, binary, "defaults");

        resetUi();                                                              // preset + no water
        box("opt_plip")->setChecked(true);
        box("opt_water")->setChecked(false);
        w.analyze(true);
        const QString c2 = w.lastCommandLine();
        CHECK(c2.contains("--no-water") && c2.contains("--hb-long 4.1") && c2.contains("--salt-centroid 5.5") && c2.contains("--hal-max 4"));
        CHECK(!c2.contains("--plip-like"));                                     // explicit options: valid even if the preset changes
        checkEquivalence(w, binary, "plip-like, no water");

        resetUi();                                                              // many settings
        box("opt_backbone")->setChecked(true);
        box("opt_merge")->setChecked(false);
        box("opt_metalDonors")->setChecked(true);
        spin("thr_hbMax")->setValue(3.8); spin("thr_hydMax")->setValue(4.5); spin("thr_piAngle")->setValue(40); spin("thr_saltMax")->setValue(4.5);
        w.analyze(true);
        const QString c3 = w.lastCommandLine();
        for (const char* opt : {"--no-merge", "--backbone-hyd", "--metal-donors", "--hb-max 3.8", "--hyd-max 4.5", "--pi-angle 40", "--salt-max 4.5"})
            CHECK(c3.contains(opt));
        checkEquivalence(w, binary, "many settings");

        resetUi();                                                              // back to the defaults
        CHECK(w.commandLine() == base);

        // the copy action puts the last command in the clipboard
        for (QAction* a : w.findChildren<QAction*>())
            if (a->text().contains("command line") && !a->text().contains("ingerprint")) { CHECK(a->isEnabled()); a->trigger(); }
        CHECK(QApplication::clipboard()->text() == w.lastCommandLine());

        // fingerprint command: same figure as the window, with the pocket definition currently loaded
        CHECK(w.fingerprintCommandLine("csv").contains("--fingerprint csv") && !w.fingerprintCommandLine("csv").contains("--pocket"));
        CHECK(w.loadPocket(pocketPath));
        const QString fpCmd = w.fingerprintCommandLine("svg");
        CHECK(fpCmd.contains("--fingerprint svg --pocket"));
        bool ok = false;
        const QString svg = runShell(withBinary(fpCmd, binary, ""), &ok);
        CHECK(ok && svg.toStdString() == ligint::fingerprintSvg(w.fingerprint()));
        w.clearPocket();

        // a template adds its path to the command; the CLI result matches the window
        CHECK(w.openStructure(dataDir + "/1efy.pdb"));
        w.selectLigand(0);
        resetUi();
        CHECK(w.addTemplate(dataDir + "/1efy.sdf"));
        w.analyze(true);
        const QString c4 = w.lastCommandLine();
        CHECK(c4.contains("--sdf") && c4.contains("1efy.sdf") && c4.contains(" BZC "));
        CHECK(w.logText().contains("template"));
        checkEquivalence(w, binary, "SDF template");
        CHECK(w.openStructure(argv[1]));                                        // templates stay for the next structure
        CHECK(w.commandLine().contains("--sdf"));

        // two residues that differ only in the insertion code: the code appears in the command and selects exactly that one
        {
            const QString tricky = dir.filePath("insertion.pdb");
            QFile f(tricky);
            CHECK(f.open(QIODevice::WriteOnly));
            auto line = [](int serial, const char* name, const char* res, char chain, int resi, char icode, double x, double y, double z, const char* el) {
                return QString::asprintf("HETATM%5d %-4s %3s %c%4d%c   %8.3f%8.3f%8.3f  1.00 20.00          %2s\n", serial, name, res, chain, resi, icode, x, y, z, el).toUtf8();
            };
            f.write(line(1, " C1 ", "LIG", 'A', 5, ' ', 0, 0, 0, "C"));
            f.write(line(2, " O1 ", "LIG", 'A', 5, ' ', 1.43, 0, 0, "O"));
            f.write(line(3, " C1 ", "LIG", 'A', 5, 'A', 20, 0, 0, "C"));
            f.write(line(4, " O1 ", "LIG", 'A', 5, 'A', 21.43, 0, 0, "O"));
            f.write("ATOM      5  OD1 ASN A  10      23.700   2.420   0.000  1.00 20.00           O\nEND\n");
            f.close();
            gui::MainWindow w2;
            w2.setInteractive(false);
            CHECK(w2.openStructure(tricky) && w2.ligandCount() == 2);
            for (int i = 0; i < w2.ligandCount(); i++) {
                w2.selectLigand(i);
                w2.analyze(true);
                const QString c = w2.lastCommandLine();
                CHECK(c.endsWith(" 5A") || c.endsWith(" 5-"));
                checkEquivalence(w2, binary, "insertion code", true);        // the residue without a code has no contacts
                if (c.endsWith(" 5A")) CHECK(!w2.result().ligands.front().hits.empty());
                CHECK(w2.result().ligands.size() == 1);
            }
        }


        // alternate conformations: the combo lists those of the file and the choice goes into the command line
        {
            const QString altPdb = dir.filePath("altloc.pdb");
            QFile f(altPdb);
            CHECK(f.open(QIODevice::WriteOnly));
            f.write("HETATM    1  C1  LIG L   1       0.000   0.000   0.000  1.00 20.00           C\n"
                    "HETATM    2  O1  LIG L   1       1.430   0.000   0.000  1.00 20.00           O\n"
                    "ATOM      3  OD1AASN A  10       9.000   0.000   0.000  0.50 20.00           O\n"
                    "ATOM      4  OD1BASN A  10       2.150   2.420   0.000  0.50 20.00           O\n"
                    "END\n");
            f.close();
            gui::MainWindow w3;
            w3.setInteractive(false);
            QComboBox* combo = w3.findChild<QComboBox*>("opt_altloc");
            CHECK(combo != nullptr && !combo->isEnabled());                       // no alternate conformations before opening a file
            CHECK(w3.openStructure(argv[1]) && !combo->isEnabled() && combo->count() == 1);   // 1AZM has none
            CHECK(w3.openStructure(altPdb) && combo->isEnabled() && combo->count() == 3);
            CHECK(combo->itemText(1) == "A" && combo->itemText(2) == "B");
            w3.selectLigand(0);
            w3.analyze(true);
            CHECK(!w3.lastCommandLine().contains("--altloc") && w3.result().ligands.front().hits.empty());
            checkEquivalence(w3, binary, "altloc default", true);
            combo->setCurrentIndex(2);                                             // B
            w3.analyze(true);
            CHECK(w3.lastCommandLine().contains("--altloc B") && !w3.result().ligands.front().hits.empty());
            checkEquivalence(w3, binary, "altloc B");
            CHECK(w3.logText().contains("Alternate conformation B"));
            CHECK(w3.openStructure(altPdb) && combo->currentText() == "B");        // the choice survives reopening a file that has it
            CHECK(w3.openStructure(argv[1]) && combo->currentIndex() == 0);        // and is reset when the file has none
        }

        // shell quoting: every argument comes back unchanged through a real shell
        const QStringList tricky = {"plain", "a b", "it's", "$HOME", "`id`", "*", "é", "a\"b", "x;y", "-5", " ", "back\\slash"};
        QString cmd = "printf '%s\\n'";
        for (const QString& t : tricky) cmd += " " + gui::shellQuote(t);
        bool q = false;
        const QStringList back = runShell(cmd, &q).split('\n');
        CHECK(q);
        for (int i = 0; i < tricky.size(); i++) CHECK(i < back.size() && back[i] == tricky[i]);
        CHECK(gui::shellQuote("") == "''" && gui::shellQuote("plain-1.2/x") == "plain-1.2/x" && gui::shellQuote("a b") == "'a b'");
    }
#endif

    if (failures) { std::cerr << failures << " checks failed\n"; return 1; }
    std::cout << "ligint_gui_smoke: ok (" << w.hitsModel()->rowCount() << " interactions, image saved)\n";
    return 0;
}
