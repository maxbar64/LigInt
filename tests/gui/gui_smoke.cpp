// GUI smoke test: opens a PDB, analyzes, checks tables, filters and rendering ("offscreen" platform).
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QClipboard>
#include <QDir>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFile>
#include <QFileInfo>
#include <QAction>
#include <QListWidget>
#include <QProcess>
#include <QRegularExpression>
#include <QTabWidget>
#include <QTextBrowser>
#include <QTreeWidget>
#include <QImage>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSet>
#include <QTableWidget>
#include <QTemporaryDir>
#include <iostream>

#include "CommandLine.h"
#include "HelpWindow.h"
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
            if (a->text().contains("command line") && !a->text().contains("ingerprint") && !a->text().contains("matrix")) { CHECK(a->isEnabled()); a->trigger(); }
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

    // the manual: embedded pages, navigation, search, links and context help
    {
        gui::HelpWindow h;
        CHECK(h.pageCount() >= 12 && h.contents()->topLevelItemCount() == h.pageCount());
        CHECK(h.currentPage() == "index.html");
        CHECK(h.showPage("gui.html#the-log") && h.currentPage() == "gui.html" && h.browser()->source().fragment() == "the-log");
        CHECK(!h.showPage("no-such-page.html"));
        // every internal link and every image of the embedded copy of every page
        const QRegularExpression hrefRe("href=\"([^\"]+)\""), imgRe("<img [^>]*src=\"([^\"]+)\"");
        int links = 0, images = 0;
        for (const QString& file : h.pageFiles()) {
            QFile f(":/manual/" + file);
            CHECK(f.open(QIODevice::ReadOnly));
            const QString html = QString::fromUtf8(f.readAll());
            for (auto it = hrefRe.globalMatch(html); it.hasNext();) {
                const QString href = it.next().captured(1);
                if (href.startsWith("http") || href.startsWith("mailto:")) continue;
                CHECK(gui::HelpWindow::topicExists(href.startsWith('#') ? file + href : href));
                links++;
            }
            for (auto im = imgRe.globalMatch(html); im.hasNext();) { CHECK(QFile::exists(":/manual/" + im.next().captured(1))); images++; }
        }
        CHECK(links > 50 && images >= 4);
        // history and links
        h.showPage("index.html");
        h.showPage("cli.html");
        h.browser()->backward();
        CHECK(h.currentPage() == "index.html");
        h.browser()->forward();
        CHECK(h.currentPage() == "cli.html");
        emit h.browser()->anchorClicked(QUrl("gui.html#exporting-results"));          // a click on a link of the page
        CHECK(h.currentPage() == "gui.html" && h.browser()->source().fragment() == "exporting-results");
        emit h.browser()->anchorClicked(QUrl("troubleshooting.html"));
        CHECK(h.currentPage() == "troubleshooting.html");
        CHECK(h.browser()->find("gzip"));                                              // find in the page
        // search: all the words must be in the page; option names and rare words can be found
        const QStringList r = h.search("fingerprint pocket definition");
        CHECK(!r.isEmpty() && h.results()->count() == r.size() && r.contains("fingerprint.html"));
        CHECK(h.search("zzzqqqxxx").isEmpty() && h.search("fingerprint zzzqqqxxx").isEmpty() && h.search("   ").isEmpty());
        CHECK(h.search("hb-long").contains("cli.html"));
        // context help: every tab has a target that exists, also an unknown tab
        for (const char* tab : {"Interactions", "Residues", "Fingerprint", "Log", "Unknown"})
            CHECK(gui::HelpWindow::topicExists(gui::HelpWindow::topicForTab(tab)));
        CHECK(!gui::HelpWindow::topicExists("gui.html#there-is-no-such-section"));
        // F1 in the main window opens the manual at the page of the current tab
        gui::MainWindow m;
        m.setInteractive(false);
        CHECK(m.helpWindow() == nullptr);
        QAction* f1 = nullptr;
        for (QAction* a : m.findChildren<QAction*>())
            if (a->shortcut() == QKeySequence(QKeySequence::HelpContents)) f1 = a;
        CHECK(f1 != nullptr);
        QTabWidget* tabs = m.findChild<QTabWidget*>();
        CHECK(tabs != nullptr);
        for (int i = 0; tabs && i < tabs->count(); i++) {
            tabs->setCurrentIndex(i);
            if (f1) f1->trigger();
            CHECK(m.helpWindow() != nullptr && m.helpWindow()->currentPage() == "gui.html");
            CHECK(m.helpWindow()->browser()->source().fragment() == gui::HelpWindow::topicForTab(tabs->tabText(i)).section('#', 1));
        }
        m.showManual("cli.html");
        CHECK(m.helpWindow()->currentPage() == "cli.html");
        bool hasMenu = false;
        for (QAction* a : m.findChildren<QAction*>()) if (a->text().contains("User manual")) hasMenu = true;
        CHECK(hasMenu);
    }

    // RCSB downloads (here from the local cache: no network) and KLIFS pocket definitions (from a local mapping), in the background
    {
        const QString dataDir = QFileInfo(argv[1]).absolutePath();
        const QString cache = dir.filePath("rcsb");
        QDir().mkpath(cache);
        CHECK(QFile::copy(dataDir + "/1AZM.pdb", cache + "/1AZM.pdb"));
        gui::MainWindow m;
        m.setInteractive(false);
        CHECK(!m.downloadRcsb("../x", "pdb", cache, true) && m.logText().contains("valid four-character"));     // invalid IDs never reach the library
        CHECK(!m.downloadRcsb("1AZM", "xyz", cache, true) && !m.downloadRcsb("", "pdb", cache, true));
        CHECK(m.downloadRcsb("1azm", "pdb", cache, true));                                                       // waits: the structure is open on return
        CHECK(m.ligandCount() >= 1 && m.logText().contains("RCSB: structure 1AZM saved to"));
        auto waitFor = [](gui::MainWindow& w, int ms) {
            QElapsedTimer t;
            t.start();
            while (w.taskRunning() && t.elapsed() < ms) QApplication::processEvents(QEventLoop::AllEvents, 20);
            QApplication::processEvents();
            return !w.taskRunning();
        };
        gui::MainWindow bg;                                                                                      // the same in a background thread
        bg.setInteractive(false);
        CHECK(bg.downloadRcsb("1AZM", "pdb", cache));
        CHECK(waitFor(bg, 10000) && bg.ligandCount() >= 1);
        gui::MainWindow fail;                                                                                    // a structure that is not in the cache needs the network, which is
        fail.setInteractive(false);                                                                              // not used here: use a format that is never cached and a proxy that does not exist
        qputenv("https_proxy", "http://127.0.0.1:9");
        qputenv("HTTPS_PROXY", "http://127.0.0.1:9");
        CHECK(fail.downloadRcsb("1AZM", "cif", cache, true));                                                    // the request is accepted, the transfer fails
        CHECK(fail.logText().contains("RCSB download failed") && fail.ligandCount() == 0 && !QFile::exists(cache + "/1AZM.cif"));
        CHECK(QDir(cache).entryList(QDir::Files).size() == 1);                                                   // no partial file left

        // KLIFS: a mapping of 8 residues of 1AZM
        const QString mapping = dir.filePath("klifs_map.tsv");
        {
            QFile f(mapping);
            CHECK(f.open(QIODevice::WriteOnly));
            f.write("residue.klifs_id\tresidue.id\tresidue.klifs_region\tresidue.aa\n1\t91\tI\tF\n2\t92\tI\tQ\n3\t94\tI\tH\n4\t96\thinge\tH\n"
                    "5\t119\thinge\tH\n6\t199\thinge\tT\n7\t200\thinge\tH\n8\t201\thinge\tP\n");
        }
        gui::MainWindow k;
        k.setInteractive(false);
        CHECK(k.openStructure(argv[1]));
        k.selectLigand(0);
        k.analyze(true);
        const QString pocketOut = dir.filePath("klifs_pocket.tsv");
        CHECK(k.generateKlifsPocket("1AZM", "A", "", mapping, pocketOut, true));
        CHECK(QFile::exists(pocketOut) && QString::fromStdString(k.fingerprint().pocketName).contains("KLIFS"));
        CHECK(k.logText().contains("KLIFS: generated 8 mapped residues") && k.fingerprint().entries.size() >= 8);
        CHECK(k.fingerprint().entries[0].position == "1" && k.fingerprint().entries[0].region == "I");
        // errors are reported in the log, nothing is written and nothing is loaded
        const QString bad = dir.filePath("klifs_bad.tsv");
        {
            QFile f(bad);
            CHECK(f.open(QIODevice::WriteOnly));
            f.write("residue.klifs_id\tresidue.id\tresidue.klifs_region\n1\t5\tI\n1\t6\tI\n");
        }
        const QString never = dir.filePath("klifs_never.tsv");
        CHECK(k.generateKlifsPocket("1AZM", "A", "", bad, never, true) && k.logText().contains("KLIFS generation failed") && !QFile::exists(never));
        CHECK(!k.generateKlifsPocket("1AZM", "A", "", mapping, QString(), true));                               // no output file
        CHECK(k.generateKlifsPocket("1AZM", "", "", mapping, never, true) && !QFile::exists(never));            // chain is required
        const QString bgOut = dir.filePath("klifs_bg.tsv");
        CHECK(k.generateKlifsPocket("1AZM", "A", "", mapping, bgOut));                                          // in the background
        CHECK(waitFor(k, 10000) && QFile::exists(bgOut));
    }

    // ---- interaction matrix: a list of ligands, the matrix in a background thread, exports, and the command line that reproduces it
    {
        const QString dataDir = QFileInfo(argv[1]).absolutePath();
        const QString root = QFileInfo(dataDir + "/../..").canonicalFilePath();
        const QString pocketFile = root + "/docs/examples/hiv-protease-pocket.tsv";
        CHECK(QFile::exists(pocketFile) && QFile::exists(dataDir + "/1HSG.pdb") && QFile::exists(dataDir + "/1HVR.pdb") && QFile::exists(dataDir + "/1HXB.pdb"));
        auto writeFile = [](const QString& path, const QString& text) {
            QFile f(path);
            const bool ok = f.open(QIODevice::WriteOnly);
            f.write(text.toUtf8());
            return ok;
        };
        auto readFile = [](const QString& path) {
            QFile f(path);
            return f.open(QIODevice::ReadOnly) ? QString::fromUtf8(f.readAll()) : QString();
        };
        const QString list = dir.filePath("matrix_list.tsv");
        CHECK(writeFile(list, QString("structure\tligand\tlabel\n%1/1HSG.pdb\tMK1\tindinavir\n%1/1HVR.pdb\tXK2\tXK263\n%1/1HXB.pdb\tROC\tsaquinavir\n").arg(dataDir)));
        gui::MainWindow m;
        m.setInteractive(false);
        CHECK(!m.runMatrix(true) && m.logText().contains("Open a ligand list first"));                            // nothing to compute yet
        CHECK(writeFile(dir.filePath("matrix_bad.tsv"), "ligand\tchain\nMK1\tA\n"));
        CHECK(!m.loadMatrixList(dir.filePath("matrix_bad.tsv")) && m.logText().contains("'structure' column"));    // a list that is not valid is refused
        CHECK(!m.loadMatrixList(dir.filePath("no_such_list.tsv")) && m.matrixTable()->rowCount() == 0);
        CHECK(m.loadMatrixList(list));
        QTabWidget* tabs = m.findChild<QTabWidget*>();
        CHECK(tabs && tabs->tabText(tabs->currentIndex()) == "Matrix");
        CHECK(m.matrixTable()->rowCount() == 3 && m.matrixTable()->item(0, 0)->text() == "indinavir" && m.matrixTable()->item(0, 3)->text() == QString::fromUtf8("—"));
        CHECK(m.matrixCommandLine().isEmpty());                                                                    // before the first computation
        CHECK(m.loadPocket(pocketFile));
        CHECK(m.runMatrix(true));
        CHECK(m.matrix().pocketMode && m.matrix().rows.size() == 3 && m.matrix().classes.size() == 7 && m.matrix().positions.size() > 20);
        int total = 0;
        for (const auto& row : m.matrix().rows) total += row.bitCount();
        CHECK(total > 30);
        CHECK(m.matrixTable()->rowCount() == 3 && m.matrixTable()->item(0, 3)->text().toInt() == m.matrix().rows[0].bitCount() && m.matrixTable()->item(2, 0)->text() == "saquinavir");
        CHECK(m.logText().contains("Matrix: 3 ligands") && m.logText().contains("Command line: ligint --matrix"));
        // the three exports
        const QString csvFile = dir.filePath("matrix_gui.csv"), jsonFile = dir.filePath("matrix_gui.json"), svgFile = dir.filePath("matrix_gui.svg");
        CHECK(m.exportMatrix(csvFile) && m.exportMatrix(jsonFile) && m.exportMatrix(svgFile) && !m.exportMatrix(dir.filePath("matrix_gui.txt")));
        const QString csv = readFile(csvFile);
        CHECK(csv.startsWith("label,structure,ligand,") && csv.split('\n', Qt::SkipEmptyParts).size() == 4 && csv.contains("indinavir,"));
        const QJsonDocument json = QJsonDocument::fromJson(readFile(jsonFile).toUtf8());
        CHECK(json.isObject() && json.object()["rows"].toArray().size() == 3 && json.object()["pocket"].toBool() &&
              json.object()["rows"].toArray()[0].toObject()["bits"].toArray().size() == static_cast<int>(m.matrix().positions.size() * 7));
        const QString svg = readFile(svgFile);
        CHECK(svg.startsWith("<svg ") && svg.contains("indinavir") && svg.contains("matrix_list"));                // the title is the name of the list
        CHECK(m.matrixCommandLine().contains(" --pocket ") && m.matrixCommandLine("svg").contains("--matrix-format svg --matrix-title matrix_list") && !m.matrixCommandLine("csv").contains("--matrix-format"));
#ifndef Q_OS_WIN
        auto matrixEquivalence = [&](const char* what) {
            if (argc <= 2) return;
            for (const char* fmt : {"csv", "json", "svg"}) {
                const QString file = dir.filePath(QString("matrix_eq.") + fmt);
                CHECK(m.exportMatrix(file));
                bool ok = false;
                const QString out = runShell(withBinary(m.matrixCommandLine(fmt), argv[2], ""), &ok);
                if (!(ok && out == readFile(file))) std::cerr << "matrix command line differs from the window: " << what << " " << fmt << "\n";
                CHECK(ok && out == readFile(file));
            }
        };
        matrixEquivalence("defaults");
#endif
        // the options of the tab, and the analysis settings of the panel
        m.findChild<QCheckBox*>("matrix_all")->setChecked(true);
        CHECK(m.runMatrix(true) && m.matrix().classes.size() == 11 && m.matrixCommandLine().contains("--matrix-classes all"));
        const std::size_t withEmpty = m.matrix().positions.size();
        m.findChild<QCheckBox*>("matrix_hide")->setChecked(true);
        CHECK(m.runMatrix(true) && m.matrix().positions.size() < withEmpty && m.matrixCommandLine().contains("--matrix-hide-empty"));
        m.findChild<QCheckBox*>("opt_water")->setChecked(false);
        m.findChild<QDoubleSpinBox*>("thr_hbMax")->setValue(3.2);
        CHECK(m.runMatrix(true) && m.matrixCommandLine().contains("--no-water") && m.matrixCommandLine().contains("--hb-max 3.2"));
#ifndef Q_OS_WIN
        matrixEquivalence("all classes, no empty positions, no water, hb-max 3.2");
#endif
        m.findChild<QCheckBox*>("opt_water")->setChecked(true);
        m.findChild<QDoubleSpinBox*>("thr_hbMax")->setValue(3.5);
        m.findChild<QCheckBox*>("matrix_all")->setChecked(false);
        m.findChild<QCheckBox*>("matrix_hide")->setChecked(false);
        // in a background thread: the result arrives through the event loop
        CHECK(m.runMatrix(false));
        QElapsedTimer timer;
        timer.start();
        while (m.matrixRunning() && timer.elapsed() < 30000) QApplication::processEvents(QEventLoop::AllEvents, 20);
        QApplication::processEvents();
        CHECK(!m.matrixRunning() && m.matrix().rows.size() == 3 && m.matrix().classes.size() == 7);
        // an entry with a ligand that the structure does not have: reported, and no matrix (nothing that a command could not reproduce)
        CHECK(writeFile(dir.filePath("matrix_wrong.tsv"), QString("structure\tligand\tlabel\n%1/1HSG.pdb\tMK1\tok\n%1/1HVR.pdb\tZZZ\twrong\n").arg(dataDir)));
        CHECK(m.loadMatrixList(dir.filePath("matrix_wrong.tsv")));
        CHECK(m.runMatrix(true) && m.logText().contains("Matrix failed") && m.logText().contains("Entry 2") && m.matrix().rows.empty() && m.matrixCommandLine().isEmpty());
        CHECK(!m.exportMatrix(dir.filePath("matrix_never.csv")) && !QFile::exists(dir.filePath("matrix_never.csv")));
        CHECK(m.matrixTable()->rowCount() == 2 && m.matrixTable()->item(1, 1)->text() == "ZZZ");                   // back to the list of the entries
    }

    // a PDBx/mmCIF file gives the same interactions as the PDB file of the same entry
    {
        const QString cifPath = QFileInfo(argv[1]).absolutePath() + "/1AZM.cif";
        CHECK(QFileInfo::exists(cifPath));
        auto keysOf = [](gui::MainWindow& win, const QString& path) {
            QStringList keys;
            if (!win.openStructure(path) || win.ligandCount() < 1) return keys;
            win.selectLigand(0);
            win.analyze(true);
            if (win.result().ligands.empty()) return keys;
            for (const auto& h : win.result().ligands.front().hits)
                keys << (QString::fromStdString(h.type) + "|" + QString::fromStdString(h.ligandAtom) + "|" + QString::fromStdString(h.proteinAtom) + "|" +
                         QString::number(h.distance, 'f', 2));
            keys.sort();
            return keys;
        };
        gui::MainWindow wp, wc;
        wp.setInteractive(false);
        wc.setInteractive(false);
        const QStringList fromPdb = keysOf(wp, argv[1]), fromCif = keysOf(wc, cifPath);
        CHECK(!fromPdb.isEmpty() && fromPdb == fromCif);
        CHECK(wc.logText().contains("Command line: ligint") && wc.lastCommandLine().contains("1AZM.cif"));
    }

    if (failures) { std::cerr << failures << " checks failed\n"; return 1; }
    std::cout << "ligint_gui_smoke: ok (" << w.hitsModel()->rowCount() << " interactions, image saved)\n";
    return 0;
}
