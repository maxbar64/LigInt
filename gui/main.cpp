// ligint-gui: Qt graphical interface.
#include <QApplication>
#include <QIcon>
#include <QCommandLineParser>

#include "MainWindow.h"
#include "ligint/ligint.hpp"

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    QApplication::setOrganizationName("ligint");
    QApplication::setApplicationName("ligint-gui");
    QApplication::setApplicationVersion(LIGINT_VERSION);
    QApplication::setWindowIcon(QIcon(":/icons/ligint-256.png"));
    QGuiApplication::setDesktopFileName("ligint-gui");      // ties the window to the menu entry (ligint-gui.desktop)

    QCommandLineParser parser;
    parser.setApplicationDescription("Ligand-protein interactions from PDB and mmCIF files");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("file", "PDB or mmCIF file to open (optional)");
    parser.process(app);

    gui::MainWindow w;
    w.show();
    const QStringList args = parser.positionalArguments();
    if (!args.isEmpty() && w.openStructure(args.first()) && w.ligandCount() > 0) w.analyze(false);
    return app.exec();
}
