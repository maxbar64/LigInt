// ligint-gui: Qt graphical interface.
#include <QApplication>
#include <QCommandLineParser>

#include "MainWindow.h"
#include "ligint/ligint.hpp"

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    QApplication::setOrganizationName("ligint");
    QApplication::setApplicationName("ligint-gui");
    QApplication::setApplicationVersion(LIGINT_VERSION);

    QCommandLineParser parser;
    parser.setApplicationDescription("Ligand-protein interactions from PDB files");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("file", "PDB file to open (optional)");
    parser.process(app);

    gui::MainWindow w;
    w.show();
    const QStringList args = parser.positionalArguments();
    if (!args.isEmpty() && w.openStructure(args.first()) && w.ligandCount() > 0) w.analyze(false);
    return app.exec();
}
