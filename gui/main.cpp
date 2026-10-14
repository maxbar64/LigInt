// ligint-gui: interfaccia grafica Qt.
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
    parser.setApplicationDescription("Interazioni ligando-proteina da file PDB");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("file", "File PDB da aprire (facoltativo)");
    parser.process(app);

    gui::MainWindow w;
    w.show();
    const QStringList args = parser.positionalArguments();
    if (!args.isEmpty() && w.openStructure(args.first()) && w.ligandCount() > 0) w.analyze(false);
    return app.exec();
}
