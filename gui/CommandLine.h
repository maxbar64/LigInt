// Command-line strings for reproducing a GUI analysis with the `ligint` program.
#pragma once

#include <QString>
#include <QStringList>

namespace gui {

// Quotes one argument for a POSIX shell (bash, zsh, sh): plain arguments are left as they are, the others are
// wrapped in single quotes ("it's" -> 'it'\''s'). Other shells (cmd, PowerShell) quote differently.
QString shellQuote(const QString& arg);

// "program arg1 arg2 ..." with each argument quoted.
QString joinCommand(const QString& program, const QStringList& args);

}  // namespace gui
