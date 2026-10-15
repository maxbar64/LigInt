#include "CommandLine.h"

namespace gui {

QString shellQuote(const QString& arg) {
    if (arg.isEmpty()) return QStringLiteral("''");
    bool plain = true;
    for (QChar c : arg) {
        const ushort u = c.unicode();
        const bool ok = (u >= 'a' && u <= 'z') || (u >= 'A' && u <= 'Z') || (u >= '0' && u <= '9') ||
                        u == '_' || u == '@' || u == '%' || u == '+' || u == '=' || u == ':' || u == ',' || u == '.' ||
                        u == '/' || u == '-';
        if (!ok) { plain = false; break; }
    }
    if (plain) return arg;
    QString q = arg;
    q.replace(QLatin1Char('\''), QStringLiteral("'\\''"));
    return QLatin1Char('\'') + q + QLatin1Char('\'');
}

QString joinCommand(const QString& program, const QStringList& args) {
    QString s = shellQuote(program);
    for (const QString& a : args) s += QLatin1Char(' ') + shellQuote(a);
    return s;
}

}  // namespace gui
