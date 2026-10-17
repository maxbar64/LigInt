// Optional 3Dmol.js-backed WebGL viewer. Built only when Qt WebEngine is available.
#pragma once

#include <QJsonObject>
#include <QWidget>
#include <vector>

#include "ligint/ligint.hpp"

class QCheckBox;
class QComboBox;
class QWebEngineView;

namespace gui {

class Molecular3DWidget : public QWidget {
public:
    explicit Molecular3DWidget(QWidget* parent = nullptr);

    void setStructure(const QString& data, const QString& format,
                      const QString& ligandName, const QString& chain,
                      int residueNumber, char insertionCode,
                      const std::vector<ligint::HitRecord>& hits);

private:
    void sendPayload();
    void runViewerCommand(const QString& command);

    QWebEngineView* webView_ = nullptr;
    QComboBox* styleCombo_ = nullptr;
    QCheckBox* contactsCheck_ = nullptr;
    QCheckBox* surfaceCheck_ = nullptr;
    QJsonObject payload_;
    bool pageReady_ = false;
};

}  // namespace gui
