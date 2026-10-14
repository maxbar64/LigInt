// Minimal widget that shows an SVG document (used for the fingerprint figure). Requires Qt Svg.
#pragma once

#include <QByteArray>
#include <QWidget>

#ifdef LIGINT_GUI_HAS_SVG
class QSvgRenderer;
#endif

namespace gui {

class SvgView : public QWidget {
    Q_OBJECT
public:
    explicit SvgView(QWidget* parent = nullptr);
    ~SvgView() override;

    // Loads the document and resizes the widget to its natural size; returns false if the SVG is invalid
    // or if the program was built without Qt Svg.
    bool load(const QByteArray& svg);
    static bool available();
    QSize sizeHint() const override;

protected:
    void paintEvent(QPaintEvent* e) override;

private:
#ifdef LIGINT_GUI_HAS_SVG
    QSvgRenderer* renderer_ = nullptr;
#endif
};

}  // namespace gui
