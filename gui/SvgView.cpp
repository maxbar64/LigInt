#include "SvgView.h"

#include <QPainter>

#ifdef LIGINT_GUI_HAS_SVG
#include <QSvgRenderer>
#endif

namespace gui {

SvgView::SvgView(QWidget* parent) : QWidget(parent) {
#ifdef LIGINT_GUI_HAS_SVG
    renderer_ = new QSvgRenderer(this);
#endif
}

SvgView::~SvgView() = default;

bool SvgView::available() {
#ifdef LIGINT_GUI_HAS_SVG
    return true;
#else
    return false;
#endif
}

bool SvgView::load(const QByteArray& svg) {
#ifdef LIGINT_GUI_HAS_SVG
    if (!renderer_->load(svg)) return false;
    resize(sizeHint());
    updateGeometry();
    update();
    return true;
#else
    Q_UNUSED(svg);
    return false;
#endif
}

QSize SvgView::sizeHint() const {
#ifdef LIGINT_GUI_HAS_SVG
    if (renderer_ && renderer_->isValid()) return renderer_->defaultSize();
#endif
    return QSize(400, 200);
}

void SvgView::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.fillRect(rect(), Qt::white);
#ifdef LIGINT_GUI_HAS_SVG
    if (renderer_ && renderer_->isValid()) renderer_->render(&p, QRectF(rect()));
#endif
}

}  // namespace gui
