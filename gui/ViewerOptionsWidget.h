// Dedicated panel for the viewer's graphics and camera controls.
#pragma once

#include <QWidget>

namespace gui {

class ViewerWidget;

class ViewerOptionsWidget final : public QWidget {
public:
    explicit ViewerOptionsWidget(ViewerWidget* viewer, QWidget* parent = nullptr);
};

}  // namespace gui
