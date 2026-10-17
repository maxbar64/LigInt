// Panel with the graphics choices of the 3D view: protein, surface, ligand, interactions, camera.
#pragma once

#include <QWidget>

namespace gui {

class ViewerWidget;

class ViewerOptionsWidget final : public QWidget {
    Q_OBJECT
public:
    explicit ViewerOptionsWidget(ViewerWidget* viewer, QWidget* parent = nullptr);

    // Enables the protein and surface controls only when the view has the protein (call it after a new analysis).
    void updateAvailability();

private:
    ViewerWidget* viewer_;
    QList<QWidget*> proteinControls_;
};

}  // namespace gui
