#include "ViewerOptionsWidget.h"

#include <QCheckBox>
#include <QGroupBox>
#include <QPushButton>
#include <QVBoxLayout>

#include "ViewerWidget.h"

namespace gui {

ViewerOptionsWidget::ViewerOptionsWidget(ViewerWidget* viewer, QWidget* parent) : QWidget(parent) {
    setObjectName("graphicsOptionsWidget");
    setMinimumWidth(185);

    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(8, 8, 8, 8);
    layout->setSpacing(8);

    auto addGroup = [&](const QString& title) {
        auto* group = new QGroupBox(title, this);
        group->setFlat(true);
        auto* groupLayout = new QVBoxLayout(group);
        groupLayout->setContentsMargins(8, 6, 8, 6);
        groupLayout->setSpacing(4);
        layout->addWidget(group);
        return groupLayout;
    };
    auto addToggle = [&](QVBoxLayout* groupLayout, const QString& title, const QString& objectName, bool checked, auto setter) {
        auto* check = new QCheckBox(title, this);
        check->setObjectName(objectName);
        check->setChecked(checked);
        groupLayout->addWidget(check);
        connect(check, &QCheckBox::toggled, viewer, setter);
        return check;
    };

    auto* scene = addGroup(tr("Structure"));
    addToggle(scene, tr("Show ligand"), "graphics_show_ligand", true, &ViewerWidget::setShowLigand);
    addToggle(scene, tr("Show contact residues"), "graphics_show_residues", true, &ViewerWidget::setShowResidues);
    auto* protein = addToggle(scene, tr("Show full protein (C-alpha trace)"), "graphics_show_protein", false, &ViewerWidget::setShowProtein);
    protein->setToolTip(tr("Show the protein backbone as a separate C-alpha trace layer."));

    auto* overlays = addGroup(tr("Interactions"));
    addToggle(overlays, tr("Residue labels"), "graphics_show_labels", true, &ViewerWidget::setShowLabels);
    addToggle(overlays, tr("Distances"), "graphics_show_distances", true, &ViewerWidget::setShowDistances);
    addToggle(overlays, tr("Legend"), "graphics_show_legend", true, &ViewerWidget::setShowLegend);

    auto* camera = addGroup(tr("Camera"));
    auto* reset = new QPushButton(tr("Reset view"), this);
    reset->setObjectName("graphics_reset_view");
    reset->setToolTip(tr("Restore the initial orientation, zoom and pan"));
    camera->addWidget(reset);
    connect(reset, &QPushButton::clicked, viewer, &ViewerWidget::resetView);

    layout->addStretch(1);
}

}  // namespace gui
