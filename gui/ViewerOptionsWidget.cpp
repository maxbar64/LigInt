#include "ViewerOptionsWidget.h"

#include <QCheckBox>
#include <QComboBox>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QSlider>
#include <QVBoxLayout>

#include "ViewerWidget.h"

namespace gui {

ViewerOptionsWidget::ViewerOptionsWidget(ViewerWidget* viewer, QWidget* parent) : QWidget(parent), viewer_(viewer) {
    setObjectName("graphicsOptionsWidget");
    setMinimumWidth(200);

    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(8, 8, 8, 8);
    layout->setSpacing(8);

    auto addGroup = [&](const QString& title) {
        auto* group = new QGroupBox(title, this);
        group->setFlat(true);
        auto* gl = new QVBoxLayout(group);
        gl->setContentsMargins(8, 6, 8, 6);
        gl->setSpacing(4);
        layout->addWidget(group);
        return gl;
    };
    // a combo box: the items' data are the enum values; `apply` receives the current style to modify
    auto addCombo = [&](QVBoxLayout* gl, const QString& label, const QString& objectName, const QStringList& items, int current, bool protein,
                        std::function<void(ViewStyle&, int)> apply) {
        auto* l = new QLabel(label, this);
        gl->addWidget(l);
        auto* box = new QComboBox(this);
        box->setObjectName(objectName);
        box->addItems(items);
        box->setCurrentIndex(current);
        gl->addWidget(box);
        connect(box, qOverload<int>(&QComboBox::currentIndexChanged), viewer_, [this, apply](int i) {
            ViewStyle s = viewer_->style();
            apply(s, i);
            viewer_->setStyle(s);
        });
        if (protein) { proteinControls_ << l << box; }
        return box;
    };
    auto addCheck = [&](QVBoxLayout* gl, const QString& text, const QString& objectName, bool checked, std::function<void(ViewStyle&, bool)> apply) {
        auto* c = new QCheckBox(text, this);
        c->setObjectName(objectName);
        c->setChecked(checked);
        gl->addWidget(c);
        connect(c, &QCheckBox::toggled, viewer_, [this, apply](bool on) {
            ViewStyle s = viewer_->style();
            apply(s, on);
            viewer_->setStyle(s);
        });
        return c;
    };
    const ViewStyle d = viewer_->style();

    // ---- protein
    auto* protein = addGroup(tr("Protein"));
    addCombo(protein, tr("Representation"), "graphics_protein_style", {tr("Cartoon"), tr("Tube"), tr("Backbone trace"), tr("Hidden")}, 0, true,
             [](ViewStyle& s, int i) { static const ProteinStyle v[] = {ProteinStyle::Cartoon, ProteinStyle::Tube, ProteinStyle::Trace, ProteinStyle::Off}; s.protein = v[i]; });
    addCombo(protein, tr("Color"), "graphics_protein_color", {tr("Secondary structure"), tr("Chain"), tr("Rainbow (N to C)"), tr("B-factor"), tr("Single color")}, 0, true,
             [](ViewStyle& s, int i) {
                 static const ProteinColor v[] = {ProteinColor::SecondaryStructure, ProteinColor::Chain, ProteinColor::Rainbow, ProteinColor::BFactor, ProteinColor::Uniform};
                 s.proteinColor = v[i];
             });

    // ---- surface
    auto* surface = addGroup(tr("Surface"));
    addCombo(surface, tr("Show"), "graphics_surface_mode", {tr("Hidden"), tr("Pocket around the ligand"), tr("Whole protein")}, 0, true,
             [](ViewStyle& s, int i) { static const SurfaceMode v[] = {SurfaceMode::Off, SurfaceMode::Pocket, SurfaceMode::Protein}; s.surface = v[i]; });
    addCombo(surface, tr("Type"), "graphics_surface_type", {tr("Solvent-excluded"), tr("Solvent-accessible")}, 0, true,
             [](ViewStyle& s, int i) { s.surfaceType = i == 0 ? SurfaceType::SolventExcluded : SurfaceType::SolventAccessible; });
    addCombo(surface, tr("Color"), "graphics_surface_color", {tr("Neutral"), tr("By atom"), tr("Hydrophobicity"), tr("Like the protein")}, 0, true,
             [](ViewStyle& s, int i) {
                 static const SurfaceColor v[] = {SurfaceColor::Neutral, SurfaceColor::Atoms, SurfaceColor::Hydrophobicity, SurfaceColor::Protein};
                 s.surfaceColor = v[i];
             });
    auto* opacityLabel = new QLabel(tr("Opacity"), this);
    surface->addWidget(opacityLabel);
    auto* opacity = new QSlider(Qt::Horizontal, this);
    opacity->setObjectName("graphics_surface_opacity");
    opacity->setRange(10, 100);
    opacity->setValue(static_cast<int>(d.surfaceOpacity * 100));
    surface->addWidget(opacity);
    connect(opacity, &QSlider::valueChanged, viewer_, [this](int v) {
        ViewStyle s = viewer_->style();
        s.surfaceOpacity = v / 100.0;
        viewer_->setStyle(s);
    });
    proteinControls_ << opacityLabel << opacity;

    // ---- ligand
    auto* ligand = addGroup(tr("Ligand and residues"));
    addCheck(ligand, tr("Show ligand"), "graphics_show_ligand", d.showLigand, [](ViewStyle& s, bool on) { s.showLigand = on; });
    addCombo(ligand, tr("Ligand style"), "graphics_ligand_style", {tr("Ball and stick"), tr("Sticks"), tr("Spheres")}, 0, false,
             [](ViewStyle& s, int i) { static const LigandStyle v[] = {LigandStyle::BallStick, LigandStyle::Sticks, LigandStyle::Spheres}; s.ligand = v[i]; });
    addCombo(ligand, tr("Aromatic bonds"), "graphics_aromatic_style",
             {tr("Aromatic dashed"), tr("Pseudo-Kekule"), tr("Simple bonds")}, 0, false,
             [](ViewStyle& s,int i) { static const AromaticDisplay modes[]={
                 AromaticDisplay::Dashed,AromaticDisplay::KekuleLike,AromaticDisplay::Single
             }; s.aromaticDisplay=modes[i]; });
    addCheck(ligand, tr("Residues in contact"), "graphics_show_residues", d.residues == ResidueStyle::Sticks,
             [](ViewStyle& s, bool on) { s.residues = on ? ResidueStyle::Sticks : ResidueStyle::Off; });
    addCheck(ligand, tr("Waters of the bridges"), "graphics_show_waters", d.showWaters, [](ViewStyle& s, bool on) { s.showWaters = on; });

    // ---- interactions
    auto* inter = addGroup(tr("Interactions"));
    addCheck(inter, tr("Residue labels"), "graphics_show_labels", d.showLabels, [](ViewStyle& s, bool on) { s.showLabels = on; });
    addCheck(inter, tr("Distances"), "graphics_show_distances", d.showDistances, [](ViewStyle& s, bool on) { s.showDistances = on; });
    addCheck(inter, tr("Legend"), "graphics_show_legend", d.showLegend, [](ViewStyle& s, bool on) { s.showLegend = on; });

    // ---- view
    auto* view = addGroup(tr("View"));
    addCheck(view, tr("Perspective"), "graphics_perspective", d.perspective, [](ViewStyle& s, bool on) { s.perspective = on; });
    addCheck(view, tr("Depth cue (fog)"), "graphics_depth_cue", d.depthCue, [](ViewStyle& s, bool on) { s.depthCue = on; });
    addCombo(view, tr("Background"), "graphics_background", {tr("White"), tr("Light gray"), tr("Dark")}, 0, false,
             [](ViewStyle& s, int i) { static const Background v[] = {Background::White, Background::Light, Background::Dark}; s.background = v[i]; });
    auto* clipLabel = new QLabel(this);
    view->addWidget(clipLabel);
    auto* clip = new QSlider(Qt::Horizontal, this);
    clip->setObjectName("graphics_clip");
    clip->setRange(0, 40);
    clip->setToolTip(tr("Cuts away the protein in front of the ligand so that it is not hidden. 0 = show everything."));
    view->addWidget(clip);
    auto setClipText = [clipLabel](int v) { clipLabel->setText(v == 0 ? tr("Clip in front of the ligand: off") : tr("Clip in front of the ligand: %1 Å").arg(v)); };
    clip->setValue(static_cast<int>(d.clipMargin));
    setClipText(clip->value());
    connect(clip, &QSlider::valueChanged, viewer_, [this, setClipText](int v) {
        setClipText(v);
        ViewStyle s = viewer_->style();
        s.clipMargin = v;
        viewer_->setStyle(s);
    });
    proteinControls_ << clipLabel << clip;

    auto* camera = addGroup(tr("Camera"));
    auto addButton = [&](const QString& text, const QString& objectName, const QString& tip, auto slot) {
        auto* b = new QPushButton(text, this);
        b->setObjectName(objectName);
        b->setToolTip(tip);
        camera->addWidget(b);
        connect(b, &QPushButton::clicked, viewer_, slot);
        return b;
    };
    addButton(tr("Reset view"), "graphics_reset_view", tr("Restore the initial orientation, zoom and pan"), &ViewerWidget::resetView);
    addButton(tr("Focus ligand"), "graphics_focus_ligand", tr("Frame the ligand and its pocket"), &ViewerWidget::focusLigand);
    proteinControls_ << addButton(tr("Fit protein"), "graphics_fit_protein", tr("Frame the whole protein"), &ViewerWidget::fitProtein);

    layout->addStretch(1);
    updateAvailability();
}

void ViewerOptionsWidget::updateAvailability() {
    const bool on = viewer_->hasProtein();
    for (QWidget* w : proteinControls_) w->setEnabled(on);
}

}  // namespace gui
