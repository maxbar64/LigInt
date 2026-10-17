// What the 3D view shows and how: the choices of the Graphics panel.
#pragma once

#include <QString>

namespace gui {

enum class ProteinStyle { Off, Cartoon, Tube, Trace };
enum class ProteinColor { SecondaryStructure, Chain, Rainbow, BFactor, Uniform };
enum class SurfaceMode { Off, Pocket, Protein };
enum class SurfaceType { SolventExcluded, SolventAccessible };
enum class SurfaceColor { Neutral, Atoms, Hydrophobicity, Protein };    // Protein = the same colors as the cartoon
enum class LigandStyle { Sticks, BallStick, Spheres };
enum class ResidueStyle { Off, Sticks };
enum class Background { White, Light, Dark };

struct ViewStyle {
    ProteinStyle protein = ProteinStyle::Cartoon;
    ProteinColor proteinColor = ProteinColor::SecondaryStructure;
    SurfaceMode surface = SurfaceMode::Off;
    SurfaceType surfaceType = SurfaceType::SolventExcluded;
    SurfaceColor surfaceColor = SurfaceColor::Neutral;
    double surfaceOpacity = 0.45;          // 1 = opaque
    LigandStyle ligand = LigandStyle::BallStick;
    ResidueStyle residues = ResidueStyle::Sticks;
    bool showLigand = true;
    bool showWaters = true;
    bool showLabels = true, showDistances = true, showLegend = true;
    bool perspective = true;
    bool depthCue = true;
    double clipMargin = 3.0;               // Angstrom of protein kept in front of the ligand; 0 = no clipping
    Background background = Background::White;
};

}  // namespace gui
