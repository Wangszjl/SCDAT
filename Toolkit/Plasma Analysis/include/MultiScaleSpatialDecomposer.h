#pragma once

#include "FluidAlgorithmConfig.h"

#include "../../Tools/Output/include/ResultTypes.h"

#include <cstddef>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

enum class MultiScaleRegionType
{
    FluidMacro,
    HybridTransition,
    LocalPIC
};

struct MultiScaleSpatialThresholds
{
    // Criterion 1: L_cell / lambda_D
    double debye_ratio_hybrid_threshold = 1.2;
    double debye_ratio_pic_threshold = 2.5;

    // Criterion 2: L_cell * |grad(n_e)| / max(n_e, eps)
    double density_gradient_hybrid_threshold = 0.20;
    double density_gradient_pic_threshold = 0.60;

    // Criterion 3: L_cell * |grad(E)| / max(|E|, eps)
    double field_gradient_hybrid_threshold = 0.20;
    double field_gradient_pic_threshold = 0.55;

    // Criterion 4: nu_collision / omega_pe (scaled surrogate)
    double collisionality_hybrid_threshold = 0.05;
    double collisionality_pic_threshold = 0.015;

    // Criterion 5: max(T_e / T_i, T_i / T_e)
    double non_equilibrium_hybrid_threshold = 1.8;
    double non_equilibrium_pic_threshold = 2.8;

    // Small deterministic smoothing to suppress one-cell islands.
    std::size_t smoothing_passes = 1;
};

struct MultiScaleCellIndicators
{
    double cell_length_m = 0.0;
    double debye_length_m = 0.0;
    double debye_ratio = 0.0;
    double density_gradient_indicator = 0.0;
    double field_gradient_indicator = 0.0;
    double collisionality_ratio = 0.0;
    double non_equilibrium_ratio = 0.0;
    MultiScaleRegionType region = MultiScaleRegionType::FluidMacro;
};

struct MultiScaleDecompositionSummary
{
    std::size_t fluid_macro_cells = 0;
    std::size_t hybrid_transition_cells = 0;
    std::size_t local_pic_cells = 0;
    std::string reproducibility_signature;
};

struct MultiScaleSpatialDecomposition
{
    std::vector<MultiScaleRegionType> regions;
    std::vector<MultiScaleCellIndicators> indicators;
    MultiScaleDecompositionSummary summary;
};

class MultiScaleSpatialDecomposer
{
  public:
    explicit MultiScaleSpatialDecomposer(const MultiScaleSpatialThresholds& thresholds = {})
        : thresholds_(thresholds)
    {
    }

    MultiScaleSpatialDecomposition
    decomposeFromProfile(const Output::ColumnarDataSet& profile,
                         const PlasmaParameters& plasma) const;

    const MultiScaleSpatialThresholds& thresholds() const { return thresholds_; }

  private:
    static double computeDebyeLengthM(double electron_density_m3, double electron_temperature_ev);
    static double computePlasmaFrequencyHz(double electron_density_m3);
    static double computeCellLengthM(const std::vector<double>& axis_values, std::size_t index);
    static std::string makeReproducibilitySignature(const std::vector<MultiScaleRegionType>& regions);

    std::vector<MultiScaleRegionType>
    smoothRegions(const std::vector<MultiScaleRegionType>& regions) const;

    MultiScaleSpatialThresholds thresholds_{};
};

std::string toString(MultiScaleRegionType region);

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
