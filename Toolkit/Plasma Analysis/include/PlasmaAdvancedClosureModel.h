#pragma once

#include "FluidAlgorithmConfig.h"
#include "PlasmaReactionCollisionLibrary.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct AdvancedClosureState
{
    double non_equilibrium_ratio = 1.0;
    double non_equilibrium_relaxation_rate_hz = 0.0;
    double electron_temperature_delta_ev = 0.0;
    double turbulence_intensity = 0.0;
    double turbulence_eddy_diffusivity_m2_per_s = 0.0;
    double turbulence_dissipation_rate_w_per_m3 = 0.0;
    double density_correction_m3 = 0.0;
    std::size_t active_terms = 0;
};

class PlasmaAdvancedClosureModel
{
  public:
    PlasmaAdvancedClosureModel() = default;

    void configure(const AdvancedClosureConfig& config);
    const AdvancedClosureConfig& getConfig() const { return config_; }

    AdvancedClosureState evaluate(const PlasmaParameters& parameters,
                                  const DensePlasmaAssessment& assessment,
                                  const ReactionCollisionState& reaction_state,
                                  double characteristic_length_m,
                                  double dt) const;

  private:
    AdvancedClosureConfig config_;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
