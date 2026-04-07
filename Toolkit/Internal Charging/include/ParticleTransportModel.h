#pragma once

#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

struct VolumeLayerState
{
    double depth_m = 0.0;
    double charge_density_c_per_m3 = 0.0;
    double deposited_energy_j_per_m3 = 0.0;
  double dose_gy = 0.0;
  double effective_conductivity_s_per_m = 0.0;
  double deposited_charge_rate_c_per_m3_s = 0.0;
  double deposited_particle_rate_per_m3_s = 0.0;
};

class ParticleTransportModel
{
  public:
    bool initialize(double thickness_m, std::size_t layers);
    void depositFlux(double current_density_a_per_m2, double mean_energy_ev,
                     double incident_charge_state_abs, double dt);
    const std::vector<VolumeLayerState>& getLayerStates() const { return layers_; }
    std::vector<VolumeLayerState>& getLayerStates() { return layers_; }
    void reset();

  private:
    double thickness_m_ = 0.0;
    double layer_thickness_m_ = 0.0;
    std::vector<VolumeLayerState> layers_;
};

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
