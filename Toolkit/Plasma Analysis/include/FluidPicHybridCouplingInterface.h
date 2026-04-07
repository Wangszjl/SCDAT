#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct FluidPicCouplingConfig
{
    double potential_jump_tolerance_v = 5.0;
    double charge_imbalance_tolerance_c = 1.0e-8;
    double energy_imbalance_tolerance_j = 1.0e-8;
    double correction_relaxation = 1.0;
    std::size_t max_correction_iterations = 4;
};

struct FluidPicInterfaceState
{
    std::string interface_id;
    std::string fluid_region_id;
    std::string pic_region_id;

    // 0 -> prefer fluid side, 1 -> prefer PIC side.
    double transition_weight = 0.5;

    double fluid_potential_v = 0.0;
    double pic_potential_v = 0.0;

    // Charge/current exchange terms.
    double fluid_charge_current_a = 0.0;
    double pic_charge_current_a = 0.0;

    // Energy exchange terms.
    double fluid_energy_power_w = 0.0;
    double pic_energy_power_w = 0.0;

    // Effective coupling conductivity on each side.
    double fluid_effective_conductivity_s_per_m = 0.0;
    double pic_effective_conductivity_s_per_m = 0.0;
};

struct FluidPicExchangeRecord
{
    std::string interface_id;
    std::string fluid_region_id;
    std::string pic_region_id;
    std::string units_tag;

    double time_s = 0.0;
    double macro_dt_s = 0.0;
    double boundary_potential_v = 0.0;
    double exchange_charge_current_a = 0.0;
    double exchange_energy_power_w = 0.0;
    double effective_conductivity_s_per_m = 0.0;
    double charge_imbalance_c = 0.0;
    double energy_imbalance_j = 0.0;
};

struct FluidPicCouplingMetrics
{
    std::size_t interface_count = 0;
    std::size_t correction_iterations = 0;
    double max_potential_jump_before_v = 0.0;
    double max_potential_jump_after_v = 0.0;
    double total_charge_imbalance_c = 0.0;
    double total_energy_imbalance_j = 0.0;
    bool continuity_within_tolerance = false;
    bool conservation_within_tolerance = false;
};

class FluidPicHybridCouplingInterface
{
  public:
    bool initialize(const FluidPicCouplingConfig& config);

    std::vector<FluidPicExchangeRecord>
    couple(double time_s, double macro_dt_s, std::vector<FluidPicInterfaceState>& interfaces);

    const FluidPicCouplingConfig& getConfig() const { return config_; }
    const FluidPicCouplingMetrics& getMetrics() const { return metrics_; }

    void reset();

  private:
    void updateImbalanceMetrics(double macro_dt_s,
                                const std::vector<FluidPicInterfaceState>& interfaces);

    FluidPicCouplingConfig config_{};
    FluidPicCouplingMetrics metrics_{};
    bool initialized_ = false;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
