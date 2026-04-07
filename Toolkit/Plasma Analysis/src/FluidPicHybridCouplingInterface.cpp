#include "FluidPicHybridCouplingInterface.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{

double clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

bool finiteAndNonNegative(double value)
{
    return std::isfinite(value) && value >= 0.0;
}

} // namespace

bool FluidPicHybridCouplingInterface::initialize(const FluidPicCouplingConfig& config)
{
    if (!std::isfinite(config.potential_jump_tolerance_v) ||
        !std::isfinite(config.charge_imbalance_tolerance_c) ||
        !std::isfinite(config.energy_imbalance_tolerance_j) ||
        !std::isfinite(config.correction_relaxation) ||
        !finiteAndNonNegative(config.potential_jump_tolerance_v) ||
        !finiteAndNonNegative(config.charge_imbalance_tolerance_c) ||
        !finiteAndNonNegative(config.energy_imbalance_tolerance_j) ||
        config.max_correction_iterations == 0 || config.correction_relaxation <= 0.0)
    {
        return false;
    }

    config_ = config;
    metrics_ = FluidPicCouplingMetrics{};
    initialized_ = true;
    return true;
}

void FluidPicHybridCouplingInterface::updateImbalanceMetrics(
    double macro_dt_s, const std::vector<FluidPicInterfaceState>& interfaces)
{
    metrics_.total_charge_imbalance_c = 0.0;
    metrics_.total_energy_imbalance_j = 0.0;

    for (const auto& state : interfaces)
    {
        const double charge_imbalance_c =
            std::abs(state.fluid_charge_current_a - state.pic_charge_current_a) * macro_dt_s;
        const double energy_imbalance_j =
            std::abs(state.fluid_energy_power_w - state.pic_energy_power_w) * macro_dt_s;
        metrics_.total_charge_imbalance_c += charge_imbalance_c;
        metrics_.total_energy_imbalance_j += energy_imbalance_j;
    }

    metrics_.conservation_within_tolerance =
        metrics_.total_charge_imbalance_c <= config_.charge_imbalance_tolerance_c &&
        metrics_.total_energy_imbalance_j <= config_.energy_imbalance_tolerance_j;
}

std::vector<FluidPicExchangeRecord>
FluidPicHybridCouplingInterface::couple(double time_s,
                                        double macro_dt_s,
                                        std::vector<FluidPicInterfaceState>& interfaces)
{
    metrics_ = FluidPicCouplingMetrics{};
    metrics_.interface_count = interfaces.size();

    if (!initialized_ || !std::isfinite(time_s) || !std::isfinite(macro_dt_s) ||
        macro_dt_s <= 0.0)
    {
        return {};
    }

    for (auto& state : interfaces)
    {
        const double jump_before = std::abs(state.fluid_potential_v - state.pic_potential_v);
        metrics_.max_potential_jump_before_v =
            std::max(metrics_.max_potential_jump_before_v, jump_before);

        const double weight = clamp01(state.transition_weight);
        const double boundary_potential_v =
            (1.0 - weight) * state.fluid_potential_v + weight * state.pic_potential_v;

        // Enforce continuity on the switching face by writing the same potential back to both sides.
        state.fluid_potential_v = boundary_potential_v;
        state.pic_potential_v = boundary_potential_v;

        const double boundary_conductivity =
            (1.0 - weight) * state.fluid_effective_conductivity_s_per_m +
            weight * state.pic_effective_conductivity_s_per_m;
        state.fluid_effective_conductivity_s_per_m = boundary_conductivity;
        state.pic_effective_conductivity_s_per_m = boundary_conductivity;

        const double jump_after = std::abs(state.fluid_potential_v - state.pic_potential_v);
        metrics_.max_potential_jump_after_v =
            std::max(metrics_.max_potential_jump_after_v, jump_after);
    }

    metrics_.continuity_within_tolerance =
        metrics_.max_potential_jump_after_v <= config_.potential_jump_tolerance_v;

    updateImbalanceMetrics(macro_dt_s, interfaces);
    while (!metrics_.conservation_within_tolerance &&
           metrics_.correction_iterations < config_.max_correction_iterations)
    {
        for (auto& state : interfaces)
        {
            const double delta_current_a = state.fluid_charge_current_a - state.pic_charge_current_a;
            const double current_correction_a = 0.5 * config_.correction_relaxation * delta_current_a;
            state.fluid_charge_current_a -= current_correction_a;
            state.pic_charge_current_a += current_correction_a;

            const double delta_power_w = state.fluid_energy_power_w - state.pic_energy_power_w;
            const double power_correction_w = 0.5 * config_.correction_relaxation * delta_power_w;
            state.fluid_energy_power_w -= power_correction_w;
            state.pic_energy_power_w += power_correction_w;
        }

        ++metrics_.correction_iterations;
        updateImbalanceMetrics(macro_dt_s, interfaces);
    }

    std::vector<FluidPicExchangeRecord> records;
    records.reserve(interfaces.size());
    for (const auto& state : interfaces)
    {
        const double exchange_charge_current_a =
            0.5 * (state.fluid_charge_current_a + state.pic_charge_current_a);
        const double exchange_energy_power_w =
            0.5 * (state.fluid_energy_power_w + state.pic_energy_power_w);

        FluidPicExchangeRecord record;
        record.interface_id = state.interface_id;
        record.fluid_region_id = state.fluid_region_id;
        record.pic_region_id = state.pic_region_id;
        record.units_tag = "SI[A,W,V,S/m,s,C,J]";
        record.time_s = time_s;
        record.macro_dt_s = macro_dt_s;
        record.boundary_potential_v = state.fluid_potential_v;
        record.exchange_charge_current_a = exchange_charge_current_a;
        record.exchange_energy_power_w = exchange_energy_power_w;
        record.effective_conductivity_s_per_m = state.fluid_effective_conductivity_s_per_m;
        record.charge_imbalance_c =
            std::abs(state.fluid_charge_current_a - state.pic_charge_current_a) * macro_dt_s;
        record.energy_imbalance_j =
            std::abs(state.fluid_energy_power_w - state.pic_energy_power_w) * macro_dt_s;
        records.push_back(record);
    }

    return records;
}

void FluidPicHybridCouplingInterface::reset()
{
    *this = FluidPicHybridCouplingInterface{};
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
