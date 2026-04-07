#pragma once

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

enum class ElectronCollectionModelKind
{
    OmlLike,
    ShiftedEnergy,
    BarrierLimited
};

class SurfaceElectronCollectionModel
{
  public:
    virtual ~SurfaceElectronCollectionModel() = default;

    virtual double populationFactor(double surface_potential_v, double temperature_ev) const = 0;

    virtual double fluxCollectionWeight(double surface_potential_v, double incident_energy_ev,
                                        double reference_energy_ev) const = 0;

    virtual double fluxEmissionWeight(double surface_potential_v, double incident_energy_ev,
                                      double reference_energy_ev) const = 0;
};

namespace detail
{

inline double positiveFloor(double value, double minimum)
{
    return std::max(value, minimum);
}

inline double clampedExp(double exponent)
{
    return std::exp(std::clamp(exponent, -200.0, 60.0));
}

} // namespace detail

class OmlLikeSurfaceElectronCollectionModel final : public SurfaceElectronCollectionModel
{
  public:
    double populationFactor(double surface_potential_v, double temperature_ev) const override
    {
        const double thermal_ev = detail::positiveFloor(temperature_ev, 1.0e-3);
        if (surface_potential_v <= 0.0)
        {
            return detail::clampedExp(surface_potential_v / thermal_ev);
        }
        return 1.0 + surface_potential_v / thermal_ev;
    }

    double fluxCollectionWeight(double surface_potential_v, double incident_energy_ev,
                                double reference_energy_ev) const override
    {
        if (surface_potential_v <= 0.0 && incident_energy_ev + surface_potential_v <= 0.0)
        {
            return 0.0;
        }
        if (surface_potential_v > 0.0)
        {
            return 1.0 + surface_potential_v / detail::positiveFloor(reference_energy_ev, 1.0e-3);
        }
        return 1.0;
    }

    double fluxEmissionWeight(double surface_potential_v, double incident_energy_ev,
                              double /*reference_energy_ev*/) const override
    {
        if (surface_potential_v <= 0.0 && incident_energy_ev + surface_potential_v <= 0.0)
        {
            return 0.0;
        }
        return 1.0;
    }
};

class ShiftedEnergySurfaceElectronCollectionModel final : public SurfaceElectronCollectionModel
{
  public:
    double populationFactor(double surface_potential_v, double temperature_ev) const override
    {
        const double thermal_ev = detail::positiveFloor(temperature_ev, 1.0e-3);
        if (surface_potential_v <= 0.0)
        {
            return detail::clampedExp(surface_potential_v / thermal_ev);
        }
        return std::sqrt(std::max(0.0, 1.0 + surface_potential_v / thermal_ev));
    }

    double fluxCollectionWeight(double surface_potential_v, double incident_energy_ev,
                                double reference_energy_ev) const override
    {
        const double incident = detail::positiveFloor(incident_energy_ev, 1.0e-6);
        const double shifted = incident + surface_potential_v;
        if (shifted <= 0.0)
        {
            return 0.0;
        }

        if (surface_potential_v <= 0.0)
        {
            return std::clamp(shifted / incident, 0.0, 1.0);
        }

        const double reference = detail::positiveFloor(reference_energy_ev, 1.0e-3);
        return 1.0 + surface_potential_v / (incident + reference);
    }

    double fluxEmissionWeight(double surface_potential_v, double incident_energy_ev,
                              double reference_energy_ev) const override
    {
        const double incident = detail::positiveFloor(incident_energy_ev, 1.0e-6);
        const double shifted = incident + surface_potential_v;
        if (shifted <= 0.0)
        {
            return 0.0;
        }

        if (surface_potential_v <= 0.0)
        {
            return std::clamp(shifted / incident, 0.0, 1.0);
        }

        const double reference = detail::positiveFloor(reference_energy_ev, 1.0e-3);
        return 1.0 + 0.5 * surface_potential_v / (incident + reference);
    }
};

class BarrierLimitedSurfaceElectronCollectionModel final : public SurfaceElectronCollectionModel
{
  public:
    double populationFactor(double surface_potential_v, double temperature_ev) const override
    {
        const double thermal_ev = detail::positiveFloor(temperature_ev, 1.0e-3);
        if (surface_potential_v <= 0.0)
        {
            return detail::clampedExp(surface_potential_v / thermal_ev);
        }
        return std::clamp(1.0 + 0.5 * surface_potential_v / thermal_ev, 1.0, 3.0);
    }

    double fluxCollectionWeight(double surface_potential_v, double incident_energy_ev,
                                double /*reference_energy_ev*/) const override
    {
        if (surface_potential_v <= 0.0 && incident_energy_ev + surface_potential_v <= 0.0)
        {
            return 0.0;
        }
        return 1.0;
    }

    double fluxEmissionWeight(double surface_potential_v, double incident_energy_ev,
                              double /*reference_energy_ev*/) const override
    {
        if (surface_potential_v <= 0.0 && incident_energy_ev + surface_potential_v <= 0.0)
        {
            return 0.0;
        }
        return 1.0;
    }
};

inline const SurfaceElectronCollectionModel&
selectSurfaceElectronCollectionModel(ElectronCollectionModelKind model)
{
    static const OmlLikeSurfaceElectronCollectionModel kOmlLikeModel;
    static const ShiftedEnergySurfaceElectronCollectionModel kShiftedEnergyModel;
    static const BarrierLimitedSurfaceElectronCollectionModel kBarrierLimitedModel;

    switch (model)
    {
    case ElectronCollectionModelKind::ShiftedEnergy:
        return kShiftedEnergyModel;
    case ElectronCollectionModelKind::BarrierLimited:
        return kBarrierLimitedModel;
    case ElectronCollectionModelKind::OmlLike:
    default:
        return kOmlLikeModel;
    }
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
