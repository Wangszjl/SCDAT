#include "FluidAlgorithmAdapter.h"

#include <algorithm>
#include <numeric>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

bool FluidAlgorithmAdapter::initialize(const FluidAlgorithmConfig& config)
{
    config_ = config;
    plasma_state_ = config.initial_plasma;
    detector_ = DensePlasmaDetector(config.dense_plasma_threshold_m3);

    FieldSolver::BoltzmannParameters boltzmann_parameters;
    boltzmann_parameters.reduced_field_td =
        std::max(1.0, plasma_state_.electric_field_v_per_m / plasma_state_.neutral_density_m3 * 1.0e21);
    boltzmann_solver_.setParameters(boltzmann_parameters);

    if (!diffusion_solver_.initialize(config_, plasma_state_.electron_density_m3))
    {
        return false;
    }

    poisson_solver_.setParameters(FieldSolver::NonlinearPoissonParameters{});
    if (!poisson_solver_.initializeGrid(config.domain_size, config.resolution))
    {
        return false;
    }

    poisson_solver_.setBoundaryPotential([config](const Geometry::Point3D& point) {
        return point.z() >= config.domain_size.z() ? config.initial_potential_v : 0.0;
    });
    poisson_solver_.setChargeDensityFunction([this](const Geometry::Point3D&, double) {
        return 1.602176634e-19 * (plasma_state_.ion_density_m3 - plasma_state_.electron_density_m3);
    });
    poisson_solver_.setPermittivityFunction([](const Geometry::Point3D&, double) { return 1.0; });
    poisson_solver_.solve();

    assessment_ = detector_.assess(plasma_state_);
    boundary_state_ = boundary_layer_.analyze(plasma_state_, assessment_, config.initial_potential_v);
    status_.average_density_m3 = plasma_state_.electron_density_m3;
    status_.average_potential_v = 0.5 * config.initial_potential_v;
    status_.debye_length_m = assessment_.debye_length_m;
    status_.sheath_thickness_m = boundary_state_.sheath_thickness_m;
    status_.dense_plasma_detected = assessment_.is_dense;
    initialized_ = true;
    return true;
}

bool FluidAlgorithmAdapter::advance(double dt)
{
    if (!initialized_)
    {
        return false;
    }

    assessment_ = detector_.assess(plasma_state_);
    const auto& boltzmann_state = boltzmann_solver_.solve(plasma_state_.electron_temperature_ev);
    plasma_state_.electric_field_v_per_m = boltzmann_state.mean_energy_ev * 500.0;

    const auto& electric_field = poisson_solver_.getElectricField();
    if (!diffusion_solver_.advance(electric_field.empty()
                                       ? std::vector<Geometry::Vector3D>(
                                             diffusion_solver_.getState().density.size())
                                       : electric_field,
                                   dt))
    {
        return false;
    }

    const auto& density = diffusion_solver_.getDensityDistribution();
    const double average_density =
        density.empty() ? 0.0
                        : std::accumulate(density.begin(), density.end(), 0.0) /
                              static_cast<double>(density.size());
    plasma_state_.electron_density_m3 = average_density;

    poisson_solver_.setChargeDensityFunction([this](const Geometry::Point3D& point, double phi) {
        const double sheath_boost =
            1.0 + 0.05 * std::abs(phi) / std::max(1.0e-3, config_.initial_potential_v + 1.0);
        const double ne = plasma_state_.electron_density_m3 / sheath_boost;
        (void)point;
        return 1.602176634e-19 * (plasma_state_.ion_density_m3 - ne);
    });
    poisson_solver_.solve();

    const auto& potential = poisson_solver_.getPotentialField();
    status_.average_potential_v =
        potential.empty()
            ? 0.0
            : std::accumulate(potential.begin(), potential.end(), 0.0) /
                  static_cast<double>(potential.size());
    status_.average_density_m3 = average_density;
    status_.time_s += dt;
    status_.steps_completed += 1;
    status_.dense_plasma_detected = assessment_.is_dense;
    status_.debye_length_m = assessment_.debye_length_m;
    boundary_state_ = boundary_layer_.analyze(plasma_state_, assessment_, status_.average_potential_v);
    status_.sheath_thickness_m = boundary_state_.sheath_thickness_m;
    return true;
}

void FluidAlgorithmAdapter::reset()
{
    *this = FluidAlgorithmAdapter{};
}

Output::ColumnarDataSet FluidAlgorithmAdapter::buildProfileDataSet() const
{
    Output::ColumnarDataSet data_set;
    data_set.axis_name = "z_m";

    const auto& diffusion_state = diffusion_solver_.getState();
    const auto& density = diffusion_state.density;
    const auto& potential = poisson_solver_.getPotentialField();
    const auto& electric_field = poisson_solver_.getElectricField();

    const std::size_t nx = static_cast<std::size_t>(diffusion_state.resolution.x());
    const std::size_t ny = static_cast<std::size_t>(diffusion_state.resolution.y());
    const std::size_t nz = static_cast<std::size_t>(diffusion_state.resolution.z());
    data_set.scalar_series["potential_v"] = std::vector<double>(nz, 0.0);
    data_set.scalar_series["electric_field_z_v_per_m"] = std::vector<double>(nz, 0.0);
    data_set.scalar_series["electron_density_m3"] = std::vector<double>(nz, 0.0);
    data_set.scalar_series["ion_density_m3"] = std::vector<double>(nz, plasma_state_.ion_density_m3);
    data_set.scalar_series["electron_temperature_ev"] =
        std::vector<double>(nz, plasma_state_.electron_temperature_ev);

    for (std::size_t k = 0; k < nz; ++k)
    {
        data_set.axis_values.push_back(static_cast<double>(k) * diffusion_state.spacing.z());
        double potential_sum = 0.0;
        double field_sum = 0.0;
        double density_sum = 0.0;

        for (std::size_t j = 0; j < ny; ++j)
        {
            for (std::size_t i = 0; i < nx; ++i)
            {
                const std::size_t index = i + nx * (j + ny * k);
                if (index < potential.size())
                {
                    potential_sum += potential[index];
                    field_sum += electric_field[index].z();
                }
                if (index < density.size())
                {
                    density_sum += density[index];
                }
            }
        }

        const double plane_scale = 1.0 / static_cast<double>(nx * ny);
        data_set.scalar_series["potential_v"][k] = potential_sum * plane_scale;
        data_set.scalar_series["electric_field_z_v_per_m"][k] = field_sum * plane_scale;
        data_set.scalar_series["electron_density_m3"][k] = density_sum * plane_scale;
    }

    data_set.metadata["module"] = "Plasma Analysis";
    return data_set;
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
