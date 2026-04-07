#include "MultiScaleSpatialDecomposer.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{
namespace
{

constexpr double kEpsilon0 = 8.8541878128e-12;
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kElectronMass = 9.1093837015e-31;
constexpr double kMinDensity = 1.0e6;
constexpr double kMinTemperatureEv = 1.0e-3;
constexpr double kMinCellLength = 1.0e-12;

std::size_t stableNeighborCount(const std::vector<MultiScaleRegionType>& regions,
                                std::size_t index,
                                MultiScaleRegionType target)
{
    std::size_t matches = 0;
    if (index > 0 && regions[index - 1] == target)
    {
        ++matches;
    }
    if (index + 1 < regions.size() && regions[index + 1] == target)
    {
        ++matches;
    }
    return matches;
}

} // namespace

std::string toString(MultiScaleRegionType region)
{
    switch (region)
    {
    case MultiScaleRegionType::FluidMacro:
        return "fluid_macro";
    case MultiScaleRegionType::HybridTransition:
        return "hybrid_transition";
    case MultiScaleRegionType::LocalPIC:
        return "local_pic";
    }
    return "unknown";
}

double MultiScaleSpatialDecomposer::computeDebyeLengthM(double electron_density_m3,
                                                        double electron_temperature_ev)
{
    const double ne = std::max(kMinDensity, electron_density_m3);
    const double te_ev = std::max(kMinTemperatureEv, electron_temperature_ev);
    return std::sqrt(kEpsilon0 * te_ev / (ne * kElementaryCharge));
}

double MultiScaleSpatialDecomposer::computePlasmaFrequencyHz(double electron_density_m3)
{
    const double ne = std::max(kMinDensity, electron_density_m3);
    return std::sqrt(ne * kElementaryCharge * kElementaryCharge / (kEpsilon0 * kElectronMass)) /
           (2.0 * 3.14159265358979323846);
}

double MultiScaleSpatialDecomposer::computeCellLengthM(const std::vector<double>& axis_values,
                                                       std::size_t index)
{
    if (axis_values.empty())
    {
        return 0.0;
    }
    if (axis_values.size() == 1)
    {
        return kMinCellLength;
    }

    if (index == 0)
    {
        return std::max(kMinCellLength, std::abs(axis_values[1] - axis_values[0]));
    }
    if (index + 1 == axis_values.size())
    {
        return std::max(kMinCellLength,
                        std::abs(axis_values[index] - axis_values[index - 1]));
    }

    const double left = std::abs(axis_values[index] - axis_values[index - 1]);
    const double right = std::abs(axis_values[index + 1] - axis_values[index]);
    return std::max(kMinCellLength, 0.5 * (left + right));
}

std::string MultiScaleSpatialDecomposer::makeReproducibilitySignature(
    const std::vector<MultiScaleRegionType>& regions)
{
    std::uint64_t hash = 1469598103934665603ULL;
    constexpr std::uint64_t kPrime = 1099511628211ULL;

    for (const auto region : regions)
    {
        hash ^= static_cast<std::uint64_t>(region);
        hash *= kPrime;
    }

    std::ostringstream stream;
    stream << "msd-v1-" << std::hex << std::setfill('0') << std::setw(16) << hash;
    return stream.str();
}

std::vector<MultiScaleRegionType>
MultiScaleSpatialDecomposer::smoothRegions(const std::vector<MultiScaleRegionType>& regions) const
{
    if (regions.size() < 3 || thresholds_.smoothing_passes == 0)
    {
        return regions;
    }

    std::vector<MultiScaleRegionType> current = regions;
    for (std::size_t pass = 0; pass < thresholds_.smoothing_passes; ++pass)
    {
        std::vector<MultiScaleRegionType> next = current;
        for (std::size_t i = 1; i + 1 < current.size(); ++i)
        {
            if (current[i - 1] == current[i + 1] && current[i] != current[i - 1])
            {
                next[i] = current[i - 1];
                continue;
            }

            // Prefer the lower-complexity region if both neighbors differ and the
            // center has no same-type neighbor, to keep partition boundaries stable.
            if (stableNeighborCount(current, i, current[i]) == 0 &&
                current[i - 1] != current[i + 1])
            {
                next[i] = std::min(current[i - 1], current[i + 1]);
            }
        }
        current.swap(next);
    }

    return current;
}

MultiScaleSpatialDecomposition
MultiScaleSpatialDecomposer::decomposeFromProfile(const Output::ColumnarDataSet& profile,
                                                  const PlasmaParameters& plasma) const
{
    const auto density_it = profile.scalar_series.find("electron_density_m3");
    const auto field_it = profile.scalar_series.find("electric_field_z_v_per_m");
    if (density_it == profile.scalar_series.end() || field_it == profile.scalar_series.end())
    {
        throw std::runtime_error(
            "MultiScaleSpatialDecomposer requires electron_density_m3 and electric_field_z_v_per_m columns");
    }

    const std::vector<double>& density = density_it->second;
    const std::vector<double>& field = field_it->second;
    const std::vector<double>& axis = profile.axis_values;

    if (density.empty() || field.empty() || axis.empty())
    {
        throw std::runtime_error("MultiScaleSpatialDecomposer received empty profile arrays");
    }
    if (density.size() != field.size() || density.size() != axis.size())
    {
        throw std::runtime_error("MultiScaleSpatialDecomposer profile column length mismatch");
    }

    const auto te_it = profile.scalar_series.find("electron_temperature_ev");
    const auto ti_it = profile.scalar_series.find("ion_temperature_ev");

    std::vector<MultiScaleCellIndicators> indicators(density.size());
    std::vector<MultiScaleRegionType> regions(density.size(), MultiScaleRegionType::FluidMacro);

    for (std::size_t i = 0; i < density.size(); ++i)
    {
        const double ne = std::max(kMinDensity, density[i]);
        const double te = (te_it != profile.scalar_series.end() && i < te_it->second.size())
                              ? std::max(kMinTemperatureEv, te_it->second[i])
                              : std::max(kMinTemperatureEv, plasma.electron_temperature_ev);
        const double ti = (ti_it != profile.scalar_series.end() && i < ti_it->second.size())
                              ? std::max(kMinTemperatureEv, ti_it->second[i])
                              : std::max(kMinTemperatureEv, plasma.ion_temperature_ev);

        const double cell_length_m = computeCellLengthM(axis, i);
        const double debye_length_m = computeDebyeLengthM(ne, te);
        const double debye_ratio = cell_length_m / std::max(debye_length_m, kMinCellLength);

        double dne_dz = 0.0;
        double dfield_dz = 0.0;
        if (i == 0)
        {
            const double dz = std::max(kMinCellLength, axis[1] - axis[0]);
            dne_dz = (density[1] - density[0]) / dz;
            dfield_dz = (field[1] - field[0]) / dz;
        }
        else if (i + 1 == density.size())
        {
            const double dz = std::max(kMinCellLength, axis[i] - axis[i - 1]);
            dne_dz = (density[i] - density[i - 1]) / dz;
            dfield_dz = (field[i] - field[i - 1]) / dz;
        }
        else
        {
            const double dz = std::max(kMinCellLength, axis[i + 1] - axis[i - 1]);
            dne_dz = (density[i + 1] - density[i - 1]) / dz;
            dfield_dz = (field[i + 1] - field[i - 1]) / dz;
        }

        const double density_gradient_indicator =
            cell_length_m * std::abs(dne_dz) / std::max(ne, kMinDensity);
        const double field_gradient_indicator =
            cell_length_m * std::abs(dfield_dz) / std::max(std::abs(field[i]), 1.0);

        // Scaled surrogate for nu_collision / omega_pe; coefficient is chosen so
        // typical plasma analysis presets span fluid/hybrid/PIC partitions.
        const double nu_collision_hz =
            std::max(0.0, plasma.neutral_density_m3) * 1.0e-9 * std::sqrt(std::max(ti, 0.0));
        const double omega_pe_hz = computePlasmaFrequencyHz(ne);
        const double collisionality_ratio = nu_collision_hz / std::max(omega_pe_hz, 1.0e-9);

        const double non_equilibrium_ratio = std::max(te / ti, ti / te);

        std::size_t pic_votes = 0;
        std::size_t hybrid_votes = 0;

        if (debye_ratio >= thresholds_.debye_ratio_pic_threshold)
        {
            ++pic_votes;
        }
        else if (debye_ratio >= thresholds_.debye_ratio_hybrid_threshold)
        {
            ++hybrid_votes;
        }

        if (density_gradient_indicator >= thresholds_.density_gradient_pic_threshold)
        {
            ++pic_votes;
        }
        else if (density_gradient_indicator >= thresholds_.density_gradient_hybrid_threshold)
        {
            ++hybrid_votes;
        }

        if (field_gradient_indicator >= thresholds_.field_gradient_pic_threshold)
        {
            ++pic_votes;
        }
        else if (field_gradient_indicator >= thresholds_.field_gradient_hybrid_threshold)
        {
            ++hybrid_votes;
        }

        if (collisionality_ratio <= thresholds_.collisionality_pic_threshold)
        {
            ++pic_votes;
        }
        else if (collisionality_ratio <= thresholds_.collisionality_hybrid_threshold)
        {
            ++hybrid_votes;
        }

        if (non_equilibrium_ratio >= thresholds_.non_equilibrium_pic_threshold)
        {
            ++pic_votes;
        }
        else if (non_equilibrium_ratio >= thresholds_.non_equilibrium_hybrid_threshold)
        {
            ++hybrid_votes;
        }

        MultiScaleRegionType region = MultiScaleRegionType::FluidMacro;
        if (pic_votes >= 2 || (pic_votes >= 1 && hybrid_votes >= 2))
        {
            region = MultiScaleRegionType::LocalPIC;
        }
        else if (hybrid_votes >= 2)
        {
            region = MultiScaleRegionType::HybridTransition;
        }

        regions[i] = region;
        indicators[i] = MultiScaleCellIndicators{
            cell_length_m,
            debye_length_m,
            debye_ratio,
            density_gradient_indicator,
            field_gradient_indicator,
            collisionality_ratio,
            non_equilibrium_ratio,
            region,
        };
    }

    const std::vector<MultiScaleRegionType> smoothed_regions = smoothRegions(regions);
    for (std::size_t i = 0; i < indicators.size(); ++i)
    {
        indicators[i].region = smoothed_regions[i];
    }

    MultiScaleDecompositionSummary summary;
    for (const auto region : smoothed_regions)
    {
        switch (region)
        {
        case MultiScaleRegionType::FluidMacro:
            ++summary.fluid_macro_cells;
            break;
        case MultiScaleRegionType::HybridTransition:
            ++summary.hybrid_transition_cells;
            break;
        case MultiScaleRegionType::LocalPIC:
            ++summary.local_pic_cells;
            break;
        }
    }
    summary.reproducibility_signature = makeReproducibilitySignature(smoothed_regions);

    return MultiScaleSpatialDecomposition{
        smoothed_regions,
        indicators,
        summary,
    };
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
