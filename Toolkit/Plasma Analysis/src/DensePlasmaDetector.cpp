#include "DensePlasmaDetector.h"

#include <cmath>

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
}

DensePlasmaAssessment DensePlasmaDetector::assess(const PlasmaParameters& parameters) const
{
    DensePlasmaAssessment assessment;
    const double ne = std::max(1.0e6, parameters.electron_density_m3);
    const double te_ev = std::max(1.0e-3, parameters.electron_temperature_ev);

    assessment.debye_length_m = std::sqrt(kEpsilon0 * te_ev / (ne * kElementaryCharge));
    assessment.plasma_frequency_hz =
        std::sqrt(ne * kElementaryCharge * kElementaryCharge / (kEpsilon0 * kElectronMass)) /
        (2.0 * 3.14159265358979323846);
    assessment.collisionality = parameters.neutral_density_m3 * 1.0e-20 *
                                std::sqrt(std::max(0.0, parameters.ion_temperature_ev));
    assessment.is_dense = ne >= dense_threshold_m3_ || assessment.debye_length_m < 5.0e-5;
    return assessment;
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
