#include "../include/FieldMonitor.h"

#include <algorithm>
#include <numeric>

namespace SCDAT
{
namespace Diagnostics
{

FieldMonitor::FieldMonitor(std::string name, VectorFieldProvider electric_field_provider,
                           ScalarFieldProvider potential_provider)
    : Monitor(std::move(name)),
      electric_field_provider_(std::move(electric_field_provider)),
      potential_provider_(std::move(potential_provider))
{
}

VoidResult FieldMonitor::sample(double time)
{
    if (!electric_field_provider_)
    {
        return VoidResult::failure("FieldMonitor has no electric field provider");
    }

    const auto electric_field = electric_field_provider_();
    std::unordered_map<std::string, double> scalars;
    if (!electric_field.empty())
    {
        double max_magnitude = 0.0;
        double mean_magnitude = 0.0;
        for (const auto& value : electric_field)
        {
            const double magnitude = value.magnitude();
            max_magnitude = std::max(max_magnitude, magnitude);
            mean_magnitude += magnitude;
        }
        mean_magnitude /= static_cast<double>(electric_field.size());
        scalars["field_max_v_per_m"] = max_magnitude;
        scalars["field_mean_v_per_m"] = mean_magnitude;
    }

    if (potential_provider_)
    {
        const auto potential = potential_provider_();
        if (!potential.empty())
        {
            const auto [min_it, max_it] = std::minmax_element(potential.begin(), potential.end());
            scalars["potential_min_v"] = *min_it;
            scalars["potential_max_v"] = *max_it;
        }
    }

    recordSample(time, std::move(scalars));
    return VoidResult::success();
}

} // namespace Diagnostics
} // namespace SCDAT
