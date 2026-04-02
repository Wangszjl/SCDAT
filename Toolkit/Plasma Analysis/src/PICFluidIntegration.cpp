#include "PICFluidIntegration.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

bool PICFluidIntegration::initialize(const FluidAlgorithmConfig& config)
{
    if (!adapter_.initialize(config))
    {
        return false;
    }

    auto field_monitor = std::make_shared<Diagnostics::FieldMonitor>(
        "plasma_field",
        [this]() {
            std::vector<Geometry::Vector3D> field_values;
            const auto data_set = adapter_.buildProfileDataSet();
            if (const auto it = data_set.scalar_series.find("electric_field_z_v_per_m");
                it != data_set.scalar_series.end())
            {
                field_values.reserve(it->second.size());
                for (const double ez : it->second)
                {
                    field_values.emplace_back(0.0, 0.0, ez);
                }
            }
            return field_values;
        },
        [this]() {
            const auto data_set = adapter_.buildProfileDataSet();
            if (const auto it = data_set.scalar_series.find("potential_v");
                it != data_set.scalar_series.end())
            {
                return it->second;
            }
            return std::vector<double>{};
        });
    monitor_manager_.addMonitor(field_monitor);

    initialized_ = true;
    return true;
}

bool PICFluidIntegration::advance(double dt)
{
    if (!initialized_ || !adapter_.advance(dt))
    {
        return false;
    }
    return static_cast<bool>(monitor_manager_.sampleAll(adapter_.getStatus().time_s));
}

void PICFluidIntegration::reset()
{
    adapter_.reset();
    monitor_manager_ = Diagnostics::MonitorManager{};
    initialized_ = false;
}

bool PICFluidIntegration::exportResults(const std::filesystem::path& csv_path) const
{
    return static_cast<bool>(exporter_.exportDataSet(csv_path, adapter_.buildProfileDataSet()));
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
