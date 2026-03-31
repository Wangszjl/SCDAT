#pragma once

#include "FluidAlgorithmAdapter.h"

#include "../../Tools/Diagnostics/include/FieldMonitor.h"
#include "../../Tools/Diagnostics/include/MonitorManager.h"
#include "../../Tools/Output/include/ResultExporter.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

class PICFluidIntegration
{
  public:
    bool initialize(const FluidAlgorithmConfig& config);
    bool advance(double dt);
    const FluidAlgorithmStatus& getStatus() const { return adapter_.getStatus(); }
    void reset();
    bool exportResults(const std::filesystem::path& csv_path) const;

  private:
    FluidAlgorithmAdapter adapter_;
    Diagnostics::MonitorManager monitor_manager_;
    Output::ResultExporter exporter_;
    bool initialized_ = false;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
