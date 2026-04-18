#pragma once

#include "DensePlasmaSurfaceCharging.h"

#include "../../Tools/Output/include/ResultTypes.h"

#include <cstddef>
#include <filesystem>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

void appendSurfaceInstrumentObservabilityMetadata(
    Output::ColumnarDataSet& data_set, const std::filesystem::path& csv_path,
    SurfaceInstrumentSetKind instrument_set_kind, std::size_t node_count);

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
