#pragma once

#include "Monitor.h"

#include "../../Geometry/include/Vector3D.h"

#include <functional>
#include <vector>

namespace SCDAT
{
namespace Diagnostics
{

class FieldMonitor : public Monitor
{
  public:
    using VectorFieldProvider = std::function<std::vector<Geometry::Vector3D>()>;
    using ScalarFieldProvider = std::function<std::vector<double>()>;

    FieldMonitor(std::string name, VectorFieldProvider electric_field_provider,
                 ScalarFieldProvider potential_provider = {});

    VoidResult sample(double time) override;

  private:
    VectorFieldProvider electric_field_provider_;
    ScalarFieldProvider potential_provider_;
};

} // namespace Diagnostics
} // namespace SCDAT
