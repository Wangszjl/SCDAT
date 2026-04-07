#pragma once

namespace SCDAT
{
namespace Toolkit
{
namespace VacuumArc
{

struct PlasmaChannelParameters
{
  double field_conductivity_gain = 1.0e-5;
  double emission_conductivity_gain = 10.0;
  double field_current_gain = 1.0e-8;
  double resistance_scale = 1.0e3;
  double min_conductivity_s_per_m = 1.0;
  double max_conductivity_s_per_m = 1.0e7;
  double min_current_density_a_per_m2 = 0.0;
  double max_current_density_a_per_m2 = 1.0e8;
  double min_resistance_ohm = 1.0e-9;
  double max_resistance_ohm = 1.0e12;
};

struct ArcChannelState
{
    double conductivity_s_per_m = 0.0;
    double current_density_a_per_m2 = 0.0;
    double resistance_ohm = 0.0;
};

class PlasmaChannelModel
{
  public:
    void setParameters(const PlasmaChannelParameters& parameters)
    {
        parameters_ = parameters;
    }

    const PlasmaChannelParameters& getParameters() const
    {
        return parameters_;
    }

    ArcChannelState update(double electric_field_v_per_m, double emission_current_density_a_per_m2,
                           double gap_distance_m) const;

  private:
    PlasmaChannelParameters parameters_{};
};

} // namespace VacuumArc
} // namespace Toolkit
} // namespace SCDAT
