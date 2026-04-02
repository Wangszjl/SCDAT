#pragma once

#include "../../Plasma Analysis/include/FluidAlgorithmConfig.h"

#include "../../Tools/Material/include/MaterialProperty.h"
#include "../../Tools/Particle/include/ParticleSource.h"

#include <cstddef>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct PicMccSurfaceSamplerConfig
{
    double surface_area_m2 = 1.0e-4;
    double gap_distance_m = 1.0e-3;
    double surface_potential_v = 0.0;
    double plasma_reference_potential_v = 0.0;
    double electron_flow_coupling = 0.0;
    double bulk_flow_velocity_m_per_s = 0.0;
    double flow_alignment_cosine = 1.0;
    double flow_angle_deg = 0.0;
    double ion_directed_velocity_m_per_s = 0.0;
    std::size_t z_layers = 8;
    std::size_t particles_per_element = 2;
    std::size_t window_steps = 12;
    bool enable_mcc = true;
    PlasmaAnalysis::PlasmaParameters plasma;
    Particle::ResolvedSpectrum electron_spectrum;
    Particle::ResolvedSpectrum ion_spectrum;
    bool has_electron_spectrum = false;
    bool has_ion_spectrum = false;
    Material::MaterialProperty material{2, Mesh::MaterialType::DIELECTRIC, "surface"};
};

struct PicMccCurrentSample
{
    bool valid = false;
    bool mcc_enabled = false;
    double sampled_surface_potential_v = 0.0;
    double electron_collection_current_density_a_per_m2 = 0.0;
    double ion_collection_current_density_a_per_m2 = 0.0;
    double emitted_electron_current_density_a_per_m2 = 0.0;
    double net_collection_current_density_a_per_m2 = 0.0;
    double current_derivative_a_per_m2_per_v = 0.0;
    int total_collisions = 0;
};

class PicMccSurfaceCurrentSampler
{
  public:
    PicMccCurrentSample sample(const PicMccSurfaceSamplerConfig& config) const;
    PicMccCurrentSample sampleWithDerivative(const PicMccSurfaceSamplerConfig& config,
                                            double probe_delta_v) const;
};

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
