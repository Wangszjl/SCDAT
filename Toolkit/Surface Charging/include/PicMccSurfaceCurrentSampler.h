#pragma once

#include "../../Plasma Analysis/include/FluidAlgorithmConfig.h"

#include "../../Tools/Interactions/Collisions/include/SurfacePicMccContracts.h"
#include "../../Tools/Material/include/MaterialProperty.h"
#include "../../Tools/Particle/include/ParticleSource.h"

#include <cstddef>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{

struct PicMccSurfaceSamplerConfig : public Collision::PicMccControlContract
{
    double surface_area_m2 = 1.0e-4;
    double gap_distance_m = 1.0e-3;
    std::string node_name;
    std::string boundary_group_id;
    std::size_t linked_boundary_face_count = 0;
    double projection_weight_sum = 1.0;
    double surface_aspect_ratio = 1.0;
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
    std::string deposition_scheme = "pic_window_cic";
    unsigned int seed = 20260408u;
    std::string sampling_policy = "deterministic";
    PlasmaAnalysis::PlasmaParameters plasma;
    Particle::ResolvedSpectrum electron_spectrum;
    Particle::ResolvedSpectrum ion_spectrum;
    bool has_electron_spectrum = false;
    bool has_ion_spectrum = false;
    std::vector<std::string> source_keys;
    Material::MaterialProperty material{2, Mesh::MaterialType::DIELECTRIC, "surface"};
};

struct PicMccSourceResolvedSample
{
    std::string source_key;
    double collected_current_density_a_per_m2 = 0.0;
    double superparticle_count = 0.0;
    // True only when this slot is backed by an explicit source attribution path.
    bool strictly_attributed = false;
};

struct PicMccCurrentSample : public Collision::PicMccTelemetryContract
{
    bool valid = false;
    std::string surface_domain_family = "boundary_topology_patch_domain_v1";
    std::string sampled_node_name;
    std::string sampled_boundary_group_id;
    double topology_signature = 0.0;
    double topology_area_scale = 1.0;
    double sampled_surface_potential_v = 0.0;
    double electron_collection_current_density_a_per_m2 = 0.0;
    double ion_collection_current_density_a_per_m2 = 0.0;
    double emitted_electron_current_density_a_per_m2 = 0.0;
    double net_collection_current_density_a_per_m2 = 0.0;
    double current_derivative_a_per_m2_per_v = 0.0;
    std::size_t deposition_segments = 0;
    std::string deposition_kernel = "linear_segment_cloud";
    unsigned int seed_used = 0u;
    std::string sampling_policy_resolved = "deterministic";
    std::vector<PicMccSourceResolvedSample> source_resolved_samples;
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
