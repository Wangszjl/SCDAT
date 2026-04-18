#pragma once

#include "SurfaceFieldVolumeContracts.h"

#include "../../Mesh/include/MeshParsing.h"
#include "../../Particle/include/ParticleDefinitions.h"

#include <string>
#include <vector>

namespace SCDAT
{
namespace FieldSolver
{

double blendFieldVolumeScalar(double base_value, double target_value, double blend_factor);

double updateFieldVolumeMismatchMetric(double current_metric, double candidate_value,
                                       double baseline_value, double tolerance);

double computeFieldVolumeExternalBlendFactor(double mismatch_metric,
                                             double linear_residual_norm,
                                             double max_delta_v,
                                             bool linear_converged,
                                             double potential_tolerance_v,
                                             double volume_linear_relaxation);

std::vector<NativeVolumeBoundaryConditionFamily> resolveNativeVolumeBoundaryConditionFamilies(
    NativeVolumeParityRoute route);
std::vector<NativeVolumeBoundaryConditionFamily> resolveActiveNativeVolumeBoundaryConditionFamilies(
    NativeVolumeParityRoute route);
std::vector<NativeVolumeBoundaryConditionFamily> normalizeNativeVolumeBoundaryConditionFamilies(
    const std::vector<NativeVolumeBoundaryConditionFamily>& families);

std::vector<NativeVolumeFieldFamily> resolveNativeVolumeFieldFamilies(
    NativeVolumeParityRoute route);
std::vector<NativeVolumeFieldFamily> resolveActiveNativeVolumeFieldFamilies(
    NativeVolumeParityRoute route);
std::vector<NativeVolumeFieldFamily> normalizeNativeVolumeFieldFamilies(
    const std::vector<NativeVolumeFieldFamily>& families);

std::vector<NativeVolumeDistributionFamily> resolveNativeVolumeDistributionFamilies(
    NativeVolumeParityRoute route);
std::vector<NativeVolumeDistributionFamily> normalizeNativeVolumeDistributionFamilies(
    const std::vector<NativeVolumeDistributionFamily>& families);
std::vector<NativeVolumeInteractionFamily> resolveActiveNativeVolumeInteractionFamilies(
    NativeVolumeParityRoute route);

std::vector<NativeVolumeInteractionFamily> resolveNativeVolumeInteractionFamilies(
    NativeVolumeParityRoute route);
std::vector<NativeVolumeInteractionFamily> normalizeNativeVolumeInteractionFamilies(
    const std::vector<NativeVolumeInteractionFamily>& families);

std::string nativeVolumeBoundaryConditionFamilySignature(NativeVolumeParityRoute route);
std::string nativeVolumeActiveBoundaryConditionFamilySignature(NativeVolumeParityRoute route);
std::string nativeVolumeBoundaryConditionFamilySignature(
    const std::vector<NativeVolumeBoundaryConditionFamily>& families);
std::string nativeVolumeFieldFamilySignature(NativeVolumeParityRoute route);
std::string nativeVolumeActiveFieldFamilySignature(NativeVolumeParityRoute route);
std::string nativeVolumeFieldFamilySignature(
    const std::vector<NativeVolumeFieldFamily>& families);
std::string nativeVolumeDistributionFamilySignature(NativeVolumeParityRoute route);
std::string nativeVolumeDistributionFamilySignature(
    const std::vector<NativeVolumeDistributionFamily>& families);
std::string nativeVolumeActiveInteractionFamilySignature(NativeVolumeParityRoute route);
std::string nativeVolumeInteractionFamilySignature(NativeVolumeParityRoute route);
std::string nativeVolumeInteractionFamilySignature(
    const std::vector<NativeVolumeInteractionFamily>& families);

NativeVolumeMeshSummary summarizeNativeVolumeMesh(const Mesh::VolMesh& mesh);

NativeVolumeFieldSummary solveNativeVolumePotential(const Mesh::VolMeshPtr& mesh,
                                                    double lower_dirichlet_v,
                                                    double upper_dirichlet_v,
                                                    NativeVolumeParityRoute route =
                                                        NativeVolumeParityRoute::NativeMinimal);

NativeVolumeDistributionSummary evaluateNativeVolumeDistributions(
    std::size_t node_count, const std::vector<Particle::ParticlePtr>& particles,
    NativeVolumeParityRoute route = NativeVolumeParityRoute::NativeMinimal);

NativeVolumeInteractionSummary executeNativeVolumeInteractions(
    std::vector<Particle::ParticlePtr>& particles, double dt,
    NativeVolumeParityRoute route = NativeVolumeParityRoute::NativeMinimal);
NativeVolumeInteractionSummary executeNativeVolumeInteractions(
    std::vector<Particle::ParticlePtr>& particles, double dt,
    const std::vector<NativeVolumeInteractionFamily>& active_families);

NativeVolumeParitySnapshot buildNativeVolumeParitySnapshot(
    const Mesh::VolMeshPtr& mesh, std::vector<Particle::ParticlePtr>& particles,
    double lower_dirichlet_v, double upper_dirichlet_v, double dt,
    NativeVolumeParityRoute route = NativeVolumeParityRoute::NativeMinimal);

} // namespace FieldSolver
} // namespace SCDAT
