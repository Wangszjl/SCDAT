#include "PICFluidIntegration.h"

#include <filesystem>
#include <fstream>
#include <iomanip>

namespace
{

const char* environmentModelName(
    SCDAT::Toolkit::PlasmaAnalysis::PlasmaEnvironmentModelKind model)
{
    using SCDAT::Toolkit::PlasmaAnalysis::PlasmaEnvironmentModelKind;

    switch (model)
    {
    case PlasmaEnvironmentModelKind::SpisCcpReference:
        return "spis_ccp_reference";
    case PlasmaEnvironmentModelKind::SpisOrbitalWake:
        return "spis_orbital_wake";
    case PlasmaEnvironmentModelKind::SpisThrusterPlume:
    default:
        return "spis_thruster_plume";
    }
}

const char* distributionModelName(
    SCDAT::Toolkit::PlasmaAnalysis::PlasmaDistributionModelKind model)
{
    using SCDAT::Toolkit::PlasmaAnalysis::PlasmaDistributionModelKind;

    switch (model)
    {
    case PlasmaDistributionModelKind::MaxwellianProjected:
        return "maxwellian_projected";
    case PlasmaDistributionModelKind::WakeAnisotropic:
        return "wake_anisotropic";
    case PlasmaDistributionModelKind::MultipleSurf:
        return "multiple_surf";
    case PlasmaDistributionModelKind::LocalModifiedPearsonIV:
        return "local_modified_pearson_iv";
    case PlasmaDistributionModelKind::LocalTabulated:
        return "local_tabulated";
    case PlasmaDistributionModelKind::TwoAxesTabulatedVelocity:
        return "two_axes_tabulated_velocity";
    case PlasmaDistributionModelKind::GlobalMaxwellBoltzmann:
        return "global_maxwell_boltzmann";
    case PlasmaDistributionModelKind::GlobalMaxwellBoltzmann2:
        return "global_maxwell_boltzmann2";
    case PlasmaDistributionModelKind::GlobalMaxwell:
        return "global_maxwell";
    case PlasmaDistributionModelKind::LocalMaxwell:
        return "local_maxwell";
    case PlasmaDistributionModelKind::LocalMaxwell2:
        return "local_maxwell2";
    case PlasmaDistributionModelKind::RecollMaxwell:
        return "recoll_maxwell";
    case PlasmaDistributionModelKind::PICSurf:
        return "pic_surf";
    case PlasmaDistributionModelKind::NonPICSurf:
        return "non_pic_surf";
    case PlasmaDistributionModelKind::GenericSurf:
        return "generic_surf";
    case PlasmaDistributionModelKind::GlobalSurf:
        return "global_surf";
    case PlasmaDistributionModelKind::LocalGenericSurf:
        return "local_generic_surf";
    case PlasmaDistributionModelKind::TestableSurf:
        return "testable_surf";
    case PlasmaDistributionModelKind::TestableForA:
        return "testable_for_a";
    case PlasmaDistributionModelKind::MaxwellianThruster:
        return "maxwellian_thruster";
    case PlasmaDistributionModelKind::UniformVelocity:
        return "uniform_velocity";
    case PlasmaDistributionModelKind::Fluid:
        return "fluid";
    case PlasmaDistributionModelKind::FowlerNordheim:
        return "fowler_nordheim";
    case PlasmaDistributionModelKind::AxisymTabulatedVelocity:
        return "axisym_tabulated_velocity";
    case PlasmaDistributionModelKind::MultiPopulationHybrid:
    default:
        return "multi_population_hybrid";
    }
}

const char* reactionRegistryName(
    SCDAT::Toolkit::PlasmaAnalysis::PlasmaReactionRegistryKind registry)
{
    using SCDAT::Toolkit::PlasmaAnalysis::PlasmaReactionRegistryKind;

    switch (registry)
    {
    case PlasmaReactionRegistryKind::SpisDensePlasma:
        return "spis_dense_plasma";
    case PlasmaReactionRegistryKind::SpisCore:
    default:
        return "spis_core";
    }
}

const char* diagnosticSetName(
    SCDAT::Toolkit::PlasmaAnalysis::PlasmaDiagnosticSetKind diagnostic_set)
{
    using SCDAT::Toolkit::PlasmaAnalysis::PlasmaDiagnosticSetKind;

    switch (diagnostic_set)
    {
    case PlasmaDiagnosticSetKind::SpisCoreDiagnostics:
        return "spis_core_diagnostics";
    case PlasmaDiagnosticSetKind::SheathMultiscaleDiagnostics:
        return "sheath_multiscale_diagnostics";
    case PlasmaDiagnosticSetKind::FullPhysicsDiagnostics:
    default:
        return "full_physics_diagnostics";
    }
}

std::filesystem::path plasmaPhysicsDiagnosticsPath(const std::filesystem::path& csv_path)
{
    auto diagnostics_path = csv_path;
    diagnostics_path.replace_extension(".physics_diagnostics.json");
    return diagnostics_path;
}

bool exportPlasmaPhysicsDiagnosticsArtifact(
    const std::filesystem::path& diagnostics_path,
    const SCDAT::Toolkit::PlasmaAnalysis::FluidAlgorithmConfig& config,
    const SCDAT::Toolkit::PlasmaAnalysis::FluidAlgorithmStatus& status)
{
    const auto parent_path = diagnostics_path.parent_path();
    if (!parent_path.empty())
    {
        std::filesystem::create_directories(parent_path);
    }

    std::ofstream output(diagnostics_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        return false;
    }

    const bool boundary_layer_resolved =
        status.sheath_thickness_m > 0.0 && status.presheath_thickness_m > 0.0;
    const bool reaction_balance_active = status.reaction_active_processes > 0u;
    const bool non_equilibrium_active =
        status.non_equilibrium_ratio > 1.0 || status.non_equilibrium_relaxation_rate_hz > 0.0;
    const bool turbulence_active =
        status.turbulence_eddy_diffusivity_m2_per_s > 0.0 ||
        status.turbulence_dissipation_rate_w_per_m3 > 0.0;
    const double net_reaction_source_m3_per_s =
        status.ionization_source_m3_per_s - status.recombination_sink_m3_per_s;

    output << std::fixed << std::setprecision(12);
    output << "{\n";
    output << "  \"schema_version\": \"scdat.plasma.physics_diagnostics.v1\",\n";
    output << "  \"contract_id\": \"plasma-physics-diagnostics-v1\",\n";
    output << "  \"reaction_contract_id\": \"plasma-reaction-balance-v1\",\n";
    output << "  \"module\": \"Plasma Analysis\",\n";
    output << "  \"organization_family\": \""
           << (config.enable_spis_style_organization ? "spis_numeric_v1" : "native") << "\",\n";
    output << "  \"environment_model\": \"" << environmentModelName(config.environment_model)
           << "\",\n";
    output << "  \"distribution_model\": \"" << distributionModelName(config.distribution_model)
           << "\",\n";
    output << "  \"reaction_registry\": \"" << reactionRegistryName(config.reaction_registry)
           << "\",\n";
    output << "  \"diagnostic_set\": \"" << diagnosticSetName(config.diagnostic_set) << "\",\n";
    output << "  \"dense_plasma_detected\": "
           << (status.dense_plasma_detected ? "true" : "false") << ",\n";
    output << "  \"time_s\": " << status.time_s << ",\n";
    output << "  \"steps_completed\": " << status.steps_completed << ",\n";
    output << "  \"transport_summary\": {\n";
    output << "    \"average_density_m3\": " << status.average_density_m3 << ",\n";
    output << "    \"average_potential_v\": " << status.average_potential_v << ",\n";
    output << "    \"debye_length_m\": " << status.debye_length_m << "\n";
    output << "  },\n";
    output << "  \"boundary_layer_summary\": {\n";
    output << "    \"sheath_thickness_m\": " << status.sheath_thickness_m << ",\n";
    output << "    \"presheath_thickness_m\": " << status.presheath_thickness_m << ",\n";
    output << "    \"sheath_potential_drop_v\": " << status.sheath_potential_drop_v << ",\n";
    output << "    \"presheath_potential_drop_v\": " << status.presheath_potential_drop_v
           << ",\n";
    output << "    \"boundary_wall_field_v_per_m\": " << status.boundary_wall_field_v_per_m
           << ",\n";
    output << "    \"ion_bohm_velocity_m_per_s\": " << status.ion_bohm_velocity_m_per_s
           << ",\n";
    output << "    \"ion_mach_at_sheath_edge\": " << status.ion_mach_at_sheath_edge << "\n";
    output << "  },\n";
    output << "  \"reaction_balance_summary\": {\n";
    output << "    \"ionization_source_m3_per_s\": " << status.ionization_source_m3_per_s
           << ",\n";
    output << "    \"recombination_sink_m3_per_s\": " << status.recombination_sink_m3_per_s
           << ",\n";
    output << "    \"net_reaction_source_m3_per_s\": " << net_reaction_source_m3_per_s
           << ",\n";
    output << "    \"effective_collision_frequency_hz\": "
           << status.effective_collision_frequency_hz << ",\n";
    output << "    \"charge_exchange_frequency_hz\": " << status.charge_exchange_frequency_hz
           << ",\n";
    output << "    \"reaction_energy_loss_ev_per_s\": " << status.reaction_energy_loss_ev_per_s
           << ",\n";
    output << "    \"reaction_momentum_transfer_ratio\": "
           << status.reaction_momentum_transfer_ratio << ",\n";
    output << "    \"reaction_active_processes\": " << status.reaction_active_processes << "\n";
    output << "  },\n";
    output << "  \"advanced_closure_summary\": {\n";
    output << "    \"advanced_closure_enabled\": "
           << (status.advanced_closure_enabled ? "true" : "false") << ",\n";
    output << "    \"non_equilibrium_ratio\": " << status.non_equilibrium_ratio << ",\n";
    output << "    \"non_equilibrium_relaxation_rate_hz\": "
           << status.non_equilibrium_relaxation_rate_hz << ",\n";
    output << "    \"electron_temperature_closure_delta_ev\": "
           << status.electron_temperature_closure_delta_ev << ",\n";
    output << "    \"turbulence_intensity\": " << status.turbulence_intensity << ",\n";
    output << "    \"turbulence_eddy_diffusivity_m2_per_s\": "
           << status.turbulence_eddy_diffusivity_m2_per_s << ",\n";
    output << "    \"turbulence_dissipation_rate_w_per_m3\": "
           << status.turbulence_dissipation_rate_w_per_m3 << ",\n";
    output << "    \"closure_density_correction_m3\": " << status.closure_density_correction_m3
           << ",\n";
    output << "    \"closure_active_terms\": " << status.closure_active_terms << "\n";
    output << "  },\n";
    output << "  \"consistency_flags\": {\n";
    output << "    \"boundary_layer_resolved\": "
           << (boundary_layer_resolved ? "true" : "false") << ",\n";
    output << "    \"reaction_balance_active\": "
           << (reaction_balance_active ? "true" : "false") << ",\n";
    output << "    \"non_equilibrium_active\": "
           << (non_equilibrium_active ? "true" : "false") << ",\n";
    output << "    \"turbulence_active\": " << (turbulence_active ? "true" : "false") << "\n";
    output << "  },\n";
    output << "  \"config_snapshot\": {\n";
    output << "    \"domain_size_xyz_m\": ["
           << config.domain_size.x() << ", " << config.domain_size.y() << ", "
           << config.domain_size.z() << "],\n";
    output << "    \"resolution_xyz\": ["
           << config.resolution.x() << ", " << config.resolution.y() << ", "
           << config.resolution.z() << "],\n";
    output << "    \"dense_plasma_threshold_m3\": " << config.dense_plasma_threshold_m3 << ",\n";
    output << "    \"initial_potential_v\": " << config.initial_potential_v << "\n";
    output << "  }\n";
    output << "}\n";

    return static_cast<bool>(output);
}

} // namespace

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
    auto data_set = adapter_.buildProfileDataSet();
    const auto diagnostics_path = plasmaPhysicsDiagnosticsPath(csv_path);
    data_set.metadata["plasma_physics_diagnostics_artifact_path"] =
        diagnostics_path.filename().string();
    data_set.metadata["plasma_physics_diagnostics_schema_version"] =
        "scdat.plasma.physics_diagnostics.v1";

    if (!static_cast<bool>(exporter_.exportDataSet(csv_path, data_set)))
    {
        return false;
    }

    return exportPlasmaPhysicsDiagnosticsArtifact(diagnostics_path, adapter_.getConfig(),
                                                  adapter_.getStatus());
}

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
