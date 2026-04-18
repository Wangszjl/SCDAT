#include "SurfaceInstrumentObservability.h"

namespace SCDAT
{
namespace Toolkit
{
namespace SurfaceCharging
{
namespace
{

std::string instrumentObserverSeriesEnabled(SurfaceInstrumentSetKind instrument_set_kind)
{
    return instrument_set_kind == SurfaceInstrumentSetKind::SurfacePicObserverSet ? "1" : "0";
}

} // namespace

void appendSurfaceInstrumentObservabilityMetadata(
    Output::ColumnarDataSet& data_set, const std::filesystem::path& csv_path,
    SurfaceInstrumentSetKind instrument_set_kind, std::size_t node_count)
{
    auto shared_runtime_observer_metadata_path = csv_path;
    shared_runtime_observer_metadata_path.replace_extension(".shared_runtime_observer.json");
    auto monitor_json_path = csv_path;
    monitor_json_path.replace_extension(".monitor.json");

    data_set.metadata["surface_instrument_contract_id"] = "surface-instrument-observer-v1";
    data_set.metadata["surface_instrument_contract_version"] = "v1";
    data_set.metadata["surface_instrument_observability_profile"] =
        "runtime_observer+metadata_catalog+aggregated_post_sensors";
    data_set.metadata["surface_instrument_runtime_class_signature"] =
        "LangmuirProbe+MonitorInstrument+NumericsMonitor+TotalCurrentOnSC+"
        "MultipleESNPotentialSensors+MultipleESNTotalCurrentMonitor";
    data_set.metadata["surface_instrument_metadata_only_class_signature"] =
        "Instrument+VirtualInstrument+InstrumentsCatalogue+TopInstrumentFactory+"
        "BundleAPI_DefaultInstruments+BundleAPI_InstrumentsCatalogue";
    data_set.metadata["surface_instrument_pending_observability_class_signature"] =
        "PotentialPS+PotentialLPS+PointPS+LinePS+SphericalPS+ParticleDetector+"
        "VirtualParticleDetector+ParticleSPS+VelocityDistFunctionPS+EnergyDistFuncPS+"
        "SurfaceFluxDistFunctionPS+VolDistribMomentPS";
    data_set.metadata["surface_instrument_runtime_class_count"] = "6";
    data_set.metadata["surface_instrument_metadata_only_class_count"] = "6";
    data_set.metadata["surface_instrument_pending_observability_class_count"] = "12";

    data_set.metadata["surface_instrument_monitor_artifact"] =
        monitor_json_path.filename().string();
    data_set.metadata["surface_instrument_spatial_post_sensor_artifact"] =
        shared_runtime_observer_metadata_path.filename().string();
    data_set.metadata["surface_instrument_potential_sensor_artifact"] =
        shared_runtime_observer_metadata_path.filename().string();
    data_set.metadata["surface_instrument_particle_detector_artifact"] =
        shared_runtime_observer_metadata_path.filename().string();
    data_set.metadata["surface_instrument_distribution_sensor_artifact"] =
        csv_path.filename().string() + ".metadata.json";
    data_set.metadata["surface_instrument_metadata_catalog_artifact"] =
        csv_path.filename().string() + ".metadata.json";

    data_set.metadata["surface_instrument_exports_node_level_pic_series"] =
        instrumentObserverSeriesEnabled(instrument_set_kind);
    data_set.metadata["surface_instrument_node_count"] = std::to_string(node_count);

    data_set.metadata["surface_instrument_monitor_class_signature"] =
        "LangmuirProbe+MonitorInstrument+NumericsMonitor+TotalCurrentOnSC+"
        "MultipleESNTotalCurrentMonitor";
    data_set.metadata["surface_instrument_monitor_class_count"] = "5";
    data_set.metadata["surface_instrument_spatial_post_sensor_class_signature"] =
        "PotentialPS+PotentialLPS+PointPS+LinePS+SphericalPS+"
        "MultipleESNPotentialSensors+MultipleESNPotentialSensorsSCFrame";
    data_set.metadata["surface_instrument_spatial_post_sensor_class_count"] = "7";
    data_set.metadata["surface_instrument_particle_detector_class_signature"] =
        "AbstractParticleDetector+ParticleDetector+ParticleDetectorInterface+"
        "VirtualParticleDetector+ParticleSPS+TotalSuperParticlePS";
    data_set.metadata["surface_instrument_particle_detector_class_count"] = "6";
    data_set.metadata["surface_instrument_distribution_sensor_class_signature"] =
        "DensityPS+DistFuncOfOneVariable+EnergyDistFuncPS+TotalEnergyDistFuncPS+"
        "VelocityDistFunctionPS+SurfaceFluxDistFunctionPS+VolDistribMomentPS+"
        "KineticEnergyPS+TotalEnergyPS";
    data_set.metadata["surface_instrument_distribution_sensor_class_count"] = "9";

    data_set.metadata["surface_instrument_monitor_observability"] =
        "monitor_json+probe_current_voltage_series";
    data_set.metadata["surface_instrument_spatial_post_sensor_observability"] =
        "shared_runtime_observer+node_potential_histories";
    data_set.metadata["surface_instrument_particle_detector_observability"] =
        "computeSurfaceCurrents+node_level_current_histories";
    data_set.metadata["surface_instrument_distribution_sensor_observability"] =
        "aggregated_sampler_distribution_metadata";
    data_set.metadata["surface_instrument_metadata_catalog_observability"] =
        "metadata_sidecar+class_signature_catalog";
}

} // namespace SurfaceCharging
} // namespace Toolkit
} // namespace SCDAT
