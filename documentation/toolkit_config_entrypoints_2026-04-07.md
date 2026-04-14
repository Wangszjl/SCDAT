# Toolkit Config Entrypoints (2026-04-07)

## Commands

```powershell
SCDAT plasma-config <config_json> [output_csv]
SCDAT surface-config <config_json> [output_csv]
SCDAT internal-config <config_json> [output_csv]
```

## Example JSON Files

- [surface_surface_pic_hybrid_example.json](/e:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/documentation/examples/surface_surface_pic_hybrid_example.json)
- [plasma_spis_thruster_example.json](/e:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/documentation/examples/plasma_spis_thruster_example.json)
- [internal_spis_geant4_example.json](/e:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/documentation/examples/internal_spis_geant4_example.json)

## Surface Charging

Use `surface_runtime_route` or `runtime_route` to choose the route.

- `scdat_unified`
- `surface_pic`
- `surface_pic_hybrid`
- `legacy_benchmark`

Use `surface_pic_strategy` to control how surface-PIC currents enter the solver.

- `surface_pic_direct`
- `surface_pic_calibrated`
- `surface_pic_hybrid_reference`

Recommended surface-PIC fields:

- `enable_live_pic_window`
- `enable_live_pic_mcc`
- `enable_body_patch_circuit`
- `live_pic_collision_cross_section_set_id`

## Plasma Analysis

The new SPIS-style organization layer is driven by:

- `enable_spis_style_organization`
- `environment_model`
- `distribution_model`
- `reaction_registry`
- `diagnostic_set`

Typical values:

- `environment_model = spis_ccp_reference | spis_orbital_wake | spis_thruster_plume`
- `distribution_model = maxwellian_projected | wake_anisotropic | multi_population_hybrid`
- `reaction_registry = spis_core | spis_dense_plasma`
- `diagnostic_set = spis_core_diagnostics | sheath_multiscale_diagnostics | full_physics_diagnostics`

Plasma exports now also write a dedicated physics diagnostics artifact alongside the
unified CSV/VTK/H5 sidecars:

- `*.physics_diagnostics.json`
- `diagnostic_contract_id = plasma-physics-diagnostics-v1`
- `reaction_contract_id = plasma-reaction-balance-v1`
- `plasma_physics_diagnostics_artifact_path`

## Internal Charging

The new SPIS-style plus Geant4-style organization layer is driven by:

- `enable_spis_style_organization`
- `source_mode`
- `material_stack_model`
- `geometry_model`
- `primary_source_model`
- `physics_process_list`
- `energy_deposition_model`
- `charge_response_model`

Typical values:

- `source_mode = preset | radiation`
- `material_stack_model = spis_layered_stack | spis_harness_bundle | spis_backsheet_stack`
- `geometry_model = layer_stack_1d | shielded_layer_stack_1d`
- `primary_source_model = preset_monoenergetic_flux | radiation_drive_coupled`
- `physics_process_list = geant4_em_standard_like | geant4_shielding_like`
- `energy_deposition_model = continuous_slab_deposition | geant4_step_recorder_like`
- `charge_response_model = spis_layered_dielectric | radiation_induced_conductivity_relaxation`

## Contract

The unified schema file is:

- [config_schema_v1.json](/e:/3-Code/1-Cplusplus/PIC-Surface-Charging-master/documentation/contracts/config_schema_v1.json)

This schema now covers:

- `plasma`
- `surface`
- `internal`
- `radiation`
- `arc`

