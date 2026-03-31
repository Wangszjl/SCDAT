# Toolkit Scenario Presets

## Command line entry

Run:

```bash
SCDAT.exe presets
SCDAT.exe presets plasma
SCDAT.exe plasma ccp_argon_10pa
SCDAT.exe surface geo_eclipse_dielectric
SCDAT.exe internal geo_electron_belt
SCDAT.exe arc triple_junction_flashover
```

Each domain command accepts:

```bash
SCDAT.exe <module> [preset] [output_csv]
```

If no preset is given, the module default is used. If no output path is given, a preset-specific file under `results/` is used.

## Plasma Analysis presets

| preset | purpose | default output |
|---|---|---|
| `ccp_argon_10pa` | low-pressure CCP argon sheath and transport reference case | `results/plasma_ccp_argon_10pa.csv` |
| `leo_wake` | low Earth orbit wake-side plasma relaxation case | `results/plasma_leo_wake.csv` |
| `hall_thruster_plume` | near-field Hall thruster plume transport case | `results/plasma_hall_thruster_plume.csv` |

## Surface Charging presets

| preset | purpose | default output |
|---|---|---|
| `leo_daylight_kapton` | sunlit LEO dielectric panel charging case | `results/surface_leo_daylight_kapton.csv` |
| `geo_eclipse_dielectric` | GEO eclipse negative charging case | `results/surface_geo_eclipse_dielectric.csv` |
| `thruster_plume_dielectric` | plume-driven surface current balance case | `results/surface_thruster_plume_dielectric.csv` |

## Internal Charging presets

| preset | purpose | default output |
|---|---|---|
| `geo_electron_belt` | outer radiation belt dielectric charging case | `results/internal_geo_electron_belt.csv` |
| `meo_dielectric_harness` | MEO harness volume charging case | `results/internal_meo_dielectric_harness.csv` |
| `solar_array_backsheet` | solar-array backsheet accumulation case | `results/internal_solar_array_backsheet.csv` |

## Vacuum Arc presets

| preset | purpose | default output |
|---|---|---|
| `microgap_flashover` | compact micro-gap flashover growth case | `results/arc_microgap_flashover.csv` |
| `triple_junction_flashover` | triple-junction initiated flashover case | `results/arc_triple_junction_flashover.csv` |
| `restrike_recovery` | post-breakdown restrike and channel recovery case | `results/arc_restrike_recovery.csv` |

## Plotting

For PIC-MCC benchmark output:

```bash
python scripts/python/plot_pic_mcc_results.py --input results/pic_mcc_refactor_round.csv --output results/pic_mcc_refactor_round.png
```

For any toolkit csv output:

```bash
python scripts/python/plot_toolkit_results.py --input results/plasma_ccp_argon_10pa.csv --output results/plasma_ccp_argon_10pa.png
```
