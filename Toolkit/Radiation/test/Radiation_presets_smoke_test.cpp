#include "RadiationCases.h"
#include "RadiationDoseAlgorithm.h"

#include <gtest/gtest.h>

#include <cmath>

using SCDAT::Toolkit::Radiation::RadiationDoseAlgorithm;
using SCDAT::Toolkit::Radiation::RadiationScenarioPreset;

TEST(RadiationPresetsSmokeTest, PresetsRunWithConservation)
{
    const auto names = SCDAT::Toolkit::Radiation::listRadiationScenarioPresetNames();
    ASSERT_FALSE(names.empty());

    for (const auto& name : names)
    {
        RadiationScenarioPreset preset;
        ASSERT_TRUE(SCDAT::Toolkit::Radiation::tryGetRadiationScenarioPreset(name, preset));

        RadiationDoseAlgorithm algorithm;
        ASSERT_TRUE(algorithm.initialize(preset.config));

        for (std::size_t i = 0; i < preset.steps; ++i)
        {
            ASSERT_TRUE(algorithm.advance(preset.time_step_s));
        }

        const auto& status = algorithm.getStatus();
        EXPECT_GT(status.time_s, 0.0);
        if (preset.config.particle_flux_m2_s > 0.0 && preset.config.mean_energy_ev > 0.0)
        {
            EXPECT_GT(status.incident_energy_j_per_m2, 0.0);
        }
        else
        {
            EXPECT_GE(status.incident_energy_j_per_m2, 0.0);
        }
        EXPECT_GE(status.deposited_energy_j_per_m2, 0.0);
        EXPECT_GE(status.escaped_energy_j_per_m2, 0.0);
        EXPECT_LE(status.deposited_energy_j_per_m2, status.incident_energy_j_per_m2 + 1.0e-12);
        EXPECT_LE(status.energy_conservation_error, 1.0e-9);

        ASSERT_EQ(status.layers.size(), preset.config.layers);
        EXPECT_TRUE(std::isfinite(status.layers.front().dose_gy));
        EXPECT_TRUE(std::isfinite(status.layers.back().dose_gy));
    }
}
