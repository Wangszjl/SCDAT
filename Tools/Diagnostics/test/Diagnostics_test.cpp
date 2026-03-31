#include "ConservationMonitor.h"
#include "FieldMonitor.h"
#include "MonitorManager.h"
#include "ParticleMonitor.h"

#include "ParticleManager.h"

#include <gtest/gtest.h>

using SCDAT::Diagnostics::ConservationMonitor;
using SCDAT::Diagnostics::FieldMonitor;
using SCDAT::Diagnostics::MonitorManager;
using SCDAT::Diagnostics::ParticleMonitor;
using SCDAT::Geometry::Point3D;
using SCDAT::Geometry::Vector3D;
using SCDAT::Particle::ParticleManager;

TEST(DiagnosticsTest, FieldMonitorCapturesMagnitudeAndPotentialRange)
{
    FieldMonitor monitor(
        "field",
        [] { return std::vector<Vector3D>{Vector3D(1.0, 0.0, 0.0), Vector3D(0.0, 2.0, 0.0)}; },
        [] { return std::vector<double>{-2.0, 3.0}; });

    ASSERT_TRUE(monitor.sample(1.0));
    const auto* sample = monitor.latestSample();
    ASSERT_NE(sample, nullptr);
    EXPECT_DOUBLE_EQ(sample->scalars.at("potential_min_v"), -2.0);
    EXPECT_DOUBLE_EQ(sample->scalars.at("potential_max_v"), 3.0);
    EXPECT_DOUBLE_EQ(sample->scalars.at("field_max_v_per_m"), 2.0);
}

TEST(DiagnosticsTest, ParticleMonitorReadsParticleManagerStatistics)
{
    ParticleManager manager;
    manager.createElectron(Point3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0), 2.0);

    ParticleMonitor monitor("particles", [&manager] { return manager.getStatistics(); });
    ASSERT_TRUE(monitor.sample(0.0));
    const auto* sample = monitor.latestSample();
    ASSERT_NE(sample, nullptr);
    EXPECT_EQ(sample->scalars.at("particles_active"), 1.0);
    EXPECT_NEAR(sample->scalars.at("total_charge_c"), -1.602176634e-19, 1.0e-30);
}

TEST(DiagnosticsTest, MonitorManagerAggregatesLatestSamples)
{
    double charge = 1.0;
    double energy = 2.0;
    auto conservation =
        std::make_shared<ConservationMonitor>("cons", [&charge] { return charge; },
                                              [&energy] { return energy; });

    MonitorManager manager;
    manager.addMonitor(conservation);

    ASSERT_TRUE(manager.sampleAll(0.0));
    charge = 0.8;
    energy = 2.4;
    ASSERT_TRUE(manager.sampleAll(1.0));

    const auto scalars = manager.latestScalars();
    EXPECT_NEAR(scalars.at("cons.charge_drift_c"), -0.2, 1.0e-12);
    EXPECT_NEAR(scalars.at("cons.energy_drift_j"), 0.4, 1.0e-12);
}
