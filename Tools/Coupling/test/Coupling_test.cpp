#include "CouplingBase.h"
#include "MultiPhysicsManager.h"

#include <gtest/gtest.h>

using SCDAT::Coupling::ConvergenceCriteria;
using SCDAT::Coupling::CouplingType;
using SCDAT::Coupling::FunctionalCoupling;
using SCDAT::Coupling::MultiPhysicsManager;

TEST(CouplingTest, IterativeManagerDrivesFunctionalCouplingToConvergence)
{
    double residual = 1.0;
    auto coupling = std::make_shared<FunctionalCoupling>(
        "test", CouplingType::ITERATIVE,
        [&residual](double) {
            residual *= 0.25;
            return SCDAT::VoidResult::success();
        },
        [&residual] { return std::pair<double, double>{residual, residual}; });

    ConvergenceCriteria criteria;
    criteria.relative_tolerance = 0.01;
    criteria.absolute_tolerance = 0.01;
    criteria.max_iterations = 20;
    criteria.min_iterations = 1;
    coupling->setConvergenceCriteria(criteria);

    MultiPhysicsManager manager;
    manager.addCoupling(coupling);

    ASSERT_TRUE(manager.executeIterativeCouplings(1.0));
    const auto stats = manager.getAllStatistics();
    ASSERT_TRUE(stats.contains("test"));
    EXPECT_TRUE(stats.at("test").converged);
    EXPECT_LT(stats.at("test").final_absolute_error, 0.01);
}
