#pragma once

#include "DensePlasmaBoundaryLayer.h"
#include "DensePlasmaDetector.h"
#include "DiffusionSolverAdapter.h"
#include "PlasmaAdvancedClosureModel.h"
#include "PlasmaReactionCollisionLibrary.h"

#include "../../Tools/FieldSolver/include/BoltzmannSolver.h"
#include "../../Tools/FieldSolver/include/NonlinearPoissonSolver.h"
#include "../../Tools/Output/include/ResultTypes.h"

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

class FluidAlgorithmAdapter
{
  public:
    bool initialize(const FluidAlgorithmConfig& config);
    bool advance(double dt);
    void reset();

    const FluidAlgorithmStatus& getStatus() const { return status_; }
    const FluidAlgorithmConfig& getConfig() const { return config_; }
    Output::ColumnarDataSet buildProfileDataSet() const;

  private:
    FluidAlgorithmConfig config_;
    FluidAlgorithmStatus status_;
    DensePlasmaDetector detector_;
    DensePlasmaBoundaryLayer boundary_layer_;
    PlasmaReactionCollisionLibrary reaction_collision_library_;
    PlasmaAdvancedClosureModel advanced_closure_model_;
    DiffusionSolverAdapter diffusion_solver_;
    FieldSolver::BoltzmannSolver boltzmann_solver_;
    FieldSolver::NonlinearPoissonSolver poisson_solver_;
    PlasmaParameters plasma_state_;
    DensePlasmaAssessment assessment_;
    BoundaryLayerState boundary_state_;
    ReactionCollisionState reaction_state_;
    AdvancedClosureState closure_state_;
    bool initialized_ = false;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
