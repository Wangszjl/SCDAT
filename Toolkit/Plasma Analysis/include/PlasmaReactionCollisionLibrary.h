#pragma once

#include "FluidAlgorithmConfig.h"

#include <memory>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace PlasmaAnalysis
{

struct ReactionCollisionState
{
    double ionization_source_m3_per_s = 0.0;
    double recombination_sink_m3_per_s = 0.0;
    double effective_collision_frequency_hz = 0.0;
    double charge_exchange_frequency_hz = 0.0;
    double electron_energy_loss_ev_per_s = 0.0;
    double momentum_transfer_ratio = 0.0;
    std::size_t active_processes = 0;
};

class IReactionCollisionProcess
{
  public:
    virtual ~IReactionCollisionProcess() = default;

    virtual void evaluate(const PlasmaParameters& parameters,
                          const DensePlasmaAssessment& assessment,
                          ReactionCollisionState& state) const = 0;
};

class PlasmaReactionCollisionLibrary
{
  public:
    PlasmaReactionCollisionLibrary();

    void configure(const ReactionCollisionConfig& config);
    const ReactionCollisionConfig& getConfig() const { return config_; }

    void registerProcess(std::shared_ptr<const IReactionCollisionProcess> process);
    void clearCustomProcesses();

    ReactionCollisionState evaluate(const PlasmaParameters& parameters,
                                    const DensePlasmaAssessment& assessment) const;

  private:
    void rebuildBuiltinProcesses();

    ReactionCollisionConfig config_;
    std::vector<std::shared_ptr<const IReactionCollisionProcess>> builtin_processes_;
    std::vector<std::shared_ptr<const IReactionCollisionProcess>> custom_processes_;
};

} // namespace PlasmaAnalysis
} // namespace Toolkit
} // namespace SCDAT
