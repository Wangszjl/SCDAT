#ifndef SCDAT_PICCORE_TIMESTEPCONTROLLER_H
#define SCDAT_PICCORE_TIMESTEPCONTROLLER_H

/**
 * @file TimeStepController.h
 * @ingroup FieldSolverModule
 */

#include "../../Geometry/include/Point3D.h"
#include "../../Geometry/include/Vector3D.h"
#include "../../Mesh/include/MeshAlgorithms.h"
#include "../../Mesh/include/MeshParsing.h"
#include "../../Mesh/include/MeshPartitioning.h"
#include "../../Particle/include/ParticleTypes.h"
#include <functional>
#include <memory>
#include <vector>

namespace SCDAT
{
namespace PICcore
{

class TimeStepController
{
  public:
    struct Parameters
    {
        double initial_dt = 1e-12;
        double min_dt = 1e-15;
        double max_dt = 1e-9;
        double cfl_factor = 0.5;
        double acceleration_factor = 0.1;
        double field_change_factor = 0.05;
        double safety_factor = 0.9;
        bool adaptive_enabled = true;
        bool verbose = false;
    };

    struct Statistics
    {
        double current_dt = 0.0;
        double average_dt = 0.0;
        double min_dt_used = 1e10;
        double max_dt_used = 0.0;
        int total_steps = 0;
        int dt_reductions = 0;
        int dt_increases = 0;
        double total_time = 0.0;
    };

    enum class LimitType
    {
        CFL_CONDITION,
        ACCELERATION,
        FIELD_CHANGE,
        USER_DEFINED,
        STABILITY
    };

    TimeStepController();
    explicit TimeStepController(const Parameters& params);
    virtual ~TimeStepController() = default;

    void setParameters(const Parameters& params)
    {
        params_ = params;
    }
    const Parameters& getParameters() const
    {
        return params_;
    }

    void initialize(double initial_dt = 0.0);
    void reset();

    double computeTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                           const std::vector<Utils::Vector3D>& electric_field,
                           const std::vector<Utils::Vector3D>& magnetic_field,
                           const std::vector<Mesh::ElementPtr>& elements);

    double computeTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                           const std::vector<Utils::Vector3D>& electric_field,
                           const std::vector<Mesh::ElementPtr>& elements);

    double computeCFLTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                              const std::vector<Mesh::ElementPtr>& elements) const;

    double computeAccelerationTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                                       const std::vector<Utils::Vector3D>& electric_field,
                                       const std::vector<Utils::Vector3D>& magnetic_field,
                                       const std::vector<Mesh::ElementPtr>& elements) const;

    double computeFieldChangeTimeStep(const std::vector<Utils::Vector3D>& electric_field,
                                      const std::vector<Utils::Vector3D>& previous_field,
                                      const std::vector<Mesh::ElementPtr>& elements) const;

    void updateTimeStep(double suggested_dt, LimitType limit_type);
    void forceTimeStep(double dt);

    const Statistics& getStatistics() const
    {
        return stats_;
    }
    void resetStatistics();
    void printStatistics() const;

    void setTimeStepCallback(std::function<void(double, LimitType)> callback);

    double getCurrentTimeStep() const
    {
        return current_dt_;
    }
    double getNextTimeStep() const
    {
        return next_dt_;
    }

    std::vector<double> getTimeStepHistory() const
    {
        return dt_history_;
    }
    void clearHistory();

  protected:
    double clampTimeStep(double dt) const;
    void updateStatistics(double dt, LimitType limit_type);

  private:
    Parameters params_;
    Statistics stats_;

    double current_dt_;
    double next_dt_;
    double previous_dt_;

    std::vector<double> dt_history_;
    std::vector<Utils::Vector3D> previous_electric_field_;

    std::function<void(double, LimitType)> time_step_callback_;

    double computeCharacteristicLength(const std::vector<Mesh::ElementPtr>& elements) const;
    double computeMaxVelocity(const std::vector<Particle::ParticlePtr>& particles) const;
    double computeMaxAcceleration(const std::vector<Particle::ParticlePtr>& particles,
                                  const std::vector<Utils::Vector3D>& electric_field,
                                  const std::vector<Utils::Vector3D>& magnetic_field) const;

    double computeFieldChangeRate(const std::vector<Utils::Vector3D>& current_field,
                                  const std::vector<Utils::Vector3D>& previous_field) const;

    double adjustTimeStep(double suggested_dt, LimitType limit_type);
    bool isTimeStepStable(double dt) const;

    LimitType getMostRestrictiveLimitType(double cfl_dt, double acc_dt, double field_dt) const;
};

using TimeStepControllerPtr = std::shared_ptr<TimeStepController>;

class AdaptiveTimeStepController : public TimeStepController
{
  public:
    struct AdaptiveParameters
    {
        double error_tolerance = 1e-6;
        double growth_factor = 1.2;
        double shrink_factor = 0.8;
        int max_consecutive_reductions = 5;
        bool use_error_estimation = true;
        bool use_richardson_extrapolation = false;
    };

    AdaptiveTimeStepController();
    AdaptiveTimeStepController(const Parameters& params, const AdaptiveParameters& adaptive_params);

    void setAdaptiveParameters(const AdaptiveParameters& params)
    {
        adaptive_params_ = params;
    }
    const AdaptiveParameters& getAdaptiveParameters() const
    {
        return adaptive_params_;
    }

    double computeAdaptiveTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                                   const std::vector<Utils::Vector3D>& electric_field,
                                   const std::vector<Mesh::ElementPtr>& elements,
                                   double estimated_error);

    double estimateLocalError(const std::vector<Particle::ParticlePtr>& particles,
                              const std::vector<Utils::Vector3D>& electric_field, double dt) const;

  private:
    AdaptiveParameters adaptive_params_;

    int consecutive_reductions_;
    std::vector<double> error_history_;

    double richardsonErrorEstimate(const std::vector<Particle::ParticlePtr>& particles,
                                   const std::vector<Utils::Vector3D>& electric_field,
                                   double dt) const;

    Utils::Vector3D
    interpolateFieldAtPosition(const Utils::Point3D& position,
                               const std::vector<Utils::Vector3D>& electric_field) const;

    double estimatePositionError(const Utils::Point3D& position,
                                 const Utils::Vector3D& electric_field, double dt) const;

    void performSingleStep(const Utils::Point3D& pos0, const Utils::Vector3D& vel0, double charge,
                           double mass, const std::vector<Utils::Vector3D>& electric_field,
                           double dt, Utils::Point3D& pos1, Utils::Vector3D& vel1) const;
};

using AdaptiveTimeStepControllerPtr = std::shared_ptr<AdaptiveTimeStepController>;

} // namespace PICcore
} // namespace SCDAT

#endif // SCDAT_PICCORE_TIMESTEPCONTROLLER_H
