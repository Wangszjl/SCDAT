#include "TimeStepController.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace SCDAT
{
namespace PICcore
{

TimeStepController::TimeStepController() : current_dt_(1e-12), next_dt_(1e-12), previous_dt_(1e-12)
{
}

TimeStepController::TimeStepController(const Parameters& params)
    : params_(params), current_dt_(params.initial_dt), next_dt_(params.initial_dt),
      previous_dt_(params.initial_dt)
{
}

void TimeStepController::initialize(double initial_dt)
{
    const double dt = (initial_dt > 0.0) ? initial_dt : params_.initial_dt;
    current_dt_ = dt;
    next_dt_ = dt;
    previous_dt_ = dt;

    resetStatistics();
    dt_history_.clear();
    previous_electric_field_.clear();
}

void TimeStepController::reset()
{
    initialize(params_.initial_dt);
}

double TimeStepController::computeTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                                           const std::vector<Utils::Vector3D>& electric_field,
                                           const std::vector<Utils::Vector3D>& magnetic_field,
                                           const std::vector<Mesh::ElementPtr>& elements)
{
    if (!params_.adaptive_enabled)
    {
        return current_dt_;
    }

    const double cfl_dt = computeCFLTimeStep(particles, elements);
    const double acc_dt =
        computeAccelerationTimeStep(particles, electric_field, magnetic_field, elements);

    double field_dt = params_.max_dt;
    if (!previous_electric_field_.empty() &&
        previous_electric_field_.size() == electric_field.size())
    {
        field_dt = computeFieldChangeTimeStep(electric_field, previous_electric_field_, elements);
    }

    double suggested_dt = std::min({cfl_dt, acc_dt, field_dt}) * params_.safety_factor;
    LimitType limit_type = getMostRestrictiveLimitType(cfl_dt, acc_dt, field_dt);

    updateTimeStep(suggested_dt, limit_type);
    previous_electric_field_ = electric_field;
    return next_dt_;
}

double TimeStepController::computeTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                                           const std::vector<Utils::Vector3D>& electric_field,
                                           const std::vector<Mesh::ElementPtr>& elements)
{
    std::vector<Utils::Vector3D> zero_b(electric_field.size(), Utils::Vector3D(0.0, 0.0, 0.0));
    return computeTimeStep(particles, electric_field, zero_b, elements);
}

double TimeStepController::computeCFLTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                                              const std::vector<Mesh::ElementPtr>& elements) const
{
    if (particles.empty() || elements.empty())
    {
        return params_.max_dt;
    }

    const double dx = computeCharacteristicLength(elements);
    const double vmax = computeMaxVelocity(particles);
    if (vmax <= 1e-30)
    {
        return params_.max_dt;
    }

    return clampTimeStep(params_.cfl_factor * dx / vmax);
}

double
TimeStepController::computeAccelerationTimeStep(const std::vector<Particle::ParticlePtr>& particles,
                                                const std::vector<Utils::Vector3D>& electric_field,
                                                const std::vector<Utils::Vector3D>& magnetic_field,
                                                const std::vector<Mesh::ElementPtr>& elements) const
{
    if (particles.empty() || elements.empty())
    {
        return params_.max_dt;
    }

    const double dx = computeCharacteristicLength(elements);
    const double amax = computeMaxAcceleration(particles, electric_field, magnetic_field);
    if (amax <= 1e-30)
    {
        return params_.max_dt;
    }

    return clampTimeStep(params_.acceleration_factor * std::sqrt(2.0 * dx / amax));
}

double
TimeStepController::computeFieldChangeTimeStep(const std::vector<Utils::Vector3D>& electric_field,
                                               const std::vector<Utils::Vector3D>& previous_field,
                                               const std::vector<Mesh::ElementPtr>& elements) const
{
    if (electric_field.size() != previous_field.size() || elements.empty())
    {
        return params_.max_dt;
    }

    const double dEdt = computeFieldChangeRate(electric_field, previous_field);
    if (dEdt <= 1e-30)
    {
        return params_.max_dt;
    }

    double eavg = 0.0;
    for (const auto& e : electric_field)
    {
        eavg += e.magnitude();
    }
    eavg /= static_cast<double>(electric_field.size());

    if (eavg <= 1e-30)
    {
        return params_.max_dt;
    }

    return clampTimeStep(params_.field_change_factor * eavg / dEdt);
}

void TimeStepController::updateTimeStep(double suggested_dt, LimitType limit_type)
{
    previous_dt_ = current_dt_;
    current_dt_ = next_dt_;

    const double max_increase = 1.5 * current_dt_;
    const double max_decrease = 0.5 * current_dt_;

    if (suggested_dt >= current_dt_)
    {
        next_dt_ = std::min(suggested_dt, max_increase);
    }
    else
    {
        next_dt_ = std::max(suggested_dt, max_decrease);
    }

    next_dt_ = clampTimeStep(next_dt_);
    updateStatistics(next_dt_, limit_type);

    dt_history_.push_back(next_dt_);
    if (dt_history_.size() > 1000)
    {
        dt_history_.erase(dt_history_.begin());
    }

    if (time_step_callback_)
    {
        time_step_callback_(next_dt_, limit_type);
    }
}

void TimeStepController::forceTimeStep(double dt)
{
    next_dt_ = clampTimeStep(dt);
    updateStatistics(next_dt_, LimitType::USER_DEFINED);
}

void TimeStepController::resetStatistics()
{
    stats_ = Statistics{};
    stats_.current_dt = current_dt_;
    stats_.average_dt = current_dt_;
    stats_.min_dt_used = current_dt_;
    stats_.max_dt_used = current_dt_;
}

void TimeStepController::printStatistics() const
{
    std::cout << std::scientific << std::setprecision(3)
              << "dt(current/avg/min/max) = " << stats_.current_dt << " / " << stats_.average_dt
              << " / " << stats_.min_dt_used << " / " << stats_.max_dt_used
              << ", steps=" << stats_.total_steps << '\n';
}

void TimeStepController::setTimeStepCallback(std::function<void(double, LimitType)> callback)
{
    time_step_callback_ = std::move(callback);
}

void TimeStepController::clearHistory()
{
    dt_history_.clear();
}

double TimeStepController::clampTimeStep(double dt) const
{
    return std::max(params_.min_dt, std::min(dt, params_.max_dt));
}

void TimeStepController::updateStatistics(double dt, LimitType /*limit_type*/)
{
    stats_.current_dt = dt;
    stats_.total_steps += 1;
    stats_.total_time += dt;
    stats_.min_dt_used = std::min(stats_.min_dt_used, dt);
    stats_.max_dt_used = std::max(stats_.max_dt_used, dt);
    stats_.average_dt = stats_.total_time / static_cast<double>(stats_.total_steps);

    if (dt < current_dt_ * 0.99)
    {
        stats_.dt_reductions += 1;
    }
    else if (dt > current_dt_ * 1.01)
    {
        stats_.dt_increases += 1;
    }
}

double
TimeStepController::computeCharacteristicLength(const std::vector<Mesh::ElementPtr>& elements) const
{
    double hmin = 1e30;
    for (const auto& e : elements)
    {
        if (!e)
        {
            continue;
        }
        hmin = std::min(hmin, std::pow(std::max(1e-30, e->getVolume()), 1.0 / 3.0));
    }
    return (hmin < 1e29) ? hmin : 1e-6;
}

double
TimeStepController::computeMaxVelocity(const std::vector<Particle::ParticlePtr>& particles) const
{
    double vmax = 0.0;
    for (const auto& p : particles)
    {
        if (p && p->isActive())
        {
            vmax = std::max(vmax, p->getVelocity().magnitude());
        }
    }
    return vmax;
}

double
TimeStepController::computeMaxAcceleration(const std::vector<Particle::ParticlePtr>& particles,
                                           const std::vector<Utils::Vector3D>& electric_field,
                                           const std::vector<Utils::Vector3D>& magnetic_field) const
{
    Utils::Vector3D avg_e(0.0, 0.0, 0.0);
    Utils::Vector3D avg_b(0.0, 0.0, 0.0);

    if (!electric_field.empty())
    {
        for (const auto& e : electric_field)
        {
            avg_e += e;
        }
        avg_e /= static_cast<double>(electric_field.size());
    }

    if (!magnetic_field.empty())
    {
        for (const auto& b : magnetic_field)
        {
            avg_b += b;
        }
        avg_b /= static_cast<double>(magnetic_field.size());
    }

    double amax = 0.0;
    for (const auto& p : particles)
    {
        if (!p || !p->isActive())
        {
            continue;
        }

        const double m = p->getMass();
        if (m <= 0.0)
        {
            continue;
        }

        const Utils::Vector3D v = p->getVelocity();
        const Utils::Vector3D f = p->getCharge() * (avg_e + v.cross(avg_b));
        amax = std::max(amax, f.magnitude() / m);
    }

    return amax;
}

double
TimeStepController::computeFieldChangeRate(const std::vector<Utils::Vector3D>& current_field,
                                           const std::vector<Utils::Vector3D>& previous_field) const
{
    if (current_field.size() != previous_field.size() || current_field.empty())
    {
        return 0.0;
    }

    const double dt = std::max(1e-30, current_dt_);
    double max_rate = 0.0;
    for (size_t i = 0; i < current_field.size(); ++i)
    {
        const double rate = (current_field[i] - previous_field[i]).magnitude() / dt;
        max_rate = std::max(max_rate, rate);
    }
    return max_rate;
}

double TimeStepController::adjustTimeStep(double suggested_dt, LimitType /*limit_type*/)
{
    return clampTimeStep(suggested_dt);
}

bool TimeStepController::isTimeStepStable(double dt) const
{
    return dt >= params_.min_dt && dt <= params_.max_dt;
}

TimeStepController::LimitType
TimeStepController::getMostRestrictiveLimitType(double cfl_dt, double acc_dt, double field_dt) const
{
    const double m = std::min({cfl_dt, acc_dt, field_dt});
    if (std::abs(m - cfl_dt) < 1e-15)
    {
        return LimitType::CFL_CONDITION;
    }
    if (std::abs(m - acc_dt) < 1e-15)
    {
        return LimitType::ACCELERATION;
    }
    return LimitType::FIELD_CHANGE;
}

AdaptiveTimeStepController::AdaptiveTimeStepController()
    : TimeStepController(), consecutive_reductions_(0)
{
}

AdaptiveTimeStepController::AdaptiveTimeStepController(const Parameters& params,
                                                       const AdaptiveParameters& adaptive_params)
    : TimeStepController(params), adaptive_params_(adaptive_params), consecutive_reductions_(0)
{
}

double AdaptiveTimeStepController::computeAdaptiveTimeStep(
    const std::vector<Particle::ParticlePtr>& particles,
    const std::vector<Utils::Vector3D>& electric_field,
    const std::vector<Mesh::ElementPtr>& elements, double estimated_error)
{
    const double base_dt = computeTimeStep(particles, electric_field, elements);
    if (!adaptive_params_.use_error_estimation)
    {
        return base_dt;
    }

    double factor = 1.0;
    if (estimated_error > adaptive_params_.error_tolerance)
    {
        factor = adaptive_params_.shrink_factor;
        consecutive_reductions_ += 1;
    }
    else if (estimated_error < 0.1 * adaptive_params_.error_tolerance)
    {
        factor = adaptive_params_.growth_factor;
        consecutive_reductions_ = 0;
    }

    if (consecutive_reductions_ > adaptive_params_.max_consecutive_reductions)
    {
        factor = std::max(factor, 0.9);
    }

    const double dt = clampTimeStep(base_dt * factor);
    error_history_.push_back(estimated_error);
    if (error_history_.size() > 100)
    {
        error_history_.erase(error_history_.begin());
    }

    return dt;
}

double
AdaptiveTimeStepController::estimateLocalError(const std::vector<Particle::ParticlePtr>& particles,
                                               const std::vector<Utils::Vector3D>& electric_field,
                                               double dt) const
{
    if (particles.empty())
    {
        return 0.0;
    }

    double max_err = 0.0;
    const size_t sample = std::min(static_cast<size_t>(100), particles.size());
    for (size_t i = 0; i < sample; ++i)
    {
        const auto& p = particles[i];
        if (!p || !p->isActive())
        {
            continue;
        }

        const Utils::Point3D pos = p->getPosition();
        const Utils::Vector3D vel = p->getVelocity();
        const double mass = p->getMass();
        if (mass <= 0.0)
        {
            continue;
        }

        const Utils::Vector3D e = interpolateFieldAtPosition(pos, electric_field);
        const Utils::Vector3D a = e * (p->getCharge() / mass);
        const double vel_mag = std::max(1e-30, vel.magnitude());
        const double err = a.magnitude() * dt * dt / vel_mag;
        max_err = std::max(max_err, err);
    }

    return max_err;
}

double AdaptiveTimeStepController::richardsonErrorEstimate(
    const std::vector<Particle::ParticlePtr>& particles,
    const std::vector<Utils::Vector3D>& electric_field, double dt) const
{
    if (particles.empty())
    {
        return 0.0;
    }

    double max_err = 0.0;
    const size_t sample = std::min(static_cast<size_t>(50), particles.size());

    for (size_t i = 0; i < sample; ++i)
    {
        const auto& p = particles[i];
        if (!p || !p->isActive())
        {
            continue;
        }

        Utils::Point3D x1;
        Utils::Vector3D v1;
        performSingleStep(p->getPosition(), p->getVelocity(), p->getCharge(), p->getMass(),
                          electric_field, dt, x1, v1);

        Utils::Point3D xh;
        Utils::Vector3D vh;
        performSingleStep(p->getPosition(), p->getVelocity(), p->getCharge(), p->getMass(),
                          electric_field, 0.5 * dt, xh, vh);

        Utils::Point3D x2;
        Utils::Vector3D v2;
        performSingleStep(xh, vh, p->getCharge(), p->getMass(), electric_field, 0.5 * dt, x2, v2);

        const Utils::Vector3D dx(x1.x() - x2.x(), x1.y() - x2.y(), x1.z() - x2.z());
        const Utils::Vector3D dv = v1 - v2;

        max_err = std::max(max_err, std::max(dx.magnitude(), dv.magnitude()) / 3.0);
    }

    return max_err;
}

Utils::Vector3D AdaptiveTimeStepController::interpolateFieldAtPosition(
    const Utils::Point3D& position, const std::vector<Utils::Vector3D>& electric_field) const
{
    if (electric_field.empty())
    {
        return Utils::Vector3D(0.0, 0.0, 0.0);
    }

    const size_t idx = static_cast<size_t>(std::abs(position.x() + position.y() + position.z())) %
                       electric_field.size();
    return electric_field[idx];
}

double AdaptiveTimeStepController::estimatePositionError(const Utils::Point3D& /*position*/,
                                                         const Utils::Vector3D& electric_field,
                                                         double dt) const
{
    const double g = electric_field.magnitude() / 1e-3;
    const double d = electric_field.magnitude() * dt * dt;
    return g * d * d;
}

void AdaptiveTimeStepController::performSingleStep(
    const Utils::Point3D& pos0, const Utils::Vector3D& vel0, double charge, double mass,
    const std::vector<Utils::Vector3D>& electric_field, double dt, Utils::Point3D& pos1,
    Utils::Vector3D& vel1) const
{
    if (mass <= 0.0)
    {
        pos1 = pos0;
        vel1 = vel0;
        return;
    }

    const Utils::Vector3D e0 = interpolateFieldAtPosition(pos0, electric_field);
    const Utils::Vector3D a0 = e0 * (charge / mass);

    const Utils::Vector3D vh = vel0 + a0 * (0.5 * dt);
    pos1 = Utils::Point3D(pos0.x() + vh.x() * dt, pos0.y() + vh.y() * dt, pos0.z() + vh.z() * dt);

    const Utils::Vector3D e1 = interpolateFieldAtPosition(pos1, electric_field);
    const Utils::Vector3D a1 = e1 * (charge / mass);
    vel1 = vh + a1 * (0.5 * dt);
}

} // namespace PICcore
} // namespace SCDAT
