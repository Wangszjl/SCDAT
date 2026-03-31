#include "../include/BoundaryConditions.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace SCDAT
{
namespace Particle
{

using Point3D = SCDAT::Geometry::Point3D;
using Vector3D = SCDAT::Geometry::Vector3D;

namespace
{
constexpr double kTolerance = 1.0e-12;

std::mt19937& globalRng()
{
    static std::mt19937 rng(std::random_device{}());
    return rng;
}
} // namespace

void BoundaryConditionBase::updateStatistics(const ParticleTypeDef& particle,
                                             const std::string& action)
{
    const double energy = 0.5 * particle.getMass() * particle.getVelocity().magnitudeSquared();
    const double charge = particle.getCharge();

    if (action == "reflected")
    {
        statistics_.particles_reflected++;
    }
    else if (action == "absorbed")
    {
        statistics_.particles_absorbed++;
    }
    else if (action == "emitted")
    {
        statistics_.particles_emitted++;
    }
    else if (action == "transmitted")
    {
        statistics_.particles_transmitted++;
    }

    statistics_.total_energy_flux += energy;
    statistics_.total_charge_flux += charge;
}

PeriodicBoundaryCondition::PeriodicBoundaryCondition(const std::vector<Vector3D>& period_vectors)
    : BoundaryConditionBase(AdvancedBoundaryType::PERIODIC_3D, BoundaryGeometry::PLANAR),
      period_vectors_(period_vectors), period_determinant_(0.0), is_orthogonal_(false)
{
    if (period_vectors_.empty() || period_vectors_.size() > 3)
    {
        throw std::invalid_argument("periodic boundary expects 1 to 3 period vectors");
    }

    updatePeriodMatrix();
}

std::vector<ParticleTypeDef>
PeriodicBoundaryCondition::processParticle(ParticleTypeDef& particle,
                                           const Point3D& /* intersection_point */,
                                           const Vector3D& /* normal */, double /* dt */)
{
    if (!enabled_)
    {
        return {};
    }

    particle.setPosition(mapPeriodicPosition(particle.getPosition()));
    updateStatistics(particle, "transmitted");
    return {};
}

void PeriodicBoundaryCondition::setPeriodVectors(const std::vector<Vector3D>& period_vectors)
{
    if (period_vectors.empty() || period_vectors.size() > 3)
    {
        throw std::invalid_argument("periodic boundary expects 1 to 3 period vectors");
    }

    period_vectors_ = period_vectors;
    updatePeriodMatrix();
}

Point3D PeriodicBoundaryCondition::mapPeriodicPosition(const Point3D& position) const
{
    if (period_vectors_.empty())
    {
        return position;
    }

    Point3D mapped_position = position;

    if (is_orthogonal_)
    {
        for (size_t i = 0; i < period_vectors_.size() && i < 3; ++i)
        {
            const double period_length = period_vectors_[i].magnitude();
            if (period_length <= kTolerance)
            {
                continue;
            }

            double coord = 0.0;
            if (i == 0)
            {
                coord = mapped_position.x();
            }
            else if (i == 1)
            {
                coord = mapped_position.y();
            }
            else
            {
                coord = mapped_position.z();
            }

            coord = std::fmod(coord, period_length);
            if (coord < 0.0)
            {
                coord += period_length;
            }

            if (i == 0)
            {
                mapped_position = Point3D(coord, mapped_position.y(), mapped_position.z());
            }
            else if (i == 1)
            {
                mapped_position = Point3D(mapped_position.x(), coord, mapped_position.z());
            }
            else
            {
                mapped_position = Point3D(mapped_position.x(), mapped_position.y(), coord);
            }
        }

        return mapped_position;
    }

    return findNearestPeriodicImage(position);
}

bool PeriodicBoundaryCondition::isValidPeriodicConfiguration() const
{
    if (period_vectors_.empty())
    {
        return false;
    }

    if (period_vectors_.size() >= 2)
    {
        const Vector3D cross = period_vectors_[0].cross(period_vectors_[1]);
        if (cross.magnitude() < kTolerance)
        {
            return false;
        }
    }

    if (period_vectors_.size() == 3)
    {
        const double det = period_vectors_[0].dot(period_vectors_[1].cross(period_vectors_[2]));
        if (std::abs(det) < kTolerance)
        {
            return false;
        }
    }

    return true;
}

void PeriodicBoundaryCondition::updatePeriodMatrix()
{
    is_orthogonal_ = true;
    for (size_t i = 0; i < period_vectors_.size(); ++i)
    {
        for (size_t j = i + 1; j < period_vectors_.size(); ++j)
        {
            if (std::abs(period_vectors_[i].dot(period_vectors_[j])) > kTolerance)
            {
                is_orthogonal_ = false;
                break;
            }
        }

        if (!is_orthogonal_)
        {
            break;
        }
    }

    period_matrix_.clear();
    inverse_period_matrix_.clear();
    reciprocal_vectors_.clear();
    period_determinant_ = 0.0;

    if (!is_orthogonal_ && period_vectors_.size() >= 2)
    {
        buildNonOrthogonalPeriodMatrix();
        calculateReciprocalLattice();
        setupLatticeTransformations();
    }
}

Point3D PeriodicBoundaryCondition::findNearestPeriodicImage(const Point3D& position) const
{
    Point3D best_position = position;
    double min_distance_sq = std::numeric_limits<double>::max();

    for (int i = -2; i <= 2; ++i)
    {
        for (int j = -2; j <= 2; ++j)
        {
            for (int k = -2; k <= 2; ++k)
            {
                Point3D candidate = position;

                if (period_vectors_.size() >= 1)
                {
                    candidate = candidate + period_vectors_[0] * static_cast<double>(i);
                }
                if (period_vectors_.size() >= 2)
                {
                    candidate = candidate + period_vectors_[1] * static_cast<double>(j);
                }
                if (period_vectors_.size() >= 3)
                {
                    candidate = candidate + period_vectors_[2] * static_cast<double>(k);
                }

                const double distance_sq = (candidate - Point3D(0.0, 0.0, 0.0)).magnitudeSquared();
                if (distance_sq < min_distance_sq)
                {
                    min_distance_sq = distance_sq;
                    best_position = candidate;
                }
            }
        }
    }

    return best_position;
}

ReflectiveBoundaryCondition::ReflectiveBoundaryCondition(double reflection_coefficient,
                                                         double energy_loss_factor)
    : BoundaryConditionBase(AdvancedBoundaryType::REFLECTIVE_CURVED, BoundaryGeometry::ARBITRARY),
      reflection_coefficient_(reflection_coefficient), energy_loss_factor_(energy_loss_factor)
{
    if (reflection_coefficient_ < 0.0 || reflection_coefficient_ > 1.0)
    {
        throw std::invalid_argument("reflection coefficient must be in [0, 1]");
    }

    if (energy_loss_factor_ < 0.0 || energy_loss_factor_ > 1.0)
    {
        throw std::invalid_argument("energy loss factor must be in [0, 1]");
    }
}

std::vector<ParticleTypeDef>
ReflectiveBoundaryCondition::processParticle(ParticleTypeDef& particle,
                                             const Point3D& intersection_point,
                                             const Vector3D& normal, double /* dt */)
{
    if (!enabled_)
    {
        return {};
    }

    static std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

    if (unit_dist(globalRng()) < reflection_coefficient_)
    {
        const Vector3D incident_velocity = particle.getVelocity();
        const Vector3D surface_normal =
            normal.magnitude() > kTolerance ? Vector3D(normal.normalized())
                                            : calculateSurfaceNormal(intersection_point);
        Vector3D reflected_velocity =
            calculateReflectedVelocity(incident_velocity, surface_normal);
        reflected_velocity = applyEnergyLoss(reflected_velocity);

        particle.setVelocity(reflected_velocity);
        particle.setPosition(intersection_point + surface_normal * 1.0e-10);
        updateStatistics(particle, "reflected");
        return {};
    }

    particle.setStatus(ParticleStatus::ABSORBED);
    updateStatistics(particle, "absorbed");
    return {};
}

void ReflectiveBoundaryCondition::setSurfaceFunction(
    std::function<double(const Point3D&)> surface_function,
    std::function<Vector3D(const Point3D&)> gradient_function)
{
    surface_function_ = std::move(surface_function);
    gradient_function_ = std::move(gradient_function);
}

Vector3D ReflectiveBoundaryCondition::calculateSurfaceNormal(const Point3D& point) const
{
    if (gradient_function_)
    {
        const Vector3D gradient = gradient_function_(point);
        const double magnitude = gradient.magnitude();
        if (magnitude > kTolerance)
        {
            return Vector3D(gradient / magnitude);
        }
    }

    return Vector3D(0.0, 0.0, 1.0);
}

Vector3D ReflectiveBoundaryCondition::calculateReflectedVelocity(
    const Vector3D& incident_velocity, const Vector3D& normal) const
{
    const double dot_product = incident_velocity.dot(normal);
    return Vector3D(incident_velocity - normal * (2.0 * dot_product));
}

Vector3D ReflectiveBoundaryCondition::applyEnergyLoss(const Vector3D& velocity) const
{
    if (energy_loss_factor_ <= 0.0)
    {
        return velocity;
    }

    if (velocity.magnitude() < kTolerance)
    {
        return velocity;
    }

    const double speed_factor = std::sqrt(1.0 - energy_loss_factor_);
    return Vector3D(velocity * speed_factor);
}

AbsorbingBoundaryCondition::AbsorbingBoundaryCondition(double absorption_probability)
    : BoundaryConditionBase(AdvancedBoundaryType::ABSORBING_STATISTICAL,
                            BoundaryGeometry::ARBITRARY),
      absorption_probability_(absorption_probability), total_absorbed_energy_(0.0),
      total_absorbed_charge_(0.0), absorption_count_(0)
{
    if (absorption_probability_ < 0.0 || absorption_probability_ > 1.0)
    {
        throw std::invalid_argument("absorption probability must be in [0, 1]");
    }
}

std::vector<ParticleTypeDef>
AbsorbingBoundaryCondition::processParticle(ParticleTypeDef& particle,
                                            const Point3D& /* intersection_point */,
                                            const Vector3D& normal, double /* dt */)
{
    if (!enabled_)
    {
        return {};
    }

    static std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

    if (unit_dist(globalRng()) < absorption_probability_)
    {
        const double energy = 0.5 * particle.getMass() * particle.getVelocity().magnitudeSquared();
        const double charge = particle.getCharge();

        total_absorbed_energy_ += energy;
        total_absorbed_charge_ += charge;
        absorption_count_++;

        const double angle = calculateIncidenceAngle(particle.getVelocity(), normal);
        statistics_.average_angle =
            (statistics_.average_angle * static_cast<double>(absorption_count_ - 1) + angle) /
            static_cast<double>(absorption_count_);

        particle.setStatus(ParticleStatus::ABSORBED);
        updateStatistics(particle, "absorbed");
        return {};
    }

    updateStatistics(particle, "transmitted");
    return {};
}

void AbsorbingBoundaryCondition::setAbsorptionProbability(double probability)
{
    if (probability < 0.0 || probability > 1.0)
    {
        throw std::invalid_argument("absorption probability must be in [0, 1]");
    }

    absorption_probability_ = probability;
}

double AbsorbingBoundaryCondition::getAverageAbsorptionEnergy() const
{
    if (absorption_count_ == 0)
    {
        return 0.0;
    }

    return total_absorbed_energy_ / static_cast<double>(absorption_count_);
}

double AbsorbingBoundaryCondition::getAverageIncidenceAngle() const
{
    return statistics_.average_angle;
}

void AbsorbingBoundaryCondition::resetAbsorptionStatistics()
{
    total_absorbed_energy_ = 0.0;
    total_absorbed_charge_ = 0.0;
    absorption_count_ = 0;
    resetStatistics();
}

double AbsorbingBoundaryCondition::calculateIncidenceAngle(const Vector3D& velocity,
                                                           const Vector3D& normal) const
{
    const double velocity_magnitude = velocity.magnitude();
    const double normal_magnitude = normal.magnitude();
    if (velocity_magnitude < kTolerance || normal_magnitude < kTolerance)
    {
        return 0.0;
    }

    const double cos_angle =
        std::abs(velocity.dot(normal)) / (velocity_magnitude * normal_magnitude);
    return std::acos(std::clamp(cos_angle, 0.0, 1.0));
}

ThermalEmissionBoundaryCondition::ThermalEmissionBoundaryCondition(double temperature,
                                                                   double work_function,
                                                                   double richardson_constant)
    : BoundaryConditionBase(AdvancedBoundaryType::EMISSION_THERMAL, BoundaryGeometry::ARBITRARY),
      temperature_(temperature), work_function_(work_function),
      richardson_constant_(richardson_constant), emission_area_(1.0e-4)
{
    if (temperature_ <= 0.0)
    {
        throw std::invalid_argument("temperature must be positive");
    }
    if (work_function_ <= 0.0)
    {
        throw std::invalid_argument("work function must be positive");
    }

    statistics_.surface_temperature = temperature_;
}

std::vector<ParticleTypeDef>
ThermalEmissionBoundaryCondition::processParticle(ParticleTypeDef& /* particle */,
                                                  const Point3D& intersection_point,
                                                  const Vector3D& normal, double dt)
{
    std::vector<ParticleTypeDef> emitted_particles;

    if (!enabled_)
    {
        return emitted_particles;
    }

    constexpr double kBoltzmann = 1.380649e-23;
    constexpr double kElementaryCharge = 1.602176634e-19;

    const double current_density = richardson_constant_ * temperature_ * temperature_ *
                                   std::exp(-work_function_ * kElementaryCharge /
                                            (kBoltzmann * temperature_));
    const double electrons_per_second = current_density * emission_area_ / kElementaryCharge;
    const double expected_electrons = std::max(0.0, electrons_per_second * dt);

    std::poisson_distribution<int> poisson_dist(expected_electrons);
    const int electrons_to_emit = poisson_dist(globalRng());

    for (int i = 0; i < electrons_to_emit; ++i)
    {
        ParticleTypeDef electron = generateThermalElectron(intersection_point, normal);
        emitted_particles.push_back(electron);
        updateStatistics(emitted_particles.back(), "emitted");
    }

    return emitted_particles;
}

void ThermalEmissionBoundaryCondition::setTemperature(double temperature)
{
    if (temperature <= 0.0)
    {
        throw std::invalid_argument("temperature must be positive");
    }

    temperature_ = temperature;
    statistics_.surface_temperature = temperature_;
}

void ThermalEmissionBoundaryCondition::setEmissionArea(double area)
{
    if (area <= 0.0)
    {
        throw std::invalid_argument("emission area must be positive");
    }

    emission_area_ = area;
}

double ThermalEmissionBoundaryCondition::getCurrentDensity() const
{
    constexpr double kBoltzmann = 1.380649e-23;
    constexpr double kElementaryCharge = 1.602176634e-19;

    return richardson_constant_ * temperature_ * temperature_ *
           std::exp(-work_function_ * kElementaryCharge / (kBoltzmann * temperature_));
}

ParticleTypeDef ThermalEmissionBoundaryCondition::generateThermalElectron(
    const Point3D& position, const Vector3D& normal)
{
    const ParticleId id = ParticleFactory::getNextId();
    ParticleTypeDef thermal_electron = ParticleFactory::createThermalElectron(
        id, position + normal * 1.0e-10, Vector3D(0.0, 0.0, 0.0), temperature_, work_function_,
        1.0);

    thermal_electron.setVelocity(generateMaxwellBoltzmannVelocity(normal));
    thermal_electron.setStatus(ParticleStatus::ACTIVE);
    return thermal_electron;
}

Vector3D
ThermalEmissionBoundaryCondition::generateMaxwellBoltzmannVelocity(const Vector3D& normal)
{
    constexpr double kBoltzmann = 1.380649e-23;
    constexpr double kElectronMass = 9.10938356e-31;

    const double thermal_velocity = std::sqrt(kBoltzmann * temperature_ / kElectronMass);
    std::normal_distribution<double> velocity_dist(0.0, thermal_velocity);

    const double vx = velocity_dist(globalRng());
    const double vy = velocity_dist(globalRng());
    const double vz = std::abs(velocity_dist(globalRng()));

    Vector3D velocity(vx, vy, vz);
    if (std::abs(normal.z() - 1.0) > 1.0e-6)
    {
        velocity = Vector3D(vx, vy, std::abs(vz));
    }

    return velocity;
}

void PeriodicBoundaryCondition::buildNonOrthogonalPeriodMatrix()
{
    period_matrix_.assign(3, std::vector<double>(3, 0.0));
    for (int i = 0; i < 3; ++i)
    {
        period_matrix_[i][i] = 1.0;
    }

    for (size_t i = 0; i < period_vectors_.size() && i < 3; ++i)
    {
        period_matrix_[0][i] = period_vectors_[i].x();
        period_matrix_[1][i] = period_vectors_[i].y();
        period_matrix_[2][i] = period_vectors_[i].z();
    }

    period_determinant_ = calculateMatrixDeterminant(period_matrix_);
    if (std::abs(period_determinant_) < kTolerance)
    {
        throw std::runtime_error("period vectors are linearly dependent");
    }
}

void PeriodicBoundaryCondition::calculateReciprocalLattice()
{
    reciprocal_vectors_.assign(period_vectors_.size(), Vector3D(0.0, 0.0, 0.0));

    if (period_vectors_.size() == 1)
    {
        const double length = period_vectors_[0].magnitude();
        if (length > kTolerance)
        {
            reciprocal_vectors_[0] = Vector3D(period_vectors_[0] / (length * length));
        }
        return;
    }

    if (period_vectors_.size() == 2)
    {
        const Vector3D a1 = period_vectors_[0];
        const Vector3D a2 = period_vectors_[1];
        const Vector3D normal = a1.cross(a2);
        const double area = normal.magnitude();
        if (area > kTolerance)
        {
            reciprocal_vectors_[0] = Vector3D(a2.cross(normal) / area);
            reciprocal_vectors_[1] = Vector3D(normal.cross(a1) / area);
        }
        return;
    }

    const Vector3D a1 = period_vectors_[0];
    const Vector3D a2 = period_vectors_[1];
    const Vector3D a3 = period_vectors_[2];
    const double volume = std::abs(a1.dot(a2.cross(a3)));
    if (volume > kTolerance)
    {
        reciprocal_vectors_[0] = Vector3D(a2.cross(a3) / volume);
        reciprocal_vectors_[1] = Vector3D(a3.cross(a1) / volume);
        reciprocal_vectors_[2] = Vector3D(a1.cross(a2) / volume);
    }
}

void PeriodicBoundaryCondition::setupLatticeTransformations()
{
    if (period_vectors_.size() < 2)
    {
        return;
    }

    inverse_period_matrix_ = calculateMatrixInverse(period_matrix_);
    const auto identity_check = multiplyMatrices(period_matrix_, inverse_period_matrix_);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            const double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(identity_check[i][j] - expected) > 1.0e-10)
            {
                throw std::runtime_error("failed to build inverse periodic matrix");
            }
        }
    }
}

double PeriodicBoundaryCondition::calculateMatrixDeterminant(
    const std::vector<std::vector<double>>& matrix) const
{
    if (matrix.size() != 3 || matrix[0].size() != 3 || matrix[1].size() != 3 ||
        matrix[2].size() != 3)
    {
        return 0.0;
    }

    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

std::vector<std::vector<double>> PeriodicBoundaryCondition::calculateMatrixInverse(
    const std::vector<std::vector<double>>& matrix) const
{
    const double det = calculateMatrixDeterminant(matrix);
    if (std::abs(det) < kTolerance)
    {
        throw std::runtime_error("matrix is singular");
    }

    std::vector<std::vector<double>> inverse(3, std::vector<double>(3, 0.0));

    inverse[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det;
    inverse[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det;
    inverse[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;

    inverse[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det;
    inverse[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
    inverse[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / det;

    inverse[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det;
    inverse[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / det;
    inverse[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det;

    return inverse;
}

std::vector<std::vector<double>>
PeriodicBoundaryCondition::multiplyMatrices(const std::vector<std::vector<double>>& A,
                                            const std::vector<std::vector<double>>& B) const
{
    std::vector<std::vector<double>> result(3, std::vector<double>(3, 0.0));
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

} // namespace Particle
} // namespace SCDAT
