/**
 * @file CollisonTypes.cpp
 * @brief 粒子-粒子相互作用模型与碰撞截面实现
 */

#include "../include/CollisonTypes.h"
#include "../../../Basic/include/Constants.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace SCDAT
{
namespace Collision
{

namespace
{
constexpr double kEps = 1e-30;

inline double kineticEnergyEv(const ParticleObject& p)
{
    const double v = p.getVelocity().magnitude();
    return 0.5 * p.getMass() * v * v / SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge;
}

inline double atomicMassToKg(double amu)
{
    return amu * SCDAT::Basic::Constants::PhysicsConstants::AtomicMassUnit;
}
} // namespace

CollisionCrossSection::CollisionCrossSection(CollisionType type, const CollisionParameters& params)
    : collision_type_(type), params_(params)
{
}

bool CollisionCrossSection::isEnergyAboveThreshold(double energy_ev) const
{
    return energy_ev >= params_.energy_threshold;
}

CollisionType CollisionCrossSection::getCollisionType() const
{
    return collision_type_;
}

const CollisionParameters& CollisionCrossSection::getParameters() const
{
    return params_;
}

ElasticCollisionCrossSection::ElasticCollisionCrossSection(const CollisionParameters& params)
    : CollisionCrossSection(CollisionType::ELASTIC, params)
{
}

double ElasticCollisionCrossSection::calculateCrossSection(double energy_ev) const
{
    if (!isEnergyAboveThreshold(energy_ev))
    {
        return 0.0;
    }

    const double reference_energy =
        std::max(1e-6, params_.additional_params.count("reference_energy")
                           ? params_.additional_params.at("reference_energy")
                           : 1.0);
    const double power_index = params_.additional_params.count("power_index")
                                   ? params_.additional_params.at("power_index")
                                   : 0.5;

    const double sigma =
        params_.cross_section_max * std::pow(reference_energy / energy_ev, power_index);
    return std::max(0.0, std::min(sigma, params_.cross_section_max));
}

IonizationCollisionCrossSection::IonizationCollisionCrossSection(const CollisionParameters& params)
    : CollisionCrossSection(CollisionType::IONIZATION, params)
{
}

double IonizationCollisionCrossSection::calculateCrossSection(double energy_ev) const
{
    if (!isEnergyAboveThreshold(energy_ev) || energy_ev <= params_.energy_threshold)
    {
        return 0.0;
    }

    const double I = std::max(1e-6, params_.energy_threshold);
    const double ln_term = std::log(energy_ev / I);
    const double threshold_term = 1.0 - I / energy_ev;
    const double sigma = params_.cross_section_max * ln_term * threshold_term / energy_ev;
    return std::max(0.0, std::min(sigma, params_.cross_section_max));
}

ExcitationCollisionCrossSection::ExcitationCollisionCrossSection(const CollisionParameters& params)
    : CollisionCrossSection(CollisionType::EXCITATION, params)
{
}

double ExcitationCollisionCrossSection::calculateCrossSection(double energy_ev) const
{
    if (!isEnergyAboveThreshold(energy_ev) || energy_ev <= params_.energy_threshold)
    {
        return 0.0;
    }

    const double Eth = std::max(1e-6, params_.energy_threshold);
    const double alpha =
        params_.additional_params.count("alpha") ? params_.additional_params.at("alpha") : 1.2;
    const double beta =
        params_.additional_params.count("beta") ? params_.additional_params.at("beta") : 0.35;

    const double x = energy_ev / Eth;
    const double rise = std::pow(std::max(0.0, 1.0 - 1.0 / x), alpha);
    const double tail = std::pow(1.0 / x, beta);
    const double sigma = params_.cross_section_max * rise * tail;
    return std::max(0.0, std::min(sigma, params_.cross_section_max));
}

RecombinationCollisionCrossSection::RecombinationCollisionCrossSection(
    const CollisionParameters& params)
    : CollisionCrossSection(CollisionType::RECOMBINATION, params)
{
}

double RecombinationCollisionCrossSection::calculateCrossSection(double energy_ev) const
{
    const double reference_energy =
        std::max(1e-6, params_.additional_params.count("reference_energy")
                           ? params_.additional_params.at("reference_energy")
                           : 1.0);
    const double power_index = params_.additional_params.count("power_index")
                                   ? params_.additional_params.at("power_index")
                                   : 0.7;

    const double E = std::max(reference_energy, energy_ev);
    const double sigma = params_.cross_section_max * std::pow(reference_energy / E, power_index);
    return std::max(0.0, std::min(sigma, params_.cross_section_max));
}

CoulombCollisionCrossSection::CoulombCollisionCrossSection(const CollisionParameters& params)
    : CollisionCrossSection(CollisionType::COULOMB, params), electron_density_(1e16),
      electron_temperature_(1000.0), coulomb_logarithm_(10.0)
{
    if (params.additional_params.count("electron_density"))
    {
        electron_density_ = params.additional_params.at("electron_density");
    }
    if (params.additional_params.count("electron_temperature"))
    {
        electron_temperature_ = params.additional_params.at("electron_temperature");
    }
    coulomb_logarithm_ = calculateCoulombLogarithm(electron_temperature_, electron_density_);
}

double CoulombCollisionCrossSection::calculateCrossSection(double energy_ev) const
{
    if (energy_ev <= 0.0)
    {
        return 0.0;
    }

    const double E = energy_ev * SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge;
    const double e = SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge;
    const double m = SCDAT::Basic::Constants::PhysicsConstants::ElectronMass;
    const double eps0 = SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity;
    const double k = 1.0 / (4.0 * M_PI * eps0);

    const double sigma = M_PI * coulomb_logarithm_ * std::pow(k * e * e, 2.0) /
                         std::max(kEps, std::pow(2.0 * E, 2.0));
    return std::max(0.0, std::min(sigma, params_.cross_section_max));
}

double CoulombCollisionCrossSection::calculateDebyeLength(double temperature_k, double density_m3)
{
    if (temperature_k <= 0.0 || density_m3 <= 0.0)
    {
        return 1e-3;
    }

    const double e = SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge;
    const double eps0 = SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity;
    const double kB = SCDAT::Basic::Constants::PhysicsConstants::BoltzmannConstant;
    return std::sqrt(eps0 * kB * temperature_k / (density_m3 * e * e));
}

double CoulombCollisionCrossSection::calculateCoulombLogarithm(double temperature_k,
                                                               double density_m3)
{
    const double lambda_d = calculateDebyeLength(temperature_k, density_m3);
    const double bmin = std::max(1e-12, lambda_d * 1e-4);
    const double ln_lambda = std::log(std::max(1.000001, lambda_d / bmin));
    return std::clamp(ln_lambda, 1.0, 30.0);
}

HardSpherePairCalculator::HardSpherePairCalculator(double radius1, double radius2)
    : radius1_(radius1), radius2_(radius2)
{
}

double HardSpherePairCalculator::calculateCrossSection(const ParticleObject&, const ParticleObject&,
                                                       double) const
{
    const double r = radius1_ + radius2_;
    return M_PI * r * r;
}

CoulombPairCalculator::CoulombPairCalculator(double impact_parameter_limit)
    : impact_parameter_limit_(impact_parameter_limit)
{
}

double CoulombPairCalculator::calculateCrossSection(const ParticleObject& p1,
                                                    const ParticleObject& p2,
                                                    double relative_velocity) const
{
    if (std::abs(p1.getCharge()) < 1e-20 || std::abs(p2.getCharge()) < 1e-20)
    {
        return 0.0;
    }
    if (relative_velocity < 1e-10)
    {
        return 0.0;
    }
    // 计算约化质量和相对动能，基于经典库仑碰撞理论计算碰撞截面。引入 impact_parameter_limit_ 作为数值上界，防止极低能情况下截面数值发散。
    const double m_red =
        (p1.getMass() * p2.getMass()) / std::max(kEps, p1.getMass() + p2.getMass());
    // 计算相对动能（eV）
    const double e_kin = 0.5 * m_red * relative_velocity * relative_velocity;
    if (e_kin <= 0.0)
    {
        return 0.0;
    }
    // 计算冲量参数b，基于经典库仑碰撞理论：b = k * |q1 * q2| / (m_red * v_rel^2)，其中 k = 1/(4πε0)。最终截面为 σ = π * min(b^2, impact_parameter_limit_^2)。
    const double eps0 = SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity;
    const double k = 1.0 / (4.0 * M_PI * eps0);
    const double b = k * std::abs(p1.getCharge() * p2.getCharge()) / e_kin;
    return std::min(M_PI * b * b, M_PI * impact_parameter_limit_ * impact_parameter_limit_);
}

CollisionModel::CollisionModel(const std::string& name)
    : name_(name), enabled_(true), total_collisions_(0)
{
}

CollisionModel::~CollisionModel() = default;

void CollisionModel::setEnabled(bool enabled)
{
    enabled_ = enabled;
}

bool CollisionModel::isEnabled() const
{
    return enabled_;
}

const std::string& CollisionModel::getName() const
{
    return name_;
}

size_t CollisionModel::getTotalCollisions() const
{
    return total_collisions_;
}

void CollisionModel::resetStatistics()
{
    total_collisions_ = 0;
}

ElasticCollisionModel::ElasticCollisionModel(const ParticleTypeInfo& type1,
                                             const ParticleTypeInfo& type2, double cross_section)
    : CollisionModel("Elastic_" + type1.name + "_" + type2.name), type1_(type1), type2_(type2),
      cross_section_(cross_section)
{
}

size_t ElasticCollisionModel::processCollisions(std::vector<ParticleObject>& particles,
                                                std::vector<ParticleObject>&, double dt,
                                                double volume)
{
    if (!enabled_ || particles.empty() || volume <= 0.0)
    {
        return 0;
    }

    std::vector<size_t> lhs;
    std::vector<size_t> rhs;
    for (size_t i = 0; i < particles.size(); ++i)
    {
        if (particles[i].getType() == type1_.type)
        {
            lhs.push_back(i);
        }
        else if (particles[i].getType() == type2_.type)
        {
            rhs.push_back(i);
        }
    }

    if (lhs.empty() || rhs.empty())
    {
        return 0;
    }

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const double n2 = static_cast<double>(rhs.size()) / volume;
    size_t count = 0;

    for (size_t i : lhs)
    {
        for (size_t j : rhs)
        {
            if (i == j)
            {
                continue;
            }

            const double vrel =
                (particles[i].getVelocity() - particles[j].getVelocity()).magnitude();
            const double p = std::clamp(cross_section_ * vrel * n2 * dt, 0.0, 1.0);
            if (u(rng) < p)
            {
                performElasticCollision(particles[i], particles[j]);
                ++count;
                ++total_collisions_;
            }
        }
    }

    return count;
}

void ElasticCollisionModel::performElasticCollision(ParticleObject& p1, ParticleObject& p2)
{
    const double m1 = atomicMassToKg(type1_.mass_amu);
    const double m2 = atomicMassToKg(type2_.mass_amu);

    const Vector3D vcm = (p1.getVelocity() * m1 + p2.getVelocity() * m2) / (m1 + m2);
    const Vector3D vrel = p1.getVelocity() - p2.getVelocity();
    const double mag = vrel.magnitude();
    if (mag < 1e-12)
    {
        return;
    }

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double cos_theta = 2.0 * u(rng) - 1.0;
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    const double phi = 2.0 * M_PI * u(rng);

    const Vector3D new_rel(mag * sin_theta * std::cos(phi), mag * sin_theta * std::sin(phi),
                           mag * cos_theta);

    p1.setVelocity(vcm + new_rel * (m2 / (m1 + m2)));
    p2.setVelocity(vcm - new_rel * (m1 / (m1 + m2)));
}

IonizationCollisionModel::IonizationCollisionModel(const ParticleTypeInfo& electron_type,
                                                   const ParticleTypeInfo& neutral_type,
                                                   const ParticleTypeInfo& ion_type,
                                                   double cross_section,
                                                   double ionization_energy_ev)
    : CollisionModel("Ionization_" + neutral_type.name), electron_type_(electron_type),
      neutral_type_(neutral_type), ion_type_(ion_type), cross_section_(cross_section),
      ionization_energy_ev_(ionization_energy_ev)
{
}

size_t IonizationCollisionModel::processCollisions(std::vector<ParticleObject>& particles,
                                                   std::vector<ParticleObject>& new_particles,
                                                   double dt, double volume)
{
    if (!enabled_ || particles.empty() || volume <= 0.0)
    {
        return 0;
    }

    std::vector<size_t> electrons;
    std::vector<size_t> neutrals;
    for (size_t i = 0; i < particles.size(); ++i)
    {
        if (particles[i].getType() == electron_type_.type)
        {
            electrons.push_back(i);
        }
        else if (particles[i].getType() == neutral_type_.type)
        {
            neutrals.push_back(i);
        }
    }

    if (electrons.empty() || neutrals.empty())
    {
        return 0;
    }

    const double n = static_cast<double>(neutrals.size()) / volume;
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::uniform_int_distribution<size_t> pick_neutral(0, neutrals.size() - 1);

    size_t count = 0;
    for (size_t ie : electrons)
    {
        ParticleObject& e = particles[ie];
        if (kineticEnergyEv(e) < ionization_energy_ev_)
        {
            continue;
        }

        const double v = e.getVelocity().magnitude();
        const double p = std::clamp(cross_section_ * v * n * dt, 0.0, 1.0);
        if (u(rng) < p && !neutrals.empty())
        {
            const size_t in = neutrals[pick_neutral(rng)];
            if (performIonization(e, particles[in], new_particles))
            {
                particles[in].setStatus(SCDAT::Particle::ParticleStatus::ABSORBED);
                ++count;
                ++total_collisions_;
            }
        }
    }

    return count;
}

bool IonizationCollisionModel::performIonization(ParticleObject& electron,
                                                 const ParticleObject& neutral,
                                                 std::vector<ParticleObject>& new_particles)
{
    const double e_mass = atomicMassToKg(electron_type_.mass_amu);
    const double E_ev = kineticEnergyEv(electron);
    if (E_ev <= ionization_energy_ev_)
    {
        return false;
    }

    const double remaining_ev = E_ev - ionization_energy_ev_;
    const double each_ev = 0.5 * remaining_ev;
    const double each_j = each_ev * SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge;
    const double speed = std::sqrt(2.0 * each_j / std::max(kEps, e_mass));

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> u(0.0, 1.0);

    auto random_direction = [&]()
    {
        const double cos_theta = 2.0 * u(rng) - 1.0;
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        const double phi = 2.0 * M_PI * u(rng);
        return Vector3D(speed * sin_theta * std::cos(phi), speed * sin_theta * std::sin(phi),
                        speed * cos_theta);
    };

    electron.setVelocity(random_direction());

    static SCDAT::Particle::ParticleId next_id = 1000000;
    new_particles.push_back(SCDAT::Particle::ParticleFactory::createElectron(
        next_id++, neutral.getPosition(), random_direction(), neutral.getWeight()));
    new_particles.push_back(SCDAT::Particle::ParticleFactory::createIon(
        next_id++, neutral.getPosition(), neutral.getVelocity(), 1, 1, neutral.getWeight()));
    return true;
}

ChargeExchangeCollisionModel::ChargeExchangeCollisionModel(const ParticleTypeInfo& ion_type,
                                                           const ParticleTypeInfo& neutral_type,
                                                           double cross_section)
    : CollisionModel("ChargeExchange_" + ion_type.name + "_" + neutral_type.name),
      ion_type_(ion_type), neutral_type_(neutral_type), cross_section_(cross_section)
{
}

size_t ChargeExchangeCollisionModel::processCollisions(std::vector<ParticleObject>& particles,
                                                       std::vector<ParticleObject>&, double dt,
                                                       double volume)
{
    if (!enabled_ || particles.empty() || volume <= 0.0)
    {
        return 0;
    }

    std::vector<size_t> ions;
    std::vector<size_t> neutrals;
    for (size_t i = 0; i < particles.size(); ++i)
    {
        if (particles[i].getType() == ion_type_.type)
        {
            ions.push_back(i);
        }
        else if (particles[i].getType() == neutral_type_.type)
        {
            neutrals.push_back(i);
        }
    }
    if (ions.empty() || neutrals.empty())
    {
        return 0;
    }

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::uniform_int_distribution<size_t> pick_neutral(0, neutrals.size() - 1);

    const double n = static_cast<double>(neutrals.size()) / volume;
    size_t count = 0;
    for (size_t ii : ions)
    {
        const double v = particles[ii].getVelocity().magnitude();
        const double p = std::clamp(cross_section_ * v * n * dt, 0.0, 0.1);
        if (u(rng) < p)
        {
            performChargeExchange(particles[ii], particles[neutrals[pick_neutral(rng)]]);
            ++count;
            ++total_collisions_;
        }
    }

    return count;
}

void ChargeExchangeCollisionModel::performChargeExchange(ParticleObject& ion,
                                                         ParticleObject& neutral)
{
    const Vector3D v_ion = ion.getVelocity();
    const Vector3D v_neutral = neutral.getVelocity();

    ion.setVelocity(v_neutral);
    neutral.setVelocity(v_ion);
}

RecombinationCollisionModel::RecombinationCollisionModel(const ParticleTypeInfo& electron_type,
                                                         const ParticleTypeInfo& ion_type,
                                                         const ParticleTypeInfo& neutral_type,
                                                         double cross_section)
    : CollisionModel("Recombination_" + ion_type.name), electron_type_(electron_type),
      ion_type_(ion_type), neutral_type_(neutral_type), cross_section_(cross_section)
{
}

size_t RecombinationCollisionModel::processCollisions(std::vector<ParticleObject>& particles,
                                                      std::vector<ParticleObject>&, double dt,
                                                      double volume)
{
    if (!enabled_ || particles.empty() || volume <= 0.0)
    {
        return 0;
    }

    std::vector<size_t> electrons;
    std::vector<size_t> ions;
    for (size_t i = 0; i < particles.size(); ++i)
    {
        if (particles[i].getType() == electron_type_.type)
        {
            electrons.push_back(i);
        }
        else if (particles[i].getType() == ion_type_.type)
        {
            ions.push_back(i);
        }
    }
    if (electrons.empty() || ions.empty())
    {
        return 0;
    }

    const double ni = static_cast<double>(ions.size()) / volume;
    const double alpha = cross_section_ * std::pow(1000.0, -1.5); // 计算复合反应速率系数，假设平均速率为1000 m/s，截面随能量的-1.5次方衰减。根据反应速率系数计算碰撞概率，并限制在合理范围内。
    const double p = std::clamp(alpha * ni * dt, 0.0, 0.05);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> u(0.0, 1.0);

    size_t count = 0;
    const size_t pairs = std::min(electrons.size(), ions.size());
    for (size_t k = 0; k < pairs; ++k)
    {
        if (u(rng) < p)
        {
            particles[electrons[k]].setStatus(SCDAT::Particle::ParticleStatus::ABSORBED);
            particles[ions[k]].setStatus(SCDAT::Particle::ParticleStatus::ABSORBED);
            ++count;
            ++total_collisions_;
        }
    }

    return count;
}

CollisionModelManager::CollisionModelManager() = default;
CollisionModelManager::~CollisionModelManager() = default;

void CollisionModelManager::addModel(std::shared_ptr<CollisionModel> model)
{
    if (!model)
    {
        throw std::invalid_argument("CollisionModelManager: model must not be null");
    }

    for (const auto& m : models_)
    {
        if (m->getName() == model->getName())
        {
            throw std::invalid_argument("CollisionModelManager: duplicate model name");
        }
    }
    models_.push_back(std::move(model));
}

void CollisionModelManager::removeModel(const std::string& name)
{
    auto it = std::remove_if(models_.begin(), models_.end(),
                             [&name](const std::shared_ptr<CollisionModel>& m)
                             { return m->getName() == name; });
    models_.erase(it, models_.end());
}

std::shared_ptr<CollisionModel> CollisionModelManager::getModel(const std::string& name)
{
    auto it = std::find_if(models_.begin(), models_.end(),
                           [&name](const std::shared_ptr<CollisionModel>& m)
                           { return m->getName() == name; });
    return it == models_.end() ? nullptr : *it;
}

const std::vector<std::shared_ptr<CollisionModel>>& CollisionModelManager::getModels() const
{
    return models_;
}

void CollisionModelManager::setModelEnabled(const std::string& name, bool enabled)
{
    auto model = getModel(name);
    if (model)
    {
        model->setEnabled(enabled);
    }
}

std::vector<std::string> CollisionModelManager::getModelNames() const
{
    std::vector<std::string> names;
    names.reserve(models_.size());
    for (const auto& m : models_)
    {
        names.push_back(m->getName());
    }
    return names;
}

size_t CollisionModelManager::getTotalCollisions() const
{
    size_t total = 0;
    for (const auto& m : models_)
    {
        total += m->getTotalCollisions();
    }
    return total;
}

void CollisionModelManager::resetStatistics()
{
    for (auto& m : models_)
    {
        m->resetStatistics();
    }
}

void CollisionModelManager::clearModels()
{
    models_.clear();
}

} // namespace Collision
} // namespace SCDAT
