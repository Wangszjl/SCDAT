/**
 * @file ParticleDefinitions.cpp
 * @brief 粒子系统实现
 *
 * @author Wang Sizhan
 * @date 2026年3月19日 8:24:51
 */

#include "../include/ParticleDefinitions.h"
#include <algorithm>
#include <cmath>
#include <execution>
#include <format>
#include <iomanip>
#include <sstream>

namespace SCDAT
{
namespace Particle
{

using SCDAT::Basic::Constants::PhysicsConstants;

ParticleId ParticleFactory::next_id_ = 0;

// ============================================================================
// ModernParticle 实现 / ModernParticle Implementation
// ============================================================================

void Particle::push(const SCDAT::Geometry::Vector3D& E_field,
                    const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    // 根据粒子类型选择推进算法 / Choose push algorithm based on particle type
    if (data_.type == ParticleType::PHOTON)
    {
        // 光子直线运动 / Photon straight-line motion
        SCDAT::Geometry::Point3D new_pos = getPosition() + getVelocity() * dt;
        setPosition(new_pos);
    }
    else if (getSpeed() > 0.1 * PhysicsConstants::SpeedOfLight)
    { // 10% 光速
        // 相对论性推进 / Relativistic push
        relativisticPush(E_field, B_field, dt);
    }
    else
    {
        // 经典Boris推进 / Classical Boris push
        borisPush(E_field, B_field, dt);
    }

    updateTime(data_.last_update_time + dt);
    updateAge(dt);
}

void Particle::borisPush(const SCDAT::Geometry::Vector3D& E_field,
                         const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    if (data_.mass <= 0.0 || data_.status != ParticleStatus::ACTIVE)
    {
        return;
    }

    const double q_over_m = data_.charge / data_.mass;
    const double half_dt = 0.5 * dt;

    // 第一步：电场加速 v- = v + (q/m) * E * dt/2
    SCDAT::Geometry::Vector3D v_minus = getVelocity() + E_field * (q_over_m * half_dt);

    // 第二步：磁场旋转
    SCDAT::Geometry::Vector3D t = B_field * (q_over_m * half_dt);
    double t_magnitude_sq = t.magnitudeSquared();

    SCDAT::Geometry::Vector3D v_plus;
    if (t_magnitude_sq > 1e-30)
    {
        SCDAT::Geometry::Vector3D s = t * (2.0 / (1.0 + t_magnitude_sq));
        SCDAT::Geometry::Vector3D v_prime = v_minus + v_minus.cross(t);
        v_plus = v_minus + v_prime.cross(s);
    }
    else
    {
        v_plus = v_minus;
    }

    // 第三步：电场加速 v+ = v+ + (q/m) * E * dt/2
    SCDAT::Geometry::Vector3D new_velocity = v_plus + E_field * (q_over_m * half_dt);

    // 位置更新：x_new = x_old + v_new * dt
    SCDAT::Geometry::Point3D new_position = getPosition() + new_velocity * dt;

    setVelocity(new_velocity);
    setPosition(new_position);
}

void Particle::relativisticPush(const SCDAT::Geometry::Vector3D& E_field,
                                const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    // 相对论性Boris算法实现
    const double c = PhysicsConstants::SpeedOfLight;
    const double q_over_m = data_.charge / data_.mass;
    const double half_dt = 0.5 * dt;

    // 当前动量
    SCDAT::Geometry::Vector3D p = getVelocity() * data_.mass;
    double gamma = std::sqrt(1.0 + p.magnitudeSquared() / (data_.mass * data_.mass * c * c));

    // 电场推进
    SCDAT::Geometry::Vector3D p_minus = p + E_field * (data_.charge * half_dt);

    // 磁场旋转
    SCDAT::Geometry::Vector3D t = B_field * (q_over_m * half_dt / gamma);
    double t_magnitude_sq = t.magnitudeSquared();

    SCDAT::Geometry::Vector3D p_plus;
    if (t_magnitude_sq > 1e-30)
    {
        SCDAT::Geometry::Vector3D s = t * (2.0 / (1.0 + t_magnitude_sq));
        SCDAT::Geometry::Vector3D p_prime = p_minus + p_minus.cross(t);
        p_plus = p_minus + p_prime.cross(s);
    }
    else
    {
        p_plus = p_minus;
    }

    // 最终电场推进
    SCDAT::Geometry::Vector3D p_new = p_plus + E_field * (data_.charge * half_dt);

    // 更新速度和位置
    double gamma_new =
        std::sqrt(1.0 + p_new.magnitudeSquared() / (data_.mass * data_.mass * c * c));
    SCDAT::Geometry::Vector3D v_new = p_new / (data_.mass * gamma_new);
    SCDAT::Geometry::Point3D pos_new = getPosition() + v_new * dt;

    setVelocity(v_new);
    setPosition(pos_new);
}

double Particle::getTemperatureEV() const noexcept
{
    constexpr double joule_to_ev = 1.0 / 1.602176634e-19;
    const double kinetic_energy_ev = getKineticEnergy() * joule_to_ev;
    return (2.0 / 3.0) * kinetic_energy_ev;
}

std::string Particle::getSpeciesName() const
{
    if (data_.type != ParticleType::ION && data_.type != ParticleType::POSITIVE_ION &&
        data_.type != ParticleType::NEGATIVE_ION && data_.type != ParticleType::MOLECULAR_ION &&
        data_.type != ParticleType::CLUSTER_ION && data_.type != ParticleType::HEAVY_ION)
    {
        return "Unknown";
    }

    if (data_.mass_number == 1 && data_.charge_number == 1)
    {
        return "H+";
    }
    if (data_.mass_number == 4 && data_.charge_number == 2)
    {
        return "He++";
    }
    if (data_.mass_number == 16 && data_.charge_number == 1)
    {
        return "O+";
    }

    return std::format("Ion(A={},Z={})", data_.mass_number, data_.charge_number);
}

double Particle::getMaxKineticEnergy() const noexcept
{
    if (data_.type == ParticleType::PHOTOELECTRON)
    {
        return std::max(0.0, data_.photon_energy - data_.work_function);
    }
    return 0.0;
}

std::string Particle::toString() const
{
    std::string type_str;
    switch (data_.type)
    {
    case ParticleType::UNKNOWN:
        type_str = "Unknown";
        break;
    case ParticleType::ELECTRON:
        type_str = "Electron";
        break;
    case ParticleType::POSITRON:
        type_str = "Positron";
        break;
    case ParticleType::MUON_MINUS:
        type_str = "Muon-";
        break;
    case ParticleType::MUON_PLUS:
        type_str = "Muon+";
        break;
    case ParticleType::TAU_MINUS:
        type_str = "Tau-";
        break;
    case ParticleType::TAU_PLUS:
        type_str = "Tau+";
        break;
    case ParticleType::ELECTRON_NEUTRINO:
        type_str = "ElectronNeutrino";
        break;
    case ParticleType::MUON_NEUTRINO:
        type_str = "MuonNeutrino";
        break;
    case ParticleType::TAU_NEUTRINO:
        type_str = "TauNeutrino";
        break;
    case ParticleType::PHOTON:
        type_str = "Photon";
        break;
    case ParticleType::GAMMA:
        type_str = "Gamma";
        break;
    case ParticleType::XRAY:
        type_str = "XRay";
        break;
    case ParticleType::OPTICAL_PHOTON:
        type_str = "OpticalPhoton";
        break;
    case ParticleType::PROTON:
        type_str = "Proton";
        break;
    case ParticleType::ANTIPROTON:
        type_str = "Antiproton";
        break;
    case ParticleType::NEUTRON:
        type_str = "Neutron";
        break;
    case ParticleType::ANTINEUTRON:
        type_str = "Antineutron";
        break;
    case ParticleType::PION_PLUS:
        type_str = "Pion+";
        break;
    case ParticleType::PION_MINUS:
        type_str = "Pion-";
        break;
    case ParticleType::PION_ZERO:
        type_str = "Pion0";
        break;
    case ParticleType::KAON_PLUS:
        type_str = "Kaon+";
        break;
    case ParticleType::KAON_MINUS:
        type_str = "Kaon-";
        break;
    case ParticleType::ALPHA:
        type_str = "Alpha";
        break;
    case ParticleType::DEUTERON:
        type_str = "Deuteron";
        break;
    case ParticleType::TRITON:
        type_str = "Triton";
        break;
    case ParticleType::HELION:
        type_str = "Helion";
        break;
    case ParticleType::HEAVY_ION:
        type_str = "HeavyIon";
        break;
    case ParticleType::ATOM:
        type_str = "Atom";
        break;
    case ParticleType::MOLECULE:
        type_str = "Molecule";
        break;
    case ParticleType::ION:
        type_str = "Ion";
        break;
    case ParticleType::RADICAL:
        type_str = "Radical";
        break;
    case ParticleType::CLUSTER:
        type_str = "Cluster";
        break;
    case ParticleType::NANOPARTICLE:
        type_str = "Nanoparticle";
        break;
    case ParticleType::DUST:
        type_str = "Dust";
        break;
    case ParticleType::DROPLET:
        type_str = "Droplet";
        break;
    case ParticleType::BUBBLE:
        type_str = "Bubble";
        break;
    case ParticleType::PHONON:
        type_str = "Phonon";
        break;
    case ParticleType::PLASMON:
        type_str = "Plasmon";
        break;
    case ParticleType::MAGNON:
        type_str = "Magnon";
        break;
    case ParticleType::EXCITON:
        type_str = "Exciton";
        break;
    case ParticleType::POLARON:
        type_str = "Polaron";
        break;
    case ParticleType::HOLE:
        type_str = "Hole";
        break;
    case ParticleType::COOPER_PAIR:
        type_str = "CooperPair";
        break;
    case ParticleType::VIRTUAL_PHOTON:
        type_str = "VirtualPhoton";
        break;
    case ParticleType::GEANTINO:
        type_str = "Geantino";
        break;
    case ParticleType::TRACER:
        type_str = "Tracer";
        break;
    case ParticleType::PHOTOELECTRON:
        type_str = "Photoelectron";
        break;
    case ParticleType::SECONDARY_ELECTRON:
        type_str = "SecondaryElectron";
        break;
    case ParticleType::BACKSCATTERED_ELECTRON:
        type_str = "BackscatteredElectron";
        break;
    case ParticleType::AUGER_ELECTRON:
        type_str = "AugerElectron";
        break;
    case ParticleType::THERMAL_ELECTRON:
        type_str = "ThermalElectron";
        break;
    case ParticleType::FIELD_EMISSION_ELECTRON:
        type_str = "FieldEmissionElectron";
        break;
    case ParticleType::POSITIVE_ION:
        type_str = "PositiveIon";
        break;
    case ParticleType::NEGATIVE_ION:
        type_str = "NegativeIon";
        break;
    case ParticleType::MOLECULAR_ION:
        type_str = "MolecularIon";
        break;
    case ParticleType::CLUSTER_ION:
        type_str = "ClusterIon";
        break;
    case ParticleType::BETA_PARTICLE:
        type_str = "BetaParticle";
        break;
    case ParticleType::USER_DEFINED_1:
        type_str = "UserDefined1";
        break;
    case ParticleType::USER_DEFINED_2:
        type_str = "UserDefined2";
        break;
    case ParticleType::USER_DEFINED_3:
        type_str = "UserDefined3";
        break;
    case ParticleType::GENERIC:
        type_str = "Generic";
        break;
    default:
        type_str = "UnknownType";
        break;
    }

    std::string status_str;
    switch (data_.status)
    {
    case ParticleStatus::UNINITIALIZED:
        status_str = "Uninitialized";
        break;
    case ParticleStatus::CREATED:
        status_str = "Created";
        break;
    case ParticleStatus::ACTIVE:
        status_str = "Active";
        break;
    case ParticleStatus::SUSPENDED:
        status_str = "Suspended";
        break;
    case ParticleStatus::WAITING:
        status_str = "Waiting";
        break;
    case ParticleStatus::ABSORBED:
        status_str = "Absorbed";
        break;
    case ParticleStatus::ESCAPED:
        status_str = "Escaped";
        break;
    case ParticleStatus::DEAD:
        status_str = "Dead";
        break;
    case ParticleStatus::LOST:
        status_str = "Lost";
        break;
    case ParticleStatus::DECAYED:
        status_str = "Decayed";
        break;
    case ParticleStatus::ANNIHILATED:
        status_str = "Annihilated";
        break;
    case ParticleStatus::CAPTURED:
        status_str = "Captured";
        break;
    case ParticleStatus::FUSED:
        status_str = "Fused";
        break;
    case ParticleStatus::FISSIONED:
        status_str = "Fissioned";
        break;
    case ParticleStatus::THERMALIZED:
        status_str = "Thermalized";
        break;
    case ParticleStatus::STUCK:
        status_str = "Stuck";
        break;
    case ParticleStatus::TIME_EXCEEDED:
        status_str = "TimeExceeded";
        break;
    case ParticleStatus::STEP_EXCEEDED:
        status_str = "StepExceeded";
        break;
    case ParticleStatus::RECOMBINED:
        status_str = "Recombined";
        break;
    case ParticleStatus::ERROR:
        status_str = "Error";
        break;
    case ParticleStatus::INVALID:
        status_str = "Invalid";
        break;
    default:
        status_str = "UnknownStatus";
        break;
    }

    std::string extra;
    if (data_.type == ParticleType::ELECTRON)
    {
        extra = std::format(", thermal_vel={:.3e}, temp_eV={:.3f}", data_.thermal_velocity,
                            getTemperatureEV());
    }
    else if (data_.type == ParticleType::ION || data_.type == ParticleType::POSITIVE_ION ||
             data_.type == ParticleType::NEGATIVE_ION ||
             data_.type == ParticleType::MOLECULAR_ION ||
             data_.type == ParticleType::CLUSTER_ION || data_.type == ParticleType::HEAVY_ION)
    {
        extra = std::format(", species={}, mass_num={}, charge_num={}", getSpeciesName(),
                            data_.mass_number, data_.charge_number);
    }
    else if (data_.type == ParticleType::PHOTOELECTRON)
    {
        extra = std::format(", photon_energy={:.3f}, work_function={:.3f}, max_KE={:.3f}",
                            data_.photon_energy, data_.work_function, getMaxKineticEnergy());
    }
    else if (data_.type == ParticleType::SECONDARY_ELECTRON)
    {
        extra = std::format(", primary_energy={:.3f}, yield={:.3f}, emission_angle={:.3f}",
                            data_.primary_energy, data_.yield, data_.emission_angle);
    }
    else if (data_.type == ParticleType::THERMAL_ELECTRON)
    {
        extra = std::format(", temperature={:.3f}, work_function={:.3f}, richardson={:.3e}",
                            data_.temperature, data_.work_function, data_.richardson_constant);
    }

    return std::format("Particle[id={}, type={}, status={}, pos=({:.3f},{:.3f},{:.3f}), "
                       "vel=({:.3e},{:.3e},{:.3e}), energy={:.3e}J, weight={:.3f}, age={:.3e}{}]",
                       data_.id, type_str, status_str, data_.position[0], data_.position[1],
                       data_.position[2], data_.velocity[0], data_.velocity[1], data_.velocity[2],
                       data_.energy, data_.weight, data_.age, extra);
}

Particle ParticleFactory::createElectron(ParticleId id,
                                         const SCDAT::Geometry::Point3D& position,
                                         const SCDAT::Geometry::Vector3D& velocity,
                                         double weight)
{
    Particle p(id, ParticleType::ELECTRON, position, velocity, PhysicsConstants::ElectronMass,
               PhysicsConstants::ElectronCharge, weight);
    return p;
}

Particle ParticleFactory::createIon(ParticleId id, const SCDAT::Geometry::Point3D& position,
                                    const SCDAT::Geometry::Vector3D& velocity, int mass_number,
                                    int charge_number, double weight)
{
    Particle p(id, ParticleType::ION, position, velocity,
               mass_number * PhysicsConstants::ProtonMass,
               charge_number * PhysicsConstants::ElementaryCharge, weight);
    p.setMassNumber(mass_number);
    p.setChargeNumber(charge_number);
    return p;
}

Particle ParticleFactory::createPhotoelectron(ParticleId id,
                                              const SCDAT::Geometry::Point3D& position,
                                              const SCDAT::Geometry::Vector3D& velocity,
                                              double photon_energy, double work_function,
                                              double weight)
{
    Particle p(id, ParticleType::PHOTOELECTRON, position, velocity,
               PhysicsConstants::ElectronMass, PhysicsConstants::ElectronCharge, weight);
    p.setPhotonEnergy(photon_energy);
    p.setWorkFunction(work_function);
    return p;
}

Particle ParticleFactory::createSecondaryElectron(ParticleId id,
                                                  const SCDAT::Geometry::Point3D& position,
                                                  const SCDAT::Geometry::Vector3D& velocity,
                                                  double primary_energy, double yield,
                                                  double weight)
{
    Particle p(id, ParticleType::SECONDARY_ELECTRON, position, velocity,
               PhysicsConstants::ElectronMass, PhysicsConstants::ElectronCharge, weight);
    p.setPrimaryEnergy(primary_energy);
    p.setYield(yield);
    return p;
}

Particle ParticleFactory::createThermalElectron(ParticleId id,
                                                const SCDAT::Geometry::Point3D& position,
                                                const SCDAT::Geometry::Vector3D& velocity,
                                                double temperature, double work_function,
                                                double weight)
{
    Particle p(id, ParticleType::THERMAL_ELECTRON, position, velocity,
               PhysicsConstants::ElectronMass, PhysicsConstants::ElectronCharge, weight);
    p.setTemperature(temperature);
    p.setWorkFunction(work_function);
    p.setRichardsonConstant(1.2e6);
    return p;
}

Particle ParticleFactory::createParticle(ParticleType type, ParticleId id,
                                         const SCDAT::Geometry::Point3D& position,
                                         const SCDAT::Geometry::Vector3D& velocity,
                                         double weight)
{
    switch (type)
    {
    case ParticleType::ELECTRON:
        return createElectron(id, position, velocity, weight);
    case ParticleType::ION:
        return createIon(id, position, velocity, 1, 1, weight);
    case ParticleType::PHOTOELECTRON:
        return createPhotoelectron(id, position, velocity, 2.5, 4.5, weight);
    case ParticleType::SECONDARY_ELECTRON:
        return createSecondaryElectron(id, position, velocity, 100.0, 1.0, weight);
    case ParticleType::THERMAL_ELECTRON:
        return createThermalElectron(id, position, velocity, 300.0, 4.5, weight);
    case ParticleType::BACKSCATTERED_ELECTRON:
    case ParticleType::AUGER_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
        return createElectron(id, position, velocity, weight);
    case ParticleType::POSITIVE_ION:
    case ParticleType::NEGATIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
        return createIon(id, position, velocity, 1, 1, weight);
    default:
        return Particle(id, type, position, velocity, PhysicsConstants::ElectronMass,
                        PhysicsConstants::ElectronCharge, weight);
    }
}

// ============================================================================
// ParticleContainer 实现
// ============================================================================

// Core::Result<size_t> ModernParticleContainer::addParticle(const ModernParticle& particle)
// {
//     return addParticle(particle.getData());
// }

// Core::Result<size_t> ModernParticleContainer::addParticle(ParticleData data)
// {
//     // 设置唯一ID / Set unique ID
//     data.id = next_particle_id_++;

//     size_t index = particles_.size();
//     particles_.push_back(std::move(data));

//     LOG_DEBUG("ModernParticleContainer", "添加粒子 ID: " + std::to_string(data.id) +
//                                              ", 总粒子数: " + std::to_string(particles_.size()));
//     return Core::Result<size_t>(index);
// }

// Core::VoidResult ModernParticleContainer::removeParticle(size_t index)
// {
//     if (index >= particles_.size())
//     {
//         return Core::ErrorHandler::makeError<void>(Core::ErrorCode::OUT_OF_RANGE);
//     }

//     // 使用swap-and-pop优化 / Use swap-and-pop optimization
//     if (index < particles_.size() - 1)
//     {
//         std::swap(particles_[index], particles_.back());
//     }
//     particles_.pop_back();

//     LOG_DEBUG("ModernParticleContainer", "移除粒子索引: " + std::to_string(index) +
//                                              ", 剩余粒子数: " +
//                                              std::to_string(particles_.size()));
//     return Core::ErrorHandler::makeSuccess();
// }

// Core::VoidResult ModernParticleContainer::removeParticleById(std::size_t id)
// {
//     auto it = std::find_if(particles_.begin(), particles_.end(),
//                            [id](const ParticleData& p) { return p.id == id; });

//     if (it == particles_.end())
//     {
//         return Core::ErrorHandler::makeError<void>(Core::ErrorCode::INVALID_ARGUMENT);
//     }

//     size_t index = static_cast<size_t>(std::distance(particles_.begin(), it));
//     return removeParticle(index);
// }

void ParticleContainer::borisPushAll(const SCDAT::Geometry::Vector3D& E_field,
                                     const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    // 使用SIMD优化的批量Boris推进
    if (particles_.empty())
        return;

// 检查是否支持SIMD
#ifdef __AVX2__
    if (particles_.size() >= 4)
    {
        simdBorisPush(std::span(particles_), E_field, B_field, dt);
        return;
    }
#endif

    // 回退到并行标量实现
    std::for_each(std::execution::par_unseq, particles_.begin(), particles_.end(),
                  [&](ParticleData& particle)
                  {
                      if (particle.status == ParticleStatus::ACTIVE && particle.mass > 0.0)
                      {
                          Particle p(particle);
                          p.borisPush(E_field, B_field, dt);
                          particle = p.getData();
                      }
                  });
}

void ParticleContainer::updatePositionsAll(double dt)
{
    std::for_each(std::execution::par_unseq, particles_.begin(), particles_.end(),
                  [dt](ParticleData& particle)
                  {
                      if (particle.status == ParticleStatus::ACTIVE)
                      {
                          particle.position[0] += particle.velocity[0] * dt;
                          particle.position[1] += particle.velocity[1] * dt;
                          particle.position[2] += particle.velocity[2] * dt;
                      }
                  });
}

void ParticleContainer::updateEnergiesAll()
{
    std::for_each(std::execution::par_unseq, particles_.begin(), particles_.end(),
                  [](ParticleData& particle) { particle.updateEnergy(); });
}

std::vector<size_t> ParticleContainer::findActiveParticles() const
{
    std::vector<size_t> active_indices;

    for (size_t i = 0; i < particles_.size(); ++i)
    {
        if (particles_[i].status == ParticleStatus::ACTIVE)
        {
            active_indices.push_back(i);
        }
    }

    return active_indices;
}

std::vector<size_t> ParticleContainer::findParticlesByType(ParticleType type) const
{
    std::vector<size_t> type_indices;

    for (size_t i = 0; i < particles_.size(); ++i)
    {
        if (particles_[i].type == type)
        {
            type_indices.push_back(i);
        }
    }

    return type_indices;
}

size_t ParticleContainer::getActiveParticleCount() const
{
    return static_cast<size_t>(std::count_if(std::execution::par_unseq, particles_.begin(),
                                             particles_.end(), [](const ParticleData& p)
                                             { return p.status == ParticleStatus::ACTIVE; }));
}

double ParticleContainer::getTotalKineticEnergy() const
{
    return std::transform_reduce(std::execution::par_unseq, particles_.begin(), particles_.end(),
                                 0.0, std::plus<>{}, [](const ParticleData& p)
                                 { return (p.status == ParticleStatus::ACTIVE) ? p.energy : 0.0; });
}

size_t ParticleContainer::getMemoryUsage() const
{
    return particles_.size() * sizeof(ParticleData) + particles_.capacity() * sizeof(ParticleData);
}

// void ParticleContainer::compactMemory()
// {
//     // 移除非活跃粒子 / Remove inactive particles
//     auto new_end = std::remove_if(particles_.begin(), particles_.end(), [](const ParticleData& p)
//                                   { return p.status != ParticleStatus::ACTIVE; });

//     size_t removed_count = static_cast<size_t>(std::distance(new_end, particles_.end()));
//     particles_.erase(new_end, particles_.end());

//     // 收缩容量 / Shrink capacity
//     particles_.shrink_to_fit();

//     LOG_INFO("ParticleContainer", "内存压缩完成，移除 " + std::to_string(removed_count) +
//                                             " 个非活跃粒子，剩余 " +
//                                             std::to_string(particles_.size()) + " 个");
// }

// Core::VoidResult ModernParticleContainer::validate() const
// {
//     LOG_INFO("ModernParticleContainer", "开始验证粒子容器");

//     for (size_t i = 0; i < particles_.size(); ++i)
//     {
//         const auto& particle = particles_[i];

//         // 验证物理量 / Validate physical quantities
//         if (particle.mass < 0.0)
//         {
//             LOG_ERROR("ModernParticleContainer", "粒子 " + std::to_string(particle.id) +
//                                                      " 质量为负: " + std::to_string(particle.mass));
//             return Core::ErrorHandler::makeError<void>(Core::ErrorCode::PARTICLE_INVALID_TYPE);
//         }

//         if (!std::isfinite(particle.energy))
//         {
//             LOG_ERROR("ModernParticleContainer",
//                       "粒子 " + std::to_string(particle.id) +
//                           " 能量无效: " + std::to_string(particle.energy));
//             return Core::ErrorHandler::makeError<void>(Core::ErrorCode::PARTICLE_ENERGY_INVALID);
//         }

//         // 验证位置和速度 / Validate position and velocity
//         for (size_t j = 0; j < 3; ++j)
//         {
//             if (!std::isfinite(particle.position[j]) || !std::isfinite(particle.velocity[j]))
//             {
//                 LOG_ERROR("ModernParticleContainer",
//                           "粒子 " + std::to_string(particle.id) + " 位置或速度无效");
//                 return Core::ErrorHandler::makeError<void>(Core::ErrorCode::PARTICLE_OUT_OF_BOUNDS);
//             }
//         }
//     }

//     LOG_INFO("ModernParticleContainer", "粒子容器验证完成，无错误");
//     return Core::ErrorHandler::makeSuccess();
// }

#ifdef __AVX2__
#include <immintrin.h>

void ParticleContainer::simdBorisPush(std::span<ParticleData> particles,
                                      const SCDAT::Geometry::Vector3D& E_field,
                                      const SCDAT::Geometry::Vector3D& B_field, double dt)
{
    // 完整的SIMD优化Boris推进实现 / Complete SIMD-optimized Boris push implementation

    const size_t simd_width = 4; // AVX2 可以处理4个double
    const size_t aligned_size = (particles.size() / simd_width) * simd_width;

    // 预计算常量 / Pre-compute constants
    const double half_dt = 0.5 * dt;

    // 广播电磁场到SIMD寄存器 / Broadcast E&M fields to SIMD registers
    __m256d E_x = _mm256_set1_pd(E_field.x());
    __m256d E_y = _mm256_set1_pd(E_field.y());
    __m256d E_z = _mm256_set1_pd(E_field.z());

    __m256d B_x = _mm256_set1_pd(B_field.x());
    __m256d B_y = _mm256_set1_pd(B_field.y());
    __m256d B_z = _mm256_set1_pd(B_field.z());

    __m256d half_dt_vec = _mm256_set1_pd(half_dt);
    __m256d dt_vec = _mm256_set1_pd(dt);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d two = _mm256_set1_pd(2.0);
    __m256d epsilon = _mm256_set1_pd(1e-30);

    // 处理对齐的部分 / Process aligned portion
    for (size_t i = 0; i < aligned_size; i += simd_width)
    {
        // 检查粒子状态并加载数据 / Check particle status and load data
        __m256i active_mask = _mm256_setzero_si256();
        __m256d mass_vec = _mm256_setzero_pd();
        __m256d charge_vec = _mm256_setzero_pd();
        __m256d q_over_m_vec = _mm256_setzero_pd();

        // 位置和速度向量 / Position and velocity vectors
        __m256d pos_x = _mm256_setzero_pd();
        __m256d pos_y = _mm256_setzero_pd();
        __m256d pos_z = _mm256_setzero_pd();

        __m256d vel_x = _mm256_setzero_pd();
        __m256d vel_y = _mm256_setzero_pd();
        __m256d vel_z = _mm256_setzero_pd();

        // 加载粒子数据 / Load particle data
        for (size_t j = 0; j < simd_width; ++j)
        {
            const auto& particle = particles[i + j];

            if (particle.status == ParticleStatus::ACTIVE && particle.mass > 0.0)
            {
                // 设置活跃掩码 / Set active mask
                reinterpret_cast<int64_t*>(&active_mask)[j] = -1;

                // 加载物理量 / Load physical quantities
                reinterpret_cast<double*>(&mass_vec)[j] = particle.mass;
                reinterpret_cast<double*>(&charge_vec)[j] = particle.charge;
                reinterpret_cast<double*>(&q_over_m_vec)[j] = particle.charge / particle.mass;

                // 加载位置和速度 / Load position and velocity
                reinterpret_cast<double*>(&pos_x)[j] = particle.position[0];
                reinterpret_cast<double*>(&pos_y)[j] = particle.position[1];
                reinterpret_cast<double*>(&pos_z)[j] = particle.position[2];

                reinterpret_cast<double*>(&vel_x)[j] = particle.velocity[0];
                reinterpret_cast<double*>(&vel_y)[j] = particle.velocity[1];
                reinterpret_cast<double*>(&vel_z)[j] = particle.velocity[2];
            }
        }

        // 执行SIMD Boris算法 / Execute SIMD Boris algorithm

        // 第一步：电场加速 v- = v + (q/m) * E * dt/2
        __m256d E_accel_x = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, E_x), half_dt_vec);
        __m256d E_accel_y = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, E_y), half_dt_vec);
        __m256d E_accel_z = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, E_z), half_dt_vec);

        __m256d v_minus_x = _mm256_add_pd(vel_x, E_accel_x);
        __m256d v_minus_y = _mm256_add_pd(vel_y, E_accel_y);
        __m256d v_minus_z = _mm256_add_pd(vel_z, E_accel_z);

        // 第二步：磁场旋转
        __m256d t_x = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, B_x), half_dt_vec);
        __m256d t_y = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, B_y), half_dt_vec);
        __m256d t_z = _mm256_mul_pd(_mm256_mul_pd(q_over_m_vec, B_z), half_dt_vec);

        // 计算 |t|² / Compute |t|²
        __m256d t_mag_sq =
            _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(t_x, t_x), _mm256_mul_pd(t_y, t_y)),
                          _mm256_mul_pd(t_z, t_z));

        // 计算 s = 2t / (1 + |t|²) / Compute s = 2t / (1 + |t|²)
        __m256d denominator = _mm256_add_pd(one, t_mag_sq);
        __m256d s_factor = _mm256_div_pd(two, denominator);

        __m256d s_x = _mm256_mul_pd(t_x, s_factor);
        __m256d s_y = _mm256_mul_pd(t_y, s_factor);
        __m256d s_z = _mm256_mul_pd(t_z, s_factor);

        // 计算 v' = v- + v- × t / Compute v' = v- + v- × t
        __m256d cross1_x =
            _mm256_sub_pd(_mm256_mul_pd(v_minus_y, t_z), _mm256_mul_pd(v_minus_z, t_y));
        __m256d cross1_y =
            _mm256_sub_pd(_mm256_mul_pd(v_minus_z, t_x), _mm256_mul_pd(v_minus_x, t_z));
        __m256d cross1_z =
            _mm256_sub_pd(_mm256_mul_pd(v_minus_x, t_y), _mm256_mul_pd(v_minus_y, t_x));

        __m256d v_prime_x = _mm256_add_pd(v_minus_x, cross1_x);
        __m256d v_prime_y = _mm256_add_pd(v_minus_y, cross1_y);
        __m256d v_prime_z = _mm256_add_pd(v_minus_z, cross1_z);

        // 计算 v+ = v- + v' × s / Compute v+ = v- + v' × s
        __m256d cross2_x =
            _mm256_sub_pd(_mm256_mul_pd(v_prime_y, s_z), _mm256_mul_pd(v_prime_z, s_y));
        __m256d cross2_y =
            _mm256_sub_pd(_mm256_mul_pd(v_prime_z, s_x), _mm256_mul_pd(v_prime_x, s_z));
        __m256d cross2_z =
            _mm256_sub_pd(_mm256_mul_pd(v_prime_x, s_y), _mm256_mul_pd(v_prime_y, s_x));

        __m256d v_plus_x = _mm256_add_pd(v_minus_x, cross2_x);
        __m256d v_plus_y = _mm256_add_pd(v_minus_y, cross2_y);
        __m256d v_plus_z = _mm256_add_pd(v_minus_z, cross2_z);

        // 第三步：最终电场加速 / Final electric field acceleration
        __m256d new_vel_x = _mm256_add_pd(v_plus_x, E_accel_x);
        __m256d new_vel_y = _mm256_add_pd(v_plus_y, E_accel_y);
        __m256d new_vel_z = _mm256_add_pd(v_plus_z, E_accel_z);

        // 位置更新 / Position update
        __m256d new_pos_x = _mm256_add_pd(pos_x, _mm256_mul_pd(new_vel_x, dt_vec));
        __m256d new_pos_y = _mm256_add_pd(pos_y, _mm256_mul_pd(new_vel_y, dt_vec));
        __m256d new_pos_z = _mm256_add_pd(pos_z, _mm256_mul_pd(new_vel_z, dt_vec));

        // 存储结果 / Store results
        for (size_t j = 0; j < simd_width; ++j)
        {
            if (reinterpret_cast<int64_t*>(&active_mask)[j] != 0)
            {
                auto& particle = particles[i + j];

                particle.velocity[0] = reinterpret_cast<double*>(&new_vel_x)[j];
                particle.velocity[1] = reinterpret_cast<double*>(&new_vel_y)[j];
                particle.velocity[2] = reinterpret_cast<double*>(&new_vel_z)[j];

                particle.position[0] = reinterpret_cast<double*>(&new_pos_x)[j];
                particle.position[1] = reinterpret_cast<double*>(&new_pos_y)[j];
                particle.position[2] = reinterpret_cast<double*>(&new_pos_z)[j];

                // 更新能量 / Update energy
                particle.updateEnergy();
            }
        }
    }

    // 处理剩余的粒子（标量实现）/ Process remaining particles (scalar implementation)
    for (size_t i = aligned_size; i < particles.size(); ++i)
    {
        if (particles[i].status == ParticleStatus::ACTIVE && particles[i].mass > 0.0)
        {
            Particle p(particles[i]);
            p.borisPush(E_field, B_field, dt);
            particles[i] = p.getData();
        }
    }
}
#endif

} // namespace Particle
} // namespace SCDAT
