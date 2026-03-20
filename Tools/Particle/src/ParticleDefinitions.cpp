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
             data_.type == ParticleType::MOLECULAR_ION || data_.type == ParticleType::CLUSTER_ION ||
             data_.type == ParticleType::HEAVY_ION)
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

Particle ParticleFactory::createElectron(ParticleId id, const SCDAT::Geometry::Point3D& position,
                                         const SCDAT::Geometry::Vector3D& velocity, double weight)
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
    Particle p(id, ParticleType::PHOTOELECTRON, position, velocity, PhysicsConstants::ElectronMass,
               PhysicsConstants::ElectronCharge, weight);
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
                                         const SCDAT::Geometry::Vector3D& velocity, double weight)
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

} // namespace Particle
} // namespace SCDAT
