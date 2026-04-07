#include "../include/FieldEmissionBoundaryCondition.h"

#include "../../Basic/include/Constants.h"

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Particle
{

namespace
{
constexpr double kNormalTolerance = 1.0e-12;
}

FieldEmissionBoundaryCondition::FieldEmissionBoundaryCondition()
    : FieldEmissionBoundaryCondition(Parameters{})
{
}

FieldEmissionBoundaryCondition::FieldEmissionBoundaryCondition(const Parameters& parameters)
    : BoundaryConditionBase(AdvancedBoundaryType::EMISSION_FIELD, BoundaryGeometry::PLANAR),
      parameters_(parameters)
{
    if (parameters_.emission_area_m2 <= 0.0)
    {
        parameters_.emission_area_m2 = 1.0e-8;
    }
    if (parameters_.emitted_electron_speed_m_per_s <= 0.0)
    {
        parameters_.emitted_electron_speed_m_per_s = 2.0e6;
    }
    if (parameters_.emitted_electron_weight <= 0.0)
    {
        parameters_.emitted_electron_weight = 1.0;
    }
    if (parameters_.max_emitted_particles_per_event == 0)
    {
        parameters_.max_emitted_particles_per_event = 1;
    }
    if (parameters_.secondary_electron_speed_m_per_s <= 0.0)
    {
        parameters_.secondary_electron_speed_m_per_s = 8.0e5;
    }
    if (parameters_.secondary_electron_weight <= 0.0)
    {
        parameters_.secondary_electron_weight = 1.0;
    }
    if (parameters_.max_secondary_particles_per_event == 0)
    {
        parameters_.max_secondary_particles_per_event = 1;
    }
}

void FieldEmissionBoundaryCondition::setEmissionCurrentDensity(
    double emission_current_density_a_per_m2)
{
    parameters_.emission_current_density_a_per_m2 =
        std::max(0.0, emission_current_density_a_per_m2);
}

void FieldEmissionBoundaryCondition::setEmissionArea(double emission_area_m2)
{
    parameters_.emission_area_m2 = std::max(1.0e-12, emission_area_m2);
}

void FieldEmissionBoundaryCondition::resetEmissionAccumulator()
{
    pending_emission_particles_ = 0.0;
    pending_secondary_particles_ = 0.0;
}

bool FieldEmissionBoundaryCondition::isIonDrivenEmissionSource(const ParticleTypeDef& particle)
{
    const auto type = particle.getType();
    return (particle.getCharge() > 0.0) ||
           (type == ParticleType::ION || type == ParticleType::POSITIVE_ION ||
            type == ParticleType::MOLECULAR_ION || type == ParticleType::CLUSTER_ION ||
            type == ParticleType::PROTON || type == ParticleType::ALPHA ||
            type == ParticleType::HEAVY_ION);
}

std::vector<ParticleTypeDef> FieldEmissionBoundaryCondition::processParticle(
    ParticleTypeDef& particle, const SCDAT::Geometry::Point3D& intersection_point,
    const SCDAT::Geometry::Vector3D& normal, double dt)
{
    std::vector<ParticleTypeDef> emitted_particles;
    if (!enabled_)
    {
        return emitted_particles;
    }

    particle.markAsAbsorbed();
    updateStatistics(particle, "absorbed");

    if (!isIonDrivenEmissionSource(particle))
    {
        return emitted_particles;
    }

    const auto normal_magnitude = normal.magnitude();
    const SCDAT::Geometry::Vector3D unit_normal =
        (normal_magnitude > kNormalTolerance)
            ? SCDAT::Geometry::Vector3D(normal / normal_magnitude)
            : SCDAT::Geometry::Vector3D(0.0, 0.0, 1.0);

    const auto emit_electrons = [&](std::size_t count, ParticleType type, double speed,
                                    double weight) {
        if (count == 0)
        {
            return;
        }

        const SCDAT::Geometry::Vector3D emitted_velocity = unit_normal * speed;
        const double emitted_mass = Basic::Constants::PhysicsConstants::ElectronMass;
        const double emitted_charge = Basic::Constants::PhysicsConstants::ElectronCharge;

        for (std::size_t i = 0; i < count; ++i)
        {
            const double offset = 1.0e-10 * static_cast<double>(emitted_particles.size() + 1);
            const SCDAT::Geometry::Point3D emitted_position =
                intersection_point + unit_normal * offset;
            ParticleTypeDef emitted_particle(next_particle_id_++, type, emitted_position,
                                             emitted_velocity, emitted_mass, emitted_charge,
                                             weight);
            emitted_particles.push_back(emitted_particle);
            updateStatistics(emitted_particles.back(), "emitted");
        }
    };

    if (dt > 0.0 && parameters_.emission_current_density_a_per_m2 > 0.0)
    {
        const double elementary_charge = Basic::Constants::PhysicsConstants::ElementaryCharge;
        const double expected_emitted =
            parameters_.emission_current_density_a_per_m2 * parameters_.emission_area_m2 * dt /
            elementary_charge;
        if (std::isfinite(expected_emitted) && expected_emitted > 0.0)
        {
            pending_emission_particles_ += expected_emitted;
            std::size_t emission_count =
                static_cast<std::size_t>(std::floor(pending_emission_particles_));
            emission_count = std::min(emission_count, parameters_.max_emitted_particles_per_event);
            pending_emission_particles_ -= static_cast<double>(emission_count);
            pending_emission_particles_ = std::max(0.0, pending_emission_particles_);

            emit_electrons(emission_count, ParticleType::FIELD_EMISSION_ELECTRON,
                           parameters_.emitted_electron_speed_m_per_s,
                           parameters_.emitted_electron_weight);
        }
    }

    if (parameters_.enable_secondary_electron_emission &&
        parameters_.secondary_electron_yield_per_ion > 0.0)
    {
        const double ion_weight = std::max(0.0, particle.getWeight());
        const double normalization_weight =
            std::max(1.0e-12, parameters_.secondary_electron_weight);
        const double expected_secondary =
            parameters_.secondary_electron_yield_per_ion * ion_weight / normalization_weight;
        if (std::isfinite(expected_secondary) && expected_secondary > 0.0)
        {
            pending_secondary_particles_ += expected_secondary;
            std::size_t secondary_count =
                static_cast<std::size_t>(std::floor(pending_secondary_particles_));
            secondary_count =
                std::min(secondary_count, parameters_.max_secondary_particles_per_event);
            pending_secondary_particles_ -= static_cast<double>(secondary_count);
            pending_secondary_particles_ = std::max(0.0, pending_secondary_particles_);

            emit_electrons(secondary_count, ParticleType::SECONDARY_ELECTRON,
                           parameters_.secondary_electron_speed_m_per_s,
                           parameters_.secondary_electron_weight);
        }
    }

    return emitted_particles;
}

} // namespace Particle
} // namespace SCDAT
