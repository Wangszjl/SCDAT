#ifndef SCDAT_PICCORE_PIC_CYCLE_H
#define SCDAT_PICCORE_PIC_CYCLE_H

#include "../../Boundary/include/BoundaryConditions.h"
#include "../../FieldSolver/include/PoissonSolver.h"
#include "../../Particle/include/PICParticleAlgorithms.h"
#include <array>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace PICcore
{

using ::SCDAT::FieldSolver::PoissonSolver;
using ::SCDAT::Particle::BoundaryConditionPtr;
using ::SCDAT::Particle::ChargeDepositor;
using ::SCDAT::Particle::ChargeDepositorPtr;
using ::SCDAT::Particle::CrossTetrahedron;
using ::SCDAT::Particle::CrossTetrahedronPtr;
using ::SCDAT::Particle::ElementId;
using ::SCDAT::Particle::FieldInterpolator;
using ::SCDAT::Particle::FieldInterpolatorPtr;
using ::SCDAT::Particle::ParticleLocation;
using ::SCDAT::Particle::ParticleManager;
using ::SCDAT::Particle::ParticlePtr;
using ::SCDAT::Particle::ParticleTypeDef;

class PICCycle
{
  public:
    enum class BoundaryFace
    {
        XMin = 0,
        XMax = 1,
        YMin = 2,
        YMax = 3,
        ZMin = 4,
        ZMax = 5,
    };

    struct Parameters
    {
        double time_step = 1e-12;
        int max_iterations = 1000;
        double convergence_tolerance = 1e-6;
        bool adaptive_time_step = false;
        bool charge_conservation = true;
        int statistics_frequency = 10;
        bool verbose = false;
    };

    PICCycle();
    explicit PICCycle(const Parameters& params);
    ~PICCycle() = default;

    void setMesh(const std::vector<Mesh::NodePtr>& nodes,
                 const std::vector<Mesh::ElementPtr>& elements);
    void setParticleManager(std::shared_ptr<ParticleManager> particle_manager);
    void setPoissonSolver(std::shared_ptr<PoissonSolver> poisson_solver);
    void setBoundaryCondition(BoundaryFace face, BoundaryConditionPtr condition);
    void clearBoundaryConditions();

    void addParticleLocation(ParticlePtr particle, ElementId element_id);
    void updateParticleLocation(ParticlePtr particle, ElementId new_element_id);
    ElementId findParticleElement(ParticlePtr particle) const;

    bool initialize();
    bool executeTimeStep();
    bool executeSimulation(double total_time);

    const std::vector<double>& getChargeDensity() const
    {
        return charge_density_;
    }
    const std::vector<double>& getPotential() const
    {
        return potential_;
    }
    const std::vector<Utils::Vector3D>& getElectricField() const
    {
        return electric_field_;
    }

    void setTimeStep(double dt)
    {
        if (dt > 0.0)
        {
            params_.time_step = dt;
        }
    }

    double getTimeStep() const
    {
        return params_.time_step;
    }

    void setPotential(const std::vector<double>& potential)
    {
        potential_ = potential;
        field_solver_solution_valid_ = false;
        if (nodes_.size() == potential_.size())
        {
            for (std::size_t i = 0; i < nodes_.size(); ++i)
            {
                if (nodes_[i])
                {
                    nodes_[i]->setPotential(potential_[i]);
                }
            }
        }
        updateElectricField();
    }

    void setNodePotential(size_t node_id, double value)
    {
        if (node_id < potential_.size())
        {
            potential_[node_id] = value;
            field_solver_solution_valid_ = false;
            if (node_id < nodes_.size() && nodes_[node_id])
            {
                nodes_[node_id]->setPotential(value);
            }
        }
    }

    struct CycleStatistics
    {
        int completed_steps = 0;
        double simulation_time = 0.0;
        double total_charge = 0.0;
        double total_energy = 0.0;
        int active_particles = 0;
        double average_step_time = 0.0;
        double surface_deposited_charge = 0.0;
    };

    struct SurfaceBoundaryMetadata
    {
        int surface_id = -1;
        int material_id = -1;
        int circuit_node_id = -1;
    };

    struct SurfaceCurrentLedgerEntry
    {
        SurfaceBoundaryMetadata metadata;
        double absorbed_electron_charge_c = 0.0;
        double absorbed_ion_charge_c = 0.0;
        double emitted_electron_charge_c = 0.0;
        double incident_energy_j = 0.0;
        std::size_t absorbed_electron_particles = 0;
        std::size_t absorbed_ion_particles = 0;
        std::size_t emitted_electron_particles = 0;
    };

    const CycleStatistics& getStatistics() const
    {
        return cycle_stats_;
    }

    void setBoundaryMetadata(BoundaryFace face, const SurfaceBoundaryMetadata& metadata)
    {
        surface_boundary_metadata_[static_cast<std::size_t>(face)] = metadata;
    }

    const std::array<SurfaceCurrentLedgerEntry, 6>& getSurfaceCurrentLedger() const
    {
        return surface_current_ledger_;
    }

    double getAndResetSurfaceCharge()
    {
        const double charge = cycle_stats_.surface_deposited_charge;
        cycle_stats_.surface_deposited_charge = 0.0;
        return charge;
    }

  private:
    struct DomainBounds
    {
        Utils::Point3D min_corner{};
        Utils::Point3D max_corner{};
        bool valid = false;
    };

    struct RuntimeBoundaryEvent
    {
        std::size_t particle_index = 0;
        ElementId hint_element = static_cast<ElementId>(-1);
    };

    Parameters params_;

    std::vector<Mesh::NodePtr> nodes_;
    std::vector<Mesh::ElementPtr> elements_;
    std::shared_ptr<ParticleManager> particle_manager_;
    std::shared_ptr<PoissonSolver> field_poisson_solver_;
    CrossTetrahedronPtr cross_tetrahedron_;
    FieldInterpolatorPtr field_interpolator_;
    ChargeDepositorPtr charge_depositor_;

    std::vector<double> charge_density_;
    std::vector<double> potential_;
    std::vector<Utils::Vector3D> electric_field_;
    std::vector<Utils::Vector3D> magnetic_field_;
    std::array<BoundaryConditionPtr, 6> boundary_conditions_{};
    std::array<SurfaceBoundaryMetadata, 6> surface_boundary_metadata_{};
    std::array<SurfaceCurrentLedgerEntry, 6> surface_current_ledger_{};
    DomainBounds domain_bounds_{};
    bool field_solver_solution_valid_ = false;

    CycleStatistics cycle_stats_;
    std::unordered_map<Particle::ParticleId, ParticleLocation> particle_locations_;

    void pushParticles();
    void computeChargeDensity();
    void solvePoissonEquation();
    void updateElectricField();
    void updateStatistics();
    bool checkConvergence();
    void outputStatistics();
    Utils::Vector3D computeGradient(size_t node_id);
    void updateDomainBounds();
    bool hasRuntimeBoundaryConditions() const;
    std::optional<BoundaryFace> detectBoundaryFace(const Utils::Point3D& position,
                                                   double outside_tolerance,
                                                   bool allow_near_match = false,
                                                   double near_tolerance = 0.0) const;
    Utils::Point3D projectToBoundaryFace(BoundaryFace face, const Utils::Point3D& position,
                                         double inward_offset = 0.0) const;
    void processRuntimeBoundaryEvents(const std::vector<RuntimeBoundaryEvent>& events,
                                      std::vector<ParticleTypeDef>& emitted_particles);
    bool applyRuntimeBoundaryCondition(ParticlePtr particle, ElementId hint_element,
                                       std::vector<ParticleTypeDef>& emitted_particles,
                                       ElementId& resolved_element);
    void recordSurfaceBoundaryInteraction(BoundaryFace face, const ParticleTypeDef& particle,
                                          const std::vector<ParticleTypeDef>& emitted_particles);
    ElementId locateParticleElement(const Utils::Point3D& position,
                                    ElementId hint_element = static_cast<ElementId>(-1)) const;
};

} // namespace PICcore
} // namespace SCDAT

#endif // SCDAT_PICCORE_PIC_CYCLE_H
