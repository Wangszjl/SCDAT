#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Coupling
{

struct SurfaceCircuitNode
{
    std::string name;
    double potential_v = 0.0;
    double capacitance_f = 0.0;
    double shunt_conductance_s = 0.0;
    bool fixed_potential = false;
    double fixed_value_v = 0.0;
};

struct SurfaceCircuitBranch
{
    std::size_t from_node = 0;
    std::size_t to_node = 0;
    double conductance_s = 0.0;
    double bias_v = 0.0;
};

struct SurfaceCircuitLinearization
{
    std::vector<double> plasma_current_a;
    std::vector<double> plasma_didv_a_per_v;
    std::vector<double> additional_rhs_a;
    std::vector<double> additional_diagonal_a_per_v;

    struct OffDiagonalEntry
    {
        std::size_t row_node = 0;
        std::size_t column_node = 0;
        double coefficient_a_per_v = 0.0;
    };
    std::vector<OffDiagonalEntry> additional_off_diagonal_entries;
};

struct SurfaceCircuitAdvanceResult
{
    std::vector<double> node_potentials_v;
    std::vector<double> plasma_currents_a;
    std::vector<double> branch_currents_a;
    double max_delta_potential_v = 0.0;
    bool converged = false;
};

class SurfaceCircuitCoupling
{
  public:
    void clear();
    std::size_t addNode(const SurfaceCircuitNode& node);
    void addBranch(const SurfaceCircuitBranch& branch);

    void setNodePotential(std::size_t node_index, double potential_v);
    void setNodeFixedPotential(std::size_t node_index, bool fixed, double fixed_value_v);
    void setNodeShuntConductance(std::size_t node_index, double conductance_s);
    void setNodeCapacitance(std::size_t node_index, double capacitance_f);
    void setBranchConductance(std::size_t branch_index, double conductance_s);

    const std::vector<SurfaceCircuitNode>& getNodes() const
    {
        return nodes_;
    }

    const std::vector<SurfaceCircuitBranch>& getBranches() const
    {
        return branches_;
    }

    SurfaceCircuitAdvanceResult advanceImplicit(
        double dt, const SurfaceCircuitLinearization& linearization,
        double max_delta_potential_v = std::numeric_limits<double>::infinity());

  private:
    static std::vector<double> solveDenseLinearSystem(std::vector<std::vector<double>> matrix,
                                                      std::vector<double> rhs);

    std::vector<SurfaceCircuitNode> nodes_;
    std::vector<SurfaceCircuitBranch> branches_;
};

} // namespace Coupling
} // namespace SCDAT
