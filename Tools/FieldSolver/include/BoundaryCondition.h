#ifndef SCDAT_FIELD_BOUNDARY_CONDITION_H
#define SCDAT_FIELD_BOUNDARY_CONDITION_H

#include "../../Boundary/include/BoundaryType.h"
#include "../../Geometry/include/Point3D.h"
#include "../../Mesh/include/MeshAlgorithms.h"
#include <functional>
#include <memory>
#include <string>

namespace SCDAT
{
namespace FieldSolver
{

class BoundaryCondition
{
  public:
    BoundaryCondition(BoundaryConditionType type, double value = 0.0)
        : type_(type), value_(value), is_time_dependent_(false)
    {
    }

    BoundaryCondition(BoundaryConditionType type,
                      std::function<double(const Utils::Point3D&, double)> function)
        : type_(type), value_(0.0), time_function_(std::move(function)), is_time_dependent_(true)
    {
    }

    virtual ~BoundaryCondition() = default;

    BoundaryConditionType getType() const { return type_; }
    double getValue() const { return value_; }
    void setValue(double value) { value_ = value; }

    bool isTimeDependent() const { return is_time_dependent_; }
    double getValue(const Utils::Point3D& position, double time) const;

    virtual void apply(Mesh::NodePtr node, double time = 0.0) const;
    virtual bool isApplicable(Mesh::NodePtr node) const;
    std::string toString() const;

  private:
    BoundaryConditionType type_;
    double value_;
    std::function<double(const Utils::Point3D&, double)> time_function_;
    bool is_time_dependent_;
};

using BoundaryConditionPtr = std::shared_ptr<BoundaryCondition>;

} // namespace FieldSolver
} // namespace SCDAT

#endif // SCDAT_FIELD_BOUNDARY_CONDITION_H
