#include "VolField.h"

#include "../Basic/include/Constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SCDAT
{
namespace Field
{

using SCDAT::Geometry::Point3D;
using SCDAT::Geometry::Vector3D;
namespace Constants = SCDAT::Basic::Constants;

namespace
{
/**
 * @brief 在规则立方网格上进行标量场三线性插值。
 * @details
 * 该函数假设采样点坐标位于归一化区域 [0,1]^3，先将物理坐标映射到网格索引空间，
 * 再在包含该点的体素单元内进行 x/y/z 三个方向的线性混合。
 *
 * 为避免边界越界，函数会将基元索引限制在 [0, g-2]，并将权重限制在 [0,1]。
 * 这样即使输入点略超出边界，也能返回稳定结果而不会访问非法内存。
 *
 * @param values 按 z-major 展平的一维节点值数组。
 * @param p 归一化空间中的目标插值位置。
 * @param g 每个坐标方向的网格节点数（规则网格尺寸）。
 * @return 插值得到的标量值。
 */
double interpolateScalarRegular(const std::vector<double>& values, const Point3D& p, size_t g)
{
    const double spacing = 1.0 / static_cast<double>(g - 1);

    double gx = p.x() / spacing;
    double gy = p.y() / spacing;
    double gz = p.z() / spacing;

    int i = static_cast<int>(std::floor(gx));
    int j = static_cast<int>(std::floor(gy));
    int k = static_cast<int>(std::floor(gz));

    i = std::max(0, std::min(i, static_cast<int>(g) - 2));
    j = std::max(0, std::min(j, static_cast<int>(g) - 2));
    k = std::max(0, std::min(k, static_cast<int>(g) - 2));

    const double wx = std::clamp(gx - static_cast<double>(i), 0.0, 1.0);
    const double wy = std::clamp(gy - static_cast<double>(j), 0.0, 1.0);
    const double wz = std::clamp(gz - static_cast<double>(k), 0.0, 1.0);

    auto at = [&](int x, int y, int z)
    {
        const size_t idx =
            static_cast<size_t>(z) * g * g + static_cast<size_t>(y) * g + static_cast<size_t>(x);
        return (idx < values.size()) ? values[idx] : 0.0;
    };

    const double v000 = at(i, j, k);
    const double v100 = at(i + 1, j, k);
    const double v010 = at(i, j + 1, k);
    const double v110 = at(i + 1, j + 1, k);
    const double v001 = at(i, j, k + 1);
    const double v101 = at(i + 1, j, k + 1);
    const double v011 = at(i, j + 1, k + 1);
    const double v111 = at(i + 1, j + 1, k + 1);

    const double c00 = v000 * (1.0 - wx) + v100 * wx;
    const double c10 = v010 * (1.0 - wx) + v110 * wx;
    const double c01 = v001 * (1.0 - wx) + v101 * wx;
    const double c11 = v011 * (1.0 - wx) + v111 * wx;

    const double c0 = c00 * (1.0 - wy) + c10 * wy;
    const double c1 = c01 * (1.0 - wy) + c11 * wy;

    /** @note 在完成 x 与 y 方向混合后，最后沿 z 方向完成最终线性混合。 */
    return c0 * (1.0 - wz) + c1 * wz;
}

/**
 * @brief 非规则场数据的逆距离加权（IDW）插值回退实现。
 * @details
 * 当数据不能构成完美规则立方网格时，使用该函数进行近邻加权估计：
 * 1) 通过节点总数近似推断逻辑立方尺寸 g；
 * 2) 计算目标点到所有候选节点的距离并排序；
 * 3) 取最近若干节点（最多 8 个）按权重 @f$w=1/d^2@f$ 融合。
 *
 * 为避免除零，距离使用最小阈值保护；若命中几乎重合点，直接返回该点值。
 *
 * @param values 标量节点值数组。
 * @param position 目标插值位置。
 * @return 逆距离加权后的标量估计值。
 */
double interpolateScalarIrregular(const std::vector<double>& values, const Point3D& position)
{
    const size_t n = values.size();
    if (n == 0)
    {
        return 0.0;
    }

    const size_t g =
        std::max<size_t>(1, static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(n)))));
    std::vector<std::pair<double, size_t>> dist;
    dist.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        const size_t ix = i % g;
        const size_t iy = (i / g) % g;
        const size_t iz = i / (g * g);

        const Point3D x(static_cast<double>(ix) / std::max<size_t>(1, g - 1),
                        static_cast<double>(iy) / std::max<size_t>(1, g - 1),
                        static_cast<double>(iz) / std::max<size_t>(1, g - 1));

        dist.push_back({(position - x).magnitude(), i});
    }

    std::sort(dist.begin(), dist.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    if (dist.front().first < 1e-12)
    {
        return values[dist.front().second];
    }

    double ws = 0.0;
    double vs = 0.0;
    const size_t kmax = std::min<size_t>(8, dist.size());
    for (size_t k = 0; k < kmax; ++k)
    {
        const double d = std::max(1e-12, dist[k].first);
        const double w = 1.0 / (d * d);
        ws += w;
        vs += w * values[dist[k].second];
    }

    return (ws > 0.0) ? (vs / ws) : values[dist.front().second];
}
} // namespace

VolField::VolField(const VolMeshPtr& mesh, FieldType type)
    : mesh_(mesh), type_(type), unit_(Unit::DIMENSIONLESS)
{
    if (!mesh_)
    {
        throw std::invalid_argument("VolField mesh cannot be null");
    }
}

VolField::~VolField() = default;

void VolField::setUnit(Unit unit)
{
    unit_ = unit;
}

Unit VolField::getUnit() const
{
    return unit_;
}

FieldType VolField::getType() const
{
    return type_;
}

ScalarVolField::ScalarVolField(const VolMeshPtr& mesh)
    : VolField(mesh, FieldType::SCALAR), values_(mesh_->getNodeCount(), 0.0)
{
}

ScalarVolField::ScalarVolField(const VolMeshPtr& mesh, const std::vector<double>& values)
    : VolField(mesh, FieldType::SCALAR), values_(values)
{
    if (values_.size() != mesh_->getNodeCount())
    {
        throw std::invalid_argument("ScalarVolField size mismatch");
    }
}

void ScalarVolField::setValues(const std::vector<double>& values)
{
    if (values.size() != mesh_->getNodeCount())
    {
        throw std::invalid_argument("ScalarVolField size mismatch");
    }
    values_ = values;
}

double ScalarVolField::getValue(size_t index) const
{
    if (index >= values_.size())
    {
        throw std::out_of_range("ScalarVolField index out of range");
    }
    return values_[index];
}

void ScalarVolField::setValue(size_t index, double value)
{
    if (index >= values_.size())
    {
        throw std::out_of_range("ScalarVolField index out of range");
    }
    values_[index] = value;
}

double ScalarVolField::interpolate(const Point3D& position) const
{
    if (values_.empty())
    {
        return 0.0;
    }

    const size_t n = values_.size();
    const size_t g = static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(n))));
    if (g > 1 && g * g * g == n)
    {
        return interpolateRegularGrid(position, g);
    }
    return interpolateIrregularGrid(position);
}

void ScalarVolField::reset()
{
    std::fill(values_.begin(), values_.end(), 0.0);
}

std::string ScalarVolField::getStatistics() const
{
    std::ostringstream oss;
    if (values_.empty())
    {
        oss << "ScalarField empty";
        return oss.str();
    }

    const double sum = computeSum();
    oss << "ScalarField n=" << values_.size() << ", min=" << computeMin()
        << ", max=" << computeMax() << ", avg=" << (sum / static_cast<double>(values_.size()));
    return oss.str();
}

VectorVolFieldPtr ScalarVolField::computeGradient() const
{
    auto out = std::make_shared<VectorVolField>(mesh_);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        out->setValue(i, mesh_->computeGradient(i, values_));
    }
    return out;
}

ScalarVolFieldPtr ScalarVolField::computeLaplacian() const
{
    auto grad = computeGradient();
    auto out = std::make_shared<ScalarVolField>(mesh_);
    const auto& gv = grad->getValues();
    for (size_t i = 0; i < values_.size(); ++i)
    {
        out->setValue(i, mesh_->computeDivergence(i, gv));
    }
    return out;
}

void ScalarVolField::add(const ScalarVolField& other)
{
    if (other.values_.size() != values_.size())
    {
        throw std::invalid_argument("ScalarVolField add size mismatch");
    }
    for (size_t i = 0; i < values_.size(); ++i)
    {
        values_[i] += other.values_[i];
    }
}

void ScalarVolField::multiply(double factor)
{
    for (double& v : values_)
    {
        v *= factor;
    }
}

void ScalarVolField::applyFunction(std::function<double(double)> func)
{
    for (double& v : values_)
    {
        v = func(v);
    }
}

double ScalarVolField::computeSum() const
{
    double s = 0.0;
    for (double v : values_)
    {
        s += v;
    }
    return s;
}

double ScalarVolField::computeMax() const
{
    if (values_.empty())
    {
        return 0.0;
    }
    return *std::max_element(values_.begin(), values_.end());
}

double ScalarVolField::computeMin() const
{
    if (values_.empty())
    {
        return 0.0;
    }
    return *std::min_element(values_.begin(), values_.end());
}

double ScalarVolField::interpolateRegularGrid(const Point3D& p, size_t g) const
{
    const double spacing = 1.0 / static_cast<double>(g - 1);

    double gx = p.x() / spacing;
    double gy = p.y() / spacing;
    double gz = p.z() / spacing;

    int i = static_cast<int>(std::floor(gx));
    int j = static_cast<int>(std::floor(gy));
    int k = static_cast<int>(std::floor(gz));

    i = std::max(0, std::min(i, static_cast<int>(g) - 2));
    j = std::max(0, std::min(j, static_cast<int>(g) - 2));
    k = std::max(0, std::min(k, static_cast<int>(g) - 2));

    const double wx = std::clamp(gx - static_cast<double>(i), 0.0, 1.0);
    const double wy = std::clamp(gy - static_cast<double>(j), 0.0, 1.0);
    const double wz = std::clamp(gz - static_cast<double>(k), 0.0, 1.0);

    auto at = [&](int x, int y, int z)
    {
        const size_t idx =
            static_cast<size_t>(z) * g * g + static_cast<size_t>(y) * g + static_cast<size_t>(x);
        return (idx < values_.size()) ? values_[idx] : 0.0;
    };

    const double v000 = at(i, j, k);
    const double v100 = at(i + 1, j, k);
    const double v010 = at(i, j + 1, k);
    const double v110 = at(i + 1, j + 1, k);
    const double v001 = at(i, j, k + 1);
    const double v101 = at(i + 1, j, k + 1);
    const double v011 = at(i, j + 1, k + 1);
    const double v111 = at(i + 1, j + 1, k + 1);

    const double c00 = v000 * (1.0 - wx) + v100 * wx;
    const double c10 = v010 * (1.0 - wx) + v110 * wx;
    const double c01 = v001 * (1.0 - wx) + v101 * wx;
    const double c11 = v011 * (1.0 - wx) + v111 * wx;

    const double c0 = c00 * (1.0 - wy) + c10 * wy;
    const double c1 = c01 * (1.0 - wy) + c11 * wy;

    return c0 * (1.0 - wz) + c1 * wz;
}

double ScalarVolField::interpolateIrregularGrid(const Point3D& position) const
{
    const size_t n = values_.size();
    if (n == 0)
    {
        return 0.0;
    }

    const size_t g =
        std::max<size_t>(1, static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(n)))));
    std::vector<std::pair<double, size_t>> dist;
    dist.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        const size_t ix = i % g;
        const size_t iy = (i / g) % g;
        const size_t iz = i / (g * g);

        const Point3D x(static_cast<double>(ix) / std::max<size_t>(1, g - 1),
                        static_cast<double>(iy) / std::max<size_t>(1, g - 1),
                        static_cast<double>(iz) / std::max<size_t>(1, g - 1));

        const double d = (position - x).magnitude();
        dist.push_back({d, i});
    }

    std::sort(dist.begin(), dist.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    if (dist.front().first < 1e-12)
    {
        return values_[dist.front().second];
    }

    double ws = 0.0;
    double vs = 0.0;
    const size_t kmax = std::min<size_t>(8, dist.size());
    for (size_t k = 0; k < kmax; ++k)
    {
        const double d = std::max(1e-12, dist[k].first);
        const double w = 1.0 / (d * d);
        ws += w;
        vs += w * values_[dist[k].second];
    }

    return (ws > 0.0) ? (vs / ws) : values_[dist.front().second];
}

VectorVolField::VectorVolField(const VolMeshPtr& mesh)
    : VolField(mesh, FieldType::VECTOR), values_(mesh_->getNodeCount(), Vector3D(0.0, 0.0, 0.0))
{
}

VectorVolField::VectorVolField(const VolMeshPtr& mesh, const std::vector<Vector3D>& values)
    : VolField(mesh, FieldType::VECTOR), values_(values)
{
    if (values_.size() != mesh_->getNodeCount())
    {
        throw std::invalid_argument("VectorVolField size mismatch");
    }
}

void VectorVolField::setValues(const std::vector<Vector3D>& values)
{
    if (values.size() != mesh_->getNodeCount())
    {
        throw std::invalid_argument("VectorVolField size mismatch");
    }
    values_ = values;
}

Vector3D VectorVolField::getValue(size_t index) const
{
    if (index >= values_.size())
    {
        throw std::out_of_range("VectorVolField index out of range");
    }
    return values_[index];
}

void VectorVolField::setValue(size_t index, const Vector3D& value)
{
    if (index >= values_.size())
    {
        throw std::out_of_range("VectorVolField index out of range");
    }
    values_[index] = value;
}

Vector3D VectorVolField::interpolate(const Point3D& position) const
{
    if (values_.empty())
    {
        return Vector3D(0.0, 0.0, 0.0);
    }

    const size_t n = values_.size();
    const size_t g = static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(n))));
    if (g > 1 && g * g * g == n)
    {
        return interpolateRegularGrid(position, g);
    }
    return interpolateIrregularGrid(position);
}

void VectorVolField::reset()
{
    std::fill(values_.begin(), values_.end(), Vector3D(0.0, 0.0, 0.0));
}

std::string VectorVolField::getStatistics() const
{
    std::ostringstream oss;
    oss << "VectorField n=" << values_.size() << ", |v|max=" << computeMagnitudeMax();
    return oss.str();
}

ScalarVolFieldPtr VectorVolField::computeDivergence() const
{
    auto out = std::make_shared<ScalarVolField>(mesh_);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        out->setValue(i, mesh_->computeDivergence(i, values_));
    }
    return out;
}

VectorVolFieldPtr VectorVolField::computeCurl() const
{
    auto out = std::make_shared<VectorVolField>(mesh_);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        out->setValue(i, mesh_->computeCurl(i, values_));
    }
    return out;
}

ScalarVolFieldPtr VectorVolField::computeMagnitude() const
{
    auto out = std::make_shared<ScalarVolField>(mesh_);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        out->setValue(i, values_[i].magnitude());
    }
    return out;
}

void VectorVolField::add(const VectorVolField& other)
{
    if (other.values_.size() != values_.size())
    {
        throw std::invalid_argument("VectorVolField add size mismatch");
    }
    for (size_t i = 0; i < values_.size(); ++i)
    {
        values_[i] += other.values_[i];
    }
}

void VectorVolField::multiply(double factor)
{
    for (auto& v : values_)
    {
        v *= factor;
    }
}

double VectorVolField::computeMagnitudeMax() const
{
    double m = 0.0;
    for (const auto& v : values_)
    {
        m = std::max(m, v.magnitude());
    }
    return m;
}

Vector3D VectorVolField::interpolateRegularGrid(const Point3D& p, size_t g) const
{
    std::vector<double> sx(values_.size(), 0.0);
    std::vector<double> sy(values_.size(), 0.0);
    std::vector<double> sz(values_.size(), 0.0);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        sx[i] = values_[i].x();
        sy[i] = values_[i].y();
        sz[i] = values_[i].z();
    }
    return Vector3D(interpolateScalarRegular(sx, p, g), interpolateScalarRegular(sy, p, g),
                    interpolateScalarRegular(sz, p, g));
}

Vector3D VectorVolField::interpolateIrregularGrid(const Point3D& p) const
{
    std::vector<double> sx(values_.size(), 0.0);
    std::vector<double> sy(values_.size(), 0.0);
    std::vector<double> sz(values_.size(), 0.0);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        sx[i] = values_[i].x();
        sy[i] = values_[i].y();
        sz[i] = values_[i].z();
    }
    return Vector3D(interpolateScalarIrregular(sx, p), interpolateScalarIrregular(sy, p),
                    interpolateScalarIrregular(sz, p));
}

ElectricField::ElectricField(const VolMeshPtr& mesh) : VectorVolField(mesh)
{
    setUnit(Unit::VOLT_PER_METER);
}

void ElectricField::computeFromPotential(const ScalarVolField& potential)
{
    if (potential.getMesh().get() != mesh_.get())
    {
        throw std::invalid_argument("ElectricField potential mesh mismatch");
    }
    /**
     * @note
     * 根据静电关系式 @f$\mathbf{E}=-\nabla\phi@f$：
     * 先计算电势梯度，再整体乘以 -1 得到电场向量。
     */
    const auto grad = potential.computeGradient();
    values_ = grad->getValues();
    for (auto& e : values_)
    {
        e *= -1.0;
    }
}

void ElectricField::addUniformField(const Vector3D& uniform_field)
{
    for (auto& e : values_)
    {
        e += uniform_field;
    }
}

ScalarVolFieldPtr ElectricField::computeEnergyDensity() const
{
    auto out = std::make_shared<ScalarVolField>(mesh_);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        const double e2 = values_[i].magnitudeSquared();
        /**
         * @brief 计算电场能量密度。
         * @details 使用公式 @f$u_E=\frac{1}{2}\varepsilon_0|\mathbf{E}|^2@f$。
         */
        out->setValue(i, 0.5 * Constants::PhysicsConstants::VacuumPermittivity * e2);
    }
    out->setUnit(Unit::JOULE_PER_CUBIC_METER);
    return out;
}

MagneticField::MagneticField(const VolMeshPtr& mesh) : VectorVolField(mesh)
{
    setUnit(Unit::TESLA);
}

void MagneticField::setUniformField(const Vector3D& uniform_field)
{
    for (auto& b : values_)
    {
        b = uniform_field;
    }
}

ScalarVolFieldPtr MagneticField::computeEnergyDensity() const
{
    auto out = std::make_shared<ScalarVolField>(mesh_);
    for (size_t i = 0; i < values_.size(); ++i)
    {
        const double b2 = values_[i].magnitudeSquared();
        /**
         * @brief 计算磁场能量密度。
         * @details 使用公式 @f$u_B=\frac{|\mathbf{B}|^2}{2\mu_0}@f$。
         */
        out->setValue(i, b2 / (2.0 * Constants::PhysicsConstants::VacuumPermeability));
    }
    out->setUnit(Unit::JOULE_PER_CUBIC_METER);
    return out;
}

VectorVolFieldPtr MagneticField::computeMagneticInduction() const
{
    auto out = std::make_shared<VectorVolField>(mesh_);
    out->setValues(values_);
    out->setUnit(Unit::TESLA);
    return out;
}

void MagneticField::addDipoleField(const Point3D& dipolePosition, const Vector3D& dipoleMoment)
{
    const size_t n = values_.size();
    if (n == 0)
    {
        return;
    }

    const size_t g =
        std::max<size_t>(1, static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(n)))));
    for (size_t i = 0; i < n; ++i)
    {
        const size_t ix = i % g;
        const size_t iy = (i / g) % g;
        const size_t iz = i / (g * g);

        Point3D x(static_cast<double>(ix) / std::max<size_t>(1, g - 1),
                  static_cast<double>(iy) / std::max<size_t>(1, g - 1),
                  static_cast<double>(iz) / std::max<size_t>(1, g - 1));

        Vector3D r = x - dipolePosition;
        const double rmag = r.magnitude();
        if (rmag < 1e-9)
        {
            continue;
        }

        /** @brief 从偶极子中心指向采样点的单位方向向量。 */
        const Vector3D rh = r / rmag;
        /**
         * @brief 磁偶极场系数。
         * @details 系数形式为 @f$\mu_0/(4\pi r^3)@f$。
         */
        const double coeff =
            Constants::PhysicsConstants::VacuumPermeability /
            (4.0 * Constants::MathConstants::Pi * rmag * rmag * rmag);
        /**
         * @brief 计算偶极子在采样点处的磁场增量。
         * @details 使用表达式
         * @f$\mathbf{B}(\mathbf{r})=\frac{\mu_0}{4\pi r^3}[3(\hat{\mathbf{r}}\cdot\mathbf{m})\hat{\mathbf{r}}-\mathbf{m}]@f$。
         */
        Vector3D b = (3.0 * rh * dipoleMoment.dot(rh) - dipoleMoment) * coeff;
        values_[i] += b;
    }
}

TemperatureField::TemperatureField(const VolMeshPtr& mesh) : ScalarVolField(mesh)
{
    setUnit(Unit::KELVIN);
}

void TemperatureField::setUniformTemperature(double temperature)
{
    for (size_t i = 0; i < values_.size(); ++i)
    {
        values_[i] = temperature;
    }
}

VectorVolFieldPtr TemperatureField::computeHeatFlux(double thermal_conductivity) const
{
    auto grad = computeGradient();
    auto out = std::make_shared<VectorVolField>(mesh_);
    out->setValues(grad->getValues());
    out->multiply(-thermal_conductivity);
    return out;
}

void TemperatureField::addGradient(const Vector3D& gradient, const Point3D& reference)
{
    const size_t n = values_.size();
    const size_t g =
        std::max<size_t>(1, static_cast<size_t>(std::llround(std::cbrt(static_cast<double>(n)))));

    for (size_t i = 0; i < n; ++i)
    {
        const size_t ix = i % g;
        const size_t iy = (i / g) % g;
        const size_t iz = i / (g * g);
        Point3D x(static_cast<double>(ix) / std::max<size_t>(1, g - 1),
                  static_cast<double>(iy) / std::max<size_t>(1, g - 1),
                  static_cast<double>(iz) / std::max<size_t>(1, g - 1));
        values_[i] += gradient.dot(x - reference);
    }
}

void TemperatureField::applyHeatSource(const ScalarVolField& heatSource, double thermalDiffusivity,
                                       double dt)
{
    if (heatSource.getValues().size() != values_.size())
    {
        throw std::invalid_argument("TemperatureField heat source size mismatch");
    }

    for (size_t i = 0; i < values_.size(); ++i)
    {
        values_[i] += dt * (thermalDiffusivity * computeLaplacian(i) + heatSource.getValue(i));
    }
}

double TemperatureField::computeLaplacian(size_t nodeIndex) const
{
    auto lap = ScalarVolField::computeLaplacian();
    return lap->getValue(nodeIndex);
}

CurrentDensityField::CurrentDensityField(const VolMeshPtr& mesh) : VectorVolField(mesh)
{
    setUnit(Unit::AMPERE_PER_SQUARE_METER);
}

void CurrentDensityField::computeFromElectricField(const ElectricField& electric_field,
                                                   const ScalarVolField& conductivity_field)
{
    if (electric_field.getValues().size() != values_.size() ||
        conductivity_field.getValues().size() != values_.size())
    {
        throw std::invalid_argument("CurrentDensityField input size mismatch");
    }

    for (size_t i = 0; i < values_.size(); ++i)
    {
        values_[i] = electric_field.getValue(i) * conductivity_field.getValue(i);
    }
}

void CurrentDensityField::addConvectionCurrent(const ScalarVolField& charge_density,
                                               const VectorVolField& velocity_field)
{
    if (charge_density.getValues().size() != values_.size() ||
        velocity_field.getValues().size() != values_.size())
    {
        throw std::invalid_argument("CurrentDensityField convection input size mismatch");
    }

    for (size_t i = 0; i < values_.size(); ++i)
    {
        values_[i] += velocity_field.getValue(i) * charge_density.getValue(i);
    }
}

double CurrentDensityField::computeTotalCurrent() const
{
    double sum = 0.0;
    for (const auto& j : values_)
    {
        sum += j.magnitude();
    }
    return sum;
}

double CurrentDensityField::checkCurrentContinuity() const
{
    double maxDiv = 0.0;
    for (size_t i = 0; i < values_.size(); ++i)
    {
        maxDiv = std::max(maxDiv, std::abs(mesh_->computeDivergence(i, values_)));
    }
    return maxDiv;
}

} // namespace Field
} // namespace SCDAT
