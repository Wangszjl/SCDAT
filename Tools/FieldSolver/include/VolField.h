/**
 * @file VolField.h
 * @brief 体积场计算模块头文件
 * @details 定义体积内电场分布计算、磁场分布处理、温度场耦合等功能
 *
 * @author Wang Sizhan
 * @date 2026.3.36 16:26:01
 * @version V0.0.1
 * @ingroup FieldSolverModule
 */

#ifndef SCDAT_FIELD_VOL_FIELD_H
#define SCDAT_FIELD_VOL_FIELD_H

#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
// Legacy mesh include removed during SCDAT namespace migration.
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Field
{

// 前向声明
class VolMesh;
class VolField;
class ScalarVolField;
class VectorVolField;
class ElectricField;
class MagneticField;
class TemperatureField;

// 类型别名
using VolMeshPtr = std::shared_ptr<VolMesh>;
using VolFieldPtr = std::shared_ptr<VolField>;
using ScalarVolFieldPtr = std::shared_ptr<ScalarVolField>;
using VectorVolFieldPtr = std::shared_ptr<VectorVolField>;
using ElectricFieldPtr = std::shared_ptr<ElectricField>;
using MagneticFieldPtr = std::shared_ptr<MagneticField>;

/**
 * @brief 场类型枚举
 */
enum class FieldType
{
    SCALAR = 0,
    VECTOR = 1,
    TENSOR = 2
};

/**
 * @brief 物理单位枚举
 */
enum class Unit
{
    DIMENSIONLESS = 0,
    VOLT = 1,
    VOLT_PER_METER = 2,
    TESLA = 3,
    KELVIN = 4,
    PASCAL = 5,
    JOULE_PER_CUBIC_METER = 6,
    AMPERE_PER_SQUARE_METER = 7
};

/**
 * @brief 体积网格抽象基类
 * @details 定义三维体积网格的基本接口，支持有限元和有限差分方法
 *          提供网格拓扑信息和微分算子计算功能
 *
 * 该类为纯虚基类，具体实现包括：
 * - 结构化网格（规则六面体网格）
 * - 非结构化网格（四面体、六面体混合网格）
 * - 自适应网格（支持网格细化和粗化）
 *
 * @note 所有派生类必须实现所有纯虚函数
 * @warning 网格索引从0开始，使用前需要检查索引有效性
 */
class VolMesh
{
  public:
    /**
     * @brief 默认构造函数
     * @details 创建空的网格对象，需要后续初始化
     */
    VolMesh() = default;

    /**
     * @brief 虚析构函数
     * @details 确保派生类对象能够正确析构
     */
    virtual ~VolMesh() = default;

    /**
     * @brief 获取网格节点总数
     * @return 节点数量
     * @details 返回网格中所有节点的数量，包括边界节点和内部节点
     * @note 节点索引范围为 [0, getNodeCount()-1]
     */
    virtual size_t getNodeCount() const = 0;

    /**
     * @brief 获取网格单元总数
     * @return 单元数量
     * @details 返回网格中所有单元的数量，单元类型可能包括四面体、六面体等
     * @note 单元索引范围为 [0, getElementCount()-1]
     */
    virtual size_t getElementCount() const = 0;

    /**
     * @brief 计算标量场在指定节点的梯度
     * @param node_index 节点索引，必须在有效范围内
     * @param values 标量场值数组，大小必须等于节点数量
     * @return 梯度向量 ∇φ
     * @details 使用有限元或有限差分方法计算梯度
     *
     * 数学表达式：∇φ = (∂φ/∂x, ∂φ/∂y, ∂φ/∂z)
     *
     * @throws std::out_of_range 当节点索引超出范围时
     * @throws std::invalid_argument 当values数组大小不匹配时
     */
    virtual SCDAT::Geometry::Vector3D computeGradient(size_t node_index,
                                                      const std::vector<double>& values) const = 0;

    /**
     * @brief 计算向量场在指定节点的散度
     * @param node_index 节点索引，必须在有效范围内
     * @param vector_field 向量场数组，大小必须等于节点数量
     * @return 散度值 ∇·F
     * @details 使用有限元或有限差分方法计算散度
     *
     * 数学表达式：∇·F = ∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z
     *
     * @throws std::out_of_range 当节点索引超出范围时
     * @throws std::invalid_argument 当vector_field数组大小不匹配时
     */
    virtual double
    computeDivergence(size_t node_index,
                      const std::vector<SCDAT::Geometry::Vector3D>& vector_field) const = 0;

    /**
     * @brief 计算向量场在指定节点的旋度
     * @param node_index 节点索引，必须在有效范围内
     * @param vector_field 向量场数组，大小必须等于节点数量
     * @return 旋度向量 ∇×F
     * @details 使用有限元或有限差分方法计算旋度
     *
     * 数学表达式：∇×F = (∂Fz/∂y - ∂Fy/∂z, ∂Fx/∂z - ∂Fz/∂x, ∂Fy/∂x - ∂Fx/∂y)
     *
     * @throws std::out_of_range 当节点索引超出范围时
     * @throws std::invalid_argument 当vector_field数组大小不匹配时
     */
    virtual SCDAT::Geometry::Vector3D
    computeCurl(size_t node_index,
                const std::vector<SCDAT::Geometry::Vector3D>& vector_field) const = 0;
};

/**
 * @brief 体积场基类
 */
class VolField
{
  public:
    /**
     * @brief 构造函数
     * @param mesh 网格指针
     * @param type 场类型
     */
    VolField(const VolMeshPtr& mesh, FieldType type);

    /**
     * @brief 虚析构函数
     */
    virtual ~VolField();

    /**
     * @brief 获取网格
     */
    VolMeshPtr getMesh() const
    {
        return mesh_;
    }

    /**
     * @brief 获取场类型
     */
    FieldType getType() const;

    /**
     * @brief 设置物理单位
     */
    void setUnit(Unit unit);

    /**
     * @brief 获取物理单位
     */
    Unit getUnit() const;

    /**
     * @brief 获取场名称
     */
    virtual std::string getName() const = 0;

    /**
     * @brief 重置场值
     */
    virtual void reset() = 0;

    /**
     * @brief 获取场的统计信息
     */
    virtual std::string getStatistics() const = 0;

  protected:
    VolMeshPtr mesh_; ///< 网格指针
    FieldType type_;  ///< 场类型
    Unit unit_;       ///< 物理单位
};

/**
 * @brief 标量体积场类
 */
class ScalarVolField : public VolField
{
  public:
    /**
     * @brief 构造函数
     * @param mesh 网格指针
     */
    explicit ScalarVolField(const VolMeshPtr& mesh);

    /**
     * @brief 构造函数（带初始值）
     * @param mesh 网格指针
     * @param values 初始场值
     */
    ScalarVolField(const VolMeshPtr& mesh, const std::vector<double>& values);

    /**
     * @brief 析构函数
     */
    ~ScalarVolField() override = default;

    /**
     * @brief 获取场值
     */
    const std::vector<double>& getValues() const
    {
        return values_;
    }

    /**
     * @brief 设置场值
     */
    void setValues(const std::vector<double>& values);

    /**
     * @brief 获取指定位置的场值
     */
    double getValue(size_t index) const;

    /**
     * @brief 设置指定位置的场值
     */
    void setValue(size_t index, double value);

    /**
     * @brief 插值获取场值
     */
    double interpolate(const SCDAT::Geometry::Point3D& position) const;

    /**
     * @brief 重置场值
     */
    void reset() override;

    /**
     * @brief 获取场名称
     */
    std::string getName() const override
    {
        return "ScalarField";
    }

    /**
     * @brief 获取场的统计信息
     */
    std::string getStatistics() const override;

    /**
     * @brief 计算梯度场
     */
    VectorVolFieldPtr computeGradient() const;

    /**
     * @brief 计算拉普拉斯算子
     */
    ScalarVolFieldPtr computeLaplacian() const;

    /**
     * @brief 场值相加
     * @param other 另一个标量场
     */
    void add(const ScalarVolField& other);

    /**
     * @brief 场值乘以标量
     * @param factor 乘数
     */
    void multiply(double factor);

    /**
     * @brief 应用函数到所有场值
     * @param func 要应用的函数
     */
    void applyFunction(std::function<double(double)> func);

    /**
     * @brief 计算场值总和
     * @return 总和
     */
    double computeSum() const;

    /**
     * @brief 计算场值最大值
     * @return 最大值
     */
    double computeMax() const;

    /**
     * @brief 计算场值最小值
     * @return 最小值
     */
    double computeMin() const;

  protected:
    std::vector<double> values_; ///< 标量场值

    /**
     * @brief 规则网格三线性插值
     * @param position 插值位置
     * @param gridSize 网格尺寸
     * @return 插值结果
     */
    double interpolateRegularGrid(const SCDAT::Geometry::Point3D& position, size_t gridSize) const;

    /**
     * @brief 非规则网格插值（反距离加权）
     * @param position 插值位置
     * @return 插值结果
     */
    double interpolateIrregularGrid(const SCDAT::Geometry::Point3D& position) const;
};

/**
 * @brief 向量体积场类
 */
class VectorVolField : public VolField
{
  public:
    /**
     * @brief 构造函数
     * @param mesh 网格指针
     */
    explicit VectorVolField(const VolMeshPtr& mesh);

    /**
     * @brief 构造函数（带初始值）
     * @param mesh 网格指针
     * @param values 初始场值
     */
    VectorVolField(const VolMeshPtr& mesh, const std::vector<SCDAT::Geometry::Vector3D>& values);

    /**
     * @brief 析构函数
     */
    ~VectorVolField() override = default;

    /**
     * @brief 获取场值
     */
    const std::vector<SCDAT::Geometry::Vector3D>& getValues() const
    {
        return values_;
    }

    /**
     * @brief 设置场值
     */
    void setValues(const std::vector<SCDAT::Geometry::Vector3D>& values);

    /**
     * @brief 获取指定位置的场值
     */
    SCDAT::Geometry::Vector3D getValue(size_t index) const;

    /**
     * @brief 设置指定位置的场值
     */
    void setValue(size_t index, const SCDAT::Geometry::Vector3D& value);

    /**
     * @brief 插值获取场值
     */
    SCDAT::Geometry::Vector3D interpolate(const SCDAT::Geometry::Point3D& position) const;

    /**
     * @brief 重置场值
     */
    void reset() override;

    /**
     * @brief 获取场名称
     */
    std::string getName() const override
    {
        return "VectorField";
    }

    /**
     * @brief 获取场的统计信息
     */
    std::string getStatistics() const override;

    /**
     * @brief 计算散度场
     */
    ScalarVolFieldPtr computeDivergence() const;

    /**
     * @brief 计算旋度场
     */
    VectorVolFieldPtr computeCurl() const;

    /**
     * @brief 计算场的模长
     */
    ScalarVolFieldPtr computeMagnitude() const;

    /**
     * @brief 场值相加
     * @param other 另一个向量场
     */
    void add(const VectorVolField& other);

    /**
     * @brief 场值乘以标量
     * @param factor 乘数
     */
    void multiply(double factor);

    /**
     * @brief 计算场模长的最大值
     * @return 最大模长
     */
    double computeMagnitudeMax() const;

  protected:
    std::vector<SCDAT::Geometry::Vector3D> values_; ///< 向量场值

    /**
     * @brief 规则网格向量三线性插值
     * @param position 插值位置
     * @param gridSize 网格尺寸
     * @return 插值结果
     */
    SCDAT::Geometry::Vector3D interpolateRegularGrid(const SCDAT::Geometry::Point3D& position, size_t gridSize) const;

    /**
     * @brief 非规则网格向量插值（反距离加权）
     * @param position 插值位置
     * @return 插值结果
     */
    SCDAT::Geometry::Vector3D interpolateIrregularGrid(const SCDAT::Geometry::Point3D& position) const;
};

/**
 * @brief 电场类
 * @details 表示三维空间中的电场分布，继承自向量体积场类
 *          支持从电势计算电场、添加均匀电场、计算能量密度等功能
 *
 * 电场的基本性质：
 * - 单位：V/m (伏特每米)
 * - 物理意义：单位正电荷受到的电力
 * - 数学关系：E = -∇φ (电场等于电势的负梯度)
 *
 * 主要应用场景：
 * - PIC仿真中的粒子推进
 * - 静电场分析
 * - 电容器设计
 * - 等离子体物理仿真
 *
 * @note 电场方向指向电势降低的方向
 * @warning 强电场区域可能导致数值不稳定，需要适当的网格细化
 */
class ElectricField : public VectorVolField
{
  public:
    /**
     * @brief 构造函数
     * @param mesh 网格指针，不能为空
     * @details 创建电场对象并初始化为零场
     * @throws std::invalid_argument 当mesh为空指针时
     */
    explicit ElectricField(const VolMeshPtr& mesh);

    /**
     * @brief 从电势场计算电场
     * @param potential 电势场对象，必须与当前网格兼容
     * @details 使用数值微分计算电场：E = -∇φ
     *
     * 计算过程：
     * 1. 对每个网格节点计算电势梯度
     * 2. 电场等于梯度的负值
     * 3. 更新内部电场数据
     *
     * @note 要求电势场与电场使用相同的网格
     * @warning 电势场的边界条件会影响电场计算精度
     * @throws std::invalid_argument 当网格不匹配时
     */
    void computeFromPotential(const ScalarVolField& potential);

    /**
     * @brief 添加均匀电场
     * @param uniform_field 均匀电场向量 (V/m)
     * @details 在现有电场基础上叠加均匀电场
     *
     * 应用场景：
     * - 模拟平行板电容器
     * - 添加外部电场
     * - 设置初始条件
     *
     * @note 该操作是累加性的，不会覆盖现有电场
     * @example
     *   ElectricField field(mesh);
     *   field.addUniformField(Vector3D(1000, 0, 0)); // 添加X方向1000V/m电场
     */
    void addUniformField(const SCDAT::Geometry::Vector3D& uniform_field);

    /**
     * @brief 获取场名称
     * @return 返回"ElectricField"
     * @details 用于调试和日志输出
     */
    std::string getName() const override
    {
        return "ElectricField";
    }

    /**
     * @brief 计算电场能量密度
     * @return 电场能量密度场指针
     * @details 计算电场的能量密度分布：u = ½ε₀|E|²
     *
     * 其中：
     * - u: 能量密度 (J/m³)
     * - ε₀: 真空介电常数 (8.854×10⁻¹² F/m)
     * - |E|: 电场强度的模长
     *
     * @note 返回的是标量场，表示空间各点的能量密度
     * @warning 在强电场区域能量密度可能很大，需要注意数值范围
     */
    ScalarVolFieldPtr computeEnergyDensity() const;
};

/**
 * @brief 磁场类
 * @details 表示三维空间中的磁场分布，继承自向量体积场类
 *          支持均匀磁场设置、偶极子磁场添加、磁感应强度计算等功能
 *
 * 磁场的基本性质：
 * - 单位：T (特斯拉) 或 A/m (安培每米)
 * - 物理意义：描述磁力线的方向和强度
 * - 数学关系：∇·B = 0 (磁场散度为零，无磁单极子)
 *
 * 主要应用场景：
 * - 带电粒子在磁场中的运动
 * - 等离子体约束
 * - 磁流体力学仿真
 * - 空间环境模拟
 *
 * @note 磁场线是闭合的，没有起点和终点
 * @warning 强磁场会显著影响粒子轨迹，需要适当的时间步长
 */
class MagneticField : public VectorVolField
{
  public:
    /**
     * @brief 构造函数
     * @param mesh 网格指针，不能为空
     * @details 创建磁场对象并初始化为零场
     * @throws std::invalid_argument 当mesh为空指针时
     */
    explicit MagneticField(const VolMeshPtr& mesh);

    /**
     * @brief 设置均匀磁场
     * @param uniform_field 均匀磁场向量 (T)
     * @details 将整个计算域设置为均匀磁场，覆盖现有磁场分布
     *
     * 应用场景：
     * - 模拟地球磁场
     * - 托卡马克装置的环向磁场
     * - 实验室均匀磁场环境
     *
     * @note 该操作会覆盖现有磁场，而不是叠加
     * @example
     *   MagneticField field(mesh);
     *   field.setUniformField(Vector3D(0, 0, 1e-4)); // 设置Z方向1高斯磁场
     */
    void setUniformField(const SCDAT::Geometry::Vector3D& uniform_field);

    /**
     * @brief 获取场名称
     * @return 返回"MagneticField"
     * @details 用于调试和日志输出
     */
    std::string getName() const override
    {
        return "MagneticField";
    }

    /**
     * @brief 计算磁场能量密度
     * @return 磁场能量密度场指针
     * @details 计算磁场的能量密度分布：u = |B|²/(2μ₀)
     *
     * 其中：
     * - u: 能量密度 (J/m³)
     * - μ₀: 真空磁导率 (4π×10⁻⁷ H/m)
     * - |B|: 磁感应强度的模长
     *
     * @note 返回的是标量场，表示空间各点的磁能密度
     * @warning 在强磁场区域能量密度可能很大
     */
    ScalarVolFieldPtr computeEnergyDensity() const;

    /**
     * @brief 计算磁感应强度
     * @return 磁感应强度场指针
     * @details 在真空中，磁感应强度B与磁场强度H的关系：B = μ₀H
     *
     * 在介质中的关系：B = μH = μ₀μᵣH
     * 其中μᵣ是相对磁导率
     *
     * @note 在真空中磁场和磁感应强度仅相差常数因子
     * @warning 在磁性材料中需要考虑磁化效应
     */
    VectorVolFieldPtr computeMagneticInduction() const;

    /**
     * @brief 添加磁偶极子产生的磁场
     * @param dipolePosition 磁偶极子位置坐标
     * @param dipoleMoment 磁偶极矩向量 (A·m²)
     * @details 在现有磁场基础上叠加磁偶极子产生的磁场
     *
     * 磁偶极子场的数学表达式：
     * B(r) = (μ₀/4π) * [3(m·r̂)r̂ - m] / r³
     *
     * 其中：
     * - m: 磁偶极矩
     * - r: 从偶极子到场点的位置矢量
     * - r̂: 单位方向矢量
     *
     * 应用场景：
     * - 地球磁场建模
     * - 航天器磁场干扰
     * - 磁性材料建模
     *
     * @note 该操作是累加性的，可以添加多个磁偶极子
     * @warning 在偶极子附近磁场变化剧烈，需要足够的网格分辨率
     */
    void addDipoleField(const SCDAT::Geometry::Point3D& dipolePosition, const SCDAT::Geometry::Vector3D& dipoleMoment);
};

/**
 * @brief 温度场类
 */
class TemperatureField : public ScalarVolField
{
  public:
    /**
     * @brief 构造函数
     * @param mesh 网格指针
     */
    explicit TemperatureField(const VolMeshPtr& mesh);

    /**
     * @brief 设置均匀温度
     * @param temperature 温度值 (K)
     */
    void setUniformTemperature(double temperature);

    /**
     * @brief 获取场名称
     */
    std::string getName() const override
    {
        return "TemperatureField";
    }

    /**
     * @brief 计算热流密度
     * @param thermal_conductivity 热导率
     */
    VectorVolFieldPtr computeHeatFlux(double thermal_conductivity) const;

    /**
     * @brief 添加温度梯度
     * @param gradient 温度梯度
     * @param reference 参考点
     */
    void addGradient(const SCDAT::Geometry::Vector3D& gradient, const SCDAT::Geometry::Point3D& reference);

    /**
     * @brief 应用热源
     * @param heatSource 热源场
     * @param thermalDiffusivity 热扩散率
     * @param dt 时间步长
     */
    void applyHeatSource(const ScalarVolField& heatSource, double thermalDiffusivity, double dt);

  private:
    /**
     * @brief 计算拉普拉斯算子
     * @param nodeIndex 节点索引
     * @return 拉普拉斯算子值
     */
    double computeLaplacian(size_t nodeIndex) const;
};

/**
 * @brief 电流密度场类
 */
class CurrentDensityField : public VectorVolField
{
  public:
    /**
     * @brief 构造函数
     * @param mesh 网格指针
     */
    explicit CurrentDensityField(const VolMeshPtr& mesh);

    /**
     * @brief 从电场和电导率计算电流密度
     * @param electric_field 电场
     * @param conductivity_field 电导率场
     */
    void computeFromElectricField(const ElectricField& electric_field,
                                  const ScalarVolField& conductivity_field);

    /**
     * @brief 添加对流电流
     * @param charge_density 电荷密度
     * @param velocity_field 速度场
     */
    void addConvectionCurrent(const ScalarVolField& charge_density,
                              const VectorVolField& velocity_field);

    /**
     * @brief 计算总电流
     * @return 总电流 (A)
     */
    double computeTotalCurrent() const;

    /**
     * @brief 检查电流连续性
     * @return 最大散度值
     */
    double checkCurrentContinuity() const;

    /**
     * @brief 获取场名称
     */
    std::string getName() const override
    {
        return "CurrentDensityField";
    }
};

} // namespace Field
} // namespace SCDAT

#endif // SCDAT_FIELD_VOL_FIELD_H
