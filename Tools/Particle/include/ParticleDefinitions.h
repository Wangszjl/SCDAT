/**
 * @file ParticleDefinitions.h
 * @brief 粒子定义模块
 *
 * 基于C++20的粒子系统，支持SIMD、并行算法、内存优化，主要内容包括粒子类型、粒子事件定义，粒子属性访问和设置，以及粒子推进和时间刷新算法
 *
 * @author Wang Sizhan
 * @date 2026年3月18日 16:23:36
 * @ingroup ParticleModule
 */

#pragma once

#include "../Basic/include/Constants.h"
#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
// #include "SCDAT/core/ErrorHandling.h"
// #include "SCDAT/core/Logger.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <ranges>
#include <string>
#include <vector>

namespace SCDAT
{
namespace Particle
{

using ParticleId = std::uint64_t;
using ParticleTypeId = std::uint32_t;

class Particle;
using ParticlePtr = Particle*;

// ============================================================================
// 粒子定义
// ============================================================================

/**
 * @brief 粒子类型
 */
template <typename T>
concept ParticleLike = requires(T particle) {
    { particle.getPosition() } -> std::convertible_to<SCDAT::Geometry::Point3D>;
    { particle.getVelocity() } -> std::convertible_to<SCDAT::Geometry::Vector3D>;
    { particle.getMass() } -> std::convertible_to<double>;
    { particle.getCharge() } -> std::convertible_to<double>;
    { particle.getId() } -> std::convertible_to<std::size_t>;
};

/**
 * @brief 可推进粒子概念
 */
template <typename T>
concept PushableParticle = ParticleLike<T> && requires(T particle) {
    { particle.getKineticEnergy() } -> std::convertible_to<double>;
};

// ============================================================================
// 粒子类型枚举
// ============================================================================

/**
 * @brief 粒子类型
 */
enum class ParticleType : std::uint8_t
{
    UNKNOWN = 0,

    // --- 轻子 ---
    ELECTRON = 1,
    POSITRON = 2,
    MUON_MINUS = 3,
    MUON_PLUS = 4,
    TAU_MINUS = 5,
    TAU_PLUS = 6,
    ELECTRON_NEUTRINO = 7,
    MUON_NEUTRINO = 8,
    TAU_NEUTRINO = 9,

    // --- 光子 ---
    PHOTON = 10,
    GAMMA = 11,          // 特指高能γ射线
    XRAY = 12,           // X射线
    OPTICAL_PHOTON = 13, // 光学光子（可见光/UV/IR）

    // --- 重子 ---
    PROTON = 20,
    ANTIPROTON = 21,
    NEUTRON = 22,
    ANTINEUTRON = 23,

    // --- 介子 ---
    PION_PLUS = 30,
    PION_MINUS = 31,
    PION_ZERO = 32,
    KAON_PLUS = 33,
    KAON_MINUS = 34,

    // --- 离子/原子核 ---
    ALPHA = 40,     // He-4 核
    DEUTERON = 41,  // 氘核
    TRITON = 42,    // 氚核
    HELION = 43,    // He-3 核
    HEAVY_ION = 44, // 重离子（通用）

    // --- 原子/分子 ---
    ATOM = 50,     // 原子
    MOLECULE = 51, // 分子
    ION = 52,      // 离子
    RADICAL = 53,  // 自由基

    // --- 复合粒子 ---
    CLUSTER = 60,      // 团簇
    NANOPARTICLE = 61, // 纳米粒子
    DUST = 62,         // 尘埃粒子
    DROPLET = 63,      // 液滴
    BUBBLE = 64,       // 气泡

    // --- 准粒子 ---
    PHONON = 70,      // 声子
    PLASMON = 71,     // 等离激元
    MAGNON = 72,      // 磁振子
    EXCITON = 73,     // 激子
    POLARON = 74,     // 极化子
    HOLE = 75,        // 空穴
    COOPER_PAIR = 76, // 库珀对

    // --- 虚拟/追踪粒子 ---
    VIRTUAL_PHOTON = 80, // 虚光子
    GEANTINO = 81,       // 几何追踪粒子（Geant4）
    TRACER = 82,         // 示踪粒子

    // --- 发射机理扩展---
    PHOTOELECTRON = 90,           // 光电子
    SECONDARY_ELECTRON = 91,      // 二次电子
    BACKSCATTERED_ELECTRON = 92,  // 背散射电子
    AUGER_ELECTRON = 93,          // 俄歇电子
    THERMAL_ELECTRON = 94,        // 热电子
    FIELD_EMISSION_ELECTRON = 95, // 场发射电子
    POSITIVE_ION = 96,            // 正离子
    NEGATIVE_ION = 97,            // 负离子
    MOLECULAR_ION = 98,           // 分子离子
    CLUSTER_ION = 99,             // 团簇离子
    BETA_PARTICLE = 100,          // β粒子

    // --- 用户扩展 ---
    USER_DEFINED_1 = 240,
    USER_DEFINED_2 = 241,
    USER_DEFINED_3 = 242,
    GENERIC = 255,
};

/**
 * @brief 粒子状态
 */
enum class ParticleStatus : std::uint8_t
{
    // --- 初始态 ---
    UNINITIALIZED = 0, // 未初始化
    CREATED = 1,       // 已创建，等待注入模拟

    // --- 活跃态 ---
    ACTIVE = 10,    // 正在被追踪
    SUSPENDED = 11, // 暂停追踪（等待后续处理）
    WAITING = 12,   // 等待队列中（多线程调度）

    // --- 终止态（粒子生命结束）---
    ABSORBED = 20,      // 被材料完全吸收
    ESCAPED = 21,       // 逃逸出模拟边界
    DEAD = 22,          // 能量低于截断阈值
    LOST = 23,          // 权重过低被俄罗斯轮盘赌淘汰
    DECAYED = 24,       // 粒子衰变
    ANNIHILATED = 25,   // 湮灭（如正负电子对湮灭）
    CAPTURED = 26,      // 被俘获（如中子俘获）
    FUSED = 27,         // 聚变
    FISSIONED = 28,     // 裂变
    THERMALIZED = 29,   // 热化（能量降至热平衡）
    STUCK = 30,         // 卡住（数值问题，无法继续）
    TIME_EXCEEDED = 31, // 超过最大模拟时间
    STEP_EXCEEDED = 32, // 超过最大步数限制
    RECOMBINED = 33,    // 复合

    // --- 错误态 ---
    ERROR = 250,   // 发生错误
    INVALID = 255, // 无效状态
};

// ============================================================================
// 现代化粒子数据结构
// ============================================================================

/**
 * @brief 对齐的粒子数据
 *
 * 使用SIMD友好的内存布局优化性能
 */
struct alignas(64) ParticleData
{
    // 位置和速度（连续存储以优化SIMD）
    std::array<double, 3> position; // x, y, z
    std::array<double, 3> velocity; // vx, vy, vz

    // 物理属性
    double mass;   // 质量
    double charge; // 电荷
    double weight; // 统计权重
    double energy; // 能量

    // 标识和状态
    std::size_t id;          // 唯一标识符
    ParticleType type;       // 粒子类型
    ParticleStatus status;   // 粒子状态
    std::uint8_t padding[6]; // 填充以保持对齐

    // 时间信息
    double birth_time;       // 产生时间
    double last_update_time; // 最后更新时间
    double age;              // 粒子年龄

    // 物理字段
    double thermal_velocity;    // 热速度
    int mass_number;            // 质量数
    int charge_number;          // 电荷数
    double photon_energy;       // 光子能量
    double work_function;       // 功函数
    double primary_energy;      // 初级粒子能量
    double yield;               // 产额
    double emission_angle;      // 发射角
    double temperature;         // 温度
    double richardson_constant; // 理查德森常数

    // 构造函数
    ParticleData() = default;

    ParticleData(std::size_t particle_id, ParticleType particle_type,
                 const SCDAT::Geometry::Point3D& pos, const SCDAT::Geometry::Vector3D& vel,
                 double particle_mass, double particle_charge, double particle_weight = 1.0)
        : position{pos.x(), pos.y(), pos.z()}, velocity{vel.x(), vel.y(), vel.z()},
          mass(particle_mass), charge(particle_charge), weight(particle_weight),
          energy(0.5 * particle_mass * vel.magnitudeSquared()), id(particle_id),
          type(particle_type), status(ParticleStatus::ACTIVE), padding{}, birth_time(0.0),
          last_update_time(0.0), age(0.0), thermal_velocity(0.0), mass_number(0), charge_number(0),
          photon_energy(0.0), work_function(0.0), primary_energy(0.0), yield(0.0),
          emission_angle(0.0), temperature(0.0), richardson_constant(0.0)
    {
    }

    // 访问器
    [[nodiscard]] SCDAT::Geometry::Point3D getPosition() const noexcept
    {
        return SCDAT::Geometry::Point3D(position[0], position[1], position[2]);
    }

    [[nodiscard]] SCDAT::Geometry::Vector3D getVelocity() const noexcept
    {
        return SCDAT::Geometry::Vector3D(velocity[0], velocity[1], velocity[2]);
    }

    void setPosition(const SCDAT::Geometry::Point3D& pos) noexcept
    {
        position[0] = pos.x();
        position[1] = pos.y();
        position[2] = pos.z();
    }

    void setVelocity(const SCDAT::Geometry::Vector3D& vel) noexcept
    {
        velocity[0] = vel.x();
        velocity[1] = vel.y();
        velocity[2] = vel.z();
        updateEnergy();
    }

    void updateEnergy() noexcept
    {
        double v_squared =
            velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
        energy = 0.5 * mass * v_squared;
    }

    [[nodiscard]] double getKineticEnergy() const noexcept
    {
        return energy;
    }

    [[nodiscard]] double getSpeed() const noexcept
    {
        return std::sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] +
                         velocity[2] * velocity[2]);
    }

    [[nodiscard]] double getMomentum() const noexcept
    {
        return mass * getSpeed();
    }

    [[nodiscard]] SCDAT::Geometry::Vector3D getMomentumVector() const noexcept
    {
        return getVelocity() * mass;
    }
};

// ============================================================================
// 现代化粒子类
// ============================================================================

/**
 * @brief 现代化粒子类
 */
class Particle
{
  public:
    using Data = ParticleData;

    // 构造函数
    explicit Particle(Data data) : data_(std::move(data)) {}

    Particle(std::size_t id, ParticleType type, const SCDAT::Geometry::Point3D& position,
             const SCDAT::Geometry::Vector3D& velocity, double mass, double charge,
             double weight = 1.0)
        : data_(id, type, position, velocity, mass, charge, weight)
    {
    }

    // 默认构造、拷贝、移动 / Default, copy, move
    Particle() = default;
    Particle(const Particle&) = default;
    Particle& operator=(const Particle&) = default;
    Particle(Particle&&) noexcept = default;
    Particle& operator=(Particle&&) noexcept = default;
    ~Particle() = default;

    // 基本访问器 / Basic accessors
    [[nodiscard]] std::size_t getId() const noexcept
    {
        return data_.id;
    }
    [[nodiscard]] ParticleType getType() const noexcept
    {
        return data_.type;
    }
    [[nodiscard]] ParticleStatus getStatus() const noexcept
    {
        return data_.status;
    }

    [[nodiscard]] SCDAT::Geometry::Point3D getPosition() const noexcept
    {
        return data_.getPosition();
    }
    [[nodiscard]] SCDAT::Geometry::Vector3D getVelocity() const noexcept
    {
        return data_.getVelocity();
    }

    [[nodiscard]] double getMass() const noexcept
    {
        return data_.mass;
    }
    [[nodiscard]] double getCharge() const noexcept
    {
        return data_.charge;
    }
    [[nodiscard]] double getWeight() const noexcept
    {
        return data_.weight;
    }
    [[nodiscard]] double getAge() const noexcept
    {
        return data_.age;
    }
    [[nodiscard]] double getKineticEnergy() const noexcept
    {
        return data_.getKineticEnergy();
    }
    [[nodiscard]] double getSpeed() const noexcept
    {
        return data_.getSpeed();
    }

    // 修改器
    void setPosition(const SCDAT::Geometry::Point3D& position)
    {
        data_.setPosition(position);
    }
    void setVelocity(const SCDAT::Geometry::Vector3D& velocity) noexcept
    {
        data_.setVelocity(velocity);
    }
    void setStatus(ParticleStatus status) noexcept
    {
        data_.status = status;
    }
    void setWeight(double weight) noexcept
    {
        data_.weight = weight;
    }
    void setCharge(double charge) noexcept
    {
        data_.charge = charge;
    }
    void setMass(double mass) noexcept
    {
        data_.mass = mass;
        data_.updateEnergy();
    }

    // 生命周期管理
    void updateAge(double dt) noexcept
    {
        data_.age += dt;
    }
    void setAge(double age) noexcept
    {
        data_.age = age;
    }
    [[nodiscard]] bool isActive() const noexcept
    {
        return data_.status == ParticleStatus::ACTIVE;
    }
    void markAsAbsorbed() noexcept
    {
        data_.status = ParticleStatus::ABSORBED;
    }
    void markAsEscaped() noexcept
    {
        data_.status = ParticleStatus::ESCAPED;
    }
    void markAsRecombined() noexcept
    {
        data_.status = ParticleStatus::RECOMBINED;
    }

    // 扩展字段访问器
    [[nodiscard]] double getThermalVelocity() const noexcept
    {
        return data_.thermal_velocity;
    }
    void setThermalVelocity(double v) noexcept
    {
        data_.thermal_velocity = v;
    }

    [[nodiscard]] int getMassNumber() const noexcept
    {
        return data_.mass_number;
    }
    [[nodiscard]] int getChargeNumber() const noexcept
    {
        return data_.charge_number;
    }
    void setMassNumber(int n) noexcept
    {
        data_.mass_number = n;
    }
    void setChargeNumber(int n) noexcept
    {
        data_.charge_number = n;
    }

    [[nodiscard]] double getPhotonEnergy() const noexcept
    {
        return data_.photon_energy;
    }
    [[nodiscard]] double getWorkFunction() const noexcept
    {
        return data_.work_function;
    }
    void setPhotonEnergy(double e) noexcept
    {
        data_.photon_energy = e;
    }
    void setWorkFunction(double w) noexcept
    {
        data_.work_function = w;
    }

    [[nodiscard]] double getPrimaryEnergy() const noexcept
    {
        return data_.primary_energy;
    }
    [[nodiscard]] double getYield() const noexcept
    {
        return data_.yield;
    }
    [[nodiscard]] double getEmissionAngle() const noexcept
    {
        return data_.emission_angle;
    }
    void setPrimaryEnergy(double e) noexcept
    {
        data_.primary_energy = e;
    }
    void setYield(double y) noexcept
    {
        data_.yield = y;
    }
    void setEmissionAngle(double a) noexcept
    {
        data_.emission_angle = a;
    }

    [[nodiscard]] double getTemperature() const noexcept
    {
        return data_.temperature;
    }
    [[nodiscard]] double getRichardsonConstant() const noexcept
    {
        return data_.richardson_constant;
    }
    void setTemperature(double t) noexcept
    {
        data_.temperature = t;
    }
    void setRichardsonConstant(double c) noexcept
    {
        data_.richardson_constant = c;
    }

    // 计算属性
    [[nodiscard]] double getMomentum() const noexcept
    {
        return data_.getMomentum();
    }
    [[nodiscard]] SCDAT::Geometry::Vector3D getMomentumVector() const noexcept
    {
        return data_.getMomentumVector();
    }

    // 辅助函数
    [[nodiscard]] double getTemperatureEV() const noexcept;
    [[nodiscard]] std::string getSpeciesName() const;
    [[nodiscard]] double getMaxKineticEnergy() const noexcept;

    // 时间更新
    void updateTime(double current_time) noexcept
    {
        data_.last_update_time = current_time;
    }

    // 数据访问
    [[nodiscard]] const Data& getData() const noexcept
    {
        return data_;
    }
    [[nodiscard]] Data& getData() noexcept
    {
        return data_;
    }

    // 比较运算符
    [[nodiscard]] auto operator<=>(const Particle& other) const noexcept
    {
        return data_.id <=> other.data_.id;
    }

    [[nodiscard]] bool operator==(const Particle& other) const noexcept
    {
        return data_.id == other.data_.id;
    }

    // 字符串表示
    [[nodiscard]] std::string toString() const;

    // 克隆
    [[nodiscard]] Particle clone() const
    {
        return *this;
    }

  private:
    Data data_;
};

// 兼容性别名
using ElectronParticle = Particle;
using IonParticle = Particle;
using PhotoelectronParticle = Particle;
using SecondaryElectronParticle = Particle;
using ThermalElectronParticle = Particle;

using ElectronParticlePtr = ParticlePtr;
using IonParticlePtr = ParticlePtr;
using PhotoelectronParticlePtr = ParticlePtr;
using SecondaryElectronParticlePtr = ParticlePtr;
using ThermalElectronParticlePtr = ParticlePtr;

class ParticleFactory
{
  public:
    static Particle createElectron(ParticleId id, const SCDAT::Geometry::Point3D& position,
                                   const SCDAT::Geometry::Vector3D& velocity, double weight = 1.0);

    static Particle createIon(ParticleId id, const SCDAT::Geometry::Point3D& position,
                              const SCDAT::Geometry::Vector3D& velocity, int mass_number,
                              int charge_number, double weight = 1.0);

    static Particle createPhotoelectron(ParticleId id, const SCDAT::Geometry::Point3D& position,
                                        const SCDAT::Geometry::Vector3D& velocity,
                                        double photon_energy, double work_function,
                                        double weight = 1.0);

    static Particle createSecondaryElectron(ParticleId id, const SCDAT::Geometry::Point3D& position,
                                            const SCDAT::Geometry::Vector3D& velocity,
                                            double primary_energy, double yield,
                                            double weight = 1.0);

    static Particle createThermalElectron(ParticleId id, const SCDAT::Geometry::Point3D& position,
                                          const SCDAT::Geometry::Vector3D& velocity,
                                          double temperature, double work_function,
                                          double weight = 1.0);

    static Particle createParticle(ParticleType type, ParticleId id,
                                   const SCDAT::Geometry::Point3D& position,
                                   const SCDAT::Geometry::Vector3D& velocity, double weight = 1.0);

  private:
    static ParticleId next_id_;

  public:
    static ParticleId getNextId()
    {
        return ++next_id_;
    }
};

} // namespace Particle
} // namespace SCDAT
