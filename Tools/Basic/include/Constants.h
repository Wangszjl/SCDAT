/**
 * @file Constants.h
 * @brief 常量定义（数学常量与物理常量）
 * @details 将数学常量与物理常量拆分为两个类，避免重复定义并提高可维护性
 *
 * @author Wang Sizhan
 * @date 2026年3月18日 15:34:33
 * @version V0.0.1
 * @ingroup ToolsModule
 */

#ifndef SCDAT_BASIC_CONSTANTS_H
#define SCDAT_BASIC_CONSTANTS_H

#include <algorithm>
#include <cmath>

namespace SCDAT
{
namespace Basic
{
namespace Constants
{

/**
 * @brief 数学常量类
 */
struct MathConstants
{
    static constexpr double Pi = 3.14159265358979323846;
    static constexpr double TwoPi = 2.0 * Pi;
    static constexpr double HalfPi = Pi / 2.0;
    static constexpr double E = 2.71828182845904523536;
    static constexpr double Sqrt2 = 1.41421356237309504880;
    static constexpr double Sqrt3 = 1.73205080756887729353;
    static constexpr double Ln2 = 0.69314718055994530942;
    static constexpr double Ln10 = 2.30258509299404568402;
};

/**
 * @brief 物理常量类
 */
struct PhysicsConstants
{
    // 基本物理常数 (SI)
    static constexpr double ElementaryCharge = 1.602176634e-19;    // C
    static constexpr double ElectronMass = 9.1093837015e-31;       // kg
    static constexpr double ProtonMass = 1.67262192369e-27;        // kg
    static constexpr double NeutronMass = 1.67492749804e-27;       // kg
    static constexpr double AtomicMassUnit = 1.66053906660e-27;    // kg
    static constexpr double BoltzmannConstant = 1.380649e-23;      // J/K
    static constexpr double BoltzmannConstantEV = 8.617333262e-5;  // eV/K
    static constexpr double VacuumPermittivity = 8.8541878128e-12; // F/m
    static constexpr double VacuumPermeability = 1.25663706212e-6; // H/m
    static constexpr double SpeedOfLight = 2.99792458e8;           // m/s
    static constexpr double PlanckConstant = 6.62607015e-34;       // J*s
    static constexpr double ReducedPlanckConstant = PlanckConstant / MathConstants::TwoPi;
    static constexpr double AvogadroConstant = 6.02214076e23; // 1/mol
    static constexpr double GasConstant = 8.314462618;        // J/(mol*K)

    // 粒子电荷与特征量
    static constexpr double ElectronCharge = -ElementaryCharge;
    static constexpr double ProtonCharge = ElementaryCharge;
    static constexpr double ElectronClassicalRadius = 2.8179403262e-15;    // m
    static constexpr double ProtonClassicalRadius = 1.5346982266e-18;      // m
    static constexpr double ElectronComptonWavelength = 2.42631023867e-12; // m
    static constexpr double ProtonComptonWavelength = 1.32140985539e-15;   // m
    static constexpr double BohrMagneton = 9.2740100783e-24;               // J/T
    static constexpr double NuclearMagneton = 5.0507837461e-27;            // J/T

    // 场发射与热发射
    static constexpr double FowlerNordheimA = 1.54e-6;
    static constexpr double FowlerNordheimB = 6.83e9;
    static constexpr double SchottkyCoefficient = 3.79e-5;
    static constexpr double RichardsonConstant = 1.20173e6; // A/(m^2*K^2)

    // 单位与量纲转换
    static constexpr double EVToJoule = ElementaryCharge;                 // J/eV
    static constexpr double JouleToEV = 1.0 / EVToJoule;                  // eV/J
    static constexpr double EVToK = ElementaryCharge / BoltzmannConstant; // K/eV
    static constexpr double KToEV = BoltzmannConstant / ElementaryCharge; // eV/K
    static constexpr double HartreeToEV = 27.211386245988;
    static constexpr double RydbergToEV = 13.605693122994;
    static constexpr double BohrRadius = 5.29177210903e-11; // m
    static constexpr double AngstromToMeter = 1e-10;
    static constexpr double MeterToAngstrom = 1e10;
    static constexpr double FemtosecondToSecond = 1e-15;
    static constexpr double PicosecondToSecond = 1e-12;
    static constexpr double NanosecondToSecond = 1e-9;
    static constexpr double MicrosecondToSecond = 1e-6;
    static constexpr double CelsiusToKelvinOffset = 273.15;
    static constexpr double AtmosphereToPascal = 101325.0;
    static constexpr double TorrToPascal = 133.322368;
    static constexpr double BarToPascal = 1e5;
    static constexpr double TorrToPa = TorrToPascal;
    static constexpr double PaToTorr = 1.0 / TorrToPascal;

    // 空间环境典型值
    static constexpr double SolarConstant = 1361.0;
    static constexpr double GEOHeight = 35786000.0;
    static constexpr double LEOHeight = 400000.0;
    static constexpr double EarthRadius = 6371000.0;

    // 等离子体相关经验常数
    static constexpr double DebyeLengthFactor = 7.43e3;
    static constexpr double PlasmaFrequencyFactor = 8.98e3;
    static constexpr double CyclotronFrequencyFactor = ElementaryCharge / ElectronMass;

    // 数值计算常量
    static constexpr double MachineEpsilon = 2.220446049250313e-16;
    static constexpr double SmallNumber = 1e-15;
    static constexpr double LargeNumber = 1e15;
    static constexpr double ZeroTolerance = 1e-12;
    static constexpr double DefaultTolerance = 1e-10;
    static constexpr double SolverDefaultTolerance = 1e-6;
    static constexpr int DefaultMaxIterations = 1000;
    static constexpr double Epsilon = 1e-30;

    // 网格与仿真控制
    static constexpr int MaxInterpolationOrder = 4;
    static constexpr int DefaultGridRefinement = 2;
    static constexpr double CflSafetyFactor = 0.5;
    static constexpr double DefaultParticleWeight = 1e12;
    static constexpr int MaxParticlesPerCell = 100;
    static constexpr double ParticleMergeThreshold = 0.1;
    static constexpr int MaxPoissonIterations = 1000;
    static constexpr double PoissonConvergenceTolerance = 1e-8;
    static constexpr double FieldSmoothingFactor = 0.1;
    static constexpr double CollisionNullThreshold = 1e-20;
    static constexpr double MaxCollisionFrequency = 1e15;
    static constexpr double ElasticCollisionEnergyThreshold = 1e-3;
    static constexpr double SecondaryElectronYieldMax = 10.0;
    static constexpr double WorkFunctionTypical = 4.5;

    // 材料属性
    static constexpr double VacuumConductivity = 0.0;
    static constexpr double PerfectConductorConductivity = 1e20;
    static constexpr double DielectricLossTangentTypical = 1e-4;
};

// 单位转换宏
#define DEGREES_TO_RADIANS(deg) ((deg) * ::SCDAT::Basic::Constants::MathConstants::Pi / 180.0)
#define RADIANS_TO_DEGREES(rad) ((rad) * 180.0 / ::SCDAT::Basic::Constants::MathConstants::Pi)
#define EV_TO_JOULES(ev) ((ev) * ::SCDAT::Basic::Constants::PhysicsConstants::EVToJoule)
#define JOULES_TO_EV(j) ((j) * ::SCDAT::Basic::Constants::PhysicsConstants::JouleToEV)
#define CELSIUS_TO_KELVIN(c)                                                                       \
    ((c) + ::SCDAT::Basic::Constants::PhysicsConstants::CelsiusToKelvinOffset)
#define KELVIN_TO_CELSIUS(k)                                                                       \
    ((k) - ::SCDAT::Basic::Constants::PhysicsConstants::CelsiusToKelvinOffset)

// 物理量计算宏
#define THERMAL_VELOCITY(T, m)                                                                     \
    (std::sqrt(3.0 * ::SCDAT::Basic::Constants::PhysicsConstants::BoltzmannConstant * (T) / (m)))
#define DEBYE_LENGTH(T, n)                                                                         \
    (std::sqrt(::SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity *                   \
               ::SCDAT::Basic::Constants::PhysicsConstants::BoltzmannConstant * (T) /              \
               ((n) * ::SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge *              \
                ::SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge)))
#define PLASMA_FREQUENCY(n, m)                                                                     \
    (std::sqrt((n) * ::SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge *               \
               ::SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge /                     \
               (::SCDAT::Basic::Constants::PhysicsConstants::VacuumPermittivity * (m))))
#define CYCLOTRON_FREQUENCY(B, m)                                                                  \
    (::SCDAT::Basic::Constants::PhysicsConstants::ElementaryCharge * (B) / (m))

// 数值计算宏
#define IS_ZERO(x) ((std::abs(x)) < ::SCDAT::Basic::Constants::PhysicsConstants::ZeroTolerance)
#define IS_EQUAL(a, b)                                                                             \
    ((std::abs((a) - (b))) < ::SCDAT::Basic::Constants::PhysicsConstants::DefaultTolerance)
#define CLAMP(x, min_val, max_val) (std::max((min_val), std::min((x), (max_val))))
#define SIGN(x) (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0))

} // namespace Constants
} // namespace Basic
} // namespace SCDAT

#endif // SCDAT_BASIC_CONSTANTS_H
