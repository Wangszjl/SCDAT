/**
 * @file ParticleSource.cpp
 * @brief 粒子源实现（SCDAT）
 */

#include "../include/ParticleSource.h"
#include "../../Basic/include/Constants.h"
#include "../../Mesh/include/MeshParsing.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

namespace SCDAT
{
namespace Particle
{

namespace
{
using SCDAT::Basic::Constants::MathConstants;
using SCDAT::Basic::Constants::PhysicsConstants;

/**
 * @brief 获取线程局部随机数引擎。
 * @details 每个线程独立维护一个 mt19937，避免多线程下共享 RNG 引起竞争与相关性问题。
 */
std::mt19937& globalRng()
{
    thread_local std::mt19937 rng(std::random_device{}());
    return rng;
}

/**
 * @brief 将数值约束为正值。
 * @param value 待检查数值。
 * @param fallback 当 value 非正时返回的回退值。
 */
double clampPositive(double value, double fallback)
{
    return value > 0.0 ? value : fallback;
}

/**
 * @brief 根据粒子类型返回质量（kg）。
 * @details 电子类粒子使用电子质量，离子类粒子使用质子质量；未知类型按质子质量处理。
 */
double getParticleTypeMass(ParticleType type)
{
    switch (type)
    {
    case ParticleType::ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
        return PhysicsConstants::ElectronMass;
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::NEGATIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
        return PhysicsConstants::ProtonMass;
    default:
        return PhysicsConstants::ProtonMass;
    }
}

double getParticleTypeMassAmu(ParticleType type)
{
    return getParticleTypeMass(type) / PhysicsConstants::AtomicMassUnit;
}

double getParticleTypeChargeNumber(ParticleType type)
{
    switch (type)
    {
    case ParticleType::ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
    case ParticleType::NEGATIVE_ION:
        return -1.0;
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
        return 1.0;
    default:
        return 1.0;
    }
}

double thermalSpeedToTemperatureEv(double thermal_speed, double mass_kg)
{
    const double speed = std::max(0.0, thermal_speed);
    return 0.5 * mass_kg * speed * speed / PhysicsConstants::ElementaryCharge;
}

std::vector<double> buildLogEnergyGrid(double minimum_energy_ev, double maximum_energy_ev,
                                       size_t count)
{
    const double e_min = std::max(1.0e-3, minimum_energy_ev);
    const double e_max = std::max(e_min * 1.001, maximum_energy_ev);
    const size_t points = std::max<size_t>(count, 8);
    std::vector<double> grid(points, e_min);
    const double log_min = std::log(e_min);
    const double log_max = std::log(e_max);
    const double denom = static_cast<double>(points - 1);
    for (size_t i = 0; i < points; ++i)
    {
        const double alpha = static_cast<double>(i) / denom;
        grid[i] = std::exp(log_min + alpha * (log_max - log_min));
    }
    return grid;
}

void normalizeFluxShape(std::vector<double>& flux, const std::vector<double>& energy_grid_ev,
                        double target_number_flux_m2_s)
{
    if (flux.empty() || flux.size() != energy_grid_ev.size())
    {
        return;
    }

    double integral = 0.0;
    for (size_t i = 1; i < energy_grid_ev.size(); ++i)
    {
        const double width = std::max(0.0, energy_grid_ev[i] - energy_grid_ev[i - 1]);
        integral += 0.5 * (std::max(0.0, flux[i - 1]) + std::max(0.0, flux[i])) * width;
    }

    if (!(integral > 0.0) || !(target_number_flux_m2_s > 0.0))
    {
        return;
    }

    const double scale = target_number_flux_m2_s / integral;
    for (double& value : flux)
    {
        value = std::max(0.0, value * scale);
    }
}

double estimateNumberFluxFromMoments(const SamplingParameters& params, double mass_kg)
{
    const double density = std::max(0.0, params.density);
    const double thermal_speed = std::max(1.0, params.thermal_speed);
    const double temperature_j =
        thermalSpeedToTemperatureEv(thermal_speed, mass_kg) * PhysicsConstants::ElementaryCharge;
    const double thermal_flux =
        density * std::sqrt(std::max(1.0e-18, temperature_j) / (2.0 * MathConstants::Pi * mass_kg));
    const double advective_flux = density * std::max(0.0, params.bulk_speed);
    return std::max(0.0, thermal_flux + 0.25 * advective_flux);
}

double sampleEnergyFromDiscreteSpectrum(const std::vector<double>& energy_ev,
                                        const std::vector<double>& intensity)
{
    if (energy_ev.empty() || energy_ev.size() != intensity.size())
    {
        return 0.0;
    }

    std::vector<double> cumulative(intensity.size(), 0.0);
    double total = 0.0;
    for (size_t i = 0; i < intensity.size(); ++i)
    {
        total += std::max(0.0, intensity[i]);
        cumulative[i] = total;
    }
    if (!(total > 0.0))
    {
        return energy_ev.front();
    }

    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, total);
    const double pick = u(rng);
    const auto it = std::lower_bound(cumulative.begin(), cumulative.end(), pick);
    const size_t index =
        static_cast<size_t>(std::distance(cumulative.begin(), it == cumulative.end() ? cumulative.end() - 1 : it));
    return std::max(1.0e-3, energy_ev[index]);
}

/**
 * @brief 在单位球面上进行均匀方向采样。
 * @return 长度为 1 的随机方向向量。
 */
Vector3D randomUnitVector()
{
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const double mu = 2.0 * u(rng) - 1.0;
    const double phi = 2.0 * MathConstants::Pi * u(rng);
    const double r_xy = std::sqrt(std::max(0.0, 1.0 - mu * mu));

    return Vector3D(r_xy * std::cos(phi), r_xy * std::sin(phi), mu);
}

/**
 * @brief 计算两组序列的均方根误差（RMSE）。
 * @param y_true 真实值序列。
 * @param y_pred 预测值序列。
 * @return 当输入非法时返回 0，否则返回 RMSE。
 */
double rmse(const std::vector<double>& y_true, const std::vector<double>& y_pred)
{
    if (y_true.empty() || y_true.size() != y_pred.size())
    {
        return 0.0;
    }

    double e2 = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i)
    {
        const double d = y_true[i] - y_pred[i];
        e2 += d * d;
    }
    return std::sqrt(e2 / static_cast<double>(y_true.size()));
}

bool isFinitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

size_t selectWeightedElementByArea(const std::vector<size_t>& element_indices,
                                   const std::shared_ptr<SCDAT::Mesh::SurfaceMesh>& mesh,
                                   double total_area)
{
    if (!mesh || element_indices.empty() || !isFinitePositive(total_area))
    {
        return 0;
    }

    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, total_area);
    const double pick = u(rng);

    double accum = 0.0;
    for (size_t idx : element_indices)
    {
        const double area = mesh->getElementArea(idx);
        if (!isFinitePositive(area))
        {
            continue;
        }

        accum += area;
        if (pick <= accum)
        {
            return idx;
        }
    }

    return element_indices.back();
}

size_t selectWeightedElementByVolume(const std::vector<size_t>& element_indices,
                                     const std::shared_ptr<SCDAT::Mesh::VolumeMesh>& mesh,
                                     double total_volume)
{
    if (!mesh || element_indices.empty() || !isFinitePositive(total_volume))
    {
        return 0;
    }

    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, total_volume);
    const double pick = u(rng);

    double accum = 0.0;
    for (size_t idx : element_indices)
    {
        const double volume = mesh->getElementVolume(idx);
        if (!isFinitePositive(volume))
        {
            continue;
        }

        accum += volume;
        if (pick <= accum)
        {
            return idx;
        }
    }

    return element_indices.back();
}

} // namespace

/**
 * @brief 构造粒子源基类并初始化默认运行参数。
 */
ParticleSource::ParticleSource(const std::string& name, const ParticleType& particleType)
    : name_(name), particleType_(particleType), enabled_(true), emissionRate_(0.0),
      totalEmitted_(0), currentTime_(0.0), startTime_(0.0), stopTime_(-1.0), timeAcceleration_(1.0),
      sampling_model_(SpatialSamplingModel::UNIFORM), sampling_params_{}, last_fit_summary_{}
{
}

/**
 * @brief 粒子源析构函数。
 */
ParticleSource::~ParticleSource() = default;

/**
 * @brief 启用或禁用粒子源。
 */
void ParticleSource::setEnabled(bool enabled)
{
    enabled_ = enabled;
}

/**
 * @brief 查询粒子源是否启用。
 */
bool ParticleSource::isEnabled() const
{
    return enabled_;
}

/**
 * @brief 设置发射率（粒子数/秒）。
 * @throws std::invalid_argument 当 rate < 0。
 */
void ParticleSource::setEmissionRate(double rate)
{
    if (rate < 0.0)
    {
        throw std::invalid_argument("ParticleSource: 发射率不能为负数");
    }
    emissionRate_ = rate;
}

/**
 * @brief 获取当前发射率。
 */
double ParticleSource::getEmissionRate() const
{
    return emissionRate_;
}

/**
 * @brief 设置时间窗 [startTime, stopTime)。
 * @details 当 stopTime < 0 时表示无结束时刻。
 */
void ParticleSource::setTimeWindow(double startTime, double stopTime)
{
    if (stopTime >= 0.0 && stopTime <= startTime)
    {
        throw std::invalid_argument("ParticleSource: 停止时间必须大于开始时间");
    }
    startTime_ = startTime;
    stopTime_ = stopTime;
}

/**
 * @brief 判断给定时刻粒子源是否处于可发射状态。
 */
bool ParticleSource::isActiveAtTime(double time) const
{
    if (!enabled_)
    {
        return false;
    }
    if (time < startTime_)
    {
        return false;
    }
    if (stopTime_ >= 0.0 && time >= stopTime_)
    {
        return false;
    }
    return true;
}

/**
 * @brief 推进粒子源内部时钟。
 */
void ParticleSource::updateTime(double dt)
{
    currentTime_ += dt;
}

/**
 * @brief 获取该粒子源累计发射总数。
 */
size_t ParticleSource::getTotalEmitted() const
{
    return totalEmitted_;
}

/**
 * @brief 获取粒子源名称。
 */
const std::string& ParticleSource::getName() const
{
    return name_;
}

/**
 * @brief 获取粒子类型。
 */
const ParticleType& ParticleSource::getParticleType() const
{
    return particleType_;
}

/**
 * @brief 设置空间/速度采样模型。
 */
void ParticleSource::setSamplingModel(SpatialSamplingModel model)
{
    sampling_model_ = model;
}

/**
 * @brief 获取当前采样模型。
 */
SpatialSamplingModel ParticleSource::getSamplingModel() const
{
    return sampling_model_;
}

/**
 * @brief 设置采样参数并进行物理合理性约束。
 * @details 对密度、热速度、尺度长度等参数做正值或区间裁剪，避免后续采样出现数值异常。
 */
void ParticleSource::setSamplingParameters(const SamplingParameters& params)
{
    sampling_params_ = params;
    sampling_params_.density = clampPositive(sampling_params_.density, 1.0e15);
    sampling_params_.bulk_speed = clampPositive(sampling_params_.bulk_speed, 2.0e5);
    sampling_params_.thermal_speed = clampPositive(sampling_params_.thermal_speed, 5.0e4);
    sampling_params_.hot_thermal_speed = clampPositive(sampling_params_.hot_thermal_speed, 2.0e5);
    sampling_params_.hot_fraction = std::clamp(sampling_params_.hot_fraction, 0.0, 1.0);
    sampling_params_.kappa = std::max(1.51, sampling_params_.kappa);
    sampling_params_.power_law_index = std::max(1.01, sampling_params_.power_law_index);
    sampling_params_.spatial_scale = clampPositive(sampling_params_.spatial_scale, 0.05);
    sampling_params_.min_spatial_scale =
        std::clamp(sampling_params_.min_spatial_scale, 1.0e-7, sampling_params_.spatial_scale);
    sampling_params_.q_parameter = std::clamp(sampling_params_.q_parameter, 1.01, 3.0);
    sampling_params_.characteristic_energy =
        clampPositive(sampling_params_.characteristic_energy, 10.0);
}

/**
 * @brief 获取当前采样参数。
 */
const SamplingParameters& ParticleSource::getSamplingParameters() const
{
    return sampling_params_;
}

/**
 * @brief 设置观测能谱数据。
 * @throws std::invalid_argument 当 energy_ev 与 intensity 长度不一致。
 */
void ParticleSource::setObservedSpectrum(const std::vector<double>& energy_ev,
                                         const std::vector<double>& intensity)
{
    if (energy_ev.size() != intensity.size())
    {
        throw std::invalid_argument("ParticleSource: 能谱横纵坐标长度不一致");
    }
    observed_energy_ev_ = energy_ev;
    observed_intensity_ = intensity;
}

/**
 * @brief 设置探测器电荷剖面。
 * @throws std::invalid_argument 当 radius 与 charge_density 长度不一致。
 */
void ParticleSource::setDetectorChargeProfile(const std::vector<double>& radius,
                                              const std::vector<double>& charge_density)
{
    if (radius.size() != charge_density.size())
    {
        throw std::invalid_argument("ParticleSource: 探测半径与电荷密度长度不一致");
    }
    detector_radius_ = radius;
    detector_charge_density_ = charge_density;
}

/**
 * @brief 基于观测能谱自动拟合采样模型。
 * @details
 * 1) 构造双麦克斯韦、Kappa、幂律三种候选模型；
 * 2) 通过 RMSE 评分选择最优模型；
 * 3) 反写热速度、kappa、幂律指数、空间尺度等参数。
 */
FitSummary ParticleSource::fitSamplingModelFromObservations()
{
    FitSummary summary;

    if (observed_energy_ev_.size() < 5 || observed_intensity_.size() < 5)
    {
        summary.message = "观测能谱数据不足，保留当前分布";
        last_fit_summary_ = summary;
        return summary;
    } // 如果观测能谱数据点数不足5个，直接返回，提示数据不足

    const size_t n = observed_energy_ev_.size();
    std::vector<double> model_double(n, 0.0);
    std::vector<double> model_kappa(n, 0.0);
    std::vector<double> model_power(n, 0.0);

    double sum_w = 0.0;
    double mean_e = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double w = std::max(0.0, observed_intensity_[i]);
        sum_w += w;
        mean_e += observed_energy_ev_[i] * w;
    }
    mean_e = sum_w > 0.0 ? mean_e / sum_w : 10.0;
    // 计算观测能谱的加权平均能量mean_e，权重为强度

    double e_split = mean_e;
    double cold_sum = 0.0;
    double hot_sum = 0.0;
    double cold_w = 0.0;
    double hot_w = 0.0;
    // 用平均能量作为冷热能区分界点，对观测数据进行简单划分，统计冷热区的能量加权和与权重总和，为后续双麦克斯韦拟合提供初始参数估计依据。

    for (size_t i = 0; i < n; ++i)
    {
        const double e = std::max(1e-9, observed_energy_ev_[i]);
        const double w = std::max(0.0, observed_intensity_[i]);
        if (e <= e_split)
        {
            cold_sum += e * w;
            cold_w += w;
        }
        else
        {
            hot_sum += e * w;
            hot_w += w;
        }
    }
    // 根据冷热区的能量加权和与权重总和，估计双麦克斯韦分布的冷、热组分温度（t1、t2）和热组分占比（hot_fraction）。同时对拟合得到的kappa分布参数和幂律指数进行合理性约束，避免数值异常。

    const double t1 =
        cold_w > 0.0 ? std::max(1e-6, cold_sum / cold_w) : std::max(1e-6, mean_e * 0.5);
    const double t2 = hot_w > 0.0 ? std::max(t1, hot_sum / hot_w) : std::max(t1, mean_e * 1.5);
    const double hot_fraction = sum_w > 0.0 ? std::clamp(hot_w / sum_w, 0.02, 0.98) : 0.2;
    // 计算冷、热区的平均能量t1、t2，并做下限保护；计算热组分占比hot_fraction，并限制在[0.02,
    // 0.98]范围内，避免极端分布。
    double sx = 0.0;
    double sy = 0.0;
    double sxx = 0.0;
    double sxy = 0.0;
    double m = 0.0;
    for (size_t i = n / 2; i < n; ++i)
    {
        const double x = std::log(std::max(1e-9, observed_energy_ev_[i]));
        const double y = std::log(std::max(1e-12, observed_intensity_[i]));
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        m += 1.0;
    }
    // 对后半段能谱做对数线性回归，拟合幂律斜率，并据此估计幂律指数。通过对拟合结果进行合理性约束，确保幂律指数在物理上有意义的范围内。

    double slope = -2.5;
    if (m > 2.0)
    {
        const double denom = m * sxx - sx * sx;
        if (std::abs(denom) > 1e-12)
        {
            slope = (m * sxy - sx * sy) / denom;
        }
    }

    const double fitted_power = std::clamp(-slope, 1.1, 8.0);
    const double fitted_kappa = std::clamp(1.5 + 0.8 * fitted_power, 1.6, 12.0);

    for (size_t i = 0; i < n; ++i)
    {
        const double e = std::max(1e-9, observed_energy_ev_[i]);
        model_double[i] =
            (1.0 - hot_fraction) * std::exp(-e / t1) + hot_fraction * std::exp(-e / t2);

        const double theta = std::max(1e-6, mean_e / std::max(2.0, fitted_kappa));
        model_kappa[i] = std::pow(1.0 + e / (fitted_kappa * theta), -(fitted_kappa + 1.0));

        model_power[i] = std::pow(e + 1e-3, -fitted_power);
    }
    // 计算三种模型在每个能量点的理论值，分别为双麦克斯韦、Kappa和幂律分布。通过对模型值进行归一化，使其总强度与观测数据一致，确保后续的RMSE评分具有可比性。

    double norm_obs = std::accumulate(observed_intensity_.begin(), observed_intensity_.end(), 0.0);
    norm_obs = norm_obs > 0.0 ? norm_obs : 1.0;

    const double norm_d = std::accumulate(model_double.begin(), model_double.end(), 0.0);
    const double norm_k = std::accumulate(model_kappa.begin(), model_kappa.end(), 0.0);
    const double norm_p = std::accumulate(model_power.begin(), model_power.end(), 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        model_double[i] = model_double[i] * norm_obs / std::max(1e-12, norm_d);
        model_kappa[i] = model_kappa[i] * norm_obs / std::max(1e-12, norm_k);
        model_power[i] = model_power[i] * norm_obs / std::max(1e-12, norm_p);
    }

    summary.score_double_maxwell = rmse(observed_intensity_, model_double);
    summary.score_kappa = rmse(observed_intensity_, model_kappa);
    summary.score_power_law = rmse(observed_intensity_, model_power);

    summary.selected_model = SpatialSamplingModel::DOUBLE_MAXWELL;
    double best = summary.score_double_maxwell;
    if (summary.score_kappa < best)
    {
        best = summary.score_kappa;
        summary.selected_model = SpatialSamplingModel::KAPPA;
    }
    if (summary.score_power_law < best)
    {
        summary.selected_model = SpatialSamplingModel::POWER_LAW;
    }

    sampling_model_ = summary.selected_model;
    sampling_params_.hot_fraction = hot_fraction;
    sampling_params_.power_law_index = fitted_power;
    sampling_params_.kappa = fitted_kappa;

    const double mass = getParticleTypeMass(particleType_);
    const double cold_j = t1 * PhysicsConstants::ElementaryCharge;
    const double hot_j = t2 * PhysicsConstants::ElementaryCharge;
    sampling_params_.thermal_speed = std::sqrt(2.0 * cold_j / mass);
    sampling_params_.hot_thermal_speed = std::sqrt(2.0 * hot_j / mass);

    if (!detector_radius_.empty() && detector_radius_.size() == detector_charge_density_.size())
    {
        double wsum = 0.0;
        double rsum = 0.0;
        double dsum = 0.0;
        for (size_t i = 0; i < detector_radius_.size(); ++i)
        {
            const double w = std::max(0.0, detector_charge_density_[i]);
            wsum += w;
            rsum += std::abs(detector_radius_[i]) * w;
            dsum += w;
        }
        if (wsum > 0.0)
        {
            sampling_params_.spatial_scale = std::max(1e-4, rsum / wsum);
            sampling_params_.density =
                std::max(1e10, dsum / static_cast<double>(detector_charge_density_.size()));
        }
    }

    summary.success = true;
    summary.message = "已根据观测数据完成分布拟合与模型选择";
    last_fit_summary_ = summary;
    return summary;
}

ResolvedSpectrum ParticleSource::resolveSpectrum(SpectrumUsage usage) const
{
    return buildResolvedSpectrum(particleType_, sampling_model_, sampling_params_, usage,
                                 observed_energy_ev_, observed_intensity_, last_fit_summary_);
}

ResolvedSpectrum ParticleSource::buildResolvedSpectrum(
    const ParticleType& particleType, SpatialSamplingModel model, const SamplingParameters& params,
    SpectrumUsage usage, const std::vector<double>& observed_energy_ev,
    const std::vector<double>& observed_intensity, const FitSummary& fit_summary)
{
    ResolvedSpectrum spectrum;
    spectrum.model = model;
    spectrum.fit_summary = fit_summary;

    const double mass_kg = getParticleTypeMass(particleType);
    const double mass_amu = getParticleTypeMassAmu(particleType);
    const double charge_number = getParticleTypeChargeNumber(particleType);
    const double density = std::max(0.0, params.density);
    const double hot_fraction = std::clamp(params.hot_fraction, 0.0, 1.0);
    const double cold_density = density * (1.0 - hot_fraction);
    const double hot_density = density * hot_fraction;
    const double cold_temperature_ev = std::max(1.0e-6, thermalSpeedToTemperatureEv(params.thermal_speed, mass_kg));
    const double hot_temperature_ev =
        std::max(cold_temperature_ev, thermalSpeedToTemperatureEv(params.hot_thermal_speed, mass_kg));

    auto add_population = [&](double population_density, double temperature_ev, double drift_speed,
                              double weight) {
        if (population_density <= 0.0 || temperature_ev <= 0.0)
        {
            return;
        }
        SpectrumPopulation population;
        population.density_m3 = population_density;
        population.temperature_ev = temperature_ev;
        population.drift_speed_m_per_s = drift_speed;
        population.mass_amu = mass_amu;
        population.charge_number = charge_number;
        population.weight = weight;
        spectrum.populations.push_back(population);
    };

    switch (model)
    {
    case SpatialSamplingModel::DOUBLE_MAXWELL:
        add_population(cold_density, cold_temperature_ev, params.bulk_speed, 1.0 - hot_fraction);
        add_population(hot_density, hot_temperature_ev, params.bulk_speed, hot_fraction);
        break;
    case SpatialSamplingModel::SINGLE_MAXWELL:
    case SpatialSamplingModel::MAXWELL_BOLTZMANN:
    case SpatialSamplingModel::Q_DISTRIBUTION:
    case SpatialSamplingModel::UNIFORM:
        add_population(density, cold_temperature_ev, params.bulk_speed, 1.0);
        break;
    case SpatialSamplingModel::KAPPA:
    case SpatialSamplingModel::POWER_LAW:
    case SpatialSamplingModel::TABULATED:
        add_population(density, cold_temperature_ev, params.bulk_speed, 1.0);
        break;
    }

    const bool need_discrete_spectrum =
        usage != SpectrumUsage::SamplingOnly || model == SpatialSamplingModel::KAPPA ||
        model == SpatialSamplingModel::POWER_LAW || model == SpatialSamplingModel::TABULATED ||
        !observed_energy_ev.empty();
    if (!need_discrete_spectrum)
    {
        return spectrum;
    }

    if (!observed_energy_ev.empty() && observed_energy_ev.size() == observed_intensity.size())
    {
        spectrum.energy_grid_ev = observed_energy_ev;
        spectrum.differential_number_flux = observed_intensity;
        if (model == SpatialSamplingModel::UNIFORM)
        {
            spectrum.model = SpatialSamplingModel::TABULATED;
        }
    }
    else
    {
        const double minimum_energy_ev = std::max(1.0e-3, std::min(cold_temperature_ev, hot_temperature_ev) * 1.0e-2);
        const double characteristic_ev =
            std::max({1.0e-3, params.characteristic_energy, cold_temperature_ev, hot_temperature_ev});
        const double maximum_energy_ev = std::max(10.0 * characteristic_ev, 200.0 * cold_temperature_ev);
        spectrum.energy_grid_ev = buildLogEnergyGrid(minimum_energy_ev, maximum_energy_ev, 96);
        spectrum.differential_number_flux.assign(spectrum.energy_grid_ev.size(), 0.0);

        for (size_t i = 0; i < spectrum.energy_grid_ev.size(); ++i)
        {
            const double energy_ev = spectrum.energy_grid_ev[i];
            double value = 0.0;
            switch (model)
            {
            case SpatialSamplingModel::DOUBLE_MAXWELL:
                value = (1.0 - hot_fraction) * std::exp(-energy_ev / cold_temperature_ev) +
                        hot_fraction * std::exp(-energy_ev / hot_temperature_ev);
                break;
            case SpatialSamplingModel::KAPPA:
            {
                const double kappa = std::max(1.51, params.kappa);
                const double theta =
                    std::max(1.0e-3, params.characteristic_energy / std::max(1.0, kappa - 1.5));
                value = std::pow(1.0 + energy_ev / (kappa * theta), -(kappa + 1.0));
                break;
            }
            case SpatialSamplingModel::POWER_LAW:
                value = std::pow(energy_ev + 1.0e-3, -std::max(1.01, params.power_law_index));
                break;
            case SpatialSamplingModel::TABULATED:
            case SpatialSamplingModel::SINGLE_MAXWELL:
            case SpatialSamplingModel::MAXWELL_BOLTZMANN:
            case SpatialSamplingModel::Q_DISTRIBUTION:
            case SpatialSamplingModel::UNIFORM:
            default:
                value = std::exp(-energy_ev / std::max(1.0e-3, characteristic_ev));
                break;
            }
            spectrum.differential_number_flux[i] = std::max(0.0, value);
        }
    }

    normalizeFluxShape(spectrum.differential_number_flux, spectrum.energy_grid_ev,
                       estimateNumberFluxFromMoments(params, mass_kg));
    return spectrum;
}

/**
 * @brief 批量采样空间位置。
 */
std::vector<Point3D> ParticleSource::sampleSpatialPositions(size_t count) const
{
    std::vector<Point3D> out;
    out.reserve(count);
    for (size_t i = 0; i < count; ++i)
    {
        out.push_back(sampleSpatialPointByModel());
    }
    return out;
}

/**
 * @brief 批量采样速度向量。
 * @details 默认首选方向为 +Z。
 */
std::vector<Vector3D> ParticleSource::sampleVelocityVectors(size_t count) const
{
    std::vector<Vector3D> out;
    out.reserve(count);
    for (size_t i = 0; i < count; ++i)
    {
        out.push_back(sampleVelocityByModel(Vector3D(0.0, 0.0, 1.0)));
    }
    return out;
}

/**
 * @brief 按当前空间模型采样单个空间点。
 * @details 先按模型生成半径，再与均匀方向向量组合得到 3D 位置偏移。
 */
Point3D ParticleSource::sampleSpatialPointByModel() const
{
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const double r_max = clampPositive(sampling_params_.spatial_scale, 0.05);
    const double r_min = std::clamp(sampling_params_.min_spatial_scale, 1e-7, r_max);

    double r = 0.0;
    switch (sampling_model_)
    {
    case SpatialSamplingModel::SINGLE_MAXWELL:
    {
        // 单麦克斯韦分布：使用高斯分布的径向投影
        std::normal_distribution<double> n(0.0, r_max * 0.35);
        const double x = n(rng);
        const double y = n(rng);
        const double z = n(rng);
        r = std::min(r_max, std::sqrt(x * x + y * y + z * z));
        break;
    }
    case SpatialSamplingModel::Q_DISTRIBUTION:
    {
        // q-分布（Tsallis统计）
        const double q = std::clamp(sampling_params_.q_parameter, 1.01, 3.0);
        const double rand_q = std::max(1e-9, 1.0 - u(rng));
        const double exponent = 1.0 / (q - 1.0);
        r = r_max * std::sqrt(q - 1.0) * std::sqrt(std::pow(rand_q, -exponent) - 1.0);
        r = std::min(r, r_max);
        break;
    }
    case SpatialSamplingModel::MAXWELL_BOLTZMANN:
    {
        // 麦克斯韦-玻尔兹曼分布：分三变量麦克斯韦分布
        std::normal_distribution<double> nb(0.0, r_max * 0.4);
        const double x = nb(rng);
        const double y = nb(rng);
        const double z = nb(rng);
        r = std::min(r_max, std::sqrt(x * x + y * y + z * z));
        break;
    }
    case SpatialSamplingModel::DOUBLE_MAXWELL:
    {
        std::normal_distribution<double> n1(0.0, r_max * 0.25);
        std::normal_distribution<double> n2(0.0, r_max * 0.60);
        const bool hot = (u(rng) < std::clamp(sampling_params_.hot_fraction, 0.0, 1.0));
        const double x = hot ? n2(rng) : n1(rng);
        const double y = hot ? n2(rng) : n1(rng);
        const double z = hot ? n2(rng) : n1(rng);
        r = std::min(r_max, std::sqrt(x * x + y * y + z * z));
        break;
    }
    case SpatialSamplingModel::KAPPA:
    {
        const double kappa = std::max(1.51, sampling_params_.kappa);
        const double q = std::max(1e-9, 1.0 - u(rng));
        r = r_max * std::sqrt(kappa * (std::pow(q, -1.0 / (kappa - 1.5)) - 1.0));
        r = std::min(r, r_max);
        break;
    }
    case SpatialSamplingModel::POWER_LAW:
    {
        const double p = std::max(1.01, sampling_params_.power_law_index);
        const double x = u(rng);
        const double c1 = std::pow(r_min, 1.0 - p);
        const double c2 = std::pow(r_max, 1.0 - p);
        r = std::pow(c1 + x * (c2 - c1), 1.0 / (1.0 - p));
        r = std::clamp(r, r_min, r_max);
        break;
    }
    case SpatialSamplingModel::UNIFORM:
    default:
        r = r_max * std::cbrt(u(rng));
        break;
    }

    const Vector3D dir = randomUnitVector();
    return Point3D(dir.x() * r, dir.y() * r, dir.z() * r);
}

/**
 * @brief 按当前模型采样单个速度向量。
 * @param preferred_direction 首选漂移方向；若为零向量则自动随机方向。
 * @details 返回值 = 体漂移速度 + 热运动随机项。
 */
Vector3D ParticleSource::sampleVelocityByModel(const Vector3D& preferred_direction) const
{
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const Vector3D base_dir = preferred_direction.magnitude() > 0.0
                                  ? preferred_direction.normalized()
                                  : randomUnitVector();

    double thermal = sampling_params_.thermal_speed;

    switch (sampling_model_)
    {
    case SpatialSamplingModel::SINGLE_MAXWELL:
        // 单一热速度
        thermal = sampling_params_.thermal_speed;
        break;
    case SpatialSamplingModel::Q_DISTRIBUTION:
    {
        // q-分布的热速度设置
        const double q = std::clamp(sampling_params_.q_parameter, 1.01, 3.0);
        thermal = sampling_params_.thermal_speed * std::sqrt(3.0 * (q - 1.0) / (q - 1.0 + 1.5));
        break;
    }
    case SpatialSamplingModel::MAXWELL_BOLTZMANN:
        // 标准玻尔兹曼分布
        thermal = sampling_params_.thermal_speed;
        break;
    case SpatialSamplingModel::DOUBLE_MAXWELL:
        thermal = (u(rng) < sampling_params_.hot_fraction) ? sampling_params_.hot_thermal_speed
                                                           : sampling_params_.thermal_speed;
        break;
    case SpatialSamplingModel::KAPPA:
    {
        const double kappa = std::max(1.51, sampling_params_.kappa);
        thermal = sampling_params_.thermal_speed * std::sqrt(kappa / (kappa - 1.5));
        break;
    }
    case SpatialSamplingModel::POWER_LAW:
        thermal = sampling_params_.thermal_speed *
                  std::pow(std::max(1e-9, 1.0 - u(rng)),
                           -1.0 / std::max(1.01, sampling_params_.power_law_index));
        break;
    case SpatialSamplingModel::TABULATED:
    {
        const double sampled_energy_ev =
            sampleEnergyFromDiscreteSpectrum(observed_energy_ev_, observed_intensity_);
        const double mass = getParticleTypeMass(particleType_);
        thermal = std::sqrt(
            2.0 * std::max(1.0e-6, sampled_energy_ev) * PhysicsConstants::ElementaryCharge /
            std::max(1.0e-32, mass));
        break;
    }
    case SpatialSamplingModel::UNIFORM:
    default:
        break;
    }

    std::normal_distribution<double> n(0.0, std::max(1.0, thermal));
    Vector3D random_part(n(rng), n(rng), n(rng));
    Vector3D bulk = base_dir * sampling_params_.bulk_speed;
    return bulk + random_part;
}

/**
 * @brief 构造表面粒子源。
 * @throws std::invalid_argument 当 mesh 为空。
 */
SurfaceParticleSource::SurfaceParticleSource(const std::string& name,
                                             const ParticleType& particleType,
                                             const std::shared_ptr<Mesh::SurfaceMesh>& mesh)
    : ParticleSource(name, particleType), mesh_(mesh), emissionElements_(), totalArea_(0.0),
      temperature_(300.0), workFunction_(4.5)
{
    if (!mesh_)
    {
        throw std::invalid_argument("SurfaceParticleSource: 表面网格不能为空");
    }
}

/**
 * @brief 设置表面温度（K）。
 */
void SurfaceParticleSource::setTemperature(double temperature)
{
    if (temperature <= 0.0)
    {
        throw std::invalid_argument("SurfaceParticleSource: 温度必须为正数");
    }
    temperature_ = temperature;
}

/**
 * @brief 设置材料功函数（eV）。
 */
void SurfaceParticleSource::setWorkFunction(double workFunction)
{
    if (workFunction < 0.0)
    {
        throw std::invalid_argument("SurfaceParticleSource: 功函数不能为负数");
    }
    workFunction_ = workFunction;
}

/**
 * @brief 设置发射面片集合并累计总面积。
 */
void SurfaceParticleSource::setEmissionArea(const std::vector<size_t>& elementIndices)
{
    emissionElements_ = elementIndices;
    totalArea_ = 0.0;

    for (size_t elementIndex : emissionElements_)
    {
        totalArea_ += calculateElementArea(elementIndex);
    }
}

/**
 * @brief 在时间步 dt 内发射表面粒子。
 * @details 使用“期望值取整 + 小数概率补偿”的泊松近似策略决定整数粒子数。
 */
size_t SurfaceParticleSource::emitParticles(std::vector<ParticleClass>& particles, double dt)
{
    if (!isActiveAtTime(currentTime_))
    {
        return 0;
    }

    const double expectedParticles = emissionRate_ * dt * timeAcceleration_;
    size_t particlesToEmit = static_cast<size_t>(expectedParticles);

    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double fractionalPart = expectedParticles - static_cast<double>(particlesToEmit);
    if (u(rng) < fractionalPart)
    {
        ++particlesToEmit;
    }

    for (size_t i = 0; i < particlesToEmit; ++i)
    {
        particles.push_back(generateParticle());
        ++totalEmitted_;
    }

    return particlesToEmit;
}

/**
 * @brief 生成一个表面发射粒子。
 * @details 从发射面随机选单元，采样面内位置并沿法向构造麦克斯韦速度。
 */
ParticleClass SurfaceParticleSource::generateParticle()
{
    auto& rng = globalRng();

    size_t elementIndex = 0;
    if (!emissionElements_.empty())
    {
        elementIndex = selectWeightedElementByArea(emissionElements_, mesh_, totalArea_);
    }
    else if (mesh_ && mesh_->getElementCount() > 0)
    {
        std::uniform_int_distribution<size_t> elementDis(0, mesh_->getElementCount() - 1);
        elementIndex = elementDis(rng);
    }

    const Point3D position = generateRandomPointInElement(elementIndex);
    const Vector3D normal = calculateElementNormal(elementIndex);
    const Vector3D velocity = generateMaxwellianVelocity(normal);

    static ParticleId next_id = 0;
    const ParticleId id = next_id++;

    return ParticleFactory::createParticle(particleType_, id, position, velocity, 1.0);
}

/**
 * @brief 生成面外半空间麦克斯韦速度。
 * @details 若速度指向表面内部，则对法向分量镜像反射，保证粒子向外发射。
 */
Vector3D SurfaceParticleSource::generateMaxwellianVelocity(const Vector3D& normal)
{
    auto& rng = globalRng();
    const double mass = getParticleTypeMass(particleType_);
    const double thermalVelocity =
        std::sqrt(PhysicsConstants::BoltzmannConstant * temperature_ / mass);

    std::normal_distribution<double> normalDis(0.0, thermalVelocity);
    Vector3D velocity(normalDis(rng), normalDis(rng), normalDis(rng));

    if (velocity.dot(normal) < 0.0)
    {
        velocity = velocity - normal * (2.0 * velocity.dot(normal));
    }

    return velocity;
}

/**
 * @brief 计算面元面积。
 */
double SurfaceParticleSource::calculateElementArea(size_t elementIndex) const
{
    if (!mesh_)
    {
        return 1.0e-6;
    }

    const double area = mesh_->getElementArea(elementIndex);
    return isFinitePositive(area) ? area : 1.0e-6;
}

/**
 * @brief 在指定表面单元内生成随机点。
 */
Point3D SurfaceParticleSource::generateRandomPointInElement(size_t elementIndex) const
{
    if (!mesh_)
    {
        return Point3D(0.0, 0.0, 0.0);
    }

    return mesh_->getRandomPointInElement(elementIndex);
}

/**
 * @brief 计算面元法向。
 */
Vector3D SurfaceParticleSource::calculateElementNormal(size_t elementIndex) const
{
    if (!mesh_)
    {
        return Vector3D(0.0, 0.0, 1.0);
    }

    const Vector3D n = mesh_->getElementNormal(elementIndex);
    return n.magnitude() > 0.0 ? n.normalized() : Vector3D(0.0, 0.0, 1.0);
}

/**
 * @brief 构造体积粒子源。
 * @throws std::invalid_argument 当 mesh 为空。
 */
VolumeParticleSource::VolumeParticleSource(const std::string& name,
                                           const ParticleType& particleType,
                                           const std::shared_ptr<Mesh::VolumeMesh>& mesh)
    : ParticleSource(name, particleType), mesh_(mesh), emissionElements_(), totalVolume_(0.0),
      density_(1.0e15), temperature_(300.0)
{
    if (!mesh_)
    {
        throw std::invalid_argument("VolumeParticleSource: 体积网格不能为空");
    }
}

/**
 * @brief 设置体源密度。
 */
void VolumeParticleSource::setDensity(double density)
{
    if (density < 0.0)
    {
        throw std::invalid_argument("VolumeParticleSource: 密度不能为负数");
    }
    density_ = density;
}

/**
 * @brief 设置体源温度（K）。
 */
void VolumeParticleSource::setTemperature(double temperature)
{
    if (temperature <= 0.0)
    {
        throw std::invalid_argument("VolumeParticleSource: 温度必须为正数");
    }
    temperature_ = temperature;
}

/**
 * @brief 设置发射体元并累计总体积。
 */
void VolumeParticleSource::setEmissionRegion(const std::vector<size_t>& elementIndices)
{
    emissionElements_ = elementIndices;
    totalVolume_ = 0.0;

    for (size_t elementIndex : emissionElements_)
    {
        totalVolume_ += calculateElementVolume(elementIndex);
    }
}

/**
 * @brief 在时间步 dt 内发射体积粒子。
 */
size_t VolumeParticleSource::emitParticles(std::vector<ParticleClass>& particles, double dt)
{
    if (!isActiveAtTime(currentTime_))
    {
        return 0;
    }

    const double expectedParticles = emissionRate_ * dt * timeAcceleration_;
    size_t particlesToEmit = static_cast<size_t>(expectedParticles);

    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double fractionalPart = expectedParticles - static_cast<double>(particlesToEmit);
    if (u(rng) < fractionalPart)
    {
        ++particlesToEmit;
    }

    for (size_t i = 0; i < particlesToEmit; ++i)
    {
        particles.push_back(generateParticle());
        ++totalEmitted_;
    }

    return particlesToEmit;
}

/**
 * @brief 生成一个体发射粒子。
 * @details 在体元内取基准点，再叠加模型空间偏移，速度由模型速度采样得到。
 */
ParticleClass VolumeParticleSource::generateParticle()
{
    auto& rng = globalRng();

    size_t elementIndex = 0;
    if (!emissionElements_.empty())
    {
        elementIndex = selectWeightedElementByVolume(emissionElements_, mesh_, totalVolume_);
    }
    else if (mesh_ && mesh_->getElementCount() > 0)
    {
        std::uniform_int_distribution<size_t> elementDis(0, mesh_->getElementCount() - 1);
        elementIndex = elementDis(rng);
    }

    const Point3D base_position = generateRandomPointInVolumeElement(elementIndex);
    const Point3D offset = sampleSpatialPointByModel();
    const Point3D position(base_position.x() + offset.x(), base_position.y() + offset.y(),
                           base_position.z() + offset.z());

    const Vector3D velocity = sampleVelocityByModel(randomUnitVector());

    static ParticleId next_id = 1000;
    const ParticleId id = next_id++;

    return ParticleFactory::createParticle(particleType_, id, position, velocity, 1.0);
}

/**
 * @brief 生成麦克斯韦热速度。
 */
Vector3D VolumeParticleSource::generateMaxwellianVelocity()
{
    auto& rng = globalRng();
    const double mass = getParticleTypeMass(particleType_);
    const double thermalVelocity =
        std::sqrt(PhysicsConstants::BoltzmannConstant * temperature_ / mass);

    std::normal_distribution<double> normalDis(0.0, thermalVelocity);
    return Vector3D(normalDis(rng), normalDis(rng), normalDis(rng));
}

/**
 * @brief 计算体元体积。
 */
double VolumeParticleSource::calculateElementVolume(size_t elementIndex) const
{
    if (!mesh_)
    {
        return 1.0e-9;
    }

    const double volume = mesh_->getElementVolume(elementIndex);
    return isFinitePositive(volume) ? volume : 1.0e-9;
}

/**
 * @brief 在指定体单元内生成随机点。
 */
Point3D VolumeParticleSource::generateRandomPointInVolumeElement(size_t elementIndex) const
{
    if (!mesh_)
    {
        return Point3D(0.0, 0.0, 0.0);
    }

    return mesh_->getRandomPointInElement(elementIndex);
}

/**
 * @brief 构造束流粒子源。
 */
BeamParticleSource::BeamParticleSource(const std::string& name, const ParticleType& particleType,
                                       const Point3D& origin, const Vector3D& direction)
    : ParticleSource(name, particleType), origin_(origin), direction_(direction.normalized()),
      beamEnergy_(1000.0), beamRadius_(0.01), divergenceAngle_(0.0)
{
}

/**
 * @brief 设置束流能量（eV）。
 */
void BeamParticleSource::setBeamEnergy(double energy)
{
    if (energy <= 0.0)
    {
        throw std::invalid_argument("BeamParticleSource: 束流能量必须为正数");
    }
    beamEnergy_ = energy;
}

/**
 * @brief 设置束斑半径（m）。
 */
void BeamParticleSource::setBeamRadius(double radius)
{
    if (radius <= 0.0)
    {
        throw std::invalid_argument("BeamParticleSource: 束流半径必须为正数");
    }
    beamRadius_ = radius;
}

/**
 * @brief 设置束流发散半角（rad）。
 */
void BeamParticleSource::setDivergenceAngle(double angle)
{
    if (angle < 0.0)
    {
        throw std::invalid_argument("BeamParticleSource: 发散角不能为负数");
    }
    divergenceAngle_ = angle;
}

/**
 * @brief 在时间步 dt 内发射束流粒子。
 */
size_t BeamParticleSource::emitParticles(std::vector<ParticleClass>& particles, double dt)
{
    if (!isActiveAtTime(currentTime_))
    {
        return 0;
    }

    const double expectedParticles = emissionRate_ * dt * timeAcceleration_;
    size_t particlesToEmit = static_cast<size_t>(expectedParticles);

    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double fractionalPart = expectedParticles - static_cast<double>(particlesToEmit);
    if (u(rng) < fractionalPart)
    {
        ++particlesToEmit;
    }

    for (size_t i = 0; i < particlesToEmit; ++i)
    {
        particles.push_back(generateParticle());
        ++totalEmitted_;
    }

    return particlesToEmit;
}

/**
 * @brief 生成一个束流粒子。
 * @details 在束斑圆盘内均匀采样起点，方向按发散角扰动，再由能量换算速度模长。
 */
ParticleClass BeamParticleSource::generateParticle()
{
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const double r = beamRadius_ * std::sqrt(u(rng));
    const double theta = 2.0 * MathConstants::Pi * u(rng);

    Vector3D perpendicular1;
    if (std::abs(direction_.x()) < 0.9)
    {
        perpendicular1 = Vector3D(1.0, 0.0, 0.0).cross(direction_).normalized();
    }
    else
    {
        perpendicular1 = Vector3D(0.0, 1.0, 0.0).cross(direction_).normalized();
    }
    const Vector3D perpendicular2 = direction_.cross(perpendicular1).normalized();

    const Point3D position =
        origin_ + perpendicular1 * (r * std::cos(theta)) + perpendicular2 * (r * std::sin(theta));

    Vector3D beam_direction = direction_;
    if (divergenceAngle_ > 0.0)
    {
        const double divergence = divergenceAngle_ * u(rng);
        const double azimuth = 2.0 * MathConstants::Pi * u(rng);

        beam_direction = (perpendicular1 * (std::sin(divergence) * std::cos(azimuth)) +
                          perpendicular2 * (std::sin(divergence) * std::sin(azimuth)) +
                          direction_ * std::cos(divergence))
                             .normalized();
    }

    const double mass = getParticleTypeMass(particleType_);
    const double speed = std::sqrt(2.0 * beamEnergy_ * PhysicsConstants::ElementaryCharge / mass);
    const Vector3D velocity = beam_direction * speed;

    static ParticleId next_id = 2000;
    const ParticleId id = next_id++;

    return ParticleFactory::createParticle(particleType_, id, position, velocity, 1.0);
}

/**
 * @brief 构造粒子源管理器。
 */
ParticleSourceManager::ParticleSourceManager() = default;

/**
 * @brief 析构粒子源管理器。
 */
ParticleSourceManager::~ParticleSourceManager() = default;

/**
 * @brief 添加粒子源并检查重名冲突。
 */
void ParticleSourceManager::addSource(std::shared_ptr<ParticleSource> source)
{
    if (!source)
    {
        throw std::invalid_argument("ParticleSourceManager: 粒子源不能为空");
    }

    auto it = std::find_if(sources_.begin(), sources_.end(),
                           [&source](const std::shared_ptr<ParticleSource>& current)
                           { return current->getName() == source->getName(); });

    if (it != sources_.end())
    {
        throw std::invalid_argument("ParticleSourceManager: 粒子源名称已存在: " +
                                    source->getName());
    }

    sources_.push_back(std::move(source));
}

/**
 * @brief 按名称移除粒子源；不存在则忽略。
 */
void ParticleSourceManager::removeSource(const std::string& name)
{
    auto it = std::find_if(sources_.begin(), sources_.end(),
                           [&name](const std::shared_ptr<ParticleSource>& source)
                           { return source->getName() == name; });
    if (it != sources_.end())
    {
        sources_.erase(it);
    }
}

/**
 * @brief 按名称获取粒子源。
 * @return 找到返回对象，否则返回 nullptr。
 */
std::shared_ptr<ParticleSource> ParticleSourceManager::getSource(const std::string& name)
{
    auto it = std::find_if(sources_.begin(), sources_.end(),
                           [&name](const std::shared_ptr<ParticleSource>& source)
                           { return source->getName() == name; });
    return it != sources_.end() ? *it : nullptr;
}

/**
 * @brief 触发所有已启用粒子源发射，并累计本步发射总数。
 */
size_t ParticleSourceManager::emitAllParticles(std::vector<ParticleClass>& particles, double dt)
{
    size_t totalEmitted = 0;
    for (auto& source : sources_)
    {
        if (source->isEnabled())
        {
            totalEmitted += source->emitParticles(particles, dt);
        }
    }
    return totalEmitted;
}

/**
 * @brief 推进所有粒子源内部时钟。
 */
void ParticleSourceManager::updateTime(double dt)
{
    for (auto& source : sources_)
    {
        source->updateTime(dt);
    }
}

/**
 * @brief 按名称启用或禁用指定粒子源。
 */
void ParticleSourceManager::setSourceEnabled(const std::string& name, bool enabled)
{
    auto source = getSource(name);
    if (source)
    {
        source->setEnabled(enabled);
    }
}

/**
 * @brief 获取所有粒子源名称列表。
 */
std::vector<std::string> ParticleSourceManager::getSourceNames() const
{
    std::vector<std::string> names;
    names.reserve(sources_.size());
    for (const auto& source : sources_)
    {
        names.push_back(source->getName());
    }
    return names;
}

/**
 * @brief 获取管理器中粒子源数量。
 */
size_t ParticleSourceManager::getSourceCount() const
{
    return sources_.size();
}

/**
 * @brief 清空所有粒子源。
 */
void ParticleSourceManager::clearSources()
{
    sources_.clear();
}

/**
 * @brief 统计所有粒子源累计发射总数。
 */
size_t ParticleSourceManager::getTotalEmitted() const
{
    size_t total = 0;
    for (const auto& source : sources_)
    {
        total += source->getTotalEmitted();
    }
    return total;
}

} // namespace Particle
} // namespace SCDAT
