/**
 * @file ParticleSource.cpp
 * @brief 粒子源实现（SCDAT）
 */

#include "../include/ParticleSource.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>

namespace SCDAT
{
namespace Particle
{

namespace
{

constexpr double kBoltzmann = 1.380649e-23;
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kPi = 3.14159265358979323846;

std::mt19937& globalRng()
{
    thread_local std::mt19937 rng(std::random_device{}());
    return rng;
}

double clampPositive(double value, double fallback)
{
    return value > 0.0 ? value : fallback;
}

double getParticleTypeMass(ParticleType type)
{
    switch (type)
    {
    case ParticleType::ELECTRON:
    case ParticleType::PHOTOELECTRON:
    case ParticleType::SECONDARY_ELECTRON:
    case ParticleType::THERMAL_ELECTRON:
    case ParticleType::FIELD_EMISSION_ELECTRON:
        return 9.10938356e-31;
    case ParticleType::ION:
    case ParticleType::POSITIVE_ION:
    case ParticleType::NEGATIVE_ION:
    case ParticleType::MOLECULAR_ION:
    case ParticleType::CLUSTER_ION:
        return 1.6726219e-27;
    default:
        return 1.6726219e-27;
    }
}

Vector3D randomUnitVector()
{
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const double mu = 2.0 * u(rng) - 1.0;
    const double phi = 2.0 * kPi * u(rng);
    const double r_xy = std::sqrt(std::max(0.0, 1.0 - mu * mu));

    return Vector3D(r_xy * std::cos(phi), r_xy * std::sin(phi), mu);
}

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

} // namespace

ParticleSource::ParticleSource(const std::string& name, const ParticleType& particleType)
    : name_(name), particleType_(particleType), enabled_(true), emissionRate_(0.0),
      totalEmitted_(0), currentTime_(0.0), startTime_(0.0), stopTime_(-1.0), timeAcceleration_(1.0),
      sampling_model_(SpatialSamplingModel::UNIFORM), sampling_params_{}
{
}

ParticleSource::~ParticleSource() = default;

void ParticleSource::setEnabled(bool enabled)
{
    enabled_ = enabled;
}

bool ParticleSource::isEnabled() const
{
    return enabled_;
}

void ParticleSource::setEmissionRate(double rate)
{
    if (rate < 0.0)
    {
        throw std::invalid_argument("ParticleSource: 发射率不能为负数");
    }
    emissionRate_ = rate;
}

double ParticleSource::getEmissionRate() const
{
    return emissionRate_;
}

void ParticleSource::setTimeWindow(double startTime, double stopTime)
{
    if (stopTime >= 0.0 && stopTime <= startTime)
    {
        throw std::invalid_argument("ParticleSource: 停止时间必须大于开始时间");
    }
    startTime_ = startTime;
    stopTime_ = stopTime;
}

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

void ParticleSource::updateTime(double dt)
{
    currentTime_ += dt;
}

size_t ParticleSource::getTotalEmitted() const
{
    return totalEmitted_;
}

const std::string& ParticleSource::getName() const
{
    return name_;
}

const ParticleType& ParticleSource::getParticleType() const
{
    return particleType_;
}

void ParticleSource::setSamplingModel(SpatialSamplingModel model)
{
    sampling_model_ = model;
}

SpatialSamplingModel ParticleSource::getSamplingModel() const
{
    return sampling_model_;
}

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
    sampling_params_.characteristic_energy = clampPositive(sampling_params_.characteristic_energy, 10.0);
}

const SamplingParameters& ParticleSource::getSamplingParameters() const
{
    return sampling_params_;
}

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

FitSummary ParticleSource::fitSamplingModelFromObservations()
{
    FitSummary summary;

    if (observed_energy_ev_.size() < 5 || observed_intensity_.size() < 5)
    {
        summary.message = "观测能谱数据不足，保留当前分布";
        return summary;
    }

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

    double e_split = mean_e;
    double cold_sum = 0.0;
    double hot_sum = 0.0;
    double cold_w = 0.0;
    double hot_w = 0.0;

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

    const double t1 =
        cold_w > 0.0 ? std::max(1e-6, cold_sum / cold_w) : std::max(1e-6, mean_e * 0.5);
    const double t2 = hot_w > 0.0 ? std::max(t1, hot_sum / hot_w) : std::max(t1, mean_e * 1.5);
    const double hot_fraction = sum_w > 0.0 ? std::clamp(hot_w / sum_w, 0.02, 0.98) : 0.2;

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
    const double cold_j = t1 * kElementaryCharge;
    const double hot_j = t2 * kElementaryCharge;
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
    return summary;
}

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
    case SpatialSamplingModel::UNIFORM:
    default:
        break;
    }

    std::normal_distribution<double> n(0.0, std::max(1.0, thermal));
    Vector3D random_part(n(rng), n(rng), n(rng));
    Vector3D bulk = base_dir * sampling_params_.bulk_speed;
    return bulk + random_part;
}

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

void SurfaceParticleSource::setTemperature(double temperature)
{
    if (temperature <= 0.0)
    {
        throw std::invalid_argument("SurfaceParticleSource: 温度必须为正数");
    }
    temperature_ = temperature;
}

void SurfaceParticleSource::setWorkFunction(double workFunction)
{
    if (workFunction < 0.0)
    {
        throw std::invalid_argument("SurfaceParticleSource: 功函数不能为负数");
    }
    workFunction_ = workFunction;
}

void SurfaceParticleSource::setEmissionArea(const std::vector<size_t>& elementIndices)
{
    emissionElements_ = elementIndices;
    totalArea_ = 0.0;

    for (size_t elementIndex : emissionElements_)
    {
        totalArea_ += calculateElementArea(elementIndex);
    }
}

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

ParticleClass SurfaceParticleSource::generateParticle()
{
    auto& rng = globalRng();

    size_t elementIndex = 0;
    if (!emissionElements_.empty())
    {
        std::uniform_int_distribution<size_t> elementDis(0, emissionElements_.size() - 1);
        elementIndex = emissionElements_[elementDis(rng)];
    }

    const Point3D position = generateRandomPointInElement(elementIndex);
    const Vector3D normal = calculateElementNormal(elementIndex);
    const Vector3D velocity = generateMaxwellianVelocity(normal);

    static ParticleId next_id = 0;
    const ParticleId id = next_id++;

    return ParticleFactory::createParticle(particleType_, id, position, velocity, 1.0);
}

Vector3D SurfaceParticleSource::generateMaxwellianVelocity(const Vector3D& normal)
{
    auto& rng = globalRng();
    const double mass = getParticleTypeMass(particleType_);
    const double thermalVelocity = std::sqrt(kBoltzmann * temperature_ / mass);

    std::normal_distribution<double> normalDis(0.0, thermalVelocity);
    Vector3D velocity(normalDis(rng), normalDis(rng), normalDis(rng));

    if (velocity.dot(normal) < 0.0)
    {
        velocity = velocity - normal * (2.0 * velocity.dot(normal));
    }

    return velocity;
}

double SurfaceParticleSource::calculateElementArea(size_t elementIndex) const
{
    (void)elementIndex;
    return 1.0e-6;
}

Point3D SurfaceParticleSource::generateRandomPointInElement(size_t elementIndex) const
{
    (void)elementIndex;
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    double a = u(rng);
    double b = u(rng);
    if (a + b > 1.0)
    {
        a = 1.0 - a;
        b = 1.0 - b;
    }

    return Point3D(a, b, 0.0);
}

Vector3D SurfaceParticleSource::calculateElementNormal(size_t elementIndex) const
{
    (void)elementIndex;
    return Vector3D(0.0, 0.0, 1.0);
}

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

void VolumeParticleSource::setDensity(double density)
{
    if (density < 0.0)
    {
        throw std::invalid_argument("VolumeParticleSource: 密度不能为负数");
    }
    density_ = density;
}

void VolumeParticleSource::setTemperature(double temperature)
{
    if (temperature <= 0.0)
    {
        throw std::invalid_argument("VolumeParticleSource: 温度必须为正数");
    }
    temperature_ = temperature;
}

void VolumeParticleSource::setEmissionRegion(const std::vector<size_t>& elementIndices)
{
    emissionElements_ = elementIndices;
    totalVolume_ = 0.0;

    for (size_t elementIndex : emissionElements_)
    {
        totalVolume_ += calculateElementVolume(elementIndex);
    }
}

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

ParticleClass VolumeParticleSource::generateParticle()
{
    auto& rng = globalRng();

    size_t elementIndex = 0;
    if (!emissionElements_.empty())
    {
        std::uniform_int_distribution<size_t> elementDis(0, emissionElements_.size() - 1);
        elementIndex = emissionElements_[elementDis(rng)];
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

Vector3D VolumeParticleSource::generateMaxwellianVelocity()
{
    auto& rng = globalRng();
    const double mass = getParticleTypeMass(particleType_);
    const double thermalVelocity = std::sqrt(kBoltzmann * temperature_ / mass);

    std::normal_distribution<double> normalDis(0.0, thermalVelocity);
    return Vector3D(normalDis(rng), normalDis(rng), normalDis(rng));
}

double VolumeParticleSource::calculateElementVolume(size_t elementIndex) const
{
    (void)elementIndex;
    return 1.0e-9;
}

Point3D VolumeParticleSource::generateRandomPointInVolumeElement(size_t elementIndex) const
{
    (void)elementIndex;
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    double a = u(rng);
    double b = u(rng);
    double c = u(rng);

    if (a + b + c > 1.0)
    {
        if (a + b > 1.0)
        {
            const double ta = a;
            a = 1.0 - b - c;
            b = 1.0 - ta - c;
        }
        else if (b + c > 1.0)
        {
            const double tc = c;
            c = 1.0 - a - b;
            b = 1.0 - a - tc;
        }
        else
        {
            const double ta = a;
            a = 1.0 - b - c;
            c = ta + b + c - 1.0;
        }
    }

    return Point3D(a, b, c);
}

BeamParticleSource::BeamParticleSource(const std::string& name, const ParticleType& particleType,
                                       const Point3D& origin, const Vector3D& direction)
    : ParticleSource(name, particleType), origin_(origin), direction_(direction.normalized()),
      beamEnergy_(1000.0), beamRadius_(0.01), divergenceAngle_(0.0)
{
}

void BeamParticleSource::setBeamEnergy(double energy)
{
    if (energy <= 0.0)
    {
        throw std::invalid_argument("BeamParticleSource: 束流能量必须为正数");
    }
    beamEnergy_ = energy;
}

void BeamParticleSource::setBeamRadius(double radius)
{
    if (radius <= 0.0)
    {
        throw std::invalid_argument("BeamParticleSource: 束流半径必须为正数");
    }
    beamRadius_ = radius;
}

void BeamParticleSource::setDivergenceAngle(double angle)
{
    if (angle < 0.0)
    {
        throw std::invalid_argument("BeamParticleSource: 发散角不能为负数");
    }
    divergenceAngle_ = angle;
}

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

ParticleClass BeamParticleSource::generateParticle()
{
    auto& rng = globalRng();
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const double r = beamRadius_ * std::sqrt(u(rng));
    const double theta = 2.0 * kPi * u(rng);

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
        const double azimuth = 2.0 * kPi * u(rng);

        beam_direction = (perpendicular1 * (std::sin(divergence) * std::cos(azimuth)) +
                          perpendicular2 * (std::sin(divergence) * std::sin(azimuth)) +
                          direction_ * std::cos(divergence))
                             .normalized();
    }

    const double mass = getParticleTypeMass(particleType_);
    const double speed = std::sqrt(2.0 * beamEnergy_ * kElementaryCharge / mass);
    const Vector3D velocity = beam_direction * speed;

    static ParticleId next_id = 2000;
    const ParticleId id = next_id++;

    return ParticleFactory::createParticle(particleType_, id, position, velocity, 1.0);
}

ParticleSourceManager::ParticleSourceManager() = default;
ParticleSourceManager::~ParticleSourceManager() = default;

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

std::shared_ptr<ParticleSource> ParticleSourceManager::getSource(const std::string& name)
{
    auto it = std::find_if(sources_.begin(), sources_.end(),
                           [&name](const std::shared_ptr<ParticleSource>& source)
                           { return source->getName() == name; });
    return it != sources_.end() ? *it : nullptr;
}

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

void ParticleSourceManager::updateTime(double dt)
{
    for (auto& source : sources_)
    {
        source->updateTime(dt);
    }
}

void ParticleSourceManager::setSourceEnabled(const std::string& name, bool enabled)
{
    auto source = getSource(name);
    if (source)
    {
        source->setEnabled(enabled);
    }
}

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

size_t ParticleSourceManager::getSourceCount() const
{
    return sources_.size();
}

void ParticleSourceManager::clearSources()
{
    sources_.clear();
}

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
