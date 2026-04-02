/**
 * @file ParticleSource.h
 * @brief Particle source definitions for SCDAT.
 */

#pragma once

#include "../Geometry/include/Point3D.h"
#include "../Geometry/include/Vector3D.h"
#include "ParticleDefinitions.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace SCDAT
{
namespace Mesh
{
class Mesh;
class VolMesh;
using SurfaceMesh = VolMesh;
using VolumeMesh = VolMesh;
} // namespace Mesh

namespace Particle
{

using Point3D = SCDAT::Geometry::Point3D;
using Vector3D = SCDAT::Geometry::Vector3D;
using ParticleClass = SCDAT::Particle::Particle;
using ParticleType = SCDAT::Particle::ParticleType;
using ParticlePtr = SCDAT::Particle::ParticlePtr;
using ParticleId = SCDAT::Particle::ParticleId;

enum class SpatialSamplingModel
{
    UNIFORM = 0,
    DOUBLE_MAXWELL = 1,
    KAPPA = 2,
    POWER_LAW = 3,
    SINGLE_MAXWELL = 4,
    Q_DISTRIBUTION = 5,
    MAXWELL_BOLTZMANN = 6,
    TABULATED = 7
};

struct SamplingParameters
{
    double density = 1.0e15;
    double bulk_speed = 2.0e5;
    double thermal_speed = 5.0e4;
    double hot_thermal_speed = 2.0e5;
    double hot_fraction = 0.2;
    double kappa = 3.5;
    double power_law_index = 2.6;
    double spatial_scale = 0.05;
    double min_spatial_scale = 1.0e-4;
    double q_parameter = 1.5;
    double characteristic_energy = 10.0;
};

struct FitSummary
{
    bool success = false;
    SpatialSamplingModel selected_model = SpatialSamplingModel::UNIFORM;
    double score_double_maxwell = 0.0;
    double score_kappa = 0.0;
    double score_power_law = 0.0;
    std::string message;
};

enum class SpectrumUsage
{
    SamplingOnly,
    CurrentBalance,
    PicCalibration
};

struct SpectrumPopulation
{
    double density_m3 = 0.0;
    double temperature_ev = 0.0;
    double drift_speed_m_per_s = 0.0;
    double mass_amu = 1.0;
    double charge_number = 1.0;
    double weight = 0.0;
};

struct ResolvedSpectrum
{
    SpatialSamplingModel model = SpatialSamplingModel::UNIFORM;
    std::vector<SpectrumPopulation> populations;
    std::vector<double> energy_grid_ev;
    std::vector<double> differential_number_flux;
    FitSummary fit_summary;
};

class ParticleSource
{
  public:
    ParticleSource(const std::string& name, const ParticleType& particleType);
    virtual ~ParticleSource();

    void setEnabled(bool enabled);
    bool isEnabled() const;

    void setEmissionRate(double rate);
    double getEmissionRate() const;

    void setTimeWindow(double startTime, double stopTime = -1.0);
    bool isActiveAtTime(double time) const;
    void updateTime(double dt);

    size_t getTotalEmitted() const;
    const std::string& getName() const;
    const ParticleType& getParticleType() const;

    virtual size_t emitParticles(std::vector<ParticleClass>& particles, double dt) = 0;

    void setSamplingModel(SpatialSamplingModel model);
    SpatialSamplingModel getSamplingModel() const;

    void setSamplingParameters(const SamplingParameters& params);
    const SamplingParameters& getSamplingParameters() const;

    void setObservedSpectrum(const std::vector<double>& energy_ev,
                             const std::vector<double>& intensity);
    void setDetectorChargeProfile(const std::vector<double>& radius,
                                  const std::vector<double>& charge_density);

    FitSummary fitSamplingModelFromObservations();
    ResolvedSpectrum resolveSpectrum(SpectrumUsage usage = SpectrumUsage::SamplingOnly) const;
    static ResolvedSpectrum buildResolvedSpectrum(
        const ParticleType& particleType, SpatialSamplingModel model,
        const SamplingParameters& params, SpectrumUsage usage = SpectrumUsage::SamplingOnly,
        const std::vector<double>& observed_energy_ev = {},
        const std::vector<double>& observed_intensity = {},
        const FitSummary& fit_summary = {});

    std::vector<Point3D> sampleSpatialPositions(size_t count) const;
    std::vector<Vector3D> sampleVelocityVectors(size_t count) const;

  protected:
    virtual ParticleClass generateParticle() = 0;

    Point3D sampleSpatialPointByModel() const;
    Vector3D sampleVelocityByModel(const Vector3D& preferred_direction) const;

    std::string name_;
    ParticleType particleType_;
    bool enabled_;
    double emissionRate_;
    size_t totalEmitted_;
    double currentTime_;
    double startTime_;
    double stopTime_;
    double timeAcceleration_;

    SpatialSamplingModel sampling_model_;
    SamplingParameters sampling_params_;

    std::vector<double> observed_energy_ev_;
    std::vector<double> observed_intensity_;
    std::vector<double> detector_radius_;
    std::vector<double> detector_charge_density_;
    FitSummary last_fit_summary_;
};

class SurfaceParticleSource : public ParticleSource
{
  public:
    SurfaceParticleSource(const std::string& name, const ParticleType& particleType,
                          const std::shared_ptr<Mesh::SurfaceMesh>& mesh);

    void setTemperature(double temperature);
    void setWorkFunction(double workFunction);
    void setEmissionArea(const std::vector<size_t>& elementIndices);

    size_t emitParticles(std::vector<ParticleClass>& particles, double dt) override;

  protected:
    ParticleClass generateParticle() override;
    Vector3D generateMaxwellianVelocity(const Vector3D& normal);

  private:
    double calculateElementArea(size_t elementIndex) const;
    Point3D generateRandomPointInElement(size_t elementIndex) const;
    Vector3D calculateElementNormal(size_t elementIndex) const;

    std::shared_ptr<Mesh::SurfaceMesh> mesh_;
    std::vector<size_t> emissionElements_;
    double totalArea_;
    double temperature_;
    double workFunction_;
};

class VolumeParticleSource : public ParticleSource
{
  public:
    VolumeParticleSource(const std::string& name, const ParticleType& particleType,
                         const std::shared_ptr<Mesh::VolumeMesh>& mesh);

    void setDensity(double density);
    void setTemperature(double temperature);
    void setEmissionRegion(const std::vector<size_t>& elementIndices);

    size_t emitParticles(std::vector<ParticleClass>& particles, double dt) override;

  protected:
    ParticleClass generateParticle() override;
    Vector3D generateMaxwellianVelocity();

  private:
    double calculateElementVolume(size_t elementIndex) const;
    Point3D generateRandomPointInVolumeElement(size_t elementIndex) const;

    std::shared_ptr<Mesh::VolumeMesh> mesh_;
    std::vector<size_t> emissionElements_;
    double totalVolume_;
    double density_;
    double temperature_;
};

class BeamParticleSource : public ParticleSource
{
  public:
    BeamParticleSource(const std::string& name, const ParticleType& particleType,
                       const Point3D& origin, const Vector3D& direction);
    void setBeamEnergy(double energy);
    void setBeamRadius(double radius);
    void setDivergenceAngle(double angle);

    size_t emitParticles(std::vector<ParticleClass>& particles, double dt) override;

  protected:
    ParticleClass generateParticle() override;

    Point3D origin_;
    Vector3D direction_;
    double beamEnergy_;
    double beamRadius_;
    double divergenceAngle_;
};

class ParticleSourceManager
{
  public:
    ParticleSourceManager();
    ~ParticleSourceManager();

    void addSource(std::shared_ptr<ParticleSource> source);
    void removeSource(const std::string& name);
    std::shared_ptr<ParticleSource> getSource(const std::string& name);

    size_t emitAllParticles(std::vector<ParticleClass>& particles, double dt);
    void updateTime(double dt);
    void setSourceEnabled(const std::string& name, bool enabled);

    std::vector<std::string> getSourceNames() const;
    size_t getSourceCount() const;
    void clearSources();
    size_t getTotalEmitted() const;

  private:
    std::vector<std::shared_ptr<ParticleSource>> sources_;
};

} // namespace Particle
} // namespace SCDAT
