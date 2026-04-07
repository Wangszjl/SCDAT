#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace SCDAT
{
namespace Toolkit
{
namespace Radiation
{

enum class ParticleSpecies;
enum class RadiationPhysicsList;

struct RadiationStoppingResolution
{
    bool valid = false;
    bool cache_hit = false;
    bool rollback_applied = false;

    std::string requested_version;
    std::string resolved_version;
    std::string rollback_reason;
    std::string validation_error;

    std::string dataset_id;
    std::string dataset_source_uri;
    std::string dataset_created_utc;
    std::string material_source_hash;
    std::string physics_source_hash;

    std::string material_name;

    double z_over_a = 0.0;
    double mean_excitation_ev = 0.0;
    double stopping_scale = 1.0;
    double scattering_scale = 1.0;
};

class RadiationStoppingDataGovernance
{
  public:
    RadiationStoppingDataGovernance();

    RadiationStoppingResolution resolve(ParticleSpecies species,
                                        RadiationPhysicsList physics_list,
                                        const std::string& material_name,
                                        const std::string& requested_version,
                                        bool allow_rollback) const;

    void clearCache() const;
    std::size_t cacheSize() const;

  private:
    struct MaterialRecord
    {
        std::string version;
        std::string material_name;
        std::string source_hash;
        std::string source_uri;
        std::string created_utc;
        double z_over_a = 0.0;
        double mean_excitation_ev = 0.0;
    };

    struct PhysicsRecord
    {
        std::string version;
        RadiationPhysicsList physics_list;
        std::string source_hash;
        std::string source_uri;
        std::string created_utc;
        double stopping_scale = 1.0;
        double scattering_scale = 1.0;
    };

    std::vector<std::string> version_precedence_;
    std::vector<MaterialRecord> material_records_;
    std::vector<PhysicsRecord> physics_records_;

    mutable std::unordered_map<std::string, RadiationStoppingResolution> cache_;
};

} // namespace Radiation
} // namespace Toolkit
} // namespace SCDAT
