#include "RadiationStoppingDataGovernance.h"

#include "RadiationDoseAlgorithm.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <optional>
#include <set>
#include <sstream>

namespace SCDAT
{
namespace Toolkit
{
namespace Radiation
{
namespace
{

std::string normalizeToken(std::string text)
{
    std::string token;
    token.reserve(text.size());
    for (const unsigned char ch : text)
    {
        if (std::isalnum(ch))
        {
            token.push_back(static_cast<char>(std::tolower(ch)));
        }
    }
    return token;
}

bool isFinitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

int versionRank(const std::vector<std::string>& precedence, const std::string& version)
{
    const auto it = std::find(precedence.begin(), precedence.end(), version);
    if (it == precedence.end())
    {
        return -1;
    }
    return static_cast<int>(std::distance(precedence.begin(), it));
}

std::optional<std::string> selectLatestVersion(const std::set<std::string>& versions,
                                               const std::vector<std::string>& precedence)
{
    if (versions.empty())
    {
        return std::nullopt;
    }

    std::optional<std::string> selected;
    int selected_rank = -1;
    for (const auto& version : versions)
    {
        const int rank = versionRank(precedence, version);
        if (!selected)
        {
            selected = version;
            selected_rank = rank;
            continue;
        }

        if (rank > selected_rank)
        {
            selected = version;
            selected_rank = rank;
            continue;
        }

        if (rank == selected_rank && version > *selected)
        {
            selected = version;
        }
    }

    return selected;
}

std::optional<std::string> selectRollbackVersion(const std::set<std::string>& versions,
                                                 const std::vector<std::string>& precedence,
                                                 const std::string& requested)
{
    if (versions.empty())
    {
        return std::nullopt;
    }

    const int requested_rank = versionRank(precedence, requested);
    if (requested_rank < 0)
    {
        return selectLatestVersion(versions, precedence);
    }

    std::optional<std::string> selected;
    int selected_rank = -1;
    for (const auto& version : versions)
    {
        const int rank = versionRank(precedence, version);
        if (rank < 0 || rank > requested_rank)
        {
            continue;
        }

        if (!selected || rank > selected_rank)
        {
            selected = version;
            selected_rank = rank;
        }
    }

    if (selected)
    {
        return selected;
    }

    return selectLatestVersion(versions, precedence);
}

std::string makeCacheKey(ParticleSpecies species, RadiationPhysicsList physics_list,
                         const std::string& material_name,
                         const std::string& requested_version,
                         bool allow_rollback)
{
    std::ostringstream builder;
    builder << static_cast<int>(species) << '|'
            << static_cast<int>(physics_list) << '|'
            << normalizeToken(material_name) << '|'
            << requested_version << '|'
            << (allow_rollback ? "rollback_on" : "rollback_off");
    return builder.str();
}

} // namespace

RadiationStoppingDataGovernance::RadiationStoppingDataGovernance()
{
    version_precedence_ = {
        "2026.04-rd002",
        "2026.04-rd003",
    };

    material_records_ = {
        {"2026.04-rd002", "kapton", "sha256:rd002-kapton-zovera", "material://kapton", "2026-04-02T00:00:00Z", 0.520, 79.0},
        {"2026.04-rd002", "ptfe", "sha256:rd002-ptfe-zovera", "material://ptfe", "2026-04-02T00:00:00Z", 0.480, 99.0},
        {"2026.04-rd002", "aluminum", "sha256:rd002-aluminum-zovera", "material://aluminum", "2026-04-02T00:00:00Z", 0.481, 166.0},

        {"2026.04-rd003", "kapton", "sha256:rd003-kapton-zovera", "material://kapton", "2026-04-03T00:00:00Z", 0.521, 79.5},
        {"2026.04-rd003", "ptfe", "sha256:rd003-ptfe-zovera", "material://ptfe", "2026-04-03T00:00:00Z", 0.482, 98.0},
        {"2026.04-rd003", "aluminum", "sha256:rd003-aluminum-zovera", "material://aluminum", "2026-04-03T00:00:00Z", 0.481, 166.0},
    };

    physics_records_ = {
        {"2026.04-rd002", RadiationPhysicsList::Geant4EmStandard, "sha256:rd002-physics-standard", "physics://geant4_em_standard", "2026-04-02T00:00:00Z", 1.00, 1.00},
        {"2026.04-rd002", RadiationPhysicsList::Geant4EmLivermore, "sha256:rd002-physics-livermore", "physics://geant4_em_livermore", "2026-04-02T00:00:00Z", 1.08, 1.15},
        {"2026.04-rd002", RadiationPhysicsList::Geant4EmPenelope, "sha256:rd002-physics-penelope", "physics://geant4_em_penelope", "2026-04-02T00:00:00Z", 0.95, 0.90},
        {"2026.04-rd002", RadiationPhysicsList::Geant4SpaceShielding, "sha256:rd002-physics-shielding", "physics://geant4_space_shielding", "2026-04-02T00:00:00Z", 1.15, 1.25},

        {"2026.04-rd003", RadiationPhysicsList::Geant4EmStandard, "sha256:rd003-physics-standard", "physics://geant4_em_standard", "2026-04-03T00:00:00Z", 1.00, 1.00},
        {"2026.04-rd003", RadiationPhysicsList::Geant4EmLivermore, "sha256:rd003-physics-livermore", "physics://geant4_em_livermore", "2026-04-03T00:00:00Z", 1.07, 1.12},
        {"2026.04-rd003", RadiationPhysicsList::Geant4EmPenelope, "sha256:rd003-physics-penelope", "physics://geant4_em_penelope", "2026-04-03T00:00:00Z", 0.96, 0.92},
        {"2026.04-rd003", RadiationPhysicsList::Geant4SpaceShielding, "sha256:rd003-physics-shielding", "physics://geant4_space_shielding", "2026-04-03T00:00:00Z", 1.13, 1.22},
    };
}

RadiationStoppingResolution RadiationStoppingDataGovernance::resolve(
    ParticleSpecies species, RadiationPhysicsList physics_list,
    const std::string& material_name, const std::string& requested_version,
    bool allow_rollback) const
{
    RadiationStoppingResolution resolution;
    resolution.requested_version = requested_version;
    resolution.material_name = normalizeToken(material_name);

    const std::string cache_key =
        makeCacheKey(species, physics_list, material_name, requested_version, allow_rollback);
    if (const auto it = cache_.find(cache_key); it != cache_.end())
    {
        resolution = it->second;
        resolution.cache_hit = true;
        return resolution;
    }

    const std::string normalized_material = normalizeToken(material_name);
    std::set<std::string> material_versions;
    for (const auto& record : material_records_)
    {
        if (normalizeToken(record.material_name) == normalized_material)
        {
            material_versions.insert(record.version);
        }
    }

    if (material_versions.empty())
    {
        resolution.validation_error = "material_not_found:" + material_name;
        cache_[cache_key] = resolution;
        return resolution;
    }

    std::set<std::string> available_versions;
    for (const auto& record : physics_records_)
    {
        if (record.physics_list == physics_list && material_versions.contains(record.version))
        {
            available_versions.insert(record.version);
        }
    }

    if (available_versions.empty())
    {
        resolution.validation_error = "physics_version_not_found";
        cache_[cache_key] = resolution;
        return resolution;
    }

    std::string selected_version;
    const std::string requested = requested_version.empty() ? "latest" : requested_version;
    if (requested == "latest")
    {
        const auto latest = selectLatestVersion(available_versions, version_precedence_);
        if (!latest)
        {
            resolution.validation_error = "version_resolution_failed";
            cache_[cache_key] = resolution;
            return resolution;
        }
        selected_version = *latest;
    }
    else if (available_versions.contains(requested))
    {
        selected_version = requested;
    }
    else if (allow_rollback)
    {
        const auto rollback_version =
            selectRollbackVersion(available_versions, version_precedence_, requested);
        if (!rollback_version)
        {
            resolution.validation_error = "rollback_resolution_failed";
            cache_[cache_key] = resolution;
            return resolution;
        }

        selected_version = *rollback_version;
        resolution.rollback_applied = true;
        resolution.rollback_reason =
            "requested_version_not_available:" + requested + "->" + selected_version;
    }
    else
    {
        resolution.validation_error = "requested_version_not_available:" + requested;
        cache_[cache_key] = resolution;
        return resolution;
    }

    const MaterialRecord* material_record = nullptr;
    for (const auto& record : material_records_)
    {
        if (record.version == selected_version &&
            normalizeToken(record.material_name) == normalized_material)
        {
            material_record = &record;
            break;
        }
    }

    const PhysicsRecord* physics_record = nullptr;
    for (const auto& record : physics_records_)
    {
        if (record.version == selected_version && record.physics_list == physics_list)
        {
            physics_record = &record;
            break;
        }
    }

    if (material_record == nullptr || physics_record == nullptr)
    {
        resolution.validation_error = "resolved_record_missing";
        cache_[cache_key] = resolution;
        return resolution;
    }

    resolution.resolved_version = selected_version;
    resolution.material_name = material_record->material_name;
    resolution.z_over_a = material_record->z_over_a;
    resolution.mean_excitation_ev = material_record->mean_excitation_ev;
    resolution.stopping_scale = physics_record->stopping_scale;
    resolution.scattering_scale = physics_record->scattering_scale;
    resolution.material_source_hash = material_record->source_hash;
    resolution.physics_source_hash = physics_record->source_hash;
    resolution.dataset_source_uri = material_record->source_uri + "|" + physics_record->source_uri;
    resolution.dataset_created_utc =
        std::max(material_record->created_utc, physics_record->created_utc);
    resolution.dataset_id = "stopping-data/" + selected_version + "/" +
                            normalizeToken(toString(species)) + "/" +
                            normalizeToken(toString(physics_list)) + "/" +
                            normalizeToken(material_record->material_name);

    if (!isFinitePositive(resolution.z_over_a) || resolution.z_over_a > 1.0)
    {
        resolution.validation_error = "invalid_z_over_a";
    }
    else if (!isFinitePositive(resolution.mean_excitation_ev))
    {
        resolution.validation_error = "invalid_mean_excitation_ev";
    }
    else if (!isFinitePositive(resolution.stopping_scale))
    {
        resolution.validation_error = "invalid_stopping_scale";
    }
    else if (!isFinitePositive(resolution.scattering_scale))
    {
        resolution.validation_error = "invalid_scattering_scale";
    }
    else
    {
        resolution.valid = true;
    }

    cache_[cache_key] = resolution;
    return resolution;
}

void RadiationStoppingDataGovernance::clearCache() const
{
    cache_.clear();
}

std::size_t RadiationStoppingDataGovernance::cacheSize() const
{
    return cache_.size();
}

} // namespace Radiation
} // namespace Toolkit
} // namespace SCDAT
