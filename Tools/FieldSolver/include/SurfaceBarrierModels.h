#pragma once

#include "../../Material/include/MaterialProperty.h"
#include "../../Particle/include/SurfaceDistributionFunction.h"

#include <string>
#include <vector>

namespace SCDAT
{
namespace FieldSolver
{

struct SurfaceBarrierState
{
    double local_potential_v = 0.0;
    double reference_potential_v = 0.0;
    double barrier_potential_v = 0.0;
    double normal_electric_field_v_per_m = 0.0;
    double emission_temperature_ev = 2.0;
};

struct SurfaceBarrierEvaluation
{
    double escaped_current_a_per_m2 = 0.0;
    double recollected_current_a_per_m2 = 0.0;
    double scaling = 1.0;
    double didv_a_per_m2_per_v = 0.0;
    bool valid = false;
};

struct SurfaceEmissionBarrierInputs
{
    double electron_collection_a_per_m2 = 0.0;
    double ion_collection_a_per_m2 = 0.0;
    double secondary_emission_a_per_m2 = 0.0;
    double ion_secondary_emission_a_per_m2 = 0.0;
    double backscatter_emission_a_per_m2 = 0.0;
    double photo_emission_a_per_m2 = 0.0;
    double photo_incidence_scale = 1.0;
    double user_secondary_scale = 1.0;
    double user_ion_secondary_scale = 1.0;
    double user_backscatter_scale = 1.0;
    double user_photo_scale = 1.0;
    double fallback_secondary_yield = 0.0;
    double fallback_ion_secondary_yield = 0.0;
    double fallback_backscatter_yield = 0.0;
};

struct SurfaceEmissionBarrierComponentInputs
{
    double electron_collection_a_per_m2 = 0.0;
    double ion_collection_a_per_m2 = 0.0;
    double secondary_emission_a_per_m2 = 0.0;
    double ion_secondary_emission_a_per_m2 = 0.0;
    double backscatter_emission_a_per_m2 = 0.0;
    double photo_emission_a_per_m2 = 0.0;
    double photo_incidence_scale = 1.0;
    double user_secondary_scale = 1.0;
    double user_ion_secondary_scale = 1.0;
    double user_backscatter_scale = 1.0;
    double user_photo_scale = 1.0;
    double fallback_secondary_yield = 0.0;
    double fallback_ion_secondary_yield = 0.0;
    double fallback_backscatter_yield = 0.0;
};

struct SurfaceEmissionBarrierOutputs
{
    double secondary_emission_scale = 0.0;
    double ion_secondary_emission_scale = 0.0;
    double backscatter_scale = 0.0;
    double photo_emission_scale = 0.0;
    SurfaceBarrierEvaluation escaped_secondary{};
    SurfaceBarrierEvaluation escaped_ion_secondary{};
    SurfaceBarrierEvaluation escaped_backscatter{};
    SurfaceBarrierEvaluation escaped_photo{};
    bool valid = false;
};

class BarrierCurrentScaler
{
  public:
    virtual ~BarrierCurrentScaler() = default;
    virtual const char* family() const = 0;
    virtual SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                              const SurfaceBarrierState& state,
                                              double emitted_current_a_per_m2) const = 0;
};

enum class SurfaceBarrierScalerKind
{
        SecondaryRecollection = 0,
        VariableBarrier = 1,
        OmlCurrent = 2,
        LteOmlCurrent = 3,
        FowlerNordheimCurrent = 4,
        AutomaticBarrierCurrent = 5,
        MultipleCurrent = 6,
        GlobalTempCurrent = 7,
        SmoothedGlobalTempCurrent = 8,
        CurrentVariation = 9,
        LocalVariation = 10,
};

class SecondaryRecollectionScaler final : public BarrierCurrentScaler
{
  public:
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;
};

class VariableBarrierScaler final : public BarrierCurrentScaler
{
  public:
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;
};

class OmlCurrentScaler final : public BarrierCurrentScaler
{
    public:
        const char* family() const override;
        SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                                                            const SurfaceBarrierState& state,
                                                                            double emitted_current_a_per_m2) const override;
};

class LteOmlCurrentScaler final : public BarrierCurrentScaler
{
    public:
        const char* family() const override;
        SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                                                            const SurfaceBarrierState& state,
                                                                            double emitted_current_a_per_m2) const override;
};

class FowlerNordheimCurrentScaler final : public BarrierCurrentScaler
{
    public:
        const char* family() const override;
        SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                                                            const SurfaceBarrierState& state,
                                                                            double emitted_current_a_per_m2) const override;
};

class AutomaticBarrierCurrentScaler final : public BarrierCurrentScaler
{
  public:
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;
};

class MultipleCurrentScaler final : public BarrierCurrentScaler
{
  public:
    explicit MultipleCurrentScaler(std::string property_scope = {});
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;

  private:
    std::string property_scope_;
};

class GlobalTempCurrentScaler final : public BarrierCurrentScaler
{
  public:
    explicit GlobalTempCurrentScaler(std::string property_scope = {});
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;

  private:
    std::string property_scope_;
};

class SmoothedGlobalTempCurrentScaler final : public BarrierCurrentScaler
{
  public:
    explicit SmoothedGlobalTempCurrentScaler(std::string property_scope = {});
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;

  private:
    std::string property_scope_;
};

class CurrentVariationScaler final : public BarrierCurrentScaler
{
  public:
    explicit CurrentVariationScaler(std::string property_scope = {});
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;

  private:
    std::string property_scope_;
};

class LocalVariationScaler final : public BarrierCurrentScaler
{
  public:
    explicit LocalVariationScaler(std::string property_scope = {});
    const char* family() const override;
    SurfaceBarrierEvaluation evaluate(const Material::MaterialProperty& material,
                                      const SurfaceBarrierState& state,
                                      double emitted_current_a_per_m2) const override;

  private:
    std::string property_scope_;
};

struct SurfaceBarrierScalerDecision
{
    SurfaceBarrierScalerKind configured_kind = SurfaceBarrierScalerKind::VariableBarrier;
    SurfaceBarrierScalerKind resolved_kind = SurfaceBarrierScalerKind::VariableBarrier;
    std::string configured_family = "spis_variable_barrier_scaler_v1";
    std::string resolved_family = "spis_variable_barrier_scaler_v1";
};

SurfaceBarrierScalerKind resolveSurfaceBarrierScalerKind(
        const Material::MaterialProperty& material, const char* property_key,
        SurfaceBarrierScalerKind fallback_kind);
SurfaceBarrierScalerDecision resolveSurfaceBarrierScalerDecision(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    const char* property_key, SurfaceBarrierScalerKind fallback_kind);
const BarrierCurrentScaler& selectSurfaceBarrierScaler(SurfaceBarrierScalerKind kind);
SurfaceBarrierEvaluation evaluateConfiguredSurfaceBarrierScaler(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    double emitted_current_a_per_m2, const char* property_key,
    SurfaceBarrierScalerKind fallback_kind);

SurfaceEmissionBarrierInputs makeSurfaceEmissionBarrierInputs(
    const SurfaceEmissionBarrierComponentInputs& components);
SurfaceEmissionBarrierOutputs evaluateSurfaceEmissionBarrierBundle(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    const SurfaceEmissionBarrierComponentInputs& components);
SurfaceEmissionBarrierOutputs evaluateSurfaceEmissionBarrierBundle(
    const Material::MaterialProperty& material, const SurfaceBarrierState& state,
    const SurfaceEmissionBarrierInputs& inputs);
double evaluateLegacyThermalIonCollectionNaPerM2(
    const std::vector<Particle::LegacyPopulation>& populations,
    const Material::MaterialProperty& material,
    double surface_potential_v,
    bool include_potential_barrier);

} // namespace FieldSolver
} // namespace SCDAT
