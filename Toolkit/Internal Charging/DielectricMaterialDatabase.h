#pragma once

#include "../../Tools/Material/include/MaterialDatabase.h"

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

class DielectricMaterialDatabase
{
  public:
    DielectricMaterialDatabase();

    const Material::MaterialProperty* findByName(const std::string& name) const;

  private:
    Material::UnifiedMaterialDatabase database_;
};

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
