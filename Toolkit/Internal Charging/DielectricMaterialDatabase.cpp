#include "DielectricMaterialDatabase.h"

namespace SCDAT
{
namespace Toolkit
{
namespace InternalCharging
{

DielectricMaterialDatabase::DielectricMaterialDatabase()
{
    database_.addAlias("polyimide", 2);
    database_.addAlias("kapton_hn", 2);
    database_.addAlias("teflon", 3);
}

const Material::MaterialProperty* DielectricMaterialDatabase::findByName(const std::string& name) const
{
    return database_.findByAliasOrName(name);
}

} // namespace InternalCharging
} // namespace Toolkit
} // namespace SCDAT
