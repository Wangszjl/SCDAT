#pragma once

#include <algorithm>
#include <cctype>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

namespace SCDAT
{
namespace Basic
{

template <typename ModelType>
class StringModelRegistry
{
  public:
    using Entry = std::pair<std::string_view, ModelType>;

    StringModelRegistry(std::initializer_list<Entry> entries)
    {
        for (const auto& entry : entries)
        {
            registerAlias(entry.first, entry.second);
        }
    }

    void registerAlias(std::string_view alias, ModelType model)
    {
        std::string key = normalizeToken(alias);
        if (!key.empty())
        {
            alias_to_model_[std::move(key)] = model;
        }
    }

    std::optional<ModelType> tryParse(std::string_view text) const
    {
        const std::string key = normalizeToken(text);
        const auto it = alias_to_model_.find(key);
        if (it == alias_to_model_.end())
        {
            return std::nullopt;
        }
        return it->second;
    }

  private:
    static std::string normalizeToken(std::string_view text)
    {
        std::string lowered;
        lowered.reserve(text.size());
        std::transform(text.begin(), text.end(), std::back_inserter(lowered),
                       [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });

        std::string normalized;
        normalized.reserve(lowered.size());
        for (const unsigned char ch : lowered)
        {
            if (std::isalnum(ch))
            {
                normalized.push_back(static_cast<char>(ch));
            }
        }
        return normalized;
    }

    std::unordered_map<std::string, ModelType> alias_to_model_;
};

} // namespace Basic
} // namespace SCDAT
