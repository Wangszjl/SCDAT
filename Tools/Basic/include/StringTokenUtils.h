#pragma once

#include <cctype>
#include <string>
#include <string_view>

namespace SCDAT
{
namespace Basic
{

inline std::string toLowerAscii(std::string_view text)
{
    std::string lowered;
    lowered.reserve(text.size());
    for (const unsigned char ch : text)
    {
        lowered.push_back(static_cast<char>(std::tolower(ch)));
    }
    return lowered;
}

inline std::string trimAscii(std::string_view text)
{
    std::size_t first = 0;
    while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first])) != 0)
    {
        ++first;
    }
    if (first == text.size())
    {
        return {};
    }

    std::size_t last = text.size();
    while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1])) != 0)
    {
        --last;
    }
    return std::string(text.substr(first, last - first));
}

inline std::string normalizeAlnumToken(std::string_view text)
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

} // namespace Basic
} // namespace SCDAT