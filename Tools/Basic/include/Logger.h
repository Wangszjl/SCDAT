#pragma once

#include <iostream>
#include <string>
#include <utility>

namespace SCDAT
{
namespace Core
{

template <typename... Args>
inline void logMessage(const char* level, const std::string& component,
                       const std::string& message, Args&&... /*args*/)
{
    std::clog << "[" << level << "][" << component << "] " << message << '\n';
}

} // namespace Core
} // namespace SCDAT

#define LOG_DEBUG(component, message, ...)                                                          \
    ::SCDAT::Core::logMessage("DEBUG", component, message, ##__VA_ARGS__)
#define LOG_INFO(component, message, ...)                                                           \
    ::SCDAT::Core::logMessage("INFO", component, message, ##__VA_ARGS__)
