#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace SCDAT
{
namespace Basic
{

template <typename T>
inline double sumValue(const std::vector<T>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0,
                           [](double acc, const T& value) {
                               return acc + static_cast<double>(value);
                           });
}

template <typename T>
inline double meanValue(const std::vector<T>& values, double fallback = 0.0)
{
    if (values.empty())
    {
        return fallback;
    }
    const double sum = sumValue(values);
    return sum / static_cast<double>(values.size());
}

template <typename T>
inline double trapezoidalIntegralUniform(const std::vector<T>& samples, double step)
{
    if (samples.size() < 2 || !std::isfinite(step) || step <= 0.0)
    {
        return 0.0;
    }

    double sum = 0.0;
    for (std::size_t i = 1; i < samples.size(); ++i)
    {
        const double left = static_cast<double>(samples[i - 1]);
        const double right = static_cast<double>(samples[i]);
        sum += 0.5 * (left + right);
    }
    return sum * step;
}

} // namespace Basic
} // namespace SCDAT