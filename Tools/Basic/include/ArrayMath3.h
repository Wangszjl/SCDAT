#pragma once

#include <array>
#include <cmath>

namespace SCDAT
{
namespace Basic
{

inline double squaredNorm3(const std::array<double, 3>& value)
{
    return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

inline double norm3(const std::array<double, 3>& value)
{
    return std::sqrt(squaredNorm3(value));
}

inline std::array<double, 3> subtract3(const std::array<double, 3>& a,
                                       const std::array<double, 3>& b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline std::array<double, 3> add3(const std::array<double, 3>& a,
                                  const std::array<double, 3>& b)
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline std::array<double, 3> scale3(const std::array<double, 3>& value, double scale)
{
    return {value[0] * scale, value[1] * scale, value[2] * scale};
}

inline std::array<double, 3> normalize3(const std::array<double, 3>& value,
                                        const std::array<double, 3>& fallback)
{
    const double magnitude = norm3(value);
    if (magnitude <= 1.0e-15)
    {
        return fallback;
    }
    return scale3(value, 1.0 / magnitude);
}

inline double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

} // namespace Basic
} // namespace SCDAT