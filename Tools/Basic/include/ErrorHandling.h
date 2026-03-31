#pragma once

#include "ErrorCode.h"
#include "VoidResult.h"
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

namespace SCDAT
{
namespace Core
{

template <typename T> class Result
{
  public:
    Result(T value) : value_(std::move(value)), error_code_(ErrorCode::SUCCESS) {}

    Result(ErrorCode error_code, std::string message = {})
        : value_(std::nullopt), error_code_(error_code), message_(std::move(message))
    {
    }

    bool has_value() const
    {
        return value_.has_value();
    }

    explicit operator bool() const
    {
        return has_value();
    }

    T& value()
    {
        return *value_;
    }

    const T& value() const
    {
        return *value_;
    }

    ErrorCode error() const
    {
        return error_code_;
    }

    const std::string& message() const
    {
        return message_;
    }

  private:
    std::optional<T> value_;
    ErrorCode error_code_;
    std::string message_;
};

class ErrorHandler
{
  public:
    static VoidResult makeSuccess()
    {
        return VoidResult::success();
    }

    template <typename T> static Result<T> makeSuccess(T value)
    {
        return Result<T>(std::move(value));
    }

    template <typename T>
    static auto makeError(ErrorCode error_code, std::string message = {})
        -> std::conditional_t<std::is_void_v<T>, VoidResult, Result<T>>
    {
        if constexpr (std::is_void_v<T>)
        {
            if (message.empty())
            {
                return VoidResult(error_code);
            }
            return VoidResult::failure(message);
        }
        else
        {
            return Result<T>(error_code, std::move(message));
        }
    }
};

} // namespace Core
} // namespace SCDAT
