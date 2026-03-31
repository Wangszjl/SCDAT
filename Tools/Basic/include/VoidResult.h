/**
 * @file VoidResult.h
 * @brief VoidResult type for operation status
 * @ingroup CoreModule
 */

#pragma once

#include "ErrorCode.h"
#include <string>

namespace SCDAT {
namespace Core {

template <typename T> class Result;

/**
 * @brief Result type for void operations
 */
class VoidResult {
  public:
    VoidResult() : success_(true), error_code_(ErrorCode::SUCCESS), message_("") {}

    // Constructor from ErrorCode for compatibility
    VoidResult(ErrorCode error)
        : success_(error == ErrorCode::SUCCESS), error_code_(error), message_("Error occurred")
    {
    }

    template <typename T>
    VoidResult(const Result<T>& result)
        : success_(result.has_value()), error_code_(ErrorCode::SUCCESS), message_("")
    {
        if (!success_) {
            error_code_ = result.error();
        }
    }

    bool isSuccess() const { return success_; }
    bool isFailure() const { return !success_; }
    bool has_value() const { return success_; }
    ErrorCode error() const { return error_code_; }
    const std::string& message() const { return message_; }
    const std::string& getErrorMessage() const { return message_; }

    static VoidResult error(const std::string& msg) { return VoidResult(false, msg); }

    static VoidResult success() { return VoidResult(true, ""); }

    static VoidResult failure(const std::string& message) { return VoidResult(false, message); }

    explicit operator bool() const { return success_; }

  private:
    VoidResult(bool success, const std::string& message)
        : success_(success),
          error_code_(success ? ErrorCode::SUCCESS : ErrorCode::UNKNOWN_ERROR),
          message_(message)
    {
    }

    bool success_;
    ErrorCode error_code_;
    std::string message_;
};

} // namespace Core

// Make VoidResult available in SCDAT namespace
using Core::VoidResult;

} // namespace SCDAT
