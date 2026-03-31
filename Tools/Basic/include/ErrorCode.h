/**
 * @file ErrorCode.h
 * @brief Error code enumeration definition
 *
 * This file defines the error code enumeration used throughout the project.
 * Separated from ErrorHandling.h to avoid circular dependencies.
 * @ingroup CoreModule
 */

#pragma once

namespace SCDAT {
namespace Core {

/**
 * @brief PIC-SPIS error code enumeration
 */
enum class ErrorCode : int {
    // General errors
    SUCCESS = 0,
    UNKNOWN_ERROR = 1,
    INVALID_ARGUMENT = 2,
    OUT_OF_RANGE = 3,
    NULL_POINTER = 4,
    MEMORY_ALLOCATION_FAILED = 5,

    // File I/O errors
    FILE_NOT_FOUND = 100,
    FILE_ACCESS_DENIED = 101,
    FILE_CORRUPTED = 102,
    FILE_FORMAT_UNSUPPORTED = 103,

    // Mesh related errors
    MESH_INVALID_TOPOLOGY = 200,
    MESH_DEGENERATE_ELEMENT = 201,
    MESH_BOUNDARY_MISMATCH = 202,
    MESH_MEMORY_INSUFFICIENT = 203,

    // Particle related errors
    PARTICLE_INVALID_TYPE = 300,
    PARTICLE_OUT_OF_BOUNDS = 301,
    PARTICLE_ENERGY_INVALID = 302,
    PARTICLE_COLLISION_FAILED = 303,

    // Solver related errors
    SOLVER_CONVERGENCE_FAILED = 400,
    SOLVER_MATRIX_SINGULAR = 401,
    SOLVER_ITERATION_LIMIT = 402,
    SOLVER_NUMERICAL_INSTABILITY = 403,

    // GPU related errors
    GPU_NOT_AVAILABLE = 500,
    GPU_MEMORY_INSUFFICIENT = 501,
    GPU_KERNEL_LAUNCH_FAILED = 502,
    GPU_SYNCHRONIZATION_FAILED = 503,

    // Parallel computing errors
    PARALLEL_THREAD_CREATION_FAILED = 600,
    PARALLEL_SYNCHRONIZATION_FAILED = 601,
    PARALLEL_LOAD_IMBALANCE = 602,
    PARALLEL_DEADLOCK_DETECTED = 603
};

} // namespace Core

// Make ErrorCode available in SCDAT namespace
using Core::ErrorCode;

} // namespace SCDAT
