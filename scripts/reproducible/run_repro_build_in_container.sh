#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:-$(pwd)}"
IMAGE_TAG="${2:-scdat-repro:eng007-v1}"
BUILD_DIR="build_repro_container"

DOCKERFILE_PATH="${PROJECT_ROOT}/scripts/reproducible/Dockerfile.scdat.repro"

if [ ! -f "${DOCKERFILE_PATH}" ]; then
  echo "Missing Dockerfile: ${DOCKERFILE_PATH}" >&2
  exit 2
fi

echo "[ENG-007] Building reproducible image: ${IMAGE_TAG}"
docker build -t "${IMAGE_TAG}" -f "${DOCKERFILE_PATH}" "${PROJECT_ROOT}"

echo "[ENG-007] Running reproducible build/test in container"
docker run --rm \
  -v "${PROJECT_ROOT}:/workspace" \
  -w /workspace \
  "${IMAGE_TAG}" \
  bash -lc "cmake -S . -B ${BUILD_DIR} -G Ninja -DCMAKE_BUILD_TYPE=Debug && cmake --build ${BUILD_DIR} --target SCDAT && ctest --test-dir ${BUILD_DIR} -R 'UnifiedConfigSchemaV1Gate|UnifiedOutputContractV1Gate|ReproducibleBuildDependencyLockGate' --output-on-failure"

echo "[ENG-007] Reproducible container run completed"
