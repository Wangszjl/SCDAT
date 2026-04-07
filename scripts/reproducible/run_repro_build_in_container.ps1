Param(
    [string]$ProjectRoot = (Get-Location).Path,
    [string]$ImageTag = "scdat-repro:eng007-v1"
)

$ErrorActionPreference = "Stop"

$dockerfile = Join-Path $ProjectRoot "scripts/reproducible/Dockerfile.scdat.repro"
if (-not (Test-Path $dockerfile)) {
    throw "Missing Dockerfile: $dockerfile"
}

Write-Host "[ENG-007] Building reproducible image: $ImageTag"
docker build -t $ImageTag -f $dockerfile $ProjectRoot

Write-Host "[ENG-007] Running reproducible build/test in container"
docker run --rm `
    -v "${ProjectRoot}:/workspace" `
    -w /workspace `
    $ImageTag `
    bash -lc "cmake -S . -B build_repro_container -G Ninja -DCMAKE_BUILD_TYPE=Debug && cmake --build build_repro_container --target SCDAT && ctest --test-dir build_repro_container -R 'UnifiedConfigSchemaV1Gate|UnifiedOutputContractV1Gate|ReproducibleBuildDependencyLockGate' --output-on-failure"

Write-Host "[ENG-007] Reproducible container run completed"
