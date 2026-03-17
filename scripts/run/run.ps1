# ==============================================================================
# SCDAT - Cross-platform Auto Build & Run Script
# Spacecraft Charging and Discharging Analysis Toolkit
#
# Usage (can be invoked from any directory):
#   powershell -File scripts/run/run.ps1 [-Target SCDAT] [-Config Debug] [-Rebuild]
#   pwsh        scripts/run/run.ps1 -List
#
# Parameters:
#   -Target  <name>   Executable target to run, default: SCDAT
#                     Also accepts SCDAT_PICcore or any auto-generated target
#   -Config  <type>   Debug | Release | RelWithDebInfo, default: Debug
#   -Jobs    <N>      Parallel build jobs, default: auto-detect CPU count
#   -Rebuild          Delete build dir and reconfigure from scratch
#   -ConfigOnly       Only run CMake configure, skip build and run
#   -BuildOnly        Configure + build, skip run
#   -List             List all executables in build/bin then exit
#   -Args    <str>    Arguments forwarded to the executable (quote the whole string)
# ==============================================================================
[CmdletBinding()]
param(
    [string] $Target     = "SCDAT",
    [ValidateSet("Debug","Release","RelWithDebInfo","MinSizeRel")]
    [string] $Config     = "Debug",
    [int]    $Jobs       = 0,           # 0 = 自动检测 CPU 核数
    [switch] $Rebuild,
    [switch] $ConfigOnly,
    [switch] $BuildOnly,
    [switch] $List,
    [string] $Args       = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ──────────────────────────────────────────────────────────────────────────────
# 辅助函数
# ──────────────────────────────────────────────────────────────────────────────
function Write-Step([string]$msg) {
    Write-Host ""
    Write-Host ">>> $msg" -ForegroundColor Cyan
}

function Write-OK([string]$msg) {
    Write-Host "    [OK] $msg" -ForegroundColor Green
}

function Write-Fail([string]$msg) {
    Write-Host "    [FAIL] $msg" -ForegroundColor Red
}

function Invoke-Cmd([string]$cmd, [string[]]$cmdArgs) {
    Write-Host "    $ $cmd $($cmdArgs -join ' ')" -ForegroundColor DarkGray
    & $cmd @cmdArgs
    if ($LASTEXITCODE -ne 0) {
        Write-Fail "Command failed with exit code $LASTEXITCODE"
        exit $LASTEXITCODE
    }
}

# ──────────────────────────────────────────────────────────────────────────────
# Locate project root (works regardless of where the script is called from)
# ──────────────────────────────────────────────────────────────────────────────
$ScriptDir  = Split-Path -Parent $MyInvocation.MyCommand.Definition
$ProjectDir = (Resolve-Path (Join-Path $ScriptDir "../..")).Path
$BuildDir   = Join-Path $ProjectDir "build"

Write-Host ""
Write-Host "============================================================" -ForegroundColor DarkCyan
Write-Host "  SCDAT - Auto Build & Run Script" -ForegroundColor DarkCyan
Write-Host "  Project : $ProjectDir" -ForegroundColor DarkCyan
Write-Host "  Build   : $BuildDir" -ForegroundColor DarkCyan
Write-Host "  Config  : $Config    Target: $Target" -ForegroundColor DarkCyan
Write-Host "============================================================" -ForegroundColor DarkCyan

# ──────────────────────────────────────────────────────────────────────────────
# Detect platform and toolchain (compatible with PS 5.1 and PS 7+)
# ──────────────────────────────────────────────────────────────────────────────
if ($PSVersionTable.PSVersion.Major -ge 6) {
    $IsWin = $IsWindows
    $IsMac = $IsMacOS
    $IsLin = $IsLinux
} else {
    # PowerShell 5.1 only runs on Windows
    $IsWin = $true
    $IsMac = $false
    $IsLin = $false
}

Write-Step "Detecting toolchain"

# CMake
if (-not (Get-Command cmake -ErrorAction SilentlyContinue)) {
    Write-Fail "cmake not found. Please install CMake >= 3.20"
    exit 1
}
$cmakeVer = (cmake --version 2>&1 | Select-Object -First 1) -replace "cmake version ",""
Write-OK "CMake $cmakeVer"

# Select generator and compiler
$Generator = $null
$CxxCompiler = $null
$ExeSuffix = if ($IsWin) { ".exe" } else { "" }

if ($IsWin) {
    # Prefer Ninja (faster), fall back to MinGW Make
    if (Get-Command ninja -ErrorAction SilentlyContinue) {
        $Generator = "Ninja"
        Write-OK "Generator: Ninja"
    } elseif (Get-Command mingw32-make -ErrorAction SilentlyContinue) {
        $Generator = "MinGW Makefiles"
        Write-OK "Generator: MinGW Makefiles"
    } elseif (Get-Command make -ErrorAction SilentlyContinue) {
        $Generator = "Unix Makefiles"
        Write-OK "Generator: Unix Makefiles (MSYS/Cygwin)"
    }
    # Compiler priority: g++ > clang++ > cl
    foreach ($cc in @("g++","clang++","cl")) {
        if (Get-Command $cc -ErrorAction SilentlyContinue) {
            $CxxCompiler = (Get-Command $cc).Source
            Write-OK "CXX: $CxxCompiler"
            break
        }
    }
} else {
    # Linux / macOS
    if (Get-Command ninja -ErrorAction SilentlyContinue) {
        $Generator = "Ninja"
        Write-OK "Generator: Ninja"
    } else {
        $Generator = "Unix Makefiles"
        Write-OK "Generator: Unix Makefiles"
    }
    foreach ($cc in @("g++","clang++","c++")) {
        if (Get-Command $cc -ErrorAction SilentlyContinue) {
            $CxxCompiler = (Get-Command $cc).Source
            Write-OK "CXX: $CxxCompiler"
            break
        }
    }
}

if (-not $Generator) {
    Write-Fail "No build backend found (Ninja / MinGW Make / make). Check PATH."
    exit 1
}
if (-not $CxxCompiler) {
    Write-Fail "No C++ compiler found (g++ / clang++ / cl). Check PATH."
    exit 1
}

# Parallel job count
if ($Jobs -le 0) {
    $Jobs = if ($IsWin) {
        (Get-CimInstance -ClassName Win32_ComputerSystem).NumberOfLogicalProcessors
    } elseif ($IsMac) {
        [int](sysctl -n hw.logicalcpu)
    } else {
        [int](nproc)
    }
    if ($Jobs -le 0) { $Jobs = 4 }
}
Write-OK "Parallel jobs: $Jobs"

# ──────────────────────────────────────────────────────────────────────────────
# Optional: list compiled executables
# ──────────────────────────────────────────────────────────────────────────────
if ($List) {
    $binDir = Join-Path $BuildDir "bin"
    if (Test-Path $binDir) {
        Write-Step "Executables in build/bin"
        Get-ChildItem -Path $binDir -File | ForEach-Object {
            Write-Host "    $($_.Name)" -ForegroundColor Yellow
        }
    } else {
        Write-Host "    build/bin does not exist. Build first." -ForegroundColor Yellow
    }
    exit 0
}

# ──────────────────────────────────────────────────────────────────────────────
# Optional: clean build directory
# ──────────────────────────────────────────────────────────────────────────────
if ($Rebuild -and (Test-Path $BuildDir)) {
    Write-Step "Removing old build directory"
    Remove-Item -Recurse -Force $BuildDir
    Write-OK "Deleted $BuildDir"
}

# ──────────────────────────────────────────────────────────────────────────────
# CMake configure (only when CMakeCache.txt is absent)
# ──────────────────────────────────────────────────────────────────────────────
$cacheFile = Join-Path $BuildDir "CMakeCache.txt"
$needConfigure = (-not (Test-Path $cacheFile))

if ($needConfigure) {
    Write-Step "CMake configure"
    $cmakeConfigArgs = @(
        "-S", $ProjectDir,
        "-B", $BuildDir,
        "-G", $Generator,
        "-DCMAKE_CXX_COMPILER=$CxxCompiler",
        "-DCMAKE_BUILD_TYPE=$Config"
    )
    Invoke-Cmd cmake $cmakeConfigArgs
    Write-OK "Configure done"
} else {
    Write-Step "CMake already configured (skip; use -Rebuild to reconfigure)"
}

if ($ConfigOnly) { Write-OK "Config-only mode, exiting."; exit 0 }

# ──────────────────────────────────────────────────────────────────────────────
# Build
# ──────────────────────────────────────────────────────────────────────────────
Write-Step "Building target: $Target"
$buildArgs = @("--build", $BuildDir, "--target", $Target, "--parallel", "$Jobs")
if ($Config -and ($Generator -notmatch "Ninja|Make")) {
    # MSVC multi-config generators need --config
    $buildArgs += @("--config", $Config)
}
Invoke-Cmd cmake $buildArgs
Write-OK "Build succeeded"

if ($BuildOnly) { Write-OK "Build-only mode, exiting."; exit 0 }

# ──────────────────────────────────────────────────────────────────────────────
# Run
# ──────────────────────────────────────────────────────────────────────────────
Write-Step "Running $Target"

# Single-config generators (Ninja/Make): build/bin/
# Multi-config generators (MSVC):        build/bin/<Config>/
$binCandidates = @(
    (Join-Path $BuildDir "bin/$Target$ExeSuffix"),
    (Join-Path $BuildDir "bin/$Config/$Target$ExeSuffix"),
    (Join-Path $BuildDir "$Target$ExeSuffix")
)
$exePath = $binCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1

if (-not $exePath) {
    Write-Fail "Executable $Target$ExeSuffix not found. Searched:"
    $binCandidates | ForEach-Object { Write-Host "    $_" -ForegroundColor DarkYellow }
    exit 1
}

Write-Host "    $ $exePath $Args" -ForegroundColor DarkGray
Write-Host ""

if ($Args -ne "") {
    $argList = $Args -split "\s+"
    & $exePath @argList
} else {
    & $exePath
}

$exitCode = $LASTEXITCODE
Write-Host ""
if ($exitCode -eq 0) {
    Write-OK "$Target exited normally (code 0)"
} else {
    Write-Fail "$Target exited with code $exitCode"
}
exit $exitCode
