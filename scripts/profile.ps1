# Build a target optimized-with-symbols and record a CPU profile with samply.
# Opens the Firefox Profiler in your browser automatically.
# Requires samply on PATH:  cargo install samply
#
#   ./scripts/profile.ps1                 # profiles shipsim_cli (default)
#   ./scripts/profile.ps1 <target>        # profiles any other executable target
#
# Sampling only catches work that runs ~1s or longer.
param([string]$Target = "shipsim_cli")
$ErrorActionPreference = "Stop"
Set-Location (Resolve-Path (Join-Path $PSScriptRoot ".."))

# samply reads every process name from a system-wide ETW trace and decodes them
# as UTF-8; a non-ASCII name (common on zh-CN Windows) makes it panic. Catch that
# here instead of wasting a whole profiling run on it.
$offenders = Get-Process | Where-Object { $_.ProcessName -match '[^\x00-\x7F]' } |
    Select-Object -ExpandProperty ProcessName -Unique
if ($offenders) {
    Write-Host "These running processes have non-ASCII names and will crash samply:" -ForegroundColor Yellow
    $offenders | ForEach-Object { Write-Host "  $_" -ForegroundColor Yellow }
    Write-Host "Close them and re-run. (samply bug, not your code.)" -ForegroundColor Yellow
    return
}

cmake --preset profile
if ($LASTEXITCODE -ne 0) { throw "configure failed" }

cmake --build --preset profile --target $Target
if ($LASTEXITCODE -ne 0) { throw "build failed" }

try {
    samply record ./build/profile/bin/$Target.exe
} finally {
    # A crashed run leaves a large intermediate .etl trace in the repo root.
    Remove-Item *.etl -ErrorAction SilentlyContinue
}
