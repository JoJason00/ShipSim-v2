# Coverage report via clang-cl + llvm-cov. Requires LLVM on PATH.
# clang-cl targets the MSVC ABI, so the x64-windows vcpkg deps are reused as-is.
$ErrorActionPreference = "Stop"
Set-Location (Resolve-Path (Join-Path $PSScriptRoot ".."))

$Build = "build/coverage"
$Raw = "$Build/profraw"
$Profdata = "$Build/coverage.profdata"
$Exe = "$Build/bin/unit_tests.exe"

function Run($label, $block) {
    Write-Host "==== $label ====" -ForegroundColor Cyan
    & $block
    if ($LASTEXITCODE -ne 0) { throw "$label failed" }
}

Run "Configure" { cmake --preset coverage }
Run "Build" { cmake --build --preset coverage }

# catch_discover_tests runs the binary at build time and drops stray .profraw files.
Get-ChildItem $Build -Recurse -Filter *.profraw | Remove-Item -Force
New-Item -ItemType Directory -Force $Raw | Out-Null

# %p is the pid, so parallel ctest jobs don't clobber each other.
$env:LLVM_PROFILE_FILE = "$PWD/$Raw/unit_tests-%p.profraw"
Run "Test" { ctest --preset coverage }
Remove-Item Env:\LLVM_PROFILE_FILE

Run "Report" {
    # Each -flag=value needs quoting: PowerShell won't expand $vars in an unquoted
    # native argument that starts with a dash.
    $ignore = "(tests|vcpkg_installed|build)[\\/]"
    llvm-profdata merge -sparse (Get-ChildItem "$Raw/*.profraw").FullName -o $Profdata
    llvm-cov report $Exe "-instr-profile=$Profdata" "-ignore-filename-regex=$ignore"
    llvm-cov show $Exe "-instr-profile=$Profdata" "-ignore-filename-regex=$ignore" `
        -format=html "-output-dir=$Build/html" -show-branches=count -show-expansions
}

Start-Process (Resolve-Path "$Build/html/index.html")
