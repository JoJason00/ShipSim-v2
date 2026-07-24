@echo off
setlocal EnableDelayedExpansion

:: ShipSim build helper.
::
::   build.bat            configure + build + test  (debug)
::   build.bat release    configure + build         (release)
::   build.bat asan       configure + build + test  (AddressSanitizer)
::   build.bat profile    configure + build         (RelWithDebInfo)
::
:: Why this exists: it locates Visual Studio and initialises the MSVC
:: environment, so the project builds from any plain terminal without
:: opening a Developer prompt. The build configuration itself lives in
:: CMakePresets.json; this script does not duplicate it.
::
:: Kept ASCII-only on purpose: cmd.exe tracks its position in the batch
:: file byte-wise, and "chcp 65001" partway through a file containing
:: multi-byte characters desyncs the parser.

set "PRESET=%~1"
if "%PRESET%"=="" set "PRESET=debug"

echo =========================================
echo    ShipSim build  [preset: %PRESET%]
echo =========================================
echo.

:: ---------- locate Visual Studio ----------
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" goto :no_vswhere

for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VS_PATH=%%i"

if not defined VS_PATH goto :no_vs

:: ---------- initialise MSVC ----------
:: Ninja only looks at PATH. Without this, it silently picks up MinGW, whose
:: ABI does not match the MSVC-built vcpkg libraries; the guard in
:: CMakeLists.txt catches that and explains the fix.
call "!VS_PATH!\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64 >nul
if errorlevel 1 goto :vs_init_failed
echo [OK] MSVC ready
echo.

echo ===== Configure =====
cmake --preset %PRESET%
if errorlevel 1 goto :configure_failed

echo.
echo ===== Build =====
cmake --build --preset %PRESET%
if errorlevel 1 goto :build_failed

:: Only presets that declare a test preset: debug / coverage / asan
if "%PRESET%"=="release" goto :done
if "%PRESET%"=="profile" goto :done

echo.
echo ===== Test =====
ctest --preset %PRESET%
if errorlevel 1 goto :test_failed

:done
echo.
echo =========================================
echo    Build succeeded
echo =========================================
echo Binaries: %cd%\build\%PRESET%\bin
echo.
endlocal
exit /b 0

:no_vswhere
echo [ERROR] vswhere.exe not found - is Visual Studio installed?
exit /b 1

:no_vs
echo [ERROR] no Visual Studio installation with the C++ toolset was found
exit /b 1

:vs_init_failed
echo [ERROR] failed to initialise the MSVC environment
exit /b 1

:configure_failed
echo [ERROR] configure failed
echo         If you ever configured from a non-Developer shell, delete
echo         build\%PRESET% and retry - the failed cache pins the compiler.
exit /b 1

:build_failed
echo [ERROR] build failed
exit /b 1

:test_failed
echo [ERROR] tests failed
exit /b 1
