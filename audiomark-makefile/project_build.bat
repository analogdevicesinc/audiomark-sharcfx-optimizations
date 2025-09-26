@echo off
setlocal enabledelayedexpansion


:: List of 7 subprojects
set PROJECTS=audiomark-cces-library audiomark-sharc-fx audiomark-test-abf audiomark-test-aec audiomark-test-anr audiomark-test-kws audiomark-test-mfcc

for %%P in (%PROJECTS%) do (
    if exist "%%P" (
        echo ----------------------------------------
        echo ▶ Building project: %%P
        cd /d "%%P"

        echo Cleaning %%P...
        make clean

        echo Building %%P...
        make -j20
        if errorlevel 1 (
            echo ❌ ERROR: Build failed for %%P!
            exit /b 1
        )

        for %%F in (*.dxe) do (
            echo ✓ Build successful: %%F
        )

        cd ..
    ) else (
        echo ❌ ERROR: Folder not found: %%P
    )
)

echo ========================================
echo ✅ All projects built successfully.
endlocal