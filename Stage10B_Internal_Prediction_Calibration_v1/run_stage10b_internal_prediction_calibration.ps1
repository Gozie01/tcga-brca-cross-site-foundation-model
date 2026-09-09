$ErrorActionPreference = "Stop"

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$pipeline = $PSScriptRoot
$python = if ($env:PYTHON_EXE) { $env:PYTHON_EXE } else { "python" }

$logDir = Join-Path $root "stage10b_internal_prediction_calibration_v1\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$log = Join-Path $logDir "stage10b_internal_prediction_calibration_$stamp.log"

Write-Host "Stage 10B: internal prediction and calibration"
Write-Host "Five endpoints x five frozen outer folds; resumable by split"
Write-Host "Log: $log"

Start-Transcript -Path $log
try {
    & $python `
        (Join-Path $pipeline "stage10b_run_internal_prediction_calibration_v1.py") `
        --root $root `
        --bootstrap-reps 1000 `
        --n-jobs -1
    if ($LASTEXITCODE -ne 0) {
        throw "Stage 10B failed (exit $LASTEXITCODE)."
    }
}
finally {
    Stop-Transcript
}
