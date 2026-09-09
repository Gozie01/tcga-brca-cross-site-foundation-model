$ErrorActionPreference = "Stop"

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$pipeline = $PSScriptRoot
$python = if ($env:PYTHON_EXE) { $env:PYTHON_EXE } else { "python" }
$env:TCGA_BRCA_PROJECT_ROOT = $root

$logDir = Join-Path $root "stage10a_cohort_split_freeze_v1\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$log = Join-Path $logDir "stage10a_cohort_split_freeze_$stamp.log"

Write-Host "Stage 10A: cohort consolidation, patient pooling, and split freeze"
Write-Host "Expected: 345,937 patches -> 347 WSIs -> 345 patients"
Write-Host "Log: $log"

Start-Transcript -Path $log
try {
    & $python (Join-Path $pipeline "stage10a_build_cohort_and_splits_v1.py")
    if ($LASTEXITCODE -ne 0) {
        throw "Stage 10A failed (exit $LASTEXITCODE)."
    }
}
finally {
    Stop-Transcript
}
