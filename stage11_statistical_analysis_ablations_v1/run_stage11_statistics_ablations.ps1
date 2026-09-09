$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$pipeline = $PSScriptRoot
$python = if ($env:PYTHON_EXE) { $env:PYTHON_EXE } else { "python" }
$logDirectory = Join-Path $root "stage11_statistical_analysis_ablations_v1\logs"
New-Item -ItemType Directory -Force -Path $logDirectory | Out-Null
$log = Join-Path $logDirectory ("stage11_" + (Get-Date -Format "yyyyMMdd_HHmmss") + ".log")
Write-Host "Stage 11: paired statistics, domain shift, and representation ablations"
Write-Host "100 resumable ablation checkpoints plus paired bootstrap analyses"
Write-Host "Log: $log"
Start-Transcript -Path $log
try {
    & $python (Join-Path $pipeline "stage11_run_statistics_ablations_v1.py") `
        --root $root `
        --bootstrap-reps 1000 `
        --n-jobs -1
    if ($LASTEXITCODE -ne 0) {
        throw "Stage 11 failed (exit $LASTEXITCODE)."
    }
}
finally {
    Stop-Transcript
}
