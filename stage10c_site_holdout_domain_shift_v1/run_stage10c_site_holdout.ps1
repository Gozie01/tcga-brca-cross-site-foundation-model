$ErrorActionPreference="Stop"
$root=(Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$pipeline=$PSScriptRoot
$python=if($env:PYTHON_EXE){$env:PYTHON_EXE}else{"python"}
$logDir=Join-Path $root "stage10c_site_holdout_domain_shift_v1\logs"; New-Item -ItemType Directory -Force -Path $logDir|Out-Null
$log=Join-Path $logDir ("stage10c_"+(Get-Date -Format "yyyyMMdd_HHmmss")+".log")
Write-Host "Stage 10C: TCGA site-held-out domain-shift evaluation"; Write-Host "Five endpoints x four TCGA tissue-source sites; resumable by experiment"; Write-Host "Log: $log"
Start-Transcript -Path $log
try{& $python (Join-Path $pipeline "stage10c_run_site_holdout_v1.py") --root $root --bootstrap-reps 1000 --n-jobs -1; if($LASTEXITCODE-ne 0){throw "Stage 10C failed (exit $LASTEXITCODE)."}}finally{Stop-Transcript}
