$ErrorActionPreference = "Stop"

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$pipeline = $PSScriptRoot
$python = if ($env:PYTHON_EXE) { $env:PYTHON_EXE } else { "python" }
$env:TCGA_BRCA_PROJECT_ROOT = $root

if (-not $env:TCGA_BRCA_WSI_ROOT) {
    throw "Set TCGA_BRCA_WSI_ROOT to the directory containing A2_remaining, C8, E2, and BH."
}
$logDir = Join-Path $root "uni_multisite_embeddings_v1\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$log = Join-Path $logDir "stage9b_multisite_uni_$stamp.log"

Write-Host "Stage 9B multisite UNI extraction"
Write-Host "Sites: A2, C8, E2, BH"
Write-Host "Expected: 253,164 embeddings from 254 WSIs"
Write-Host "This run is resumable at the slide level."
Write-Host "Log: $log"

Start-Transcript -Path $log
try {
    foreach ($site in @("A2", "C8", "E2", "BH")) {
        Write-Host "`n=== STAGE 9B SITE: $site ==="
        & $python (Join-Path $pipeline "stage9b_extract_multisite_uni_v1.py") `
            --site $site `
            --batch-size 32
        if ($LASTEXITCODE -ne 0) {
            throw "Stage 9B extraction failed for $site (exit $LASTEXITCODE)."
        }
    }

    Write-Host "`n=== FINAL MULTISITE AUDIT/FREEZE ==="
    & $python (Join-Path $pipeline "stage9b_audit_freeze_multisite_uni_v1.py")
    if ($LASTEXITCODE -ne 0) {
        throw "Stage 9B final audit/freeze failed (exit $LASTEXITCODE)."
    }
}
finally {
    Stop-Transcript
}
