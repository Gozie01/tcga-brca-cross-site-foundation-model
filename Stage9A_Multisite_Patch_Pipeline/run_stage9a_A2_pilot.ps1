$ErrorActionPreference = "Stop"

$project = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$pipeline = $PSScriptRoot
$python = if ($env:PYTHON_EXE) { $env:PYTHON_EXE } else { "python" }
$env:TCGA_BRCA_PROJECT_ROOT = $project

if (-not $env:TCGA_BRCA_WSI_ROOT) {
    throw "Set TCGA_BRCA_WSI_ROOT to the directory containing A2_remaining, C8, E2, and BH."
}

Write-Host "Stage 9A: A2 full-site pilot (11 WSIs)"
Write-Host "This run is resumable. Repeating it skips verified slide checkpoints."

& $python (Join-Path $pipeline "stage9a_extract_multisite_coordinates_v1.py") --site A2
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& $python (Join-Path $pipeline "stage9a_audit_freeze_site_v1.py") --site A2
exit $LASTEXITCODE
