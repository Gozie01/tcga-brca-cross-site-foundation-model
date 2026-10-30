param(
    [string]$ProjectRoot = ".",
    [string]$AuditDirectory = "",
    [string]$OutputDirectory = "",
    [string]$PythonPath = ""
)
$ErrorActionPreference = "Stop"
$project = (Resolve-Path -LiteralPath $ProjectRoot).Path
$codeRoot = Join-Path $project "BRCA_Revision2_Pipeline\code"
if (-not $PythonPath) { $PythonPath = Join-Path $project ".venv-brca-revision\Scripts\python.exe" }
if (-not $AuditDirectory) { $AuditDirectory = Join-Path $project "brca_revision2_asset_audit_20261019_165032" }
if (-not $OutputDirectory) { $OutputDirectory = Join-Path $project "brca_revision2_hibou_complete" }
if (-not (Test-Path -LiteralPath $PythonPath -PathType Leaf)) { throw "Python environment was not found: $PythonPath" }
if (-not (Test-Path -LiteralPath (Join-Path $AuditDirectory "asset_audit.json") -PathType Leaf)) { throw "The ready asset audit was not found." }
$original = Join-Path $codeRoot "extract_hibou.py"
if (-not (Test-Path -LiteralPath $original -PathType Leaf)) { throw "Original BRCA extraction package was not found under ProjectRoot." }
$originalHash = (Get-FileHash -LiteralPath $original -Algorithm SHA256).Hash.ToLowerInvariant()
if ($originalHash -ne "be1ab441766c1e9e6fd06108a0655bfd7b73019ebaa5ff5ffd631df62331af78") {
    throw "Original extractor differs from the tested version. No checkpoint metadata was changed."
}
$lock = Get-Content -LiteralPath (Join-Path $PSScriptRoot "repair_lock.json") -Raw | ConvertFrom-Json
foreach ($entry in $lock.files) {
    $source = Join-Path $PSScriptRoot $entry.path
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) { throw "Repair package file missing: $($entry.path)" }
    if ((Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash.ToLowerInvariant() -ne $entry.sha256) {
        throw "Repair package checksum mismatch: $($entry.path)"
    }
    $destination = Join-Path $codeRoot (Split-Path $entry.path -Leaf)
    if ((Test-Path -LiteralPath $destination -PathType Leaf) -and
        (Get-FileHash -LiteralPath $destination -Algorithm SHA256).Hash.ToLowerInvariant() -ne $entry.sha256) {
        throw "A different repair already exists: $destination"
    }
}
foreach ($entry in $lock.files) {
    Copy-Item -LiteralPath (Join-Path $PSScriptRoot $entry.path) `
        -Destination (Join-Path $codeRoot (Split-Path $entry.path -Leaf)) -Force
}
& $PythonPath (Join-Path $codeRoot "test_edgefix.py")
if ($LASTEXITCODE -ne 0) { throw "Repair regression tests failed. Checkpoint metadata was not changed." }
Remove-Item Env:HF_TOKEN -ErrorAction SilentlyContinue
Remove-Item Env:HUGGING_FACE_HUB_TOKEN -ErrorAction SilentlyContinue
$env:USE_TF = "0"
$env:USE_FLAX = "0"
Write-Host "Validating existing checkpoints; this reads the source slide files and can take time."
& $PythonPath (Join-Path $codeRoot "extract_hibou_edgefix.py") `
    --audit $AuditDirectory `
    --out $OutputDirectory `
    --model-revision "7fdc8cce365ebfd9e7e60ee6ee9a8d4d59388f3c" `
    --device cuda --batch-size 32 --migrate-existing
if ($LASTEXITCODE -ne 0) { throw "Extraction stopped. Preserve the entire output folder and return the error." }
$freeze = Get-Content -LiteralPath (Join-Path $OutputDirectory "feature_freeze.json") -Raw | ConvertFrom-Json
if ($freeze.status -ne "LOCKED_COMPLETE" -or $freeze.patients -ne 345 -or $freeze.slides -ne 347 -or $freeze.patches -ne 345937) {
    throw "Complete-feature validation failed."
}
Compress-Archive -Path (Join-Path $OutputDirectory "*") -DestinationPath "$OutputDirectory.zip" -Force
Write-Host "Upload this file: $OutputDirectory.zip"
