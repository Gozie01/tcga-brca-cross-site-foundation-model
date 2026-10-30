param(
    [ValidateSet("Audit", "Pilot", "Extract", "Evaluate", "ExternalMetadata", "SelfTest")]
    [string]$Mode = "Audit",
    [string]$ProjectRoot = ".",
    [string[]]$AdditionalWsiRoots = @(),
    [string]$PythonPath = "python",
    [string]$AuditDirectory = "",
    [string]$FeatureDirectory = "",
    [string]$OutputDirectory = "",
    [int]$BatchSize = 32,
    [ValidateSet("auto", "cuda", "cpu")]
    [string]$Device = "auto"
)
$ErrorActionPreference = "Stop"
$packageRoot = $PSScriptRoot
$resolvedProject = (Resolve-Path -LiteralPath $ProjectRoot).Path
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
function Invoke-CheckedPython {
    param([string[]]$Arguments)
    & $PythonPath @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Python step failed; original data and results were preserved. See the error above." }
}
if ($Mode -eq "Audit") {
    Invoke-CheckedPython -Arguments @((Join-Path $packageRoot "code\test_asset_audit.py"))
    if (-not $OutputDirectory) { $OutputDirectory = Join-Path $resolvedProject "brca_revision2_asset_audit_$stamp" }
    $commandArgs = @((Join-Path $packageRoot "code\audit_assets.py"), "--root", $resolvedProject, "--out", $OutputDirectory)
    foreach ($wsiRoot in $AdditionalWsiRoots) { $commandArgs += @("--root", $wsiRoot) }
    Invoke-CheckedPython -Arguments $commandArgs
} elseif ($Mode -eq "SelfTest") {
    Invoke-CheckedPython -Arguments @((Join-Path $packageRoot "code\test_asset_audit.py"))
    Invoke-CheckedPython -Arguments @((Join-Path $packageRoot "code\test_pipeline.py"))
    Invoke-CheckedPython -Arguments @((Join-Path $packageRoot "code\test_manuscript_patch.py"))
} elseif ($Mode -eq "Pilot" -or $Mode -eq "Extract") {
    if (-not $AuditDirectory) { throw "Supply -AuditDirectory from the completed asset audit." }
    if (-not $OutputDirectory) { $OutputDirectory = Join-Path $resolvedProject "brca_revision2_hibou_$($Mode.ToLower())_$stamp" }
    $commandArgs = @((Join-Path $packageRoot "code\extract_hibou.py"), "--audit", $AuditDirectory, "--out", $OutputDirectory,
                     "--batch-size", "$BatchSize", "--device", $Device)
    if ($Mode -eq "Pilot") { $commandArgs += @("--pilot-patches", "128") }
    Invoke-CheckedPython -Arguments $commandArgs
} elseif ($Mode -eq "Evaluate") {
    if (-not $FeatureDirectory) { throw "Supply -FeatureDirectory containing the complete Hibou feature freeze." }
    if (-not $OutputDirectory) { $OutputDirectory = Join-Path $resolvedProject "brca_revision2_encoder_benchmark_$stamp" }
    Invoke-CheckedPython -Arguments @((Join-Path $packageRoot "code\run_encoder_benchmark.py"), "--features", $FeatureDirectory,
                                    "--out", $OutputDirectory, "--bootstrap", "1000")
} elseif ($Mode -eq "ExternalMetadata") {
    if (-not $OutputDirectory) { $OutputDirectory = Join-Path $resolvedProject "brca_revision2_external_metadata_$stamp" }
    Invoke-CheckedPython -Arguments @((Join-Path $packageRoot "code\audit_external_cohort.py"), "--out", $OutputDirectory)
}
