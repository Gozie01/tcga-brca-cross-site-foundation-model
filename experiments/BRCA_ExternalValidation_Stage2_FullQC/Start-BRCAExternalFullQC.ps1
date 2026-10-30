param(
    [string]$ProjectRoot = "$env:USERPROFILE\Desktop\BioMedical",
    [string]$WSIRoot = 'D:\CPTAC_BRCA_WSI'
)
$ErrorActionPreference = 'Stop'
$pythonExe = Join-Path $ProjectRoot '.venv-brca-revision\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonExe -PathType Leaf)) {
    throw 'Existing BRCA environment not found.'
}
& $pythonExe (Join-Path $PSScriptRoot 'code\run_qc.py') --project-root $ProjectRoot --wsi-root $WSIRoot
if ($LASTEXITCODE -ne 0) { throw 'Full QC stopped; rerun the same command to resume. Send the error output for review.' }
