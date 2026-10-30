param([string]$ProjectRoot = "$env:USERPROFILE\Desktop\BioMedical", [string]$WSIRoot = 'D:\TCGA_BRCA_WSI')
$ErrorActionPreference = 'Stop'
$pythonExe = Join-Path $ProjectRoot '.venv-brca-revision\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonExe)) { throw "BRCA environment missing: $pythonExe" }
& $pythonExe (Join-Path $PSScriptRoot 'check.py') --project-root $ProjectRoot --wsi-root $WSIRoot
if ($LASTEXITCODE -ne 0) { throw 'Compatibility check incomplete. Upload the printed report ZIP and error output.' }
