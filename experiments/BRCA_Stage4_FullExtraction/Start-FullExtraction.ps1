param([string]$ProjectRoot = "$env:USERPROFILE\Desktop\BioMedical", [string]$QCReport = '')
$ErrorActionPreference = 'Stop'
$pythonExe = Join-Path $ProjectRoot '.venv-brca-revision\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonExe)) { throw "BRCA environment missing: $pythonExe" }
$arguments = @((Join-Path $PSScriptRoot 'run.py'), '--project-root', $ProjectRoot)
if ($QCReport) { $arguments += @('--qc-report', $QCReport) }
& $pythonExe @arguments
if ($LASTEXITCODE -ne 0) { throw 'Extraction incomplete. Upload the printed report ZIP. Rerun the same command to resume unchanged checkpoints.' }
