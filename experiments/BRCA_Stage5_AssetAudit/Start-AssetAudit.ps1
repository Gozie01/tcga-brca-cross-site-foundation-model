param([string]$ProjectRoot = "$env:USERPROFILE\Desktop\BioMedical")
$ErrorActionPreference = 'Stop'
$pythonExe = Join-Path $ProjectRoot '.venv-brca-revision\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonExe)) { throw "BRCA environment missing: $pythonExe" }
& $pythonExe (Join-Path $PSScriptRoot 'collect.py') --project-root $ProjectRoot
if ($LASTEXITCODE -ne 0) { throw 'Asset collection stopped. Send the error output.' }
