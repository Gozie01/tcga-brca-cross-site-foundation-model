param(
 [string]$ProjectRoot = "$env:USERPROFILE\Desktop\BioMedical",
 [string]$Report = ''
)
$ErrorActionPreference = 'Stop'
$pythonExe = Join-Path $ProjectRoot '.venv-brca-revision\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonExe)) { throw "BRCA Python environment not found: $pythonExe" }
$arguments = @((Join-Path $PSScriptRoot 'check.py'), '--project-root', $ProjectRoot)
if ($Report) { $arguments += @('--report', $Report) }
& $pythonExe @arguments
if ($LASTEXITCODE -ne 0) { throw 'Preflight stopped. Send the error output.' }
