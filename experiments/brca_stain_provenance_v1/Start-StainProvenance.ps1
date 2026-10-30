$ErrorActionPreference = 'Stop'
& py -3.11 (Join-Path $PSScriptRoot 'check_stain_provenance.py') @args
if ($LASTEXITCODE -ne 0) { throw 'Provenance collection failed. See the message above.' }
