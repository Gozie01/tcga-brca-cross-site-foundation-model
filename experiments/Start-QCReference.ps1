$ErrorActionPreference = 'Stop'
$taskScript = Join-Path $PSScriptRoot 'collect_tcga_qc_reference.py'
& py -3.11 $taskScript
if ($LASTEXITCODE -ne 0) { throw 'QC reference collection did not complete. Read the message above; no existing models or features were changed.' }
