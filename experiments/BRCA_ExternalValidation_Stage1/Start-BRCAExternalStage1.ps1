param(
    [ValidateSet('Audit', 'Pilot', 'DownloadAll')]
    [string]$Mode = 'Audit',
    [string]$ProjectRoot = "$env:USERPROFILE\Desktop\BioMedical",
    [string]$WSIRoot = 'D:\CPTAC_BRCA_WSI'
)
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path -LiteralPath $ProjectRoot -ErrorAction Stop).Path
$slidePython = Join-Path $project '.venv-brca-revision\Scripts\python.exe'
$runner = Join-Path $PSScriptRoot 'code\download_cptac.py'

if (Test-Path -LiteralPath $slidePython -PathType Leaf) {
    $pythonExe = $slidePython
    $pythonPrefix = @()
} else {
    $pythonExe = 'py'
    $pythonPrefix = @('-3.11')
}

if ($Mode -ne 'Audit') {
    if (-not (Test-Path -LiteralPath $slidePython -PathType Leaf)) {
        throw 'The existing BRCA Python environment was not found. Run Audit and send its output before downloading.'
    }
    & $slidePython -c 'import openslide; from PIL import Image; print("OpenSlide environment available.")'
    if ($LASTEXITCODE -ne 0) { throw 'OpenSlide could not load in the existing environment. No slide download was started.' }
}

$arguments = @($runner, '--mode', $Mode, '--project-root', $project, '--wsi-root', $WSIRoot)
if ($Mode -ne 'Audit') { $arguments += @('--slide-python', $slidePython) }
& $pythonExe @pythonPrefix @arguments
if ($LASTEXITCODE -ne 0) { throw 'External download/preflight stopped. Existing TCGA files and partial CPTAC downloads were preserved.' }
