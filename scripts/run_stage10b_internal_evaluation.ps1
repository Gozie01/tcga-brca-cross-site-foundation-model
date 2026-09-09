param([string]$Root = (Get-Location).Path)
python "$PSScriptRoot/../src/prediction/run_internal_evaluation.py" --root $Root
