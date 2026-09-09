param([Parameter(Mandatory=$true)][string]$Root)
python "$PSScriptRoot/../src/statistical_analysis/run_statistics_and_ablations.py" --root $Root
