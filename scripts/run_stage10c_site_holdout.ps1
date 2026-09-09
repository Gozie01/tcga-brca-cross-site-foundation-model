param([Parameter(Mandatory=$true)][string]$Root)
python "$PSScriptRoot/../src/domain_shift/run_site_holdout.py" --root $Root
