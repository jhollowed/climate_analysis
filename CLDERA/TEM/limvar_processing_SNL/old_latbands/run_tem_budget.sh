#!/bin/bash

massmag=$1
overwrite=1

# ---- tem budget for lat band 10-daily data
python ./limvar_TEM_budget_latbands_10daily.py 1 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 2 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 3 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 4 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 5 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 6 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 7 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 8 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 9 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands_10daily.py 10 $massmag 0 36 $overwrite 0
exit

# ---- tem budget for lat bands
python ./limvar_TEM_budget_latbands.py 1 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 2 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 3 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 4 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 5 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 6 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 7 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 8 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 9 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_latbands.py 10 $massmag 0 36 $overwrite 0

# ---- tem budget for monthly data
python ./limvar_TEM_budget_monthly.py 1 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 2 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 3 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 4 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 5 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 6 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 7 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 8 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 9 $massmag 0 36 $overwrite 0
python ./limvar_TEM_budget_monthly.py 10 $massmag 0 36 $overwrite 0


