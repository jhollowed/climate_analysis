#!/bin/bash

massmag=$1

python ./limvar_TEM_budget_latbands.py 1 $massmag 0 36 1 0
python ./limvar_TEM_budget_latbands.py 2 $massmag 0 36 1 0
python ./limvar_TEM_budget_latbands.py 3 $massmag 0 36 1 0
python ./limvar_TEM_budget_latbands.py 4 $massmag 0 36 1 0
python ./limvar_TEM_budget_latbands.py 5 $massmag 0 36 1 0
python ./limvar_TEM_budget_latbands.py 6 $massmag 0 36 1 0
python ./limvar_TEM_budget_latbands.py 7 $massmag 0 36 1 0
python ./limvar_TEM_budget_latbands.py 8 $massmag 0 36 1 0

python ./limvar_TEM_budget_monthly.py 1 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 2 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 3 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 4 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 5 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 6 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 7 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 8 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 9 $massmag 0 36 1 0
python ./limvar_TEM_budget_monthly.py 9 $massmag 0 36 1 0
