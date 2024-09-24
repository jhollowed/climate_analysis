#!/bin/bash

massmag=$1
overwrite=1

python ./limvar_monthly_seasonal_means.py 1 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 2 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 3 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 4 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 5 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 6 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 7 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 8 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 9 $massmag 0 36 $overwrite 0
python ./limvar_monthly_seasonal_means.py 10 $massmag 0 36 $overwrite 0
