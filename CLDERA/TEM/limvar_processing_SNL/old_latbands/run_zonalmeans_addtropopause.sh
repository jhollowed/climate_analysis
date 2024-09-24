#!/bin/bash

massmag=$1

# ---- daily data
python ./limvar_zonalmeans.py 1 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 2 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 3 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 4 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 5 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 6 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 7 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 8 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 9 1 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 10 1 $massmag 0 36 0 0 1

exit

# ---- monthly data
python ./limvar_zonalmeans.py 1 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 2 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 3 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 4 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 5 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 6 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 7 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 8 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 9 0 $massmag 0 36 0 0 1
python ./limvar_zonalmeans.py 10 0 $massmag 0 0 36 0 0 1
