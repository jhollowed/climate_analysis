#!/bin/bash

massmag=$1
overwrite=1

# ---- monthly data
python ./limvar_timeConcat.py 1 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 2 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 3 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 4 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 5 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 6 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 7 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 8 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 9 0 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 10 0 $massmag 0 0 36 $overwrite 0

# ---- daily data
python ./limvar_timeConcat.py 1 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 2 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 3 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 4 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 5 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 6 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 7 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 8 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 9 1 $massmag 0 36 $overwrite 0
python ./limvar_timeConcat.py 10 1 $massmag 0 36 $overwrite 0
