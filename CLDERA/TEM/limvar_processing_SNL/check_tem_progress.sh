#!/bin/bash

cd /ascldap/users/jphollo/data/limvar_TEM 
ens1num=$(ls | grep ens1.eam | wc -l)
ens1cfnum=$(ls | grep ens1.cf.eam | wc -l)
ens2num=$(ls | grep ens2.eam | wc -l)
ens2cfnum=$(ls | grep ens2.cf.eam | wc -l)
ens3num=$(ls | grep ens3.eam | wc -l)
ens3cfnum=$(ls | grep ens3.cf.eam | wc -l)
ens4num=$(ls | grep ens4.eam | wc -l)
ens4cfnum=$(ls | grep ens4.cf.eam | wc -l)
ens5num=$(ls | grep ens5.eam | wc -l)
ens5cfnum=$(ls | grep ens5.cf.eam | wc -l)

echo ens1: $ens1num
echo ens1cf: $ens1cfnum
echo ens2: $ens2num
echo ens2cf: $ens2cfnum
echo ens3: $ens3num
echo ens3cf: $ens3cfnum
echo ens4: $ens4num
echo ens4cf: $ens4cfnum
echo ens5: $ens5num
echo ens5cf: $ens5cfnum

