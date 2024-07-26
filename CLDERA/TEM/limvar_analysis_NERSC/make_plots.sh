#!/bin/bash

# tracer budgets october
python ./plot_monthly_impact.py E90TEND 1991 1 10 e90
python ./plot_monthly_impact.py QTRESVEL 1991 1 10 e90
python ./plot_monthly_impact.py qtendetfd 1991 1 10 e90
python ./plot_monthly_impact.py QTLOSS 1991 1 10 e90
python ./plot_monthly_impact.py QTDIFF 1991 1 10 e90
python ./plot_monthly_impact.py QTTOTAL 1991 1 10 e90

python ./plot_monthly_impact.py AOATEND 1991 1 10 aoa
python ./plot_monthly_impact.py UTRESVEL 1991 1 10 aoa
python ./plot_monthly_impact.py qtendetfd 1991 1 10 aoa
python ./plot_monthly_impact.py QTDIFF 1991 1 10 aoa
python ./plot_monthly_impact.py QTTOTAL 1991 1 10 aoa
exit

# calendars
python ./plot_monthly_impact.py U 1991 1
python ./plot_monthly_impact.py T 1991 1
python ./plot_monthly_impact.py AOA 1991 1
python ./plot_monthly_impact.py E90j 1991 1
python ./plot_monthly_impact.py psitem 1991 1

python ./plot_monthly_impact.py UTEND 1991 1
python ./plot_monthly_impact.py UTRESVEL 1991 1
python ./plot_monthly_impact.py utendepfd 1991 1
python ./plot_monthly_impact.py UTGWTOTAL 1991 1
python ./plot_monthly_impact.py UTDIFF 1991 1
python ./plot_monthly_impact.py UTTOTAL 1991 1
exit

# October
python ./plot_monthly_impact.py U 1991 1 10
python ./plot_monthly_impact.py T 1991 1 10
python ./plot_monthly_impact.py AOA 1991 1 10
python ./plot_monthly_impact.py E90j 1991 1 10
python ./plot_monthly_impact.py psitem 1991 1 10

python ./plot_monthly_impact.py UTEND 1991 1 10
python ./plot_monthly_impact.py UTRESVEL 1991 1 10
python ./plot_monthly_impact.py utendepfd 1991 1 10
python ./plot_monthly_impact.py UTGWTOTAL 1991 1 10
python ./plot_monthly_impact.py UTDIFF 1991 1 10
python ./plot_monthly_impact.py UTTOTAL 1991 1 10


