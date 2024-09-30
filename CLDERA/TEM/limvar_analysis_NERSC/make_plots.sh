#!/bin/bash

overlay=0
    
# 5-year time series 
python ./plot_latband_impact.py U 5 $overlay
python ./plot_latband_impact.py T 5 $overlay
python ./plot_latband_impact.py AOA 5 $overlay
python ./plot_latband_impact.py E90j 5 $overlay
python ./plot_latband_impact.py psitem 5 $overlay

python ./plot_latband_impact.py UTEND 5 $overlay
python ./plot_latband_impact.py UTRESVEL 5 $overlay
python ./plot_latband_impact.py utendepfd 5 $overlay
python ./plot_latband_impact.py UTDIFF_NOGW 5 $overlay
python ./plot_latband_impact.py UTTOTAL_NOGW 5 $overlay
exit


python ./plot_monthly_impact.py psitem 1991 $overlay
python ./plot_monthly_impact.py epdiv 1991 $overlay
python ./plot_monthly_impact.py U 1991 $overlay
python ./plot_monthly_impact.py wtem 1991 $overlay
python ./plot_monthly_impact.py utendepfd 1991 $overlay
exit

if [ $1 == 1 ]; then
    # 1-year time series
    python ./plot_latband_impact.py U 1 $overlay
    python ./plot_latband_impact.py T 1 $overlay
    python ./plot_latband_impact.py AOA 1 $overlay
    python ./plot_latband_impact.py E90j 1 $overlay
    python ./plot_latband_impact.py psitem 1 $overlay

    python ./plot_latband_impact.py UTEND 1 $overlay
    python ./plot_latband_impact.py UTRESVEL 1 $overlay
    python ./plot_latband_impact.py utendepfd 1 $overlay
    python ./plot_latband_impact.py UTDIFF_NOGW 1 $overlay
    python ./plot_latband_impact.py UTTOTAL_NOGW 1 $overlay
    
    # 5-year time series 
    python ./plot_latband_impact.py U 5 $overlay
    python ./plot_latband_impact.py T 5 $overlay
    python ./plot_latband_impact.py AOA 5 $overlay
    python ./plot_latband_impact.py E90j 5 $overlay
    python ./plot_latband_impact.py psitem 5 $overlay

    python ./plot_latband_impact.py UTEND 5 $overlay
    python ./plot_latband_impact.py UTRESVEL 5 $overlay
    python ./plot_latband_impact.py utendepfd 5 $overlay
    python ./plot_latband_impact.py UTDIFF_NOGW 5 $overlay
    python ./plot_latband_impact.py UTTOTAL_NOGW 5 $overlay
    
fi

if [ $1 == 2 ]; then
    # tracer budgets time series
    python ./plot_latband_impact.py E90TEND 1 $overlay
    python ./plot_latband_impact.py QTRESVEL 1 $overlay e90
    python ./plot_latband_impact.py qtendetfd 1 $overlay e90
    python ./plot_latband_impact.py QTSRCSNK 1 $overlay e90
    python ./plot_latband_impact.py QTDIFF 1 $overlay e90
    python ./plot_latband_impact.py QTTOTAL 1 $overlay e90

    python ./plot_latband_impact.py AOATEND 1 $overlay
    python ./plot_latband_impact.py QTRESVEL 1 $overlay aoa
    python ./plot_latband_impact.py qtendetfd 1 $overlay aoa
    python ./plot_latband_impact.py QTSRCSNK 1 $overlay aoa
    python ./plot_latband_impact.py QTDIFF 1 $overlay aoa
    python ./plot_latband_impact.py QTTOTAL 1 $overlay aoa
    
    # tracer budgets time series
    python ./plot_latband_impact.py E90TEND 5 $overlay
    python ./plot_latband_impact.py QTRESVEL 5 $overlay e90
    python ./plot_latband_impact.py qtendetfd 5 $overlay e90
    python ./plot_latband_impact.py QTSRCSNK 5 $overlay e90
    python ./plot_latband_impact.py QTDIFF 5 $overlay e90
    python ./plot_latband_impact.py QTTOTAL 5 $overlay e90

    python ./plot_latband_impact.py AOATEND 5 $overlay
    python ./plot_latband_impact.py QTRESVEL 5 $overlay aoa
    python ./plot_latband_impact.py qtendetfd 5 $overlay aoa
    python ./plot_latband_impact.py QTSRCSNK 5 $overlay aoa
    python ./plot_latband_impact.py QTDIFF 5 $overlay aoa
    python ./plot_latband_impact.py QTTOTAL 5 $overlay aoa
fi

if [ $1 == 3 ]; then
    # tracer budget calendars
    python ./plot_monthly_impact.py E90TEND 1991 $overlay
    python ./plot_monthly_impact.py QTRESVEL 1991 $overlay none e90
    python ./plot_monthly_impact.py qtendetfd 1991 $overlay none e90
    python ./plot_monthly_impact.py QTSRCSNK 1991 $overlay none e90
    python ./plot_monthly_impact.py QTDIFF 1991 $overlay none e90
    python ./plot_monthly_impact.py QTTOTAL 1991 $overlay none e90

    python ./plot_monthly_impact.py AOATEND 1991 $overlay
    python ./plot_monthly_impact.py QTRESVEL 1991 $overlay none aoa
    python ./plot_monthly_impact.py qtendetfd 1991 $overlay none aoa
    python ./plot_monthly_impact.py QTSRCSNK 1991 $overlay none aoa
    python ./plot_monthly_impact.py QTDIFF 1991 $overlay none aoa
    python ./plot_monthly_impact.py QTTOTAL 1991 $overlay none aoa

    python ./plot_monthly_impact.py E90TEND 1992 $overlay
    python ./plot_monthly_impact.py QTRESVEL 1992 $overlay none e90
    python ./plot_monthly_impact.py qtendetfd 1992 $overlay none e90
    python ./plot_monthly_impact.py QTSRCSNK 1992 $overlay none e90
    python ./plot_monthly_impact.py QTDIFF 1992 $overlay none e90
    python ./plot_monthly_impact.py QTTOTAL 1992 $overlay none e90

    python ./plot_monthly_impact.py AOATEND 1992 $overlay
    python ./plot_monthly_impact.py QTRESVEL 1992 $overlay none aoa
    python ./plot_monthly_impact.py qtendetfd 1992 $overlay none aoa
    python ./plot_monthly_impact.py QTSRCSNK 1992 $overlay none aoa
    python ./plot_monthly_impact.py QTDIFF 1992 $overlay none aoa
    python ./plot_monthly_impact.py QTTOTAL 1992 $overlay none aoa
fi

if [ $1 == 4 ]; then
    # calendars
    python ./plot_monthly_impact.py U 1991 $overlay
    python ./plot_monthly_impact.py T 1991 $overlay
    python ./plot_monthly_impact.py AOA 1991 $overlay
    python ./plot_monthly_impact.py E90j 1991 $overlay
    python ./plot_monthly_impact.py psitem 1991 $overlay

    python ./plot_monthly_impact.py UTEND 1991 $overlay
    python ./plot_monthly_impact.py UTRESVEL 1991 $overlay
    python ./plot_monthly_impact.py utendepfd 1991 $overlay
    python ./plot_monthly_impact.py UTGWTOTAL 1991 $overlay
    python ./plot_monthly_impact.py UTDIFF 1991 $overlay
    python ./plot_monthly_impact.py UTTOTAL 1991 $overlay

    python ./plot_monthly_impact.py U 1992 $overlay
    python ./plot_monthly_impact.py T 1992 $overlay
    python ./plot_monthly_impact.py AOA 1992 $overlay
    python ./plot_monthly_impact.py E90j 1992 $overlay
    python ./plot_monthly_impact.py psitem 1992 $overlay

    python ./plot_monthly_impact.py UTEND 1992 $overlay
    python ./plot_monthly_impact.py UTRESVEL 1992 $overlay
    python ./plot_monthly_impact.py utendepfd 1992 $overlay
    python ./plot_monthly_impact.py UTGWTOTAL 1992 $overlay
    python ./plot_monthly_impact.py UTDIFF 1992 $overlay
    python ./plot_monthly_impact.py UTTOTAL 1992 $overlay
fi
exit

# October
python ./plot_monthly_impact.py U 1991 $overlay 10
python ./plot_monthly_impact.py T 1991 $overlay 10
python ./plot_monthly_impact.py AOA 1991 $overlay 10
python ./plot_monthly_impact.py E90j 1991 $overlay 10
python ./plot_monthly_impact.py psitem 1991 $overlay 10

python ./plot_monthly_impact.py UTEND 1991 $overlay 10
python ./plot_monthly_impact.py UTRESVEL 1991 $overlay 10
python ./plot_monthly_impact.py utendepfd 1991 $overlay 10
python ./plot_monthly_impact.py UTGWTOTAL 1991 $overlay 10
python ./plot_monthly_impact.py UTDIFF 1991 $overlay 10
python ./plot_monthly_impact.py UTTOTAL 1991 $overlay 10


# tracer budgets october
python ./plot_monthly_impact.py E90TEND 1991 $overlay 10 e90
python ./plot_monthly_impact.py QTRESVEL 1991 $overlay 10 e90
python ./plot_monthly_impact.py qtendetfd 1991 $overlay 10 e90
python ./plot_monthly_impact.py QTSRCSNK 1991 $overlay 10 e90
python ./plot_monthly_impact.py QTDIFF 1991 $overlay 10 e90
python ./plot_monthly_impact.py QTTOTAL 1991 $overlay 10 e90

python ./plot_monthly_impact.py AOATEND 1991 $overlay 10 aoa
python ./plot_monthly_impact.py UTRESVEL 1991 $overlay 10 aoa
python ./plot_monthly_impact.py qtendetfd 1991 $overlay 10 aoa
python ./plot_monthly_impact.py QTSRCSNK 1991 $overlay 10 e90
python ./plot_monthly_impact.py QTDIFF 1991 $overlay 10 aoa
python ./plot_monthly_impact.py QTTOTAL 1991 $overlay 10 aoa
