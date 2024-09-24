#!/bin/bash

overlay=0

python ./plot_EPFlux_monthly.py UTRESVEL 1992 0 2 none U none vtem wtem
#python ./plot_EPFlux_monthly.py utendepfd 1992 0 2 none U none epfy epfz

# 5-year time series in U, T
#python ./plot_latband_impact.py U 5 $overlay
#python ./plot_latband_impact.py T 5 $overlay

