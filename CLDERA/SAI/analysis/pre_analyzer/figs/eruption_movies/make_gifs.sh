#!/bin/bash

dir=$1
cd $dir
convert -delay 7 -loop 0 *.png anim.gif
gifsicle -i anim.gif -O3 --colors 64 -o anim-opt.gif
