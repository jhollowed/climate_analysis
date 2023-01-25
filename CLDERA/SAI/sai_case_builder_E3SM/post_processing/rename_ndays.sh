#!/bin/bash

# $1 = location of nc files
# make sure to run for directories first...

for file in $1/*; do mv -v "$file" "${file/900day/1200day}"; done
