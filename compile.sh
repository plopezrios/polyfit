#!/bin/bash
if [ "$1" = -d ] ; then
  opts="-Wall -Wextra -fimplicit-none -O0 -fbounds-check -g -pg -pedantic\
 -fbacktrace"
else
  opts="-O3 -ffast-math -fpeel-loops"
fi
gfortran $opts -o polyfit polyfit.f90 -llapack
