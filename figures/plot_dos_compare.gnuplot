#!/usr/bin/env gnuplot

set terminal unknown

plot '<(gunzip --stdout ../HAT-Cu-AgSlab5/2b.dos/DOS.SUM.HAT.gz)' u 1:5 w l t 'On Ag(111) 5 slab',\
     '<(gunzip --stdout ../HAT-Cu-AgSlab/3b.dos/DOS.SUM.HAT.gz)' u 1:5 w l t 'On Ag(111) 3 slab',\
     '<(gunzip --stdout ../HAT-Cu_free/3c.dos/DOS.SUM.HAT.gz)' u 1:5 w l t 'Freestanding'

set xlabel 'Energy (eV)'
set ylabel 'HAT DOS'
set xrange [-1:2.5]
set yrange [0:30]
set xzeroaxis

set terminal png
set output 'HAT_DOS_slab_comparison.png'
replot
unset output
