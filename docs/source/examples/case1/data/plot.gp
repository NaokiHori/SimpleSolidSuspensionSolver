reset

set terminal epslatex standalone color size 5.,3.5 font ',17.28'
set output 'result.tex'

set xlabel '$y$'
set ylabel '$x$'

set xrange [0:30]
set yrange [0.26:0.42]

set xtics 0,10,30
set ytics 0.26,0.04,0.42

set format x '$% .0f$'
set format y '$% .2f$'

set style line 1 lc rgb '#FF0000'
set style line 2 lc rgb '#000000'

array filenames[2] = [ \
  'artifacts/particle.dat', \
  'PanGlowinski2002.dat', \
]

array titles[2] = [ \
  'result', \
  'reference', \
]

set style line 3 lc rgb '#AAAAAA' lw 3
set key right top spacing 1.2 box ls 3

plot \
  filenames[1] u 2:1 t titles[1] ls 1 lw 5           w l, \
  filenames[2] u 2:1 t titles[2] ls 2 lw 3 pt 6 ps 2 w p

