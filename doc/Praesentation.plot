#!/usr/bin/gnuplot -persist
set terminal postscript eps size 5,3 enhanced color \
             font 'Helvetica,20' linewidth 2
             
# Fehler
set xlabel "Anzahl der Iterationen"
set ylabel "Maximaler Fehler"
set logscale y 2

#set title "Maximaler Fehler nach Iteration"
set output "plots/fehlerpres.eps"
plot "fehlerJakobi.csv"      using 1:2 with lines title "Jakobi", \
     "fehlerGaussSeidel.csv" using 1:2 with lines title "Gauss-Seidel"

# Benchmarks
set xlabel "Aufloesung des Gitters"
unset ylabel
unset logscale y
set logscale x 2
set yrange [0:1]

set output "plots/effizienzjakobi.eps"
plot "benchmark-jakobi.csv" using 1:5 with lines title "Effizienz"

set output "plots/effizienzgaussseidel.eps"
plot "benchmark-gauss-seidel.csv" using 1:5 with lines title "Effizienz"
