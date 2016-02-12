#!/usr/bin/gnuplot -persist
set terminal postscript eps size 7,2.5 enhanced color \
             font 'Helvetica,20' linewidth 2
             
set xlabel "Anzahl der Iterationen"
set ylabel "Maximaler Fehler"
set logscale y 2

#set title "Maximaler Fehler nach Iteration"
set output "plots/fehlerprot.eps"
plot "fehlerJakobi.csv"      using 1:2 with lines title "Jakobi", \
     "fehlerGaussSeidel.csv" using 1:2 with lines title "Gauss-Seidel"
