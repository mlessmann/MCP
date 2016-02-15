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
set ylabel "Effizienz"
unset logscale y
set logscale x 2
set yrange [0:1]

set output "plots/effizienzjakobi.eps"
plot "benchmark-jakobi.csv" using 1:5 with lines title "Effizienz"

set output "plots/effizienzgaussseidel.eps"
plot "benchmark-gauss-seidel.csv" using 1:5 with lines title "Effizienz"

set output "plots/effizienzmehrgitter.eps"
plot "benchmark-mehrgitter.csv" every 96::1  using 1:9 with lines title "z=4,  hmax=4 ", \
     "benchmark-mehrgitter.csv" every 96::3  using 1:9 with lines title "z=4,  hmax=16", \
     "benchmark-mehrgitter.csv" every 96::31 using 1:9 with lines title "z=16, hmax=4 ", \
     "benchmark-mehrgitter.csv" every 96::33 using 1:9 with lines title "z=16, hmax=16"
#     "benchmark-mehrgitter.csv" every 96::49 using 1:9 with lines title "alpha=2, z=4,  hmax=4 ", \
#     "benchmark-mehrgitter.csv" every 96::81 using 1:9 with lines title "alpha=2, z=16, hmax=16"

set output "plots/effizienz.eps"
plot "benchmark-jakobi.csv"       using 1:5 with lines title "Jakobi", \
     "benchmark-gauss-seidel.csv" using 1:5 with lines title "Gauss-Seidel", \
     "benchmark-mehrgitter.csv"   every 96::31 using 1:9 with lines title "Mehrgitter"

# Laufzeiten
set ylabel "Laufzeit in Sekunden"
set yrange [*:*]
set logscale y 2

set output "plots/laufzeitensequenziell.eps"
plot "benchmark-jakobi.csv"       using 1:2 with lines title "Jakobi, sequenziell", \
     "benchmark-gauss-seidel.csv" using 1:2 with lines title "Gauss-Seidel, sequenziell"

set output "plots/laufzeitenparallel.eps"
plot "benchmark-jakobi.csv"       using 1:3 with lines title "Jakobi, parallel", \
     "benchmark-gauss-seidel.csv" using 1:3 with lines title "Gauss-Seidel, parallel"

# Mehrgitterfehler
set ylabel "Fehler"
unset logscale y

set output "plots/fehlermehrgitter.eps"
plot "benchmark-mehrgitter.csv" every 96::3  using 1:10 with lines title "Mittlerer Fehler, alpha=1", \
     "benchmark-mehrgitter.csv" every 96::51 using 1:10 with lines title "Mittlerer Fehler, alpha=2", \
     "benchmark-mehrgitter.csv" every 96::3  using 1:11 with lines title "Maximaler Fehler, alpha=1", \
     "benchmark-mehrgitter.csv" every 96::51 using 1:11 with lines title "Maximaler Fehler, alpha=2"
