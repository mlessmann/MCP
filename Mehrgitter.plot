# Nach Datei
set terminal postscript eps size 5,3 enhanced color \
             font 'Helvetica,20' linewidth 2

# Plot
set grid
set zrange [-0.5:1.5]
#set output "Mehrgitter.eps"
#splot "*.matrix" matrix using 1:2:3 with lines

set output "000.eps"
set xrange [0:33]
set yrange [0:33]
splot "000 - Mehrgitter Down, vor Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Down, vor Gauss-Seidel"

set output "001.eps"
set xrange [0:33]
set yrange [0:33]
splot "001 - Mehrgitter Down, nach Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Down, nach Gauss-Seidel"

set output "002.eps"
set xrange [0:17]
set yrange [0:17]
splot "002 - Mehrgitter Down, vor Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Down, vor Gauss-Seidel"

set output "003.eps"
set xrange [0:17]
set yrange [0:17]
splot "003 - Mehrgitter Down, nach Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Down, nach Gauss-Seidel"

set output "004.eps"
set xrange [0:9]
set yrange [0:9]
splot "004 - Mehrgitter Main, vor Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Main, vor Gauss-Seidel"

set output "005.eps"
set xrange [0:9]
set yrange [0:9]
splot "005 - Mehrgitter Main, nach Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Main, nach Gauss-Seidel"

set output "006.eps"
set xrange [0:17]
set yrange [0:17]
splot "006 - Mehrgitter Up, vor Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Up, vor Gauss-Seidel"

set output "007.eps"
set xrange [0:17]
set yrange [0:17]
splot "007 - Mehrgitter Up, nach Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Up, nach Gauss-Seidel"

set output "008.eps"
set xrange [0:33]
set yrange [0:33]
splot "008 - Mehrgitter Up, vor Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Up, vor Gauss-Seidel"

set output "009.eps"
set xrange [0:33]
set yrange [0:33]
splot "009 - Mehrgitter Up, nach Gauss-Seidel.matrix" matrix using 1:2:3 with lines title "Up, nach Gauss-Seidel"
