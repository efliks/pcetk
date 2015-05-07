# Gnuplot input file
set   terminal pdf dashed size 8.3, 11.7 font "Verdana, 10" linewidth 1.5
set   output "page01.pdf"

set   xrange [0:14]
set   yrange [0:1]
set   xtics  2.0 scale 0.8 offset 0,   0.3
set   mxtics
set   ytics  0.2 scale 0.8 offset 0.3, 0
set   mytics
set   xlabel "pH"           offset 0.0, 0.8  font "Verdana, 10"
set   ylabel "Probability"  offset 1.8, 0.0  font "Verdana, 10"

set   style line  1  lc rgb "black"    lw 2     lt 1
set   style line  2  lc rgb "blue"     lw 2     lt 1
set   style line  3  lc rgb "green"    lw 2     lt 1
set   style line  4  lc rgb "red"      lw 2     lt 1
set   style line  5  lc rgb "cyan"     lw 2     lt 1
set   style line  6  lc rgb "magenta"  lw 2     lt 1
set   style line  7  lc rgb "yellow"   lw 2     lt 1

set   style line  8  lc rgb "black"    lw 4     lt 2
set   style line  9  lc rgb "blue"     lw 4     lt 2
set   style line 10  lc rgb "green"    lw 4     lt 2
set   style line 11  lc rgb "red"      lw 4     lt 2
set   style line 12  lc rgb "cyan"     lw 4     lt 2
set   style line 13  lc rgb "magenta"  lw 4     lt 2
set   style line 14  lc rgb "yellow"   lw 4     lt 2

set   multiplot layout 4, 2

unset key

set   key
set   size 0.25, 0.20
set   origin 0.000000, 0.750000
set   title "Asp52 (0)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves_analytic/PRTA_ASP52_p.dat" w l ls 1 t "analytic", \
"../curves_mc/PRTA_ASP52_p.dat" w l ls 2 t "MC" 


unset key



set   size 0.25, 0.20
set   origin 0.250000, 0.750000
set   title "Asp52 (−)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves_analytic/PRTA_ASP52_d.dat" w l ls 1 t "analytic", \
"../curves_mc/PRTA_ASP52_d.dat" w l ls 2 t "MC" 




set   size 0.25, 0.20
set   origin 0.500000, 0.750000
set   title "Glu35 (0)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves_analytic/PRTA_GLU35_p.dat" w l ls 1 t "analytic", \
"../curves_mc/PRTA_GLU35_p.dat" w l ls 2 t "MC" 




set   size 0.25, 0.20
set   origin 0.750000, 0.750000
set   title "Glu35 (−)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves_analytic/PRTA_GLU35_d.dat" w l ls 1 t "analytic", \
"../curves_mc/PRTA_GLU35_d.dat" w l ls 2 t "MC" 



unset multiplot
