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

set   style line  1  lc rgb "blue"     lw 2     lt 1
set   style line  2  lc rgb "black"    lw 2     lt 1
set   style line  3  lc rgb "green"    lw 2     lt 1
set   style line  4  lc rgb "red"      lw 2     lt 1
set   style line  5  lc rgb "cyan"     lw 2     lt 1
set   style line  6  lc rgb "magenta"  lw 2     lt 1
set   style line  7  lc rgb "yellow"   lw 2     lt 1

set   style line  8  lc rgb "blue"     lw 4     lt 2
set   style line  9  lc rgb "black"    lw 4     lt 2
set   style line 10  lc rgb "green"    lw 4     lt 2
set   style line 11  lc rgb "red"      lw 4     lt 2
set   style line 12  lc rgb "cyan"     lw 4     lt 2
set   style line 13  lc rgb "magenta"  lw 4     lt 2
set   style line 14  lc rgb "yellow"   lw 4     lt 2

set   multiplot layout 4, 4

unset key

set   key
set   size 0.25, 0.20
set   origin 0.000000, 0.750000
set   title "His260 (δ,ε-protonated+)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/PRTA_HIS260_HSP.dat" w l ls 1 t "MC" 


unset key



set   size 0.25, 0.20
set   origin 0.250000, 0.750000
set   title "His260 (ε-protonated)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/PRTA_HIS260_HSE.dat" w l ls 1 t "MC" 




set   size 0.25, 0.20
set   origin 0.500000, 0.750000
set   title "His260 (δ-protonated)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/PRTA_HIS260_HSD.dat" w l ls 1 t "MC" 




set   size 0.25, 0.20
set   origin 0.750000, 0.750000
set   title "BLV1 (full+)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/CHRO_BLF1_full.dat" w l ls 1 t "MC" 




set   size 0.25, 0.20
set   origin 0.000000, 0.550000
set   title "BLV1 (dd)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/CHRO_BLF1_ddepr.dat" w l ls 1 t "MC" 




set   size 0.25, 0.20
set   origin 0.250000, 0.550000
set   title "BLV1 (cd)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/CHRO_BLF1_cdepr.dat" w l ls 1 t "MC" 




set   size 0.25, 0.20
set   origin 0.500000, 0.550000
set   title "BLV1 (bd)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/CHRO_BLF1_bdepr.dat" w l ls 1 t "MC" 




set   size 0.25, 0.20
set   origin 0.750000, 0.550000
set   title "BLV1 (ad)" font "Verdana Bold, 10" offset 0, -0.5
plot \
"../curves/CHRO_BLF1_adepr.dat" w l ls 1 t "MC" 



unset multiplot
