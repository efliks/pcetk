#               with MP2/CBS correction                                      without MP2/CBS correction
#   pKa_exp  pKa_calc  pKa_lfit     error                                   pKa_calc  pKa_lfit     error
#       1.5      -8.8       1.5      -0.0    "2,5-dichloro-anilinium"         -8.5       1.3      -0.2
#       1.6      -8.0       2.0       0.4    "4-cyano-anilinium"              -8.4       1.4      -0.2
#       3.9      -5.7       3.5      -0.4    "4-bromo-anilinium"              -5.2       3.3      -0.6
#       4.6      -4.2       4.3      -0.3    "anilinium"                      -3.8       4.1      -0.5
#       5.4      -2.9       5.2      -0.2    "p-anisidinium"                  -1.3       5.6       0.3
#       6.7      -0.1       7.0       0.3    "2,6-dimethyl-pyridinium"         1.9       7.5       0.8
#       7.3       0.6       7.4       0.0    "2,4,6-collidinium"               3.6       8.5       1.2
#       9.3       3.6       9.3      -0.0    "benzylammonium"                  4.7       9.2      -0.1
#       9.6       4.3       9.7       0.1    "DMAP"                            6.8      10.4       0.8
#      10.4       5.3      10.3      -0.1    "lysine"                          6.1      10.1      -0.3
#      10.7       6.6      11.1       0.4    "triethylammonium"                6.5      10.2      -0.5
#      11.3       6.8      11.3      -0.0    "pyrrolidinium"                   7.4      10.8      -0.4
#      13.6      10.3      13.4      -0.2    "guanidinium"                    11.6      13.3      -0.3
set   terminal pdf size 5.0, 3.0 dashed font "Verdana, 8" linewidth 1.0
set   termopt enhanced
set   output "combined_fits_two_plots.pdf"

set   style line  1  lc rgb "blue"    lw 3.0   lt 1    #  lw 2.0   lt 1
set   style line  2  lc rgb "blue"    lw 2.5   lt 2    #  lw 1.5   lt 2
set   style line  3  lc rgb "green"   lw 3.0   lt 1    #  lw 2.0   lt 1
set   style line  4  lc rgb "green"   lw 2.5   lt 2    #  lw 1.5   lt 2

set   style line  5  lc rgb "blue"    lw 2.0   lt 8    #  lw 1.0   lt 8
set   style line  6  lc rgb "blue"    lw 2.0   lt 8    #  lw 1.0   lt 8
set   style line  7  lc rgb "green"   lw 2.0   lt 8    #  lw 1.0   lt 8
set   style line  8  lc rgb "green"   lw 2.0   lt 8    #  lw 1.0   lt 8


set   xrange[0:14]
set   yrange[-10:15]
set   xlabel "pK^{aq}_{a,exp}"  offset  0.0, 0.5
set   ylabel "pK^{aq}_{a,calc}" offset  1.5, 0.0
set   mxtics
set   mytics
set   key left top
set   multiplot layout 2, 1


#========================================================
# set fit logfile "/dev/null"
# f(x) = d + c*x
# fit f(x) "combined_fits.dat" u 1:2 via c, d
#========================================================
f(x) = 1.60196*x -11.2175
set   size 0.5, 1.0
set   origin 0.0, 0.0
set   title "a)" font "Verdana, 10"
plot  \
      "combined_fits.dat"    u 1:2 w p ls 1 t "pK^{aq}_{a,calc}" ,  \
                                    f(x) w l ls 2 notitle, \
      "combined_fits.dat"    u 1:3 w p ls 3 t "pK^{aq}_{a,fit}"  ,  \
                                    x    w l ls 4 notitle, \
      "combined_fits.dat"    u 1:2:5 w labels font "Verdana Italic, 8" rotate by 0 offset 0.75,-0.25 left notitle


#========================================================
# set fit logfile "/dev/null"
# g(x) = b + a*x
# fit g(x) "combined_fits.dat" u 1:6 via a, b
#========================================================
g(x) = 1.67712*x -10.7258
set   size 0.5, 1.0
set   origin 0.5, 0.0
unset ylabe
set   title "b)" font "Verdana, 10"
plot  \
      "combined_fits.dat"    u 1:6 w p ls 5 t "pK^{aq}_{a,calc}" ,  \
                                    g(x) w l ls 2 notitle, \
      "combined_fits.dat"    u 1:7 w p ls 7 t "pK^{aq}_{a,fit}"  ,  \
                                    x    w l ls 4 notitle, \
      "combined_fits.dat"    u 1:6:5 w labels font "Verdana Italic, 8" rotate by 0 offset 0.75,-0.25 left notitle

unset multiplot
