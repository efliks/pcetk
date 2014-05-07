#!/bin/sh
/home/mikolaj/local/bin/pqr2SolvAccVol epsin2_region > pqr2SolvAccVol_epsin2_region.out
echo -n 0%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS1\-0_p back > LYS1\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS1\-0_p LYS1-0_model_p_back2d > LYS1\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS1\-0_d back > LYS1\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS1\-0_d LYS1-0_model_d_back2d > LYS1\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site GLU7\-0_p back > GLU7\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 GLU7\-0_p GLU7-0_model_p_back2d > GLU7\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site GLU7\-0_d back > GLU7\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 GLU7\-0_d GLU7-0_model_d_back2d > GLU7\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS13\-0_p back > LYS13\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS13\-0_p LYS13-0_model_p_back2d > LYS13\-0_p_2d.out
echo -n 10%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS13\-0_d back > LYS13\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS13\-0_d LYS13-0_model_d_back2d > LYS13\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS15\-0_HSP back > HIS15\-0_HSP.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS15\-0_HSP HIS15-0_model_HSP_back2d > HIS15\-0_HSP_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS15\-0_HSD back > HIS15\-0_HSD.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS15\-0_HSD HIS15-0_model_HSD_back2d > HIS15\-0_HSD_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS15\-0_HSE back > HIS15\-0_HSE.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS15\-0_HSE HIS15-0_model_HSE_back2d > HIS15\-0_HSE_2d.out
echo -n 20%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS15\-0_full_deprot back > HIS15\-0_full_deprot.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS15\-0_full_deprot HIS15-0_model_full_deprot_back2d > HIS15\-0_full_deprot_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP18\-0_p back > ASP18\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP18\-0_p ASP18-0_model_p_back2d > ASP18\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP18\-0_d back > ASP18\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP18\-0_d ASP18-0_model_d_back2d > ASP18\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site TYR20\-0_p back > TYR20\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 TYR20\-0_p TYR20-0_model_p_back2d > TYR20\-0_p_2d.out
echo -n 30%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site TYR20\-0_d back > TYR20\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 TYR20\-0_d TYR20-0_model_d_back2d > TYR20\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site TYR23\-0_p back > TYR23\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 TYR23\-0_p TYR23-0_model_p_back2d > TYR23\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site TYR23\-0_d back > TYR23\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 TYR23\-0_d TYR23-0_model_d_back2d > TYR23\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS33\-0_p back > LYS33\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS33\-0_p LYS33-0_model_p_back2d > LYS33\-0_p_2d.out
echo -n 40%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS33\-0_d back > LYS33\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS33\-0_d LYS33-0_model_d_back2d > LYS33\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site GLU35\-0_p back > GLU35\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 GLU35\-0_p GLU35-0_model_p_back2d > GLU35\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site GLU35\-0_d back > GLU35\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 GLU35\-0_d GLU35-0_model_d_back2d > GLU35\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP48\-0_p back > ASP48\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP48\-0_p ASP48-0_model_p_back2d > ASP48\-0_p_2d.out
echo -n 50%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP48\-0_d back > ASP48\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP48\-0_d ASP48-0_model_d_back2d > ASP48\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP52\-0_p back > ASP52\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP52\-0_p ASP52-0_model_p_back2d > ASP52\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP52\-0_d back > ASP52\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP52\-0_d ASP52-0_model_d_back2d > ASP52\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site TYR53\-0_p back > TYR53\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 TYR53\-0_p TYR53-0_model_p_back2d > TYR53\-0_p_2d.out
echo -n 60%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site TYR53\-0_d back > TYR53\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 TYR53\-0_d TYR53-0_model_d_back2d > TYR53\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP66\-0_p back > ASP66\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP66\-0_p ASP66-0_model_p_back2d > ASP66\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP66\-0_d back > ASP66\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP66\-0_d ASP66-0_model_d_back2d > ASP66\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP87\-0_p back > ASP87\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP87\-0_p ASP87-0_model_p_back2d > ASP87\-0_p_2d.out
echo -n 70%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP87\-0_d back > ASP87\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP87\-0_d ASP87-0_model_d_back2d > ASP87\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS96\-0_p back > LYS96\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS96\-0_p LYS96-0_model_p_back2d > LYS96\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS96\-0_d back > LYS96\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS96\-0_d LYS96-0_model_d_back2d > LYS96\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS97\-0_p back > LYS97\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS97\-0_p LYS97-0_model_p_back2d > LYS97\-0_p_2d.out
echo -n 80%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS97\-0_d back > LYS97\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS97\-0_d LYS97-0_model_d_back2d > LYS97\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP101\-0_p back > ASP101\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP101\-0_p ASP101-0_model_p_back2d > ASP101\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP101\-0_d back > ASP101\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP101\-0_d ASP101-0_model_d_back2d > ASP101\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS116\-0_p back > LYS116\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS116\-0_p LYS116-0_model_p_back2d > LYS116\-0_p_2d.out
echo -n 90%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site LYS116\-0_d back > LYS116\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 LYS116\-0_d LYS116-0_model_d_back2d > LYS116\-0_d_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP119\-0_p back > ASP119\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP119\-0_p ASP119-0_model_p_back2d > ASP119\-0_p_2d.out
echo -n .
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site ASP119\-0_d back > ASP119\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 ASP119\-0_d ASP119-0_model_d_back2d > ASP119\-0_d_2d.out
echo -n .
echo DONE!
