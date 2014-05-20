#!/bin/sh
/home/mikolaj/local/bin/pqr2SolvAccVol epsin2_region > pqr2SolvAccVol_epsin2_region.out
echo -n 0%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS5\-0_HSP back > HIS5\-0_HSP.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS5\-0_HSP HIS5-0_model_HSP_back2d > HIS5\-0_HSP_2d.out
echo -n 10%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS5\-0_HSD back > HIS5\-0_HSD.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS5\-0_HSD HIS5-0_model_HSD_back2d > HIS5\-0_HSD_2d.out
echo -n 20%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS5\-0_HSE back > HIS5\-0_HSE.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS5\-0_HSE HIS5-0_model_HSE_back2d > HIS5\-0_HSE_2d.out
echo -n 30%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site HIS5\-0_full_deprot back > HIS5\-0_full_deprot.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 HIS5\-0_full_deprot HIS5-0_model_full_deprot_back2d > HIS5\-0_full_deprot_2d.out
echo -n 40%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site GLU8\-0_p back > GLU8\-0_p.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 GLU8\-0_p GLU8-0_model_p_back2d > GLU8\-0_p_2d.out
echo -n 50%
/home/mikolaj/local/bin/my_3diel_solver -eps2set epsin2_region -T 300 -ionicstr 0.1 -epsin1 1 -epsin2 4 -fpt site GLU8\-0_d back > GLU8\-0_d.out
/home/mikolaj/local/bin/my_2diel_solver -epsin 4 -T 300 -ionicstr 0.1 GLU8\-0_d GLU8-0_model_d_back2d > GLU8\-0_d_2d.out
echo -n 60%
echo DONE!
