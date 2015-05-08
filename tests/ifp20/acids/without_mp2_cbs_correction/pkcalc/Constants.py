#-------------------------------------------------------------------------------
# . File      : Constants.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
CONSTANT_MOLAR_GAS_KCAL_MOL =    0.001987165392 # Gas constant (R) expressed in kcal/(K*mol)
CONSTANT_LN10               =    2.302585092994
G_PROTON_AQ                 = -272.2
TEMPERATURE                 =  298.15

UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE              = 2625.5
UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE = 4.184

UNITS_ENERGY_HARTREES_TO_KILOCALORIES_PER_MOLE            = UNITS_ENERGY_HARTREES_TO_KILOJOULES_PER_MOLE / UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE

DEFAULT_A = 0.623742
DEFAULT_B = 6.990225


pKasExperimental = {
"2_chloro_pyridinium"    :    0.49   ,  # <--OUT
"4_cyano_pyridinium"     :    1.86   ,  # <--OUT
"4_bromo_anilinium"      :    3.89   ,
"anilinium"              :    4.62   ,
"pyridinium"             :    5.24   ,  # <--OUT
"246_collidinium"        :    7.33   ,
"benzylammonium"         :    9.30   ,
"triethylammonium"       :   10.72   ,
"pyrrolidinium"          :   11.27   ,
"guanidinium"            :   13.60   ,
# Previously excluded
"25_dichloro_anilinium"  :    1.53   ,
"p-anisidinium"          :    5.36   ,
"26_dimethyl_pyridinium" :    6.70   ,
"DMAP"                   :    9.60   ,
# Own compounds
"4_cyano_anilinium"      :    1.57   ,
"benzoquinoline"         :    5.05   ,  #<--To do
"acridine"               :    5.60   ,  #<--To do
"lysine"                 :   10.40   ,
"arginine"               :   12.00   ,  #<--OUT
}
