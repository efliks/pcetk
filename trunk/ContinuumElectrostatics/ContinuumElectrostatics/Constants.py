#-------------------------------------------------------------------------------
# . File      : Constants.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Constants."""


CONSTANT_MOLAR_GAS_KCAL_MOL = 1.9871653920000e-03  #1.9872041 / 1000.0
CONSTANT_LN10               = 2.302585092994046

# Only valid for T=300 K
UNITS_ENERGY_PKA_UNITS_TO_KILOCALORIES_PER_MOL = 1.0 / (CONSTANT_MOLAR_GAS_KCAL_MOL * 300.0 * CONSTANT_LN10)
