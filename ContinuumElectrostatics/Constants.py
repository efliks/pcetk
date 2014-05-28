#-------------------------------------------------------------------------------
# . File      : Constants.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Constants."""

import  os


CONSTANT_MOLAR_GAS_KCAL_MOL = 0.001987165392 # Gas constant expressed in kcal/(K*mol)
CONSTANT_LN10               = 2.302585092994

# Only valid for T=300 K
UNITS_ENERGY_PKA_UNITS_TO_KILOCALORIES_PER_MOL = 1.0 / (CONSTANT_MOLAR_GAS_KCAL_MOL * 300.0 * CONSTANT_LN10)


YAMLPATHIN       = os.path.join (os.getenv ("PDYNAMO_CONTINUUMELECTROSTATICS"), "parameters")

MAX_SITES        = 20

PREV_RESIDUE     = ("C", "O")

NEXT_RESIDUE     = ("N", "H",  "CA", "HA")

NEXT_RESIDUE_PRO = ("N", "CA", "HA", "CD",  "HD1", "HD2")

NEXT_RESIDUE_GLY = ("N", "H",  "CA", "HA1", "HA2")

PROTEIN_RESIDUES = (
   "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
   "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
   "HSP", "HSE", "HSD", "HIE", "HID"
)

REMOVE_RESIDUES  = (
    "WAT", "HOH", "TIP", "TIP3", "TP3M", "SOD",
)
