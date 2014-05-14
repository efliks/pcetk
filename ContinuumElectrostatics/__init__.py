#-------------------------------------------------------------------------------
# . File      : __init__.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""
This is a module that provides interface to the extended MEAD program and,
in the future, to the GMCT program.

It is similar in behaviour to the multiflex2qmpb.pl program by Timm Essigke.

For the calculated thermodynamic cycle, see Fig. 3.8, p. 82 in Timm's thesis.

For this module to be used, one should in the first step prepare CHARMM
topology (psf) and coordinate (crd) files for the protein of interest.

A pqr file is not needed because the radii are assigned to atoms at runtime.

All titratable residues in the protein should be set to their standard
protonation states at pH = 7, i.e. aspartates and glutamates deprotonated,
histidines doubly protonated, other residues protonated.
"""



from CEModelMEAD              import CEModelMEAD
                              
from StateVector              import StateVector

from GMCTOutputFileReader     import GMCTOutputFileReader
