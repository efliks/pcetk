#-------------------------------------------------------------------------------
# . File      : __init__.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""ContinuumElectrostatics is a module implementing a Poisson-Boltzmann continuum electrostatic model.

The model allows for predictions of protonation states of titratable residues in proteins."""

__lastchanged__ = "$Id$"

from Model                 import MEADModel
from StateVector           import StateVector
from Toolbox               import StateVector_FromProbabilities
