#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.MCModelDefault.pxd
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                  cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status                        cimport Status, Status_Continue, Status_IndexOutOfRange, Status_ValueError
from pCore.Real1DArray                   cimport Real1DArray_Set, Real1DArray_Scale
from pCore.RandomNumberGenerator         cimport CRandomNumberGenerator as CGenerator
from ContinuumElectrostatics.StateVector cimport CStateVector, StateVector, StateVector_GetPair, StateVector_Randomize
from ContinuumElectrostatics.EnergyModel cimport CEnergyModel, EnergyModel

__lastchanged__ = "$Id: $"


cdef extern from "MCModelDefault.h":
    ctypedef struct CMCModelDefault "MCModelDefault":
        Real           limit
        Integer        nprod
        Integer        nequil
        CEnergyModel  *energyModel
        CStateVector  *vector
        CGenerator    *generator


    cdef CMCModelDefault *MCModelDefault_Allocate               (Real limit, Integer nequil, Integer nprod, Integer randomSeed, Status *status)
    cdef void             MCModelDefault_Deallocate             (CMCModelDefault *self)
    cdef void             MCModelDefault_LinkToEnergyModel      (CMCModelDefault *self, CEnergyModel *energyModel, Status *status)
    cdef Real             MCModelDefault_MCScan                 (CMCModelDefault *self, Real pH, Integer nmoves, Integer *movesDone, Integer *movesAccepted, Integer *flipsDone, Integer *flipsAccepted)
    cdef Integer          MCModelDefault_FindPairs              (CMCModelDefault *self, Integer npairs, Status *status)
    cdef void             MCModelDefault_UpdateProbabilities    (CMCModelDefault *self)
    cdef void             MCModelDefault_CalculateProbabilities (CMCModelDefault *self, Real pH, Boolean equil)

#-------------------------------------------------------------------------------
cdef class MCModelDefault:
    cdef CMCModelDefault  *cObject
    cdef public object  isOwner
    cdef public object  owner
