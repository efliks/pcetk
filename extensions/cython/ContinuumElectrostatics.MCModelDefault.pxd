#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.MCModelDefault.pxd
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                  cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status                        cimport Status, Status_Continue, Status_IndexOutOfRange, Status_ValueError
from ContinuumElectrostatics.StateVector cimport CStateVector, StateVector, StateVector_GetPair
from ContinuumElectrostatics.EnergyModel cimport CEnergyModel, EnergyModel

__lastchanged__ = "$Id: $"


cdef extern from "MCModelDefault.h":
    ctypedef struct CMCModelDefault "MCModelDefault":
        Real           limit
        Integer        nprod
        Integer        nequil
        CEnergyModel  *energyModel
        CStateVector  *vector

    cdef CMCModelDefault *MCModelDefault_Allocate               (Real limit, Integer nequil, Integer nprod, Integer randomSeed, Status *status)
    cdef void             MCModelDefault_Deallocate             (CMCModelDefault *self)
    cdef void             MCModelDefault_LinkToEnergyModel      (CMCModelDefault *self, CEnergyModel *energyModel, Status *status)
    cdef Real             MCModelDefault_MCScan                 (CMCModelDefault *self, Real pH, Integer nmoves)
    cdef Integer          MCModelDefault_FindPairs              (CMCModelDefault *self, Integer npairs, Status *status)
    cdef void             MCModelDefault_CalculateProbabilities (CMCModelDefault *self, Real pH, Boolean equil)

#-------------------------------------------------------------------------------
cdef class MCModelDefault:
    cdef CMCModelDefault  *cObject
    cdef public object  isOwner
    cdef public object  owner
