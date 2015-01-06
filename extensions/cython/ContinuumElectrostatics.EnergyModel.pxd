#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.EnergyModel.pxd
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions                  cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status                        cimport Status, Status_Continue, Status_IndexOutOfRange, Status_ValueError
from pCore.Integer1DArray                cimport CInteger1DArray, Integer1DArray
from pCore.Real1DArray                   cimport CReal1DArray, Real1DArray
from pCore.Real2DArray                   cimport CReal2DArray, Real2DArray
from pCore.SymmetricMatrix               cimport CSymmetricMatrix, SymmetricMatrix
from ContinuumElectrostatics.StateVector cimport CStateVector, StateVector

__lastchanged__ = "$Id: $"


# Include EnergyModel.h in the generated C code
cdef extern from "EnergyModel.h":
  ctypedef struct CEnergyModel "EnergyModel":
    CInteger1DArray    *protons
    CReal1DArray       *intrinsic
    CReal2DArray       *interactions
    CReal1DArray       *probabilities
    CSymmetricMatrix   *symmetricmatrix
    Integer             nstates
    Integer             ninstances

  cdef CEnergyModel *EnergyModel_Allocate   (Integer ninstances, Status *status)
  cdef void          EnergyModel_Deallocate (CEnergyModel *self)
  cdef Boolean       EnergyModel_CheckInteractionsSymmetric (CEnergyModel *self, Real threshold, Real *maxDeviate)
  cdef Boolean       EnergyModel_SymmetrizeInteractions (CEnergyModel *self, Status *status)
  cdef Real          EnergyModel_CalculateMicrostateEnergy (CEnergyModel *self, CStateVector *vector, Real pH, Real temperature)
  cdef Boolean       EnergyModel_CalculateProbabilitiesAnalytically (CEnergyModel *self, CStateVector *vector, Real pH, Real temperature, Status *status)

  # Functions for accessing items
  cdef Real          EnergyModel_GetGintr                (CEnergyModel *self, Integer instIndexGlobal)
  cdef Integer       EnergyModel_GetProtons              (CEnergyModel *self, Integer instIndexGlobal)
  cdef Real          EnergyModel_GetProbability          (CEnergyModel *self, Integer instIndexGlobal)
  cdef Real          EnergyModel_GetInteraction          (CEnergyModel *self, Integer instIndexGlobalA, Integer instIndexGlobalB)
  cdef Real          EnergyModel_GetInteractionSymmetric (CEnergyModel *self, Integer instIndexGlobalA, Integer instIndexGlobalB)
  cdef Real          EnergyModel_GetDeviation            (CEnergyModel *self, Integer instIndexGlobalA, Integer instIndexGlobalB)

  cdef void          EnergyModel_SetGintr       (CEnergyModel *self, Integer instIndexGlobal, Real value)
  cdef void          EnergyModel_SetProtons     (CEnergyModel *self, Integer instIndexGlobal, Integer value)
  cdef void          EnergyModel_SetProbability (CEnergyModel *self, Integer instIndexGlobal, Real value)
  cdef void          EnergyModel_SetInteraction (CEnergyModel *self, Integer instIndexGlobalA, Integer instIndexGlobalB, Real value)


#-------------------------------------------------------------------------------
cdef class EnergyModel:
  cdef CEnergyModel  *cObject
  cdef public object  isOwner
  cdef public object  owner
