#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pxd
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions     cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Integer1DArray   cimport CInteger1DArray, Integer1DArray
from pCore.Real1DArray      cimport CReal1DArray, Real1DArray
from pCore.Real2DArray      cimport CReal2DArray, Real2DArray
from pCore.SymmetricMatrix  cimport CSymmetricMatrix, SymmetricMatrix
from pCore.Status           cimport Status, Status_Continue, Status_IndexOutOfRange, Status_ValueError

__lastchanged__ = "$Id$"


# Include StateVector.h in the generated C code
cdef extern from "StateVector.h":
  ctypedef struct CStateVector "StateVector":
    Integer *vector
    Integer *minvector
    Integer *maxvector
    Integer *substate
    Integer length
    Integer slength


  cdef CStateVector *StateVector_Allocate          (Integer size, Status *status)
  cdef CStateVector *StateVector_Clone             (CStateVector *self, Status *status)
  cdef Boolean       StateVector_CopyTo            (CStateVector *self, CStateVector *other, Status *status)
  cdef void          StateVector_Deallocate        (CStateVector *self)
  cdef void          StateVector_Reset             (CStateVector *self)
  cdef void          StateVector_ResetToMaximum    (CStateVector *self)
  cdef Integer       StateVector_GetItem           (CStateVector *self, Integer index, Status *status)
  cdef Boolean       StateVector_SetItem           (CStateVector *self, Integer index, Integer value, Status *status)
  cdef Integer       StateVector_GetActualItem     (CStateVector *self, Integer index, Status *status)
  cdef Boolean       StateVector_SetActualItem     (CStateVector *self, Integer index, Integer value, Status *status)
  cdef Boolean       StateVector_Increment         (CStateVector *self)

  # Substate-related functions
  cdef void         StateVector_ResetSubstate      (CStateVector *self)
  cdef Boolean      StateVector_AllocateSubstate   (CStateVector *self, Integer nsites, Status *status)
  cdef Boolean      StateVector_IncrementSubstate  (CStateVector *self)
  cdef Boolean      StateVector_SetSubstateItem    (CStateVector *self, Integer selectedSiteIndex, Integer index, Status *status)
  cdef Integer      StateVector_GetSubstateItem    (CStateVector *self, Integer index, Status *status)

  # Calculating microstate energy
  cdef Real         StateVector_CalculateMicrostateEnergy (CStateVector *self, CInteger1DArray *protons, CReal1DArray *intrinsic, CSymmetricMatrix *symmetricmatrix, Real pH, Real temperature)

  # Calculating probabilities of protonation states analytically
  cdef Boolean      StateVector_CalculateProbabilitiesAnalytically (CStateVector *self, CInteger1DArray *protons, CReal1DArray *intrinsic, CSymmetricMatrix *symmetricmatrix, Real pH, Real temperature, Integer nstates, CReal1DArray *probabilities, Status *status)


#-------------------------------------------------------------------------------
cdef class StateVector:
  cdef CStateVector     *cObject
