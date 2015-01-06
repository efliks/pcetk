#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pxd
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions     cimport Boolean, CFalse, CTrue, Integer, Real
from pCore.Status           cimport Status, Status_Continue, Status_IndexOutOfRange, Status_ValueError

__lastchanged__ = "$Id$"


# Include StateVector.h in the generated C code
cdef extern from "StateVector.h":
  ctypedef struct CStateVector "StateVector":
    Integer   length
    Integer   slength
    Integer  *vector
    Integer  *minvector
    Integer  *maxvector
    Integer  *substate

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
  cdef void          StateVector_ResetSubstate     (CStateVector *self)
  cdef Boolean       StateVector_AllocateSubstate  (CStateVector *self, Integer nsites, Status *status)
  cdef Boolean       StateVector_IncrementSubstate (CStateVector *self)
  cdef Boolean       StateVector_SetSubstateItem   (CStateVector *self, Integer selectedSiteIndex, Integer index, Status *status)
  cdef Integer       StateVector_GetSubstateItem   (CStateVector *self, Integer index, Status *status)


#-------------------------------------------------------------------------------
cdef class StateVector:
  cdef CStateVector     *cObject
  cdef public object  isOwner
  cdef public object  owner
