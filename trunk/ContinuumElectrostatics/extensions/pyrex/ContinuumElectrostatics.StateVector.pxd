#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pxd
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore.cDefinitions    cimport Boolean, CFalse, CTrue, Integer


# Include the file StateVector.h in the generated C code
cdef extern from "StateVector.h":

  ctypedef struct CStateVector "StateVector":
    Integer *vector
    Integer *maxvector
    Integer *substate
    Integer length
    Integer slength


  cdef CStateVector *StateVector_Allocate          (Integer size)
  cdef void          StateVector_Deallocate        (CStateVector *self)
  cdef void          StateVector_Reset             (CStateVector *self)
  cdef void          StateVector_ResetToMaximum    (CStateVector *self)
  cdef Boolean       StateVector_GetItem           (CStateVector *self, Integer index)
  cdef Boolean       StateVector_SetItem           (CStateVector *self, Integer index, Integer value)
  cdef Boolean       StateVector_Increment         (CStateVector *self)

  # Substate-related functions
  cdef void         StateVector_ResetSubstate     (CStateVector *self)
  cdef Boolean      StateVector_AllocateSubstate  (CStateVector *self, Integer nsites)
  cdef Boolean      StateVector_IncrementSubstate (CStateVector *self)
  cdef Boolean      StateVector_SetSubstateItem   (CStateVector *self, Integer selectedSiteIndex, Integer index)
  cdef Integer      StateVector_GetSubstateItem   (CStateVector *self, Integer index)

#-------------------------------------------------------------------------------
cdef class StateVector:
  cdef CStateVector     *cObject

#  cdef public object length
#  cdef public object vector
#  cdef public object maxvector
