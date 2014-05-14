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
    Integer  length


  cdef CStateVector *StateVector_Allocate       (Integer size)
  cdef void          StateVector_Deallocate     (CStateVector *self)
  cdef void          StateVector_Reset          (CStateVector *self)
  cdef void          StateVector_ResetToMaximum (CStateVector *self)
  cdef Integer       StateVector_GetItem        (CStateVector *self, Integer index)
  cdef Integer       StateVector_SetItem        (CStateVector *self, Integer index, Integer value)
  cdef Integer       StateVector_Increment      (CStateVector *self)

#-------------------------------------------------------------------------------
cdef class StateVector:
  cdef CStateVector     *cObject

#  cdef public object length
#  cdef public object vector
#  cdef public object maxvector
