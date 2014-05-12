#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pxd
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------

__revision__ = "$Revision: 6 $"


cdef extern from "StateVector.h":

  ctypedef struct CStateVector "StateVector":
    pass


  cdef CStateVector     *StateVector_Allocate (const Integer size)


cdef class StateVector:
  cdef CStateVector     *cObject
