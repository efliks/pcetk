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
  ctypedef struct CTitrSite "TitrSite":
    Integer  indexActive
    Integer  indexFirst
    Integer  indexLast
    Integer  indexSite

  ctypedef struct CStateVector "StateVector":
    CTitrSite  *sites
    CTitrSite **substateSites
    Integer     nsites
    Integer     nssites


  # Allocation and deallocation
  cdef CStateVector *StateVector_Allocate          (Integer nsites, Status *status)
  cdef void          StateVector_AllocateSubstate  (CStateVector *self, Integer nSubstateSites, Status *status)
  cdef void          StateVector_Deallocate        (CStateVector *self)

  # Copying and cloning
  cdef CStateVector *StateVector_Clone             (CStateVector *self, Status *status)
  cdef void          StateVector_CopyTo            (CStateVector *self, CStateVector *other, Status *status)

  # Setting all items at once
  cdef void          StateVector_Reset             (CStateVector *self)
  cdef void          StateVector_ResetSubstate     (CStateVector *self)
  cdef void          StateVector_ResetToMaximum    (CStateVector *self)
  cdef void          StateVector_Randomize         (CStateVector *self)

  # Accessing items
  cdef void          StateVector_SetSite           (CStateVector *self, Integer indexSite, Integer indexFirst, Integer indexLast, Status *status)
  cdef Integer       StateVector_GetItem           (CStateVector *self, Integer siteIndex, Status *status)
  cdef void          StateVector_SetItem           (CStateVector *self, Integer siteIndex, Integer instanceLocalIndex, Status *status)
  cdef Integer       StateVector_GetActualItem     (CStateVector *self, Integer siteIndex, Status *status)
  cdef void          StateVector_SetActualItem     (CStateVector *self, Integer siteIndex, Integer instanceGlobalIndex, Status *status)
  cdef Integer       StateVector_GetSubstateItem   (CStateVector *self, Integer index, Status *status)
  cdef void          StateVector_SetSubstateItem   (CStateVector *self, Integer selectedSiteIndex, Integer index, Status *status)

  # Incrementation
  cdef Boolean       StateVector_Increment         (CStateVector *self)
  cdef Boolean       StateVector_IncrementSubstate (CStateVector *self)

  # Monte Carlo-related functions
  cdef void          StateVector_Move              (CStateVector *self, Integer *siteIndex, Integer *oldIndexActive)

#-------------------------------------------------------------------------------
cdef class StateVector:
  cdef CStateVector     *cObject
  cdef public object  isOwner
  cdef public object  owner
