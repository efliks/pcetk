#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore      import logFile, LogFileActive, CLibraryError

import random

DEF ANALYTIC_STATES = 67108864
__lastchanged__ = "$Id$"


cdef class StateVector:
  """A class defining the state vector."""

  def Randomize (self):
    """Generate a random state."""
    cdef Integer index, rand
    for index from 0 <= index < self.cObject.length:
      rand = random.randint (self.cObject.minvector[index], self.cObject.maxvector[index])
      self.cObject.vector[index] = rand
 
 
  def RandomChange (self):
    """Choose a random site and set it to a random instance."""
    cdef Integer index, rand
    index = random.randint (0, self.cObject.length - 1)
    rand  = random.randint (self.cObject.minvector[index], self.cObject.maxvector[index])

    if rand == self.cObject.vector[index]:
      rand += 1
      if rand > self.cObject.maxvector[index]:
        rand = self.cObject.minvector[index]
    self.cObject.vector[index] = rand


  def __len__ (self):
    """Return the size of the vector."""
    return self.cObject.length


  def __getmodule__ (self):
    """Return the module name."""
    return "ContinuumElectrostatics.StateVector"


  def __dealloc__ (self):
    """Deallocate."""
    StateVector_Deallocate (self.cObject)


#  def __copy__ (self):
#    """Copying."""
#    return self.__deepcopy__ ()
#
#
#  def __deepcopy__ (self):
#    """Copying."""
#    cdef StateVector clone
#    clone.cObject = StateVector_Clone (self.cObject)
#    if clone.cObject == NULL:
#      raise CLibraryError ("Cannot copy vector.")
#    return clone


  def __getitem__ (self, Integer index):
    """Get an item."""
    cdef Integer item
    item = StateVector_GetItem (self.cObject, index)
    if item < 0:
      raise CLibraryError ("Index out of range.")
    return item


  def __setitem__ (self, Integer index, Integer value):
    """Set an item."""
    cdef Boolean status
    status = StateVector_SetItem (self.cObject, index, value)
    if status == CFalse:
      raise CLibraryError ("Index out of range or wrong value.")


  def __init__ (self, meadModel):
    """Constructor."""
    cdef Integer numberOfSites
    cdef Integer indexSite, indexDown, indexUp, index

    numberOfSites = len (meadModel.meadSites)
    self.cObject = StateVector_Allocate (numberOfSites)
    if (self.cObject == NULL):
      raise CLibraryError ("Cannot allocate state vector.")

    for indexSite from 0 <= indexSite < numberOfSites:
      indexDown = 9999
      indexUp   = 0
      site      = meadModel.meadSites[indexSite]

      for instance in site.instances:
        index = instance.instIndexGlobal
        if index < indexDown:
          indexDown = index
        if index > indexUp:
          indexUp = index
      self.cObject.minvector[indexSite] = indexDown
      self.cObject.maxvector[indexSite] = indexUp

    # Always reset the state vector after initialization
    StateVector_Reset (self.cObject)


  def CopyTo (StateVector self, StateVector other):
    """In-place copy of one state vector to another."""
    StateVector_CopyTo (self.cObject, other.cObject)


  def GetActualItem (self, Integer index):
    """Get an actual item."""
    cdef Integer item
    item = StateVector_GetActualItem (self.cObject, index)
    if item < 0:
      raise CLibraryError ("Index out of range.")
    return item


  def SetActualItem (self, Integer index, Integer value):
    """Set an actual item."""
    cdef Boolean status
    status = StateVector_SetActualItem (self.cObject, index, value)
    if status == CFalse:
      raise CLibraryError ("Index out of range or wrong value.")


  def Increment (self):
    """Generate a new state incrementing the state vector.

    Return False after the incrementation has finished."""
    cdef Boolean status
    status = StateVector_Increment (self.cObject)
    if status == CFalse:
      return False
    else:
      return True


  def Reset (self):
    """Set all components of the vector to their minimum values (formerly zeros)."""
    StateVector_Reset (self.cObject)


  def ResetToMaximum (self):
    """Set all components of the vector to their maximum values."""
    StateVector_ResetToMaximum (self.cObject)


# *** Try to further improve this method with data types from C ***
  def Print (self, meadModel=None, title=None, log=logFile):
    """Print the state vector."""
    cdef Integer siteIndex, selectedSiteIndex, instanceIndex
    cdef Integer ninstances, j

    if LogFileActive (log):
      if meadModel:
        tab = log.GetTable (columns = [7, 7, 7, 8, 8, 2])
        tab.Start ()
        if title is None:
          tab.Title ("State vector")
        else:
          tab.Title (title)

        tab.Heading ("Site", columnSpan = 3)
        tab.Heading ("Instance", columnSpan = 3)

        for siteIndex from 0 <= siteIndex < self.cObject.length:
          substate = "  "

          if self.cObject.substate != NULL:
            for j from 0 <= j < self.cObject.slength:
              selectedSiteIndex = StateVector_GetSubstateItem (self.cObject, j)
              if selectedSiteIndex == siteIndex:
                substate = " @"

          site = meadModel.meadSites[siteIndex]
          ninstances = len (site.instances)

          for instanceIndex from 0 <= instanceIndex < ninstances:
            instance = site.instances[instanceIndex]
            if instance.instIndexGlobal == self.cObject.vector[siteIndex]:
              break
          instance = site.instances[instanceIndex]

          tab.Entry (site.segName)
          tab.Entry (site.resName)
          tab.Entry ("%d" % site.resSerial)
          tab.Entry (instance.label)
          tab.Entry ("%d" % instanceIndex)
          tab.Entry (substate)
        tab.Stop ()


# *** Try to further improve this method with data types from C ***
  def DefineSubstate (self, meadModel, selectedSites):
    """Define a substate.

    |selectedSites| is a sequence of two-element sequences (segmentName, residueSerial)"""
    cdef Boolean status, foundSite
    cdef Integer siteIndex, substateSiteIndex
    cdef Integer nselected  = len (selectedSites)
    cdef Integer nsites     = len (meadModel.meadSites)

    status = StateVector_AllocateSubstate (self.cObject, nselected)
    if status == CFalse:
      raise CLibraryError ("Memory allocation failure.")

    for substateSiteIndex from 0 <= substateSiteIndex < nselected:
      selectedSegment, selectedSerial = selectedSites[substateSiteIndex]
      foundSite = CFalse

      for siteIndex from 0 <= siteIndex < nsites:
        site = meadModel.meadSites[siteIndex]

        if site.segName == selectedSegment and site.resSerial == selectedSerial:
          StateVector_SetSubstateItem (self.cObject, siteIndex, substateSiteIndex)
          foundSite = CTrue
          break

      if foundSite == CFalse:
        raise CLibraryError ("Site %s %d not found." % (selectedSegment, selectedSerial))


  def IncrementSubstate (self):
    """Generate a new substate.

    Return False after the incrementation has finished."""
    cdef Boolean status
    status = StateVector_IncrementSubstate (self.cObject)
    if status == CFalse:
      return False
    else:
      return True


  def ResetSubstate (self):
    """Set all components of the substate to their minimum values (formerly zeros)."""
    StateVector_ResetSubstate (self.cObject)


  def CalculateMicrostateEnergy (self, meadModel, Real pH=7.0):
    """Calculate energy of a protonation state (=microstate)."""
    cdef Real             Gmicro           =  0.
    cdef Real             temperature      =  meadModel.temperature
    cdef Integer1DArray   protons          =  meadModel._protons
    cdef Real1DArray      intrinsic        =  meadModel._intrinsic
    cdef SymmetricMatrix  symmetricmatrix  =  meadModel._symmetricmatrix

    if meadModel.isCalculated:
      Gmicro = StateVector_CalculateMicrostateEnergy (self.cObject, protons.cObject, intrinsic.cObject, symmetricmatrix.cObject, pH, temperature)
    return Gmicro


  def GetNumberOfStates (self, meadModel):
    """Calculate total number of possible protonation states."""
    cdef Integer nstates = 1, ninstances

    for meadSite in meadModel.meadSites:
      ninstances = len (meadSite.instances)
      nstates = nstates * ninstances
      if nstates > ANALYTIC_STATES:
        raise CLibraryError ("Maximum number of states (%d) exceeded." % ANALYTIC_STATES)
    return nstates


  def CalculateProbabilitiesAnalytically (self, meadModel, Real pH=7.0):
    """Calculate probabilities of protonation states analytically."""
    cdef Boolean          status
    cdef Integer          nstates          =  self.GetNumberOfStates (meadModel)
    cdef Real             temperature      =  meadModel.temperature
    cdef Integer1DArray   protons          =  meadModel._protons
    cdef Real1DArray      intrinsic        =  meadModel._intrinsic
    cdef Real1DArray      probabilities    =  meadModel._probabilities
    cdef SymmetricMatrix  symmetricmatrix  =  meadModel._symmetricmatrix

    status = StateVector_CalculateProbabilitiesAnalytically (self.cObject, protons.cObject, intrinsic.cObject, symmetricmatrix.cObject, pH, temperature, nstates, probabilities.cObject)
    if status == CFalse:
      raise CLibraryError ("Cannot allocate Boltzmann factors.")
    return nstates
