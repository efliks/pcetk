#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore    import logFile, LogFileActive, CLibraryError

__lastchanged__ = "$Id$"


cdef class StateVector:
  """A class defining the state vector."""

  def __len__ (self):
    """Return the size of the vector."""
    return self.cObject.length


  def __getitem__ (self, index):
    """Get an item."""
    item = StateVector_GetItem (self.cObject, index)
    if item < 0:
      raise CLibraryError ("Vector index out of range.")
    return item


  def __setitem__ (self, index, value):
    """Set an item."""
    ok = StateVector_SetItem (self.cObject, index, value)
    if not ok:
      raise CLibraryError ("Vector index out of range or wrong value.")


  def __dealloc__ (self):
    """Deallocate."""
    StateVector_Deallocate (self.cObject)


  def Increment (self):
    """Generate a new state incrementing the state vector.

    Return False after the incrementation has finished."""
    return StateVector_Increment (self.cObject)


  def Reset (self):
    """Set all components of the vector to their minimum values (formerly zeros)."""
    StateVector_Reset (self.cObject)


  def ResetToMaximum (self):
    """Set all components of the vector to their maximum values."""
    StateVector_ResetToMaximum (self.cObject)


  def __init__ (self, meadModel):
    """Constructor."""
    nsites = len (meadModel.meadSites)
    self.cObject = StateVector_Allocate (nsites)

    if (self.cObject == NULL): 
      raise CLibraryError ("Memory allocation failure.")

    for siteIndex, site in enumerate (meadModel.meadSites):
      indices = []
      for instance in site.instances:
        indices.append (instance.instIndexGlobal)
      self.cObject.minvector[siteIndex] = min (indices)
      self.cObject.maxvector[siteIndex] = max (indices)


  def Print (self, meadModel = None, title = None, log = logFile):
    """Print the state vector."""
    if LogFileActive (log):
      if meadModel:
        table = log.GetTable (columns = [7, 7, 7, 8, 8, 2])
        table.Start ()
        if title is None:
          table.Title ("State vector")
        else:
          table.Title (title)

        table.Heading ("Site", columnSpan = 3)
        table.Heading ("Instance", columnSpan = 3)

        for 0 <= siteIndex < self.cObject.length:
          substate = "  "

          if self.cObject.substate != NULL:
            for 0 <= j < self.cObject.slength:
              selectedSiteIndex = StateVector_GetSubstateItem (self.cObject, j)
              if selectedSiteIndex == siteIndex:
                substate = " @"

          site = meadModel.meadSites [siteIndex]
          for instanceIndex, instance in enumerate (site.instances):
            if instance.instIndexGlobal == self.cObject.vector[siteIndex]:
              break
          instance = site.instances[instanceIndex]

          table.Entry (site.segName)
          table.Entry (site.resName)
          table.Entry ("%d" % site.resSerial)
          table.Entry (instance.label)
          table.Entry ("%d" % instanceIndex)
          table.Entry (substate)
        table.Stop ()


  def DefineSubstate (self, meadModel, selectedSites):
    """Define a substate.

    |selectedSites| is a sequence of two-element sequences (segmentName, residueSerial)"""
    ok = StateVector_AllocateSubstate (self.cObject, len (selectedSites))
    if not ok:
      raise CLibraryError ("Memory allocation failure.")

    for substateSiteIndex, (selectedSegment, selectedSerial) in enumerate (selectedSites):
      foundSite = False

      for siteIndex, site in enumerate (meadModel.meadSites):
        if site.segName == selectedSegment and site.resSerial == selectedSerial:
          StateVector_SetSubstateItem (self.cObject, siteIndex, substateSiteIndex)
          foundSite = True
          break

      if not foundSite:
        raise CLibraryError ("Site %s %d not found." % (selectedSegment, selectedSerial))


  def IncrementSubstate (self):
    """Generate a new substate. 

    Return False after the incrementation has finished."""
    return StateVector_IncrementSubstate (self.cObject)


  def ResetSubstate (self):
    """Set all components of the substate to their minimum values (formerly zeros)."""
    StateVector_ResetSubstate (self.cObject)


  def CalculateMicrostateEnergy (self, meadModel, pH = 7.0):
    """Calculate energy of a protonation state (=microstate)."""
    cdef Integer1DArray  arrayProtons
    cdef Real1DArray     arrayIntrinsic
    cdef Real2DArray     arrayInteractions
    cdef Real            Gmicro

    if meadModel.isCalculated:
      arrayProtons      = meadModel.arrayProtons
      arrayIntrinsic    = meadModel.arrayIntrinsic
      arrayInteractions = meadModel.arrayInteractions
  
      Gmicro = StateVector_CalculateMicrostateEnergy (self.cObject, arrayProtons.cObject, arrayIntrinsic.cObject, arrayInteractions.cObject, pH, meadModel.temperature)
    return Gmicro
