#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore    import logFile, LogFileActive, CLibraryError


cdef class StateVector:
  """A class defining the state vector."""

  def __len__ (self):
    """Return the size of the vector."""
    return self.cObject.length


  def __getitem__ (self, index):
    """Get an item."""
    return StateVector_GetItem (self.cObject, index)


  def __setitem__ (self, index, value):
    """Set an item."""
    StateVector_SetItem (self.cObject, index, value)


  def __dealloc__ (self):
    """Deallocate."""
    StateVector_Deallocate (self.cObject)


  def Increment (self):
    """Generate a new state incrementing the state vector.

    Return False after the incrementation has finished."""
    return StateVector_Increment (self.cObject)


  def Reset (self):
    """Set all components of the vector to zero."""
    StateVector_Reset (self.cObject)


  def ResetToMaximum (self):
    """Set all components of the vector to their maximum values."""
    StateVector_ResetToMaximum (self.cObject)


  def __init__ (self, meadModel):
    """Constructor."""
    nsites       = len (meadModel.meadSites)
    self.cObject = StateVector_Allocate (nsites)
    if (self.cObject == NULL): 
      raise CLibraryError ("Memory allocation failure.")

    for isite, site in enumerate (meadModel.meadSites):
      ninstances = len (site.instances)
      self.cObject.maxvector[isite] = ninstances - 1


  def Print (self, meadModel = None, title = None, log = logFile):
    """Print the state vector."""
    if LogFileActive (log):
      if meadModel:
        table = log.GetTable (columns = [7, 7, 7, 8, 8])
        table.Start ()
        if title is None:
          table.Title ("State vector")
        else:
          table.Title (title)

        table.Heading ("Site", columnSpan = 3)
        table.Heading ("Instance", columnSpan = 2)

        # Faster than: for i in range (0, self.cObject.length):
        for 0 <= i < self.cObject.length:
          site        = meadModel.meadSites[i]
          protonation = self.cObject.vector[i]
  
          instance = site.instances[protonation]
          table.Entry (site.segName)
          table.Entry (site.resName)
          table.Entry (site.resSerial)
          table.Entry ("%s" % instance.label)
          table.Entry ("%d" % protonation)
        table.Stop ()


  def DefineSubstate (self, meadModel, selectedSites):
    """Define a substate.

    |selectedSites| is a sequence of two-element sequences (segmentName, residueSerial)"""
    nsites = len (selectedSites)
    if not (StateVector_AllocateSubstate (self.cObject, nsites)):
      raise CLibraryError ("Memory allocation failure.")

    for selectedSiteIndex, selectedSite in enumerate (selectedSites):
      selectedSegment, selectedSerial = selectedSite

      # FIXME
      selectedSerial = str (selectedSerial)

      foundSite = False
      for siteIndex, site in enumerate (meadModel.meadSites):
        if site.segName == selectedSegment and site.resSerial == selectedSerial:
          foundSite = True
          break

      if not foundSite:
        raise CLibraryError ("Site %s %d not found." % (selectedSegment, selectedSerial))

      StateVector_SetSubstateItem (self.cObject, selectedSiteIndex, siteIndex)


  def IncrementSubstate (self):
    """Generate a new substate. 

    Return False after the incrementation has finished."""
    return StateVector_IncrementSubstate (self.cObject)


  def ResetSubstate (self):
    """Set all components of the substate to zero."""
    StateVector_ResetSubstate (self.cObject)
