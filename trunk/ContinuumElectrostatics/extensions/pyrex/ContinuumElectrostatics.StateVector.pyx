#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------


from pCore         import logFile, LogFileActive, CLibraryError


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


  def __init__ (self, continuumElectrostaticModel):
    """Constructor."""
    #cdef Integer nsites, isite
    nsites       = len (continuumElectrostaticModel.meadSites)
    self.cObject = StateVector_Allocate (nsites)
    if (self.cObject == NULL): 
      raise CLibraryError ( "Memory allocation failure." )

    for isite, site in enumerate (continuumElectrostaticModel.meadSites):
      ninstances = len (site.instances)
      self.cObject.maxvector[isite] = ninstances - 1


  def Print (self, continuumElectrostaticModel = None, title = None, log = logFile):
    """Print the state vector."""
    if LogFileActive (log):
      if continuumElectrostaticModel:
        table = log.GetTable (columns = [7, 7, 7, 8, 8])
        table.Start ()
        if title is None:
          table.Title ("State vector")
        else:
          table.Title (title)

        table.Heading ("Site", columnSpan = 3)
        table.Heading ("Instance", columnSpan = 2)

        for i in range (0, self.cObject.length):
          site        = continuumElectrostaticModel.meadSites[i]
          protonation = self.cObject.vector[i]
  
          instance = site.instances[protonation]
          table.Entry (site.segName)
          table.Entry (site.resName)
          table.Entry (site.resNum)
          table.Entry ("%s" % instance.label)
          table.Entry ("%d" % protonation)
        table.Stop ()
