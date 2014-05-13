#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.StateVector.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------

__revision__ = "$Revision: 6 $"

#from CoreObjects   import CLibraryError
from pCore         import logFile, LogFileActive


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
    nsites       = len (continuumElectrostaticModel.meadSites)
    self.cObject = StateVector_Allocate (nsites)
#    if (self.cObject == NULL): 
#      raise CLibraryError ( "Memory allocation failure." )

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


#  #-------------------------------------------------------------------------------
#  # . File      : StateVector.py
#  # . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
#  # . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#  #                          Mikolaj J. Feliks (2014)
#  # . License   : CeCILL French Free Software License     (http://www.cecill.info)
#  #-------------------------------------------------------------------------------
#  
#  __revision__ = "$Revision: 6 $"
#  
#  import exceptions
#  
#  from pCore           import logFile, LogFileActive, Integer1DArray
#  
#  
#  class StateVectorError (exceptions.StandardError):
#    """A class for handling errors in StateVector."""
#    pass
#  
#  
#  #-------------------------------------------------------------------------------
#  class StateVector (object):
#    """A class defining the state vector."""
#  
#    def __len__ (self):
#      """Return the size of the vector."""
#      return len (self.vector)
#  
#  
#    def __getitem__ (self, i):
#      """Get an item."""
#      return self.vector[i]
#  
#  
#    def __setitem__ (self, i, value):
#      """Set an item."""
#      if not isinstance (i, int):
#        raise TypeError ( "Expecting integer not %s." % type (i))
#  
#      if (i < 0) or (i > len (self.vector) - 1):
#        raise IndexError ("Index %d out of range." % i)
#  
#      if not isinstance (value, int):
#        raise TypeError ( "Expecting integer not %s." % type (value))
#  
#      if value < 0 or value > self.maxvector[i]:
#        raise TypeError ( "Expecting a value in the range of 0..%d" % self.maxvector[i])
#  
#      self.vector[i] = value
#  
#  
#    def __init__ (self, continuumElectrostaticModel):
#      """Constructor."""
#  #    if not isinstance (continuumElectrostaticModel, CEModelMEAD):
#  #      raise StateVectorError ("An instance of CEModelMEAD is required.")
#  #
#  #    if not continuumElectrostaticModel.isInitialized:
#  #      raise StateVectorError ("First initialize the CEModelMEAD model.")
#  
#      nsites         = len (continuumElectrostaticModel.meadSites)
#      self.vector    = Integer1DArray.WithExtent (nsites)
#      self.maxvector = Integer1DArray.WithExtent (nsites)
#      self.length    = nsites
#      self.vector.Set (0)
#  
#      for isite, site in enumerate (continuumElectrostaticModel.meadSites):
#        ninstances = len (site.instances)
#        self.maxvector[isite] = ninstances - 1
#  
#  
#    def Print (self, continuumElectrostaticModel = None, itemsPerRow = 25, itemWidth = 4, title = None, log = logFile):
#      """Print the state vector."""
#      if LogFileActive (log):
#        if not continuumElectrostaticModel is None:
#          table = log.GetTable (columns = [7, 7, 7, 8, 8])
#          table.Start ()
#          if title is None:
#            table.Title ("State vector")
#          else:
#            table.Title (title)
#  
#          table.Heading ("Site", columnSpan = 3)
#          table.Heading ("Instance", columnSpan = 2)
#    
#          for site, protonation in zip (continuumElectrostaticModel.meadSites, self.vector):
#            instance = site.instances[protonation]
#            table.Entry (site.segName)
#            table.Entry (site.resName)
#            table.Entry (site.resNum)
#            table.Entry ("%s" % instance.label)
#            table.Entry ("%d" % protonation)
#          table.Stop ()
#        else:
#          self.vector.Print (itemFormat = ("%%%dd" % itemWidth), itemsPerRow = itemsPerRow, itemWidth = itemWidth)
#  
#  
#    def Reset (self):
#      """Set all components of the vector to zero."""
#      self.vector.Set (0)
#  
#  
#    def Increment (self):
#      """Generate a new state incrementing the state vector.
#  
#      Return False after the incrementation has finished."""
#  # This prevents zeroing the vector after the last incrementation
#  #    increment = False
#  #    for value, maxvalue in zip (self.vector, self.maxvector):
#  #      if value < maxvalue:
#  #        increment = True
#  #        break
#  #    if not increment:
#  #      return False
#      for index, (value, maxvalue) in enumerate (zip (self.vector, self.maxvector)):
#        if value < maxvalue:
#          value = value + 1
#          self.vector[index] = value
#          return True
#        else:
#          value = 0
#          self.vector[index] = value
#      return False
#  
#  
#  #  def ResetToRandom (self):
#  #    """Generate a random state vector."""
#  #    pass
#  #
#  #  def Randomize (self):
#  #    """Set a random component in the vector to a random value."""
#  #    pass
#  
#  
#  #===============================================================================
#  # Testing
#  #===============================================================================
#  if __name__ == "__main__": pass
