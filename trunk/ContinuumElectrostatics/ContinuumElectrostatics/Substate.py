#-------------------------------------------------------------------------------
# . File      : Substate.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Substate is a class for calculating energies of substates of protonation states."""

__lastchanged__ = "$Id$"

from pCore          import logFile, LogFileActive
from Error          import ContinuumElectrostaticsError
from StateVector    import StateVector


class MEADSubstate (object):
  """Substate of a protonation state."""

  def __init__ (self, meadModel, selectedSites, pH = 7.0, log = logFile):
    """Construct a substate defined by |selectedSites|.

    |selectedSites| is a sequence of three-element sequences (segmentName, residueName, residueSerial)"""
    pairs = []
    indicesOfSites = []
  
    for selectedSegment, selectedResidueName, selectedResidueSerial in selectedSites:
      foundSite = False
      for siteIndex, site in enumerate (meadModel.meadSites):
        if site.segName == selectedSegment and site.resSerial == selectedResidueSerial:
          indicesOfSites.append (siteIndex)
          foundSite = True
          break

      if not foundSite:
        raise ContinuumElectrostaticsError ("Site %s %s %d not found." % (selectedSegment, selectedResidueName, selectedResidueSerial))
      pairs.append ([selectedSegment, selectedResidueSerial])

    vector = StateVector_FromProbabilities (meadModel)
    vector.DefineSubstate (meadModel, pairs)

    self.indicesOfSites = indicesOfSites
    self.isCalculated   = False
    self.substates      = None
    self.vector         = vector
    self.model          = meadModel
    self.pH             = pH

    if LogFileActive (log):
      log.Text ("\nSubstate is initialized with %d sites.\n" % len (indicesOfSites))
 

  def CalculateSubstateEnergies (self, log = logFile): 
    """Calculate microstate energies for a substate."""
    if not self.isCalculated:
      indicesOfSites = self.indicesOfSites
      increment      = True
      substates      = []
      vector         = self.vector
      vector.ResetSubstate ()
    
      while increment:
        Gmicro = self.model.CalculateMicrostateEnergy (vector, pH = self.pH)
        indicesOfInstances = []
        for siteIndex in indicesOfSites:
          indicesOfInstances.append (vector[siteIndex])
        substates.append ([Gmicro, indicesOfInstances])
    
        increment = vector.IncrementSubstate ()
    
      substates.sort ()
      lowestSubstate = substates[0]
      lowestEnergy   = lowestSubstate[0]

      self.substates    = substates
      self.zeroEnergy   = lowestEnergy
      self.isCalculated = True

      if LogFileActive (log):
        log.Text ("\nCalculating substate energies at pH=%.1f complete.\n" % self.pH)


  def Summary (self, relativeEnergy = True, roundCharge = True, log = logFile):
    """Summarize calculated substate energies in a table."""
    if self.isCalculated:
      indicesOfSites = self.indicesOfSites
      zeroEnergy     = self.zeroEnergy
      substates      = self.substates
      model          = self.model
      nsites         = len (indicesOfSites)

      if LogFileActive (log):
        tab = log.GetTable (columns = [6, 9, 8, 8] + [14] * nsites)
        tab.Start ()
        tab.Heading ("State")
        tab.Heading ("Gmicro")
        tab.Heading ("Charge")
        tab.Heading ("Protons")
    
        for siteIndex in indicesOfSites:
          site = model.meadSites[siteIndex]
          tab.Heading ("%s %s %d" % (site.segName, site.resName, site.resSerial))
    
        for stateIndex, (energy, indicesOfInstances) in enumerate (substates):
          tab.Entry ("%6d"   % (stateIndex + 1))
          if relativeEnergy:
            energy = energy - zeroEnergy
          tab.Entry ("%9.2f" % energy)

          nprotons = 0
          charge   = 0.
          labels   = []
          for siteIndex, instanceIndex in zip (indicesOfSites, indicesOfInstances):
            site     = model.meadSites     [siteIndex]
            instance = site.instances  [instanceIndex]
            nprotons = nprotons + instance.protons
            charge   = charge + sum (instance.charges) 
            labels.append (instance.label)

          # Charges should ALWAYS sum up to integer values
          if roundCharge:
            tab.Entry ("%d" % round (charge))
          else:
            tab.Entry ("%.1f" % charge)
          tab.Entry ("%d"   % nprotons)

          for label in labels:
            tab.Entry (label.center (14))
        tab.Stop ()


#===============================================================================
# Helper functions
#===============================================================================
def StateVector_FromProbabilities (meadModel):
  """Create a state vector from the previously calculated probabilities."""
  if meadModel.isProbability:
    vector = StateVector (meadModel)

    for siteIndex, site in enumerate (meadModel.meadSites):
      pairs = []
      for instanceIndex, instance in enumerate (site.instances):
        pair = (instance.probability, instanceIndex)
        pairs.append (pair)
      maxProbPair = max (pairs)
      probability, instanceIndex = maxProbPair
      vector[siteIndex] = instanceIndex
    return vector
  else:
    raise ContinuumElectrostaticsError ("First calculate probabilities.")


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
