#-------------------------------------------------------------------------------
# . File      : Substate.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Helper functions for handling substates."""

__lastchanged__ = "$Id$"

from pCore          import logFile, LogFileActive
from Error          import ContinuumElectrostaticsError
from StateVector    import StateVector


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


def SubstateEnergies (meadModel, selectedSites, pH = 7.0, relativeEnergy = True, log = logFile):
  """Calculate and print out microstate energies of substates defined by |selectedSites|.

  |selectedSites| is a sequence of three-element sequences (segmentName, residueName, residueSerial)"""
  pairs = []
  indicesOfSites = []

  for selectedSegment, selectedResidueName, selectedResidueSerial in selectedSites:
    for siteIndex, site in enumerate (meadModel.meadSites):
      if site.segName == selectedSegment and site.resSerial == selectedResidueSerial:
        indicesOfSites.append (siteIndex)
        break
    pairs.append ([selectedSegment, selectedResidueSerial])

  vector = StateVector_FromProbabilities (meadModel)
  vector.DefineSubstate (meadModel, pairs)
  vector.ResetSubstate ()

  increment = True
  substates = []

  while increment:
    Gmicro = meadModel.CalculateMicrostateEnergy (vector, pH = pH)
    indicesOfInstances = []
    for siteIndex in indicesOfSites:
      indicesOfInstances.append (vector[siteIndex])
    substates.append ([Gmicro, indicesOfInstances])

    increment = vector.IncrementSubstate ()

  substates.sort ()
  lowestSubstate = substates[0]
  lowestEnergy   = lowestSubstate[0]


  if LogFileActive (log):
    nsites = len (indicesOfSites)
    tab = log.GetTable (columns = [6, 9] + [14] * nsites)
    tab.Start ()
    tab.Heading ("State")
    tab.Heading ("Gmicro")

    for siteIndex in indicesOfSites:
      site = meadModel.meadSites[siteIndex]
      tab.Heading ("%s %s %d" % (site.segName, site.resName, site.resSerial))

    for stateIndex, (energy, indicesOfInstances) in enumerate (substates):
      tab.Entry ("%6d"   % (stateIndex + 1))
      if relativeEnergy:
        tab.Entry ("%9.2f" % (energy - lowestEnergy))
      else:
        tab.Entry ("%9.2f" % energy)

      for siteIndex, instanceIndex in zip (indicesOfSites, indicesOfInstances):
        site     = meadModel.meadSites     [siteIndex]
        instance = site.instances      [instanceIndex]
        tab.Entry (instance.label.center (14))
    tab.Stop ()


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
