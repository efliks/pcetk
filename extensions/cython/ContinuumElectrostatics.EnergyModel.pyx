#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.EnergyModel.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore        import logFile, LogFileActive, CLibraryError
from StateVector  import StateVector

DEF ANALYTIC_STATES = 67108864
__lastchanged__ = "$Id: $"


cdef class EnergyModel:
  """A class defining the energy model."""

  def __getmodule__ (self):
    """Return the module name."""
    return "ContinuumElectrostatics.EnergyModel"

  def __dealloc__ (self):
    """Deallocate."""
    EnergyModel_Deallocate (self.cObject)

  # Setters
  def SetGintr (self, Integer instIndexGlobal, Real value):
    EnergyModel_SetGintr (self.cObject, instIndexGlobal, value)

  def SetProtons (self, Integer instIndexGlobal, Integer value):
    EnergyModel_SetProtons (self.cObject, instIndexGlobal, value)

  def SetProbability (self, Integer instIndexGlobal, Real value):
    EnergyModel_SetProbability (self.cObject, instIndexGlobal, value)

  def SetInteraction (self, Integer instIndexGlobalA, Integer instIndexGlobalB, Real value):
    EnergyModel_SetInteraction (self.cObject, instIndexGlobalA, instIndexGlobalB, value)

  # Getters
  def GetGintr (self, Integer instIndexGlobal):
    cdef Real Gintr = EnergyModel_GetGintr (self.cObject, instIndexGlobal)
    return Gintr

  def GetProtons (self, Integer instIndexGlobal):
    cdef Integer protons = EnergyModel_GetProtons (self.cObject, instIndexGlobal)
    return protons

  def GetProbability (self, Integer instIndexGlobal):
    cdef Real prob = EnergyModel_GetProbability (self.cObject, instIndexGlobal)
    return prob

  def GetInteraction (self, Integer instIndexGlobalA, Integer instIndexGlobalB):
    cdef Real interac = EnergyModel_GetInteraction (self.cObject, instIndexGlobalA, instIndexGlobalB)
    return interac

  def GetInteractionSymmetric (self, Integer instIndexGlobalA, Integer instIndexGlobalB):
    cdef Real interac = EnergyModel_GetInterSymmetric (self.cObject, instIndexGlobalA, instIndexGlobalB)
    return interac

  def GetDeviation (self, Integer instIndexGlobalA, Integer instIndexGlobalB):
    cdef Real deviate = EnergyModel_GetDeviation (self.cObject, instIndexGlobalA, instIndexGlobalB)
    return deviate


  def __init__ (self, meadModel):
    """Constructor."""
    cdef Integer nstates, ninstances, totalSites, totalInstances
    cdef Integer indexSite, indexDown, indexUp, index
    cdef Real    temperature
    cdef Status  status

    status         = Status_Continue
    totalSites     = meadModel.nsites
    totalInstances = meadModel.ninstances
    temperature    = meadModel.temperature

    self.cObject   = EnergyModel_Allocate (totalSites, totalInstances, &status)
    if status != Status_Continue:
      raise CLibraryError ("Cannot allocate energy model.")

    nstates = 1
    for indexSite from 0 <= indexSite < totalSites:
      ninstances =  0
      indexUp    =  0
      indexDown  =  9999
      site       =  meadModel.meadSites[indexSite]
      for instance in site.instances:
        ninstances += 1
        index = instance._instIndexGlobal
        if index < indexDown : indexDown = index
        if index > indexUp   : indexUp   = index
      StateVector_SetSite (self.cObject.vector, indexSite, indexDown, indexUp, &status)

      if nstates <= ANALYTIC_STATES:
        nstates = nstates * ninstances

    self.cObject.nstates     = nstates
    self.cObject.ninstances  = totalInstances
    self.cObject.temperature = temperature
    self.isOwner = False
    self.owner   = meadModel


  def CheckIfSymmetric (self, Real tolerance=0.05):
    """Check if the matrix of interactions is symmetric within the given threshold."""
    cdef Real    maxDeviation
    cdef Boolean isSymmetric
    isSymmetric = EnergyModel_CheckInteractionsSymmetric (self.cObject, tolerance, &maxDeviation)
    return (isSymmetric, maxDeviation)


  def SymmetrizeInteractions (self, log=logFile):
    """Symmetrize the matrix of interactions."""
    cdef Status status = Status_Continue
    EnergyModel_SymmetrizeInteractions (self.cObject, &status)
    if LogFileActive (log):
      log.Text ("\nSymmetrizing interactions complete.\n")


  def CalculateMicrostateEnergy (self, StateVector vector, Real pH=7.0):
    """Calculate energy of a protonation state (=microstate)."""
    cdef Real Gmicro
    meadModel = self.owner

    if not meadModel.isCalculated:
      raise CLibraryError ("First calculate electrostatic energies.")
    Gmicro = EnergyModel_CalculateMicrostateEnergy (self.cObject, vector.cObject, pH)
    return Gmicro


  def CalculateProbabilitiesAnalytically (self, Real pH=7.0):
    """Calculate probabilities of protonation states analytically."""
    cdef Status status
    status     = Status_Continue
    meadModel  = self.owner

    if self.cObject.nstates > ANALYTIC_STATES:
      raise CLibraryError ("Maximum number of states for analytic treatment (%d) exceeded." % ANALYTIC_STATES)
    if not meadModel.isCalculated:
      raise CLibraryError ("First calculate electrostatic energies.")

    EnergyModel_CalculateProbabilitiesAnalytically (self.cObject, pH, &status)
    if status != Status_Continue:
      raise CLibraryError ("Cannot allocate Boltzmann factors.")

    return self.cObject.nstates
  

  def CalculateProbabilitiesMonteCarlo (self, Real pH=7.0, Integer nequi=500, Integer nprod=20000, log=logFile):
    """Calculate probabilities of protonation states by using Metropolis Monte Carlo."""
    cdef Status  status
    status     = Status_Continue
    meadModel  = self.owner

    if not meadModel.isCalculated:
      raise CLibraryError ("First calculate electrostatic energies.")

    # Equilibration
    EnergyModel_CalculateProbabilitiesMonteCarlo (self.cObject, pH, CTrue,  nequi, &status)
    if LogFileActive (log):
      log.Text ("\nCompleted %d equilibration scans.\n" % nequi)

    # Production
    EnergyModel_CalculateProbabilitiesMonteCarlo (self.cObject, pH, CFalse, nprod, &status)
    if LogFileActive (log):
      log.Text ("\nCompleted %d production scans.\n" % nprod)


  def FindPairs (self, Real limit=2.0, log=logFile):
    """Identify pairs of sites for double moves """
    cdef Status  status = Status_Continue
    cdef Integer indexPair, indexSiteA, indexSiteB
    cdef Integer npairs
    cdef Real    Wmax
    npairs = EnergyModel_FindPairs (self.cObject, limit, -1, &status)

    if npairs > 0:
      EnergyModel_FindPairs (self.cObject, limit, npairs, &status)
      if status != Status_Continue:
        raise CLibraryError ("Cannot allocate pairs.")
  
      if LogFileActive (log):
        log.Text ("\nFound %d pair%s of strongly interacting sites:\n" % (npairs, "s" if npairs != 1 else ""))
        meadModel = self.owner
        sites     = meadModel.meadSites

        for indexPair from 0 <= indexPair < npairs:
          StateVector_GetPair (self.cObject.vector, indexPair, &indexSiteA, &indexSiteB, &Wmax, &status)
          siteFirst  = sites[indexSiteA]
          siteSecond = sites[indexSiteB]
          log.Text ("%4s %4s %4d -- %4s %4s %4d : %f\n" % (siteFirst.segName, siteFirst.resName, siteFirst.resSerial, siteSecond.segName, siteSecond.resName, siteSecond.resSerial, Wmax))
