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
    cdef Real interac = EnergyModel_GetInteractionSymmetric (self.cObject, instIndexGlobalA, instIndexGlobalB)
    return interac

  def GetDeviation (self, Integer instIndexGlobalA, Integer instIndexGlobalB):
    cdef Real deviate = EnergyModel_GetDeviation (self.cObject, instIndexGlobalA, instIndexGlobalB)
    return deviate


  def __init__ (self, meadModel):
    """Constructor."""
    cdef Status  status  = Status_Continue
    cdef Integer nstates = 1, ninstances, totalInstances

    totalInstances = meadModel.ninstances
    self.cObject = EnergyModel_Allocate (totalInstances, &status)
    if status != Status_Continue:
      raise CLibraryError ("Cannot allocate energy model.")

    for site in meadModel.meadSites:
      ninstances = site.ninstances
      nstates = nstates * ninstances
      if nstates > ANALYTIC_STATES:
          break
    self.cObject.nstates    = nstates
    self.cObject.ninstances = totalInstances
    self.isOwner = False
    self.owner   = meadModel


  def CheckIfSymmetric (self, Real tolerance=0.05):
    """Check if the matrix of interactions is symmetric within the given threshold."""
    cdef Real    maxDeviation
    cdef Boolean isSymmetric
    isSymmetric = EnergyModel_CheckInteractionsSymmetric (self.cObject, tolerance, &maxDeviation)
    return (isSymmetric, maxDeviation)


  def SymmetrizeInteractions (self):
    """Symmetrize the matrix of interactions."""
    cdef Status status = Status_Continue
    EnergyModel_SymmetrizeInteractions (self.cObject, &status)

    # The following should never happen
    if status != Status_Continue:
      raise CLibraryError ("Cannot symmetrize interactions.")


  def CalculateMicrostateEnergy (self, StateVector vector, Real pH=7.0):
    """Calculate energy of a protonation state (=microstate)."""
    cdef Real Gmicro
    cdef Real temperature
    meadModel   = self.owner
    temperature = meadModel.temperature

    if not meadModel.isCalculated:
      raise CLibraryError ("First calculate electrostatic energies.")
    Gmicro = EnergyModel_CalculateMicrostateEnergy (self.cObject, vector.cObject, pH, temperature)
    return Gmicro


  def CalculateProbabilitiesAnalytically (self, Real pH=7.0):
    """Calculate probabilities of protonation states analytically."""
    cdef Status status
    cdef Real   temperature
    status      = Status_Continue
    meadModel   = self.owner
    temperature = meadModel.temperature

    if self.cObject.nstates > ANALYTIC_STATES:
      raise CLibraryError ("Maximum number of states for analytic treatment (%d) exceeded." % ANALYTIC_STATES)
    if not meadModel.isCalculated:
      raise CLibraryError ("First calculate electrostatic energies.")

    vector = StateVector (meadModel)
    EnergyModel_CalculateProbabilitiesAnalytically (self.cObject, vector.cObject, pH, temperature, &status)

    if status != Status_Continue:
      raise CLibraryError ("Cannot allocate Boltzmann factors.")
    return self.cObject.nstates
  

  def CalculateProbabilitiesMonteCarlo (self, Real pH=7.0, Integer nequi=500, Integer nprod=20000, log=logFile):
    """Calculate probabilities of protonation states by using Metropolis Monte Carlo."""
    cdef Real    temperature, scale, Gfinal
    cdef Integer nmoves, scan
    cdef Status  status
    cdef Boolean isLogActive
    status      = Status_Continue
    meadModel   = self.owner
    temperature = meadModel.temperature
    isLogActive = CFalse

    if not meadModel.isCalculated:
      raise CLibraryError ("First calculate electrostatic energies.")

    if LogFileActive (log):
      isLogActive = CTrue

    # Initialization
    vector = StateVector (meadModel)
    nmoves = vector.cObject.nsites
    scale  = 1. / nprod

    # Equilibration
    StateVector_Randomize (vector.cObject)
    for scan from 0 <= scan < nequi:
      Gfinal = EnergyModel_MCScan (self.cObject, vector.cObject, pH, temperature, nmoves)

    if isLogActive:
      log.Text ("\nCompleted %d equilibration scans.\n" % nequi)

    # Production
    Real1DArray_Set (self.cObject.probabilities, 0.)
    for scan from 0 <= scan < nprod:
      Gfinal = EnergyModel_MCScan (self.cObject, vector.cObject, pH, temperature, nmoves)
      EnergyModel_UpdateProbabilities (self.cObject, vector.cObject)
    Real1DArray_Scale (self.cObject.probabilities, scale)

    if isLogActive:
      log.Text ("\nCompleted %d production scans.\n" % nprod)


  def FindPairs (self, Real limit=2.0, log=logFile):
    """Identify pairs of sites for double moves """
    cdef Real W, Wmax
    cdef Integer npairs
    cdef Integer indexSiteA, indexSiteB, nsites
    cdef Integer indexInstA, indexInstB, indexGlobalInstA, indexGlobalInstB, ninsta, ninstb
    meadModel = self.owner
    nsites    = meadModel.nsites
    sites     = meadModel.meadSites
    npairs    = 0

    for indexSiteA from 0 <= indexSiteA < nsites: 
      sitea  = sites[indexSiteA]
      ninsta = sitea.ninstances
 
      for indexSiteB from 0 <= indexSiteB < nsites:
        siteb  = sites[indexSiteB]
        ninstb = siteb.ninstances
        Wmax   = 0.
      
        for indexInstA from 0 <= indexInstA < ninsta: 
          insta = sitea.instances[indexInstA]
          indexGlobalInstA = insta._instIndexGlobal
      
          for indexInstB from 0 <= indexInstB < ninstb: 
            instb = siteb.instances[indexInstB]
            indexGlobalInstB = instb._instIndexGlobal
      
            W = EnergyModel_GetInteractionSymmetric (self.cObject, indexGlobalInstA, indexGlobalInstB)
            if W < 0.   : W = -1. * W
            if W > Wmax : Wmax = W
      
        if Wmax >= limit:
          npairs += 1

    if LogFileActive (log):
      log.Text ("\nFound %d pair(s) of strongly interacting sites.\n" % npairs)
