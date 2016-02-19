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
    def SetGmodel (self, Integer instIndexGlobal, Real value):
        EnergyModel_SetGmodel (self.cObject, instIndexGlobal, value)

    def SetGintr (self, Integer instIndexGlobal, Real value):
        EnergyModel_SetGintr (self.cObject, instIndexGlobal, value)

    def SetProtons (self, Integer instIndexGlobal, Integer value):
        EnergyModel_SetProtons (self.cObject, instIndexGlobal, value)

    def SetProbability (self, Integer instIndexGlobal, Real value):
        EnergyModel_SetProbability (self.cObject, instIndexGlobal, value)

    def SetInteraction (self, Integer instIndexGlobalA, Integer instIndexGlobalB, Real value):
        EnergyModel_SetInteraction (self.cObject, instIndexGlobalA, instIndexGlobalB, value)

    # Getters
    def GetGmodel (self, Integer instIndexGlobal):
        cdef Real Gmodel = EnergyModel_GetGmodel (self.cObject, instIndexGlobal)
        return Gmodel

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


    def __init__ (self, ceModel, Integer totalSites, Integer totalInstances):
        """Constructor."""
        cdef Status  status
        status        = Status_Continue
        self.cObject  = EnergyModel_Allocate (totalSites, totalInstances, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot allocate energy model.")

        self.cObject.temperature = ceModel.temperature
        self.cObject.ninstances  = totalInstances
        self.isOwner             = False
        self.owner               = ceModel


    def Initialize (self):
        """Set up sites and calculate the number of possible states."""
        cdef Integer nstates, ninstances, totalSites
        cdef Integer indexSite, indexDown, indexUp, index
        cdef Status  status
        totalSites = self.cObject.vector.nsites
        status     = Status_Continue
        ceModel    = self.owner
        nstates    = 1

        for indexSite from 0 <= indexSite < totalSites:
            ninstances =  0
            indexUp    =  0
            indexDown  =  9999
            site       =  ceModel.sites[indexSite]
            for instance in site.instances:
                ninstances += 1
                index = instance._instIndexGlobal
                if index < indexDown : indexDown = index
                if index > indexUp   : indexUp   = index
            StateVector_SetSite (self.cObject.vector, indexSite, indexDown, indexUp, &status)

            if nstates <= ANALYTIC_STATES:
                nstates = nstates * ninstances
        self.cObject.nstates = nstates


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


    def ResetInteractions (self):
        """Set all interactions to zero."""
        EnergyModel_ResetInteractions (self.cObject)


    def CalculateMicrostateEnergy (self, StateVector vector, Real pH=7.0):
        """Calculate energy of a protonation state (=microstate)."""
        cdef Real Gmicro
        ceModel = self.owner

        if not ceModel.isCalculated:
            raise CLibraryError ("First calculate electrostatic energies.")
        Gmicro = EnergyModel_CalculateMicrostateEnergy (self.cObject, vector.cObject, pH)
        return Gmicro


    def CalculateProbabilitiesAnalytically (self, Real pH=7.0):
        """Calculate probabilities of protonation states analytically."""
        cdef Status status
        status   = Status_Continue
        ceModel  = self.owner

        if self.cObject.nstates > ANALYTIC_STATES:
            raise CLibraryError ("Maximum number of states for analytic treatment (%d) exceeded." % ANALYTIC_STATES)
        if not ceModel.isCalculated:
            raise CLibraryError ("First calculate electrostatic energies.")

        EnergyModel_CalculateProbabilitiesAnalytically (self.cObject, pH, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot allocate Boltzmann factors.")
        return self.cObject.nstates


    def CalculateProbabilitiesAnalyticallyUnfolded (self, Real pH=7.0):
        """Calculate probabilities of protonation states analytically (unfolded protein)."""
        cdef Status status
        status   = Status_Continue
        ceModel  = self.owner

        if self.cObject.nstates > ANALYTIC_STATES:
            raise CLibraryError ("Maximum number of states for analytic treatment (%d) exceeded." % ANALYTIC_STATES)
        if not ceModel.isCalculated:
            raise CLibraryError ("First calculate electrostatic energies.")

        EnergyModel_CalculateProbabilitiesAnalyticallyUnfolded (self.cObject, pH, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot allocate Boltzmann factors.")
        return self.cObject.nstates


    def CalculateMicrostateEnergyUnfolded (self, StateVector vector, Real pH=7.0):
        """Calculate energy of a protonation state (=microstate) in an unfolded protein."""
        cdef Real Gmicro
        ceModel = self.owner
        if not ceModel.isCalculated:
            raise CLibraryError ("First calculate electrostatic energies.")
        Gmicro = EnergyModel_CalculateMicrostateEnergyUnfolded (self.cObject, vector.cObject, pH)
        return Gmicro


    def CalculateZfolded (self, Real Gzero=0.0, Real pH=7.0):
        """Calculate partition function of a folded protein."""
        cdef Status  status = Status_Continue
        cdef Real    Zfolded
        if self.cObject.nstates > ANALYTIC_STATES:
            raise CLibraryError ("Maximum number of states for analytic treatment (%d) exceeded." % ANALYTIC_STATES)

        Zfolded = EnergyModel_CalculateZfolded (self.cObject, pH, Gzero, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot allocate Boltzmann factors.")
        return Zfolded


    def CalculateZunfolded (self, Real Gzero=0.0, Real pH=7.0):
        """Calculate partition function of an unfolded protein."""
        cdef Status  status = Status_Continue
        cdef Real    Zunfolded
        if self.cObject.nstates > ANALYTIC_STATES:
            raise CLibraryError ("Maximum number of states for analytic treatment (%d) exceeded." % ANALYTIC_STATES)

        Zunfolded = EnergyModel_CalculateZunfolded (self.cObject, pH, Gzero, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot allocate Boltzmann factors.")
        return Zunfolded
