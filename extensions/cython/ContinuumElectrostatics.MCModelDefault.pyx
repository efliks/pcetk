#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.MCModelDefault.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore  import logFile, LogFileActive, CLibraryError

_DefaultDoubleFlip         = 2.
_DefaultProductionScans    = 20000
_DefaultEquilibrationScans = 500


cdef class MCModelDefault:
    """A class defining the default Monte Carlo model."""

    def __getmodule__ (self):
        """Return the module name."""
        return "ContinuumElectrostatics.MCModelDefault"

    def __dealloc__ (self):
        """Deallocate."""
        MCModelDefault_Deallocate (self.cObject)


    def __init__ (self, Real limit=_DefaultDoubleFlip, Integer nequil=_DefaultEquilibrationScans, Integer nprod=_DefaultProductionScans, Integer randomSeed=-1):
        """Constructor."""
        cdef Status  status
        status        = Status_Continue
        self.isOwner  = True
        self.cObject  = MCModelDefault_Allocate (limit, nequil, nprod, randomSeed, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot allocate Monte Carlo model.")


    def Initialize (self, meadModel):
        """Link Monte Carlo and energy models."""
        cdef EnergyModel energyModel = meadModel.energyModel
        cdef Status      status      = Status_Continue
        cdef Integer     npairs

        MCModelDefault_LinkToEnergyModel (self.cObject, energyModel.cObject, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot initialize Monte Carlo model.")

        npairs = MCModelDefault_FindPairs (self.cObject, -1, &status)
        if npairs > 0:
            MCModelDefault_FindPairs (self.cObject, npairs, &status)
            if status != Status_Continue:
                raise CLibraryError ("Cannot allocate pairs.")
        self.isOwner = False
        self.owner   = meadModel


    def CalculateOwnerProbabilities (self, Real pH=7.0, log=logFile):
        """Calculate probabilities of the owner"""
        owner = self.owner
        if not owner.isCalculated:
            raise CLibraryError ("First calculate electrostatic energies.")

        MCModelDefault_CalculateProbabilities (self.cObject, pH, CTrue)
        if LogFileActive (log):
            log.Text ("\nCompleted %d equilibration scans.\n" % self.cObject.nequil)
        
        MCModelDefault_CalculateProbabilities (self.cObject, pH, CFalse)
        if LogFileActive (log):
            log.Text ("\nCompleted %d production scans.\n" % self.cObject.nprod)


    def Summary (self, log=logFile):
        """Summary."""
        cdef Integer nequil = self.cObject.nequil
        cdef Integer nprod  = self.cObject.nprod
        cdef Real    limit  = self.cObject.limit

        if LogFileActive (log):
            summary = log.GetSummary ()
            summary.Start ("Default Monte Carlo sampling model")
            summary.Entry ("Equilibration scans"    , "%d"   % nequil)
            summary.Entry ("Production scans"       , "%d"   % nprod)
            summary.Entry ("Limit for double moves" , "%.1f" % limit)
            summary.Stop ()


    def PrintPairs (self, log=logFile):
        """Summary of strongly interacting pairs."""
        cdef Integer  indexPair, indexSiteA, indexSiteB
        cdef Status   status
        cdef Integer  npairs
        cdef Real     Wmax
        status = Status_Continue
        npairs = self.cObject.vector.npairs

        if LogFileActive (log):
            log.Text ("\nFound %d pair%s of strongly interacting sites:\n" % (npairs, "s" if npairs != 1 else ""))
            meadModel = self.owner
            sites     = meadModel.meadSites

            for indexPair from 0 <= indexPair < npairs:
                StateVector_GetPair (self.cObject.vector, indexPair, &indexSiteA, &indexSiteB, &Wmax, &status)
                siteFirst  = sites[indexSiteA]
                siteSecond = sites[indexSiteB]
                log.Text ("%4s %4s %4d -- %4s %4s %4d : %f\n" % (siteFirst.segName, siteFirst.resName, siteFirst.resSerial, siteSecond.segName, siteSecond.resName, siteSecond.resSerial, Wmax))
