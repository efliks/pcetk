#-------------------------------------------------------------------------------
# . File      : ContinuumElectrostatics.MCModelDefault.pyx
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2018)
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


    def __init__ (self, Real doubleFlip=_DefaultDoubleFlip, Integer nequil=_DefaultEquilibrationScans, Integer nprod=_DefaultProductionScans, Integer randomSeed=-1):
        """Constructor."""
        cdef Status  status
        status        = Status_Continue
        self.isOwner  = True
        self.cObject  = MCModelDefault_Allocate (doubleFlip, nequil, nprod, randomSeed, &status)
        if status != Status_Continue:
            raise CLibraryError ("Cannot allocate Monte Carlo model.")


    def Initialize (self, ceModel):
        """Link Monte Carlo and energy models."""
        cdef EnergyModel energyModel = ceModel.energyModel
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
        self.owner   = ceModel


    def CalculateOwnerProbabilities (self, Real pH=7.0, Integer logFrequency=0, trajectoryFilename="", log=logFile):
        """Calculate probabilities of the owner."""
        cdef CMCModelDefault   *mcModel
        cdef Integer  moves, movesAcc, flips, flipsAcc
        cdef Integer  nmoves, i, j
        cdef Real     scale, Gmicro
        cdef Boolean  active, writing

        owner = self.owner
        if (not owner.isCalculated):
            raise CLibraryError ("First calculate electrostatic energies.")

        if (logFrequency <= 0) and (trajectoryFilename == ""):
            # . Do a quiet run
            active = CTrue if (LogFileActive (log)) else CFalse

            MCModelDefault_Equilibration (self.cObject, pH)
            if active:
                log.Text ("\nCompleted %d equilibration scans.\n" % self.cObject.nequil)

            MCModelDefault_Production (self.cObject, pH)
            if active:
                log.Text ("\nCompleted %d production scans.\n" % self.cObject.nprod)
        else:
            # . Do logging or trajectory writing
            mcModel = self.cObject
            nmoves  = mcModel.vector.nsites + mcModel.vector.npairs
            active  = CTrue if (LogFileActive (log) and (logFrequency > 0)) else CFalse
    
            if active:
                j         =  0
                moves     =  0
                flips     =  0
                movesAcc  =  0
                flipsAcc  =  0
                table  = log.GetTable (columns=[8, 10, 10, 10, 10])
                table.Start ()
                table.Heading ("%8s"   % "Scan"  )
                table.Heading ("%10s"  % "Moves" )
                table.Heading ("%10s"  % "M-Acc" )
                table.Heading ("%10s"  % "Flips" )
                table.Heading ("%10s"  % "F-Acc" )
    
            # . Do equilibration phase
            StateVector_Randomize (mcModel.vector, mcModel.generator)

            for i from 0 <= i < mcModel.nequil:
                MCModelDefault_MCScan (mcModel, pH, nmoves, &moves, &movesAcc, &flips, &flipsAcc)
                if active:
                    j += 1
                    if j >= logFrequency:
                        table.Entry ("%8d"    % (i + 1)   )
                        table.Entry ("%10d"   % moves     )
                        table.Entry ("%10d"   % movesAcc  )
                        table.Entry ("%10d"   % flips     )
                        table.Entry ("%10d"   % flipsAcc  )
                        j         =  0
                        moves     =  0
                        flips     =  0
                        movesAcc  =  0
                        flipsAcc  =  0
            if active:
                j        =  0
                moves    =  0
                flips    =  0
                movesAcc =  0
                flipsAcc =  0
                table.Entry ("** Equilibration phase done **".center (8 + 10 * 4), columnSpan=5)
   
            # . Do production phase
            writing = CFalse
            if (trajectoryFilename != ""):
                writing = CTrue
                output = open (trajectoryFilename, "w")
            Real1DArray_Set (mcModel.energyModel.probabilities, 0.)
 
            for i from 0 <= i < mcModel.nprod:
                Gmicro = MCModelDefault_MCScan (mcModel, pH, nmoves, &moves, &movesAcc, &flips, &flipsAcc)
                if writing:
                    output.write ("%f\n" % Gmicro)

                MCModelDefault_UpdateProbabilities (mcModel)
                if active:
                    j += 1
                    if j >= logFrequency:
                        table.Entry ("%8d"    % (i + 1)   )
                        table.Entry ("%10d"   % moves     )
                        table.Entry ("%10d"   % movesAcc  )
                        table.Entry ("%10d"   % flips     )
                        table.Entry ("%10d"   % flipsAcc  )
                        j         =  0
                        moves     =  0
                        flips     =  0
                        movesAcc  =  0
                        flipsAcc  =  0
            if active:
                table.Entry ("** Production phase done **".center (8 + 10 * 4), columnSpan=5)
                table.Stop ()

            # . Finalize
            if writing:
                output.close ()
            scale = 1. / mcModel.nprod
            Real1DArray_Scale (mcModel.energyModel.probabilities, scale)


    def Summary (self, log=logFile):
        """Summary."""
        cdef Integer nequil     = self.cObject.nequil
        cdef Integer nprod      = self.cObject.nprod
        cdef Real    doubleFlip = self.cObject.limit
        if LogFileActive (log):
            summary = log.GetSummary ()
            summary.Start ("Default Monte Carlo sampling model")
            summary.Entry ("Equilibration scans"    , "%d"   % nequil)
            summary.Entry ("Production scans"       , "%d"   % nprod)
            summary.Entry ("Limit for double moves" , "%.1f" % doubleFlip)
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
            ceModel = self.owner
            sites   = ceModel.sites

            for indexPair from 0 <= indexPair < npairs:
                StateVector_GetPair (self.cObject.vector, indexPair, &indexSiteA, &indexSiteB, &Wmax, &status)
                siteFirst  = sites[indexSiteA]
                siteSecond = sites[indexSiteB]
                log.Text ("%4s %4s %4d -- %4s %4s %4d : %f\n" % (siteFirst.segName, siteFirst.resName, siteFirst.resSerial, siteSecond.segName, siteSecond.resName, siteSecond.resSerial, Wmax))
