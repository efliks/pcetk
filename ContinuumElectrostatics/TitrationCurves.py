#-------------------------------------------------------------------------------
# . File      : TitrationCurves.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""TitrationCurves is a class for calculating titration curves.

CurveThread is a class for running parallel calculations of titration curves."""

__lastchanged__ = "$Id$"


from pCore           import logFile, LogFileActive
from Error           import ContinuumElectrostaticsError
from InputFileWriter import WriteInputFile

import os, threading


_DefaultDirectory          = "curves"
_DefaultMethod             = "MonteCarlo"
_DefaultSampling           =   .5
_DefaultStart              =  0.
_DefaultStop               = 14.
_DefaultDoubleFlip         =  2.
_DefaultTripleFlip         =  3.
_DefaultEquilibrationScans = 500
_DefaultProductionScans    = 20000


class CurveThread (threading.Thread):
    """Calculate each pH-step in a separate thread."""

    def __init__ (self, curves, pH):
        threading.Thread.__init__ (self)
        self.curves = curves
        self.sites  = None
        self.pH     = pH

    def run (self):
        """The method that runs the calculations."""
        curves = self.curves
        model  = curves.owner
        method = curves.method
        if method == "GMCT":
            self.sites = model.CalculateProbabilitiesGMCT         (pH=self.pH, log=None, nequi=curves.mcEquilibrationScans, nprod=curves.mcProductionScans)
        elif method == "MonteCarlo":
            self.sites = model.CalculateProbabilitiesMonteCarlo   (pH=self.pH, log=None, nequi=curves.mcEquilibrationScans, nprod=curves.mcProductionScans)
        elif method == "analytically":
            self.sites = model.CalculateProbabilitiesAnalytically (pH=self.pH, log=None)
        else:
            self.sites = model.CalculateProbabilitiesAnalyticallyUnfolded (pH=self.pH, log=None)


#-------------------------------------------------------------------------------
class TitrationCurves (object):
    """Titration curves."""

    defaultAttributes = {
        "method"               :  _DefaultMethod              ,
        "curveSampling"        :  _DefaultSampling            ,
        "curveStart"           :  _DefaultStart               ,
        "curveStop"            :  _DefaultStop                ,
        "mcEquilibrationScans" :  _DefaultEquilibrationScans  ,
        "mcProductionScans"    :  _DefaultProductionScans     ,
        "mcDoubleLimit"        :  _DefaultDoubleFlip          ,
        "mcTripleLimit"        :  _DefaultTripleFlip          ,
                        }

    def __init__ (self, meadModel, log=logFile, *arguments, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)

        if not meadModel.isCalculated:
            raise ContinuumElectrostaticsError ("First calculate electrostatic energies.")
        self.owner         =  meadModel
        self.nsteps        =  int ((self.curveStop - self.curveStart) / self.curveSampling + 1)
        self.steps         =  None
        self.isCalculated  =  False


    #===============================================================================
    def _ProbabilitiesSave (self):
        """Backup the probabilities existing in the model."""
        meadModel     = self.owner
        energyModel   = meadModel.energyModel
        probabilities = [0.] * meadModel.ninstances

        for site in meadModel.meadSites:
            for instance in site.instances:
                probabilities[instance._instIndexGlobal] = energyModel.GetProbability (instance._instIndexGlobal)
        return probabilities


    def _ProbabilitiesRestore (self, probabilities):
        """Restore the original probabilities to the model."""
        meadModel   = self.owner
        energyModel = meadModel.energyModel

        for site in meadModel.meadSites:
            for instance in site.instances:
                energyModel.SetProbability (instance._instIndexGlobal, probabilities[instance._instIndexGlobal])


    #===============================================================================
    def CalculateCurves (self, forceSerial=True, printTable=False, log=logFile):
        """Calculate titration curves."""
        if not self.isCalculated:
            probabilities  =  self._ProbabilitiesSave ()
            meadModel      =  self.owner
            nthreads       =  meadModel.nthreads
            steps          =  []
            tab            =  None

            if LogFileActive (log):
                if nthreads < 2 or forceSerial:
                    log.Text ("\nStarting serial run.\n")
                else:
                    log.Text ("\nStarting parallel run on %d CPUs.\n" % nthreads)

                if printTable:
                    tab = log.GetTable (columns=[10, 10])
                    tab.Start ()
                    tab.Heading ("Step")
                    tab.Heading ("pH")


            # Serial run?
            if nthreads < 2 or forceSerial:
                for step in range (self.nsteps):
                    pH = self.curveStart + step * self.curveSampling
                    if self.method == "GMCT":
                        result = meadModel.CalculateProbabilitiesGMCT         (pH=pH, log=None, nequi=self.mcEquilibrationScans, nprod=self.mcProductionScans)
                    elif self.method == "MonteCarlo":
                        result = meadModel.CalculateProbabilitiesMonteCarlo   (pH=pH, log=None, nequi=self.mcEquilibrationScans, nprod=self.mcProductionScans)
                    elif self.method == "analytically":
                        result = meadModel.CalculateProbabilitiesAnalytically (pH=pH, log=None)
                    else:
                        result = meadModel.CalculateProbabilitiesAnalyticallyUnfolded (pH=pH, log=None)

                    steps.append (result)
                    if tab:
                        tab.Entry ("%10d"   % step)
                        tab.Entry ("%10.2f" % pH)

            # Parallel run?
            else:
                limit   = nthreads - 1
                batches = []
                batch   = []
                for step in range (self.nsteps):
                    batch.append (CurveThread (self, self.curveStart + step * self.curveSampling))
                    if len (batch) > limit:
                        batches.append (batch)
                        batch = []
                if batch:
                    batches.append (batch)

                # If GMCT is to be used, first perform a dry run in serial mode to create directories and files
                if method == "GMCT":
                    for step in range (self.nsteps):
                        meadModel.CalculateProbabilitiesGMCT (pH=(self.curveStart + step * self.curveSampling), dryRun=True, log=None)

                step = 0
                for batch in batches:
                    for thread in batch:
                        thread.start ()
                    for thread in batch:
                        thread.join ()

                    for thread in batch:
                        steps.append (thread.sites)
                        if tab:
                            tab.Entry ("%10d"   % step)
                            tab.Entry ("%10.2f" % (self.curveStart + step * self.curveSampling))
                        step = step + 1

            if LogFileActive (log):
                if tab:
                    tab.Stop ()
                log.Text ("\nCalculating titration curves complete.\n")

            # Finalize
            self.steps = steps
            self.isCalculated = True
            self._ProbabilitiesRestore (probabilities)


    #===============================================================================
    def WriteCurves (self, directory=_DefaultDirectory, log=logFile):
        """Write calculated curves to a directory."""
        if self.isCalculated:
            # Initialize output directory
            if not os.path.exists (directory):
                try:
                    os.mkdir (directory)
                except:
                    raise ContinuumElectrostaticsError ("Cannot create directory %s" % directory)

            # For each instance of each site, write a curve file
            meadModel = self.owner
            phdata    = self.steps

            for site in meadModel.meadSites:
                for instance in site.instances:
                    # Collect instance data
                    lines = []
                    for step in range (self.nsteps):
                        lines.append ("%f %f\n" % (self.curveStart + step * self.curveSampling, phdata[step][site.siteIndex][instance.instIndex]))
                    # Write instance data
                    filename = os.path.join (directory, "%s_%s.dat" % (site.label, instance.label))
                    WriteInputFile (filename, lines)

            if LogFileActive (log):
                log.Text ("\nWriting curve files complete.\n")


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
