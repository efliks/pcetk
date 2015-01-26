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


_DefaultDirectory = "curves"
_DefaultMethod    = "GMCT"
_DefaultSampling  =   .5
_DefaultStart     =  0.
_DefaultStop      = 14. 


class CurveThread (threading.Thread):
  """Calculate each pH-step in a separate thread."""

  def __init__ (self, curves, method, pH):
    model = curves.meadModel
    calculationMethods = { "GMCT"       : model.CalculateProbabilitiesGMCT          ,
                           "MonteCarlo" : model.CalculateProbabilitiesMonteCarlo    ,
                           "analytic"   : model.CalculateProbabilitiesAnalytically  }
    threading.Thread.__init__ (self)
    self.Calculate = calculationMethods[method]
    self.model     = model
    self.sites     = None
    self.pH        = pH

  def run (self):
    """The method that runs the calculations."""
    self.sites = self.Calculate (self.pH, log=None)


#-------------------------------------------------------------------------------
class TitrationCurves (object):
  """Titration curves."""

  def __init__ (self, meadModel, method=_DefaultMethod, curveSampling=_DefaultSampling, curveStart=_DefaultStart, curveStop=_DefaultStop):
    """Constructor."""
    if not meadModel.isCalculated:
      raise ContinuumElectrostaticsError ("First calculate electrostatic energies.")
    self.stop          =  curveStop
    self.start         =  curveStart
    self.sampling      =  curveSampling
    self.nsteps        =  int ((curveStop - curveStart) / curveSampling + 1)
    self.owner         =  meadModel
    self.method        =  method
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
  def CalculateCurves (self, forceSerial=False, printTable=False, log=logFile):
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
        calculationMethods = { "GMCT"       : meadModel.CalculateProbabilitiesGMCT          ,
                               "MonteCarlo" : meadModel.CalculateProbabilitiesMonteCarlo    ,
                               "analytic"   : meadModel.CalculateProbabilitiesAnalytically  }
        Calculate = calculationMethods[self.method]

        for step in range (self.nsteps):
          pH = self.start + step * self.sampling
          result = Calculate (pH=pH, log=None)
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
          batch.append (CurveThread (self, method, self.start + step * self.sampling))
          if len (batch) > limit:
            batches.append (batch)
            batch = []
        if batch:
          batches.append (batch)

        # If GMCT is to be used, first perform a dry run in serial mode to create directories and files
        if method == "GMCT":
          for step in range (self.nsteps):
            meadModel.CalculateProbabilitiesGMCT (pH=(self.start + step * self.sampling), dryRun=True, log=None)

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
              tab.Entry ("%10.2f" % (self.start + step * self.sampling))
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
            lines.append ("%f %f\n" % (self.start + step * self.sampling, phdata[step][site.siteIndex][instance.instIndex]))
          # Write instance data
          filename = os.path.join (directory, "%s_%s.dat" % (site.label, instance.label))
          WriteInputFile (filename, lines)

      if LogFileActive (log):
        log.Text ("\nWriting curve files complete.\n")


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
