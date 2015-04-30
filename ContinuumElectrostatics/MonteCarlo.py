#-------------------------------------------------------------------------------
# . File      : MonteCarlo.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MonteCarlo is a class for the handling of Monte Carlo sampling."""

from   pCore                 import logFile, LogFileActive
from   Error                 import ContinuumElectrostaticsError
from   InputFileWriter       import WriteInputFile
from   GMCTOutputFileReader  import GMCTOutputFileReader
from   Constants             import *
import os, subprocess


_DefaultPathGMCT           = "/usr/local/bin"
_DefaultDoubleFlip         = 2.
_DefaultTripleFlip         = 3.
_DefaultProductionScans    = 20000
_DefaultEquilibrationScans = 500

_DefaultSetupGMCT = """
blab        1
nconfflip   10
tlimit      %f
itraj       0
nmcfull     %d
temp        %f
icorr       0
limit       %f
nmcequi     %d
nmu         1
mu          %f  %f  0.0  0  0
"""


class MCModel (object):
    """A parent class."""

    defaultAttributes = {
        "owner"               :  None                       ,
        "label"               :  "MCModel"                  ,
        "doubleFlip"          :  _DefaultDoubleFlip         ,
        "productionScans"     :  _DefaultProductionScans    ,
        "equilibrationScans"  :  _DefaultEquilibrationScans ,
            }

    def __init__ (self, *arguments, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


    def Summary (self, log=logFile):
        """Summary."""
        if LogFileActive (log):
            summary = log.GetSummary ()
            summary.Start ("%s Monte Carlo sampling model"   % ("Default" if self.label == "MCModelDefault" else "GMCT"))
            summary.Entry ("Equilibration scans"    , "%d"   % self.equilibrationScans)
            summary.Entry ("Production scans"       , "%d"   % self.productionScans)
            summary.Entry ("Limit for double moves" , "%.1f" % self.doubleFlip)
            if hasattr (self, "tripleFlip"):
                summary.Entry ("Limit for triple moves" , "%.1f" % self.tripleFlip)
            summary.Stop ()


#-------------------------------------------------------------------------------
class MCModelDefault (MCModel):
    """The default, in-house Monte Carlo sampler."""

    def __init__ (self, *arguments, **keywordArguments):
        """Constructor."""
        super (MCModelDefault, self).__init__ (*arguments, **keywordArguments)
        self.label = "MCModelDefault"


    def CalculateOwnerProbabilities (self, pH=7.0, log=logFile):
        """Calculate probabilities of the owner."""
        owner = self.owner
        owner.energyModel.CalculateProbabilitiesMonteCarlo (pH=pH, nequi=self.equilibrationScans, nprod=self.productionScans, log=log)


#-------------------------------------------------------------------------------
class MCModelGMCT (MCModel):
    """Sampling by GMCT."""

    defaultAttributesChild = {
        "dryRun"      :  False               ,
        "pathGMCT"    :  _DefaultPathGMCT    ,
        "tripleFlip"  :  _DefaultTripleFlip  ,
            }

    def __init__ (self, *arguments, **keywordArguments):
        """Constructor."""
        super (MCModelGMCT, self).__init__ (*arguments, **keywordArguments)
        for (key, value) in self.__class__.defaultAttributesChild.iteritems (): setattr (self, key, value)
        for (key, value) in keywordArguments.iteritems (): setattr (self, key, value)

        self.pathGMCT = os.path.abspath (self.pathGMCT)
        self.label    = "MCmodelGMCT"


    def CalculateOwnerProbabilities (self, pH=7.0, log=logFile):
        """Calculate probabilities of the owner."""
        owner = self.owner

        if owner.isCalculated:
            sites       = None
            project     = "job"
            potential   = -CONSTANT_MOLAR_GAS_KCAL_MOL * owner.temperature * CONSTANT_LN10 * pH
            fileContent = _DefaultSetupGMCT % (self.tripleFlip, self.productionScans, owner.temperature, self.doubleFlip, self.equilibrationScans, potential, potential)

            # Prepare input files and directories for GMCT
            dirConf   = os.path.join ( owner.pathScratch , "gmct" , "conf"       )
            dirCalc   = os.path.join ( owner.pathScratch , "gmct" , "%s" % pH    )
            fileGint  = os.path.join ( dirConf ,           "%s.gint"  % project  )
            fileInter = os.path.join ( dirConf ,           "%s.inter" % project  )
            fileConf  = os.path.join ( dirCalc ,           "%s.conf"  % project  )
            fileSetup = os.path.join ( dirCalc ,           "%s.setup" % project  )
            linkname  = os.path.join ( dirCalc ,           "conf"                )
            if not os.path.exists ( dirConf   ): os.makedirs      ( dirConf                  )
            if not os.path.exists ( dirCalc   ): os.makedirs      ( dirCalc                  )
            if not os.path.exists ( fileGint  ): owner.WriteGintr ( fileGint  , precision=8  )
            if not os.path.exists ( fileInter ): owner.WriteW     ( fileInter , precision=8  )
            if not os.path.exists ( fileConf  ): WriteInputFile   ( fileConf  , ["conf  0.0  0.0  0.0\n"] )
            if not os.path.exists ( fileSetup ): WriteInputFile   ( fileSetup , fileContent  )
            if not os.path.exists ( linkname  ): os.symlink       ( "../conf" , linkname     )

            if not self.dryRun:
                output = os.path.join (dirCalc, "%s.gmct-out" % project)
                error  = os.path.join (dirCalc, "%s.gmct-err" % project)

                if os.path.exists (os.path.join (dirCalc, output)):
                    pass
                else:
                    command = [os.path.join (self.pathGMCT, "gmct"), project]
                    try:
                        out = open (output, "w")
                        err = open (error,  "w")
                        subprocess.check_call (command, stderr=err, stdout=out, cwd=dirCalc)
                        out.close ()
                        err.close ()
                    except:
                        raise ContinuumElectrostaticsError ("Failed running command: %s" % " ".join (command))

                # Read probabilities from the output file
                reader = GMCTOutputFileReader (output)
                reader.Parse (temperature=owner.temperature)

                # Copy probabilities from the reader to the owner
                for site in owner.meadSites:
                    for instance in site.instances:
                        key                   =  "conf_%s_%s%d_%s" % (site.segName, site.resName, site.resSerial, instance.label)
                        probability           =  reader.probabilities[key][0]
                        instance.probability  =  probability


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
