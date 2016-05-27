#-------------------------------------------------------------------------------
# . File      : InstanceMEAD.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore                  import logFile, LogFileActive
from Error                  import ContinuumElectrostaticsError
from Instance               import Instance
from MEADOutputFileReader   import MEADOutputFileReader

import os, subprocess


class InstanceMEAD (Instance):
    """A class to represent a MEAD type of instance."""

    defaultAttributes = {
        }
    defaultAttributes.update (Instance.defaultAttributes)
    # modelPqr  modelLog  modelGrid  sitePqr  siteLog  siteGrid

    def __init__ (self, **keywordArguments):
        """Constructor."""
        super (InstanceMEAD, self).__init__ (**keywordArguments)


    #-------------------------------------------------------------------------------
    def CalculateModelCompound (self, log=logFile):
        """Calculate Gborn and Gback of a site in a model compound."""
        site  = self.parent
        model = site.parent

        if model.isFilesWritten:
            if os.path.exists (self.modelLog):
                pass
            else:
                instancePqr        , ext = os.path.splitext (self.sitePqr)
                modelBackgroundPqr , ext = os.path.splitext (self.modelPqr)
                command = [
                    os.path.join (model.pathMEAD, "my_2diel_solver"), 
                    "-T", "%f" % model.temperature, 
                    "-ionicstr", "%f" % model.ionicStrength, 
                    "-epsin", "%f" % model.epsilonProtein, 
                    "-epsext", "%f" % model.epsilonWater, 
                    instancePqr, 
                    modelBackgroundPqr
                    ]
                try:
                    outFile = open (self.modelLog, "w")
                    subprocess.check_call (command, stderr=outFile, stdout=outFile)
                    outFile.close ()
                except:
                    raise ContinuumElectrostaticsError ("Failed running command: %s" % " ".join (command))
            reader = MEADOutputFileReader (self.modelLog)
            reader.Parse ()

            checks = (hasattr (reader, "born"), hasattr (reader, "back"), )
            if not all (checks):
                raise ContinuumElectrostaticsError ("Output file %s empty or corrupted. Empty the scratch directory and start anew." % self.modelLog)
            self.Gborn_model = reader.born
            self.Gback_model = reader.back


    #-------------------------------------------------------------------------------
    def CalculateProtein (self, log=logFile):
        """Calculate Gborn, Gback and Wij of a site in protein environment."""
        site  = self.parent
        model = site.parent

        if model.isFilesWritten:
            if os.path.exists (self.siteLog):
                pass
            else:
                # . Assign removing extensions, otherwise MEAD does not work
                sitesFpt             , ext = os.path.splitext (model.pathFptSites)
                proteinPqr           , ext = os.path.splitext (model.pathPqrProtein)
                proteinBackgroundPqr , ext = os.path.splitext (model.pathPqrBack)
                instancePqr          , ext = os.path.splitext (self.sitePqr)

                # . epsin1 is never used but must be given

                # . eps2set defines the whole protein
                command = [
                    os.path.join (model.pathMEAD, "my_3diel_solver"), 
                    "-T", "%f" % model.temperature, 
                    "-ionicstr", "%f" % model.ionicStrength, 
                    "-epsin1", "%f" % 1.0, 
                    "-epsin2", "%f" % model.epsilonProtein, 
                    "-epsext", "%f" % model.epsilonWater, 
                    "-eps2set", "%s" % proteinPqr, 
                    "-fpt", "%s" % sitesFpt, 
                    instancePqr, 
                    proteinBackgroundPqr
                    ]
                try:
                    outFile = open (self.siteLog, "w")
                    subprocess.check_call (command, stderr=outFile, stdout=outFile)
                    outFile.close ()
                except:
                    raise ContinuumElectrostaticsError ("Failed running command: %s" % " ".join (command))
            reader = MEADOutputFileReader (self.siteLog)
            reader.Parse ()

            checks = (hasattr (reader, "born"), hasattr (reader, "back"), hasattr (reader, "interactions"), )
            if not all (checks):
                raise ContinuumElectrostaticsError ("Output file %s empty or corrupted. Empty the scratch directory and start anew." % self.modelLog)
            self.Gborn_protein = reader.born
            self.Gback_protein = reader.back

            # . Create a list of interactions
            interactions    = []
            instances       = []
            siteIndexOld    = 99999
            parentSiteIndex = self.parent.siteIndex

            for siteIndex, instanceIndex, energy in reader.interactions:
                if siteIndex > siteIndexOld:
                    interactions.append (instances)
                    instances = []
                # . Set the interaction energy to zero if the site is interacting with itself
                if siteIndex == parentSiteIndex:
                    energy = 0.
                instances.append (energy)
                siteIndexOld = siteIndex

            if instances:
                interactions.append (instances)

            # . Copy the interactions to the centralized array
            indexGlobal = 0
            for site in interactions:
                for instance in site:
                    energy = instance
                    model.energyModel.SetInteraction (self._instIndexGlobal, indexGlobal, energy)
                    indexGlobal = indexGlobal + 1


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
