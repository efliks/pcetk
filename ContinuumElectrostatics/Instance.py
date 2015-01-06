#-------------------------------------------------------------------------------
# . File      : Instance.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADInstance is a class representing a particular instance of a site.

InstanceThread is a class used for running parallel calculations in MEAD."""

__lastchanged__ = "$Id$"


from pCore                  import logFile, LogFileActive
from Error                  import ContinuumElectrostaticsError
from MEADOutputFileReader   import MEADOutputFileReader

import os, threading, subprocess, time


#-------------------------------------------------------------------------------
class InstanceThread (threading.Thread):
  """A class for running parallel calculations for sites both in model compounds and in the protein.

  It also calculates Gintr."""

  def __init__ (self, instance, log=logFile):
    """Constructor."""
    threading.Thread.__init__ (self)
    self.instance = instance
    self.log      = log
    self.time     = 0

  def run (self):
    """The method that runs the calculations."""
    time0 = time.time ()
    self.instance.CalculateSiteInModelCompound (log=self.log)
    self.instance.CalculateSiteInProtein       (log=self.log)
    self.instance.CalculateGintr               (log=self.log)

    self.time = (time.time () - time0)


#-------------------------------------------------------------------------------
class MEADInstance (object):
  """Instance of a titratable site.

  For the time being, the instances only differ in charges. There are no rotameric instances."""
  defaultAttributes = {
                  "parent"            :  None  ,
                  "instIndex"         :  None  ,
                  "_instIndexGlobal"  :  None  ,
                  "label"             :  None  ,
                  "charges"           :  None  ,
                  "modelPqr"          :  None  ,
                  "modelLog"          :  None  ,
                  "modelGrid"         :  None  ,
                  "sitePqr"           :  None  ,
                  "siteLog"           :  None  ,
                  "siteGrid"          :  None  ,
                  "Gmodel"            :  None  ,
                  "Gborn_model"       :  None  ,
                  "Gback_model"       :  None  ,
                  "Gborn_protein"     :  None  ,
                  "Gback_protein"     :  None  ,
                      }

  @property
  def Gintr (self):
    site        = self.parent
    meadModel   = site.parent
    energyModel = meadModel.energyModel
    if energyModel:
      return energyModel.GetGintr (self._instIndexGlobal)
    else:
      return None

  @Gintr.setter
  def Gintr (self, value):
    site        = self.parent
    meadModel   = site.parent
    energyModel = meadModel.energyModel
    if energyModel:
      energyModel.SetGintr (self._instIndexGlobal, value)
    else:
      pass

  @property
  def protons (self):
    site        = self.parent
    meadModel   = site.parent
    energyModel = meadModel.energyModel
    if energyModel:
      return energyModel.GetProtons (self._instIndexGlobal)
    else:
      return None

  @protons.setter
  def protons (self, value):
    site        = self.parent
    meadModel   = site.parent
    energyModel = meadModel.energyModel
    if energyModel:
      energyModel.SetProtons (self._instIndexGlobal, value)
    else:
      pass

  @property
  def interactions (self):
    # Return a list of interactions of this particular instance with other instances
    site        = self.parent
    meadModel   = site.parent
    energyModel = meadModel.energyModel
    if energyModel:
      ninstances = model.ninstances
      energies   = [0.] * ninstances
      for _instIndexGlobal in xrange (ninstances):
        energies[_instIndexGlobal] = energyModel.GetInteractionSymmetric (self._instIndexGlobal, _instIndexGlobal)
      return energies
    else:
      return None

  @property
  def probability (self):
    site        = self.parent
    meadModel   = site.parent
    energyModel = meadModel.energyModel
    if energyModel:
      return energyModel.GetProbability (self._instIndexGlobal)
    else:
      return None

  @probability.setter
  def probability (self, value):
    site        = self.parent
    meadModel   = site.parent
    energyModel = meadModel.energyModel
    if energyModel:
      energyModel.SetProbability (self._instIndexGlobal, value)
    else:
      pass


  #===============================================================================
  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


  #===============================================================================
  def CalculateSiteInModelCompound (self, log=logFile):
    """Calculate Gborn and Gback of a site in the model compound."""
    site  = self.parent
    model = site.parent

    if model.isFilesWritten:
      if os.path.exists (self.modelLog):
        pass
      else:
        instancePqr        = self.sitePqr  [:-4]
        modelBackgroundPqr = self.modelPqr [:-4]
        command = [os.path.join (model.pathMEAD, "my_2diel_solver"), "-T", "%f" % model.temperature, "-ionicstr", "%f" % model.ionicStrength, "-epsin", "%f" % 4.0, instancePqr, modelBackgroundPqr]

        try:
          outFile = open (self.modelLog, "w")
          subprocess.check_call (command, stderr = outFile, stdout = outFile)
          outFile.close ()
        except:
          raise ContinuumElectrostaticsError ("Failed running command: %s" % " ".join (command))
      reader = MEADOutputFileReader (self.modelLog)
      reader.Parse ()

      checks = [hasattr (reader, "born"), hasattr (reader, "back")]
      if not all (checks):
        raise ContinuumElectrostaticsError ("Output file %s empty or corrupted. Empty the scratch directory and start anew." % self.modelLog)
      self.Gborn_model = reader.born
      self.Gback_model = reader.back


  #===============================================================================
  def CalculateSiteInProtein (self, log=logFile):
    """Calculate Gborn and Gback of a site in protein environment.

    Also, calculate Wij between the site and the other sites.

    Use fpt-file for other instances of sites."""
    site  = self.parent
    model = site.parent

    if model.isFilesWritten:
      if os.path.exists (self.siteLog):
        pass
      else:
        # Assign removing extensions, otherwise MEAD does not work
        sitesFpt             = model.pathFptSites   [:-4]
        proteinPqr           = model.pathPqrProtein [:-4]
        proteinBackgroundPqr = model.pathPqrBack    [:-4]
        instancePqr          = self.sitePqr         [:-4]

        # epsin1 is never used but must be given

        # eps2set defines the whole protein

        command = [os.path.join (model.pathMEAD, "my_3diel_solver"), "-T", "%f" % model.temperature, "-ionicstr", "%f" % model.ionicStrength, "-epsin1", "%f" % 1.0, "-epsin2", "%f" % 4.0, "-eps2set", "%s" % proteinPqr, "-fpt", "%s" % sitesFpt, instancePqr, proteinBackgroundPqr]

        try:
          outFile = open (self.siteLog, "w")
          subprocess.check_call (command, stderr = outFile, stdout = outFile)
          outFile.close ()
        except:
          raise ContinuumElectrostaticsError ("Failed running command: %s" % " ".join (command))
      reader = MEADOutputFileReader (self.siteLog)
      reader.Parse ()

      checks = [hasattr (reader, "born"), hasattr (reader, "back"), hasattr (reader, "interactions")]
      if not all (checks):
        raise ContinuumElectrostaticsError ("Output file %s empty or corrupted. Empty the scratch directory and start anew." % self.modelLog)
      self.Gborn_protein = reader.born
      self.Gback_protein = reader.back

      # Create a list of interactions
      interactions    = []
      instances       = []
      siteIndexOld    = 99999
      parentSiteIndex = self.parent.siteIndex

      for siteIndex, instanceIndex, energy in reader.interactions:
        if siteIndex > siteIndexOld:
          interactions.append (instances)
          instances = []
        # Set the interaction energy to zero if the site is interacting with itself
        if siteIndex == parentSiteIndex:
          energy = 0.
        instances.append (energy)
        siteIndexOld = siteIndex

      if instances:
        interactions.append (instances)

      # Copy the interactions to the centralized array
      indexGlobal = 0
      for site in interactions:
        for instance in site:
          energy = instance
          model.energyModel.SetInteraction (self._instIndexGlobal, indexGlobal, energy)
          indexGlobal = indexGlobal + 1


  #===============================================================================
  def CalculateGintr (self, log=logFile):
    """Calculate Gintr of an instance of a site in the protein."""
    # For checking, do not include the background energy of the model compound because it can sometimes be zero (for example in ligands)
    checks = [self.Gborn_protein, self.Gback_protein, self.Gborn_model]
    if all (checks):
      self.Gintr = self.Gmodel + (self.Gborn_protein - self.Gborn_model) + (self.Gback_protein - self.Gback_model)


  #===============================================================================
  def PrintInteractions (self, sort=False, log=logFile):
    """Print interactions of an instance of a site with other instances of other sites."""
    if LogFileActive (log):
      site         = self.parent
      model        = site.parent
      interactions = self.interactions

      if model.isCalculated:
        instances = []
        for site in model.meadSites:
          for instance in site.instances:
            wij = interactions[instance._instIndexGlobal]
            instances.append ([wij, site.segName, site.resName, site.resSerial, instance.label])
        if sort:
          instances.sort ()
        tab = log.GetTable (columns = [6, 6, 6, 6, 16])
        tab.Start ()
        tab.Heading ("Instance of a site", columnSpan = 4)
        tab.Heading ("Wij")

        for wij, segName, resName, resSerial, label in instances:
          entries = (
               ( "%s"     % segName   ),
               ( "%s"     % resName   ),
               ( "%d"     % resSerial ),
               ( "%s"     % label     ),
               ( "%16.4f" % wij       ),
                    )
          for entry in entries:
            tab.Entry (entry)
        tab.Stop ()


  #===============================================================================
  def _TableEntry (self, tab=None, secondsToCompletion=None):
    """Report calculated energies in a table.

    Optionally, include Estimated Time for Accomplishment (ETA).

    ETA has to be calculated outside of this method."""
    if tab:
      site = self.parent
      entries = (
            ( "%s"     % site.segName       ),
            ( "%s"     % site.resName       ),
            ( "%d"     % site.resSerial     ),
            ( "%s"     % self.label         ),
            ( "%16.4f" % self.Gborn_model   ),
            ( "%16.4f" % self.Gback_model   ),
            ( "%16.4f" % self.Gborn_protein ),
            ( "%16.4f" % self.Gback_protein ),
            ( "%16.4f" % self.Gmodel        ),
            ( "%16.4f" % self.Gintr         ),
                )
      for entry in entries:
        tab.Entry (entry)

      if isinstance (secondsToCompletion, float):
        minutes, seconds = divmod (secondsToCompletion, 60)
        hours, minutes   = divmod (minutes, 60)
        tab.Entry ("%16s" % ("%d:%02d:%02d" % (hours, minutes, seconds)))


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
