#-------------------------------------------------------------------------------
# . File      : Instance.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADInstance is a class representing a particular instance of a site.

InstanceThread is a class used for running parallel calculations in MEAD."""

__lastchanged__ = "$Id$"


from pCore                  import logFile, LogFileActive

from MEADOutputFileReader   import MEADOutputFileReader
from Error                  import ContinuumElectrostaticsError

import os, threading, subprocess, time


#-------------------------------------------------------------------------------
class InstanceThread (threading.Thread):
  """A class for running parallel calculations for sites both in model compounds and in the protein.

  It also calculates Gintr."""

  def __init__ (self, instance, log = logFile):
    """Constructor."""
    threading.Thread.__init__ (self)
    self.instance = instance
    self.log      = log
    self.time     = 0


  def run (self):
    """The method that runs the calculations."""
    time0 = time.time ()
    self.instance.CalculateSiteInModelCompound (log = self.log)
    self.instance.CalculateSiteInProtein       (log = self.log)
    self.instance.CalculateGintr               (log = self.log)

    self.time = (time.time () - time0)


#-------------------------------------------------------------------------------
class MEADInstance (object):
  """Instance of a titratable site.

  For the time being, the instances only differ in charges. There are no rotameric instances."""

  defaultAttributes = {
                  "parent"          : None ,
                  "instIndex"       : None ,
                  "instIndexGlobal" : None ,
                  "label"           : None ,
                  "charges"         : None ,
                  "modelPqr"        : None ,
                  "modelLog"        : None ,
                  "modelGrid"       : None ,
                  "sitePqr"         : None ,
                  "siteLog"         : None ,
                  "siteGrid"        : None ,
                  "Gmodel"          : None ,
                  "Gborn_model"     : None ,
                  "Gback_model"     : None ,
                  "Gborn_protein"   : None ,
                  "Gback_protein"   : None ,
                  "probability"     : None ,
                      }

  @property
  def Gintr (self):
    model = self.parent.parent
    if model._intrinsic:
      return model._intrinsic[self.instIndexGlobal]
    else:
      return None

  @Gintr.setter
  def Gintr (self, value):
    model = self.parent.parent
    if model._intrinsic:
      model._intrinsic[self.instIndexGlobal] = value
    else:
      pass

  @property
  def protons (self):
    model = self.parent.parent
    if model._protons:
      return model._protons[self.instIndexGlobal]
    else:
      return None

  @protons.setter
  def protons (self, value):
    model = self.parent.parent
    if model._protons:
      model._protons[self.instIndexGlobal] = value
    else:
      pass

  @property
  def interactions (self):
    model = self.parent.parent
    if model._interactions:
      sites = []
      for site in model.meadSites:
        instances = []
        for instance in site.instances:
          instances.append (model._interactions[self.instIndexGlobal, instance.instIndexGlobal])
        sites.append (instances)
      return sites
    else:
      return None


  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


  def CalculateSiteInModelCompound (self, log = logFile):
    """Calculate Gborn and Gback of a site in the model compound."""
    site  = self.parent
    model = site.parent

    if model.isFilesWritten:
      if os.path.exists (self.modelLog):
        pass
      else:
        instancePqr        = self.sitePqr  [:-4]
        modelBackgroundPqr = self.modelPqr [:-4]
        command = [os.path.join (model.meadPath, "my_2diel_solver"), "-T", "%f" % model.temperature, "-ionicstr", "%f" % model.ionicStrength, "-epsin", "%f" % 4.0, instancePqr, modelBackgroundPqr]
  
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


  def CalculateSiteInProtein (self, log = logFile):
    """Calculate Gborn and Gback of a site in protein environment.

    Also, calculate Wij between the site and the other sites.

    Use fpt-file for other instances of sites."""
    site  = self.parent
    model = site.parent

    if model.isFilesWritten:
      if os.path.exists (self.siteLog):
        pass
      else:
        sitesFpt             = model.sitesFpt   [:-4]
        proteinPqr           = model.proteinPqr [:-4]
        proteinBackgroundPqr = model.backPqr    [:-4]
        instancePqr          = self.sitePqr     [:-4]
  
        # epsin1 is never used but must be given
  
        # eps2set defines the whole protein
  
        command = [os.path.join (model.meadPath, "my_3diel_solver"), "-T", "%f" % model.temperature, "-ionicstr", "%f" % model.ionicStrength, "-epsin1", "%f" % 1.0, "-epsin2", "%f" % 4.0, "-eps2set", "%s" % proteinPqr, "-fpt", "%s" % sitesFpt, instancePqr, proteinBackgroundPqr]
  
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
          model._interactions[self.instIndexGlobal, indexGlobal] = energy
          indexGlobal = indexGlobal + 1


  def CalculateGintr (self, log = logFile):
    """Calculate Gintr of an instance of a site in the protein."""
    # For checking, do not include the background energy of the model compound because it can sometimes be zero (for example in ligands)
    checks = [self.Gborn_protein, self.Gback_protein, self.Gborn_model]
    if all (checks):
      self.Gintr = self.Gmodel + (self.Gborn_protein - self.Gborn_model) + (self.Gback_protein - self.Gback_model)


  def PrintInteractions (self, sort = False, log = logFile):
    """Print interactions of an instance of a site with other instances of other sites."""
    if LogFileActive (log):
      site         = self.parent
      model        = site.parent
      interactions = self.interactions

      if model.isCalculated:
        instances = []

        for site in model.meadSites:
          for instance in site.instances:
            wij = interactions[site.siteIndex][instance.instIndex]
            instances.append ([wij, site.segName, site.resName, site.resSerial, instance.label])
        if sort: 
          instances.sort ()

        tab = log.GetTable (columns = [6, 6, 6, 6, 16])
        tab.Start ()
        tab.Heading ("Instance of a site", columnSpan = 4)
        tab.Heading ("Wij")

        for wij, segName, resName, resSerial, label in instances:
          tab.Entry (segName)
          tab.Entry (resName)
          tab.Entry ("%d" % resSerial)
          tab.Entry (label)
          tab.Entry ("%16.4f" % wij)
        tab.Stop ()


  def TableEntry (self, tab = None, secondsToCompletion = None):
    """Report calculated energies in a table.

    Optionally, include Estimated Time for Accomplishment (ETA).

    ETA has to be calculated outside of this method."""
    if tab:
      site = self.parent
      tab.Entry (site.segName)
      tab.Entry (site.resName)
      tab.Entry ("%d" % site.resSerial)
      tab.Entry (self.label)
      tab.Entry ("%16.4f" % self.Gborn_model)
      tab.Entry ("%16.4f" % self.Gback_model)
      tab.Entry ("%16.4f" % self.Gborn_protein)
      tab.Entry ("%16.4f" % self.Gback_protein)
      tab.Entry ("%16.4f" % self.Gmodel)
      tab.Entry ("%16.4f" % self.Gintr)

      if secondsToCompletion:
        minutes, seconds = divmod (secondsToCompletion, 60)
        hours, minutes   = divmod (minutes, 60)
        tab.Entry ("%16s" % ("%d:%02d:%02d" % (hours, minutes, seconds)))


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
