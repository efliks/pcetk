#-------------------------------------------------------------------------------
# . File      : Instance.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADInstance is a class representing a particular instance of a site.

InstanceThread is a class used for running parallel calculations in MEAD."""


from pCore                  import logFile, LogFileActive

from MEADOutputFileReader   import MEADOutputFileReader
from Error                  import ContinuumElectrostaticsError

import  os, threading, subprocess


#-------------------------------------------------------------------------------
class InstanceThread (threading.Thread):
  """A class for running parallel calculations for sites both in model compounds and in the protein.

  It also calculates Gintr."""

  def __init__ (self, instance = None, log = logFile):
    """Constructor."""
    threading.Thread.__init__ (self)
    self.instance = instance
    self.log      = log


  def run (self):
    """The method that runs the calculations."""
    self.instance.CalculateSiteInModelCompound (log = self.log)
    self.instance.CalculateSiteInProtein       (log = self.log)
    self.instance.CalculateGintr               (log = self.log)


#-------------------------------------------------------------------------------
class MEADInstance (object):
  """Instance of a titratable site.

  For the time being, the instances only differ in charges. There are no rotameric instances."""

  defaultAttributes = {
                  "parent"          : None , # <--This should point to the instance's site
                  "instID"          : None ,
                  "instIndexGlobal" : None ,
                  "label"           : None ,
                  "protons"         : None ,
                  "charges"         : None ,
                  "interactions"    : None ,
                  "modelPqr"        : None ,
                  "modelLog"        : None ,
                  "modelGrid"       : None ,
                  "sitePqr"         : None ,
                  "siteLog"         : None ,
                  "siteGrid"        : None ,
                  "Gmodel"          : None ,
                  "Gintr"           : None ,
                  "Gborn_model"     : None ,
                  "Gback_model"     : None ,
                  "Gborn_protein"   : None ,
                  "Gback_protein"   : None ,
                      }

  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


  def CalculateSiteInModelCompound (self, log = logFile):
    """Calculate Gborn and Gback of a site in the model compound."""
    if os.path.exists (self.modelLog):
      pass
    else:
      instancePqr        = self.sitePqr  [:-4]
      modelBackgroundPqr = self.modelPqr [:-4]

      model   = self.parent.parent
      command = [os.path.join (model.meadPath, "my_2diel_solver"), "-T", "%f" % model.temperature, "-ionicstr", "%f" % model.ionicStrength, "-epsin", "%f" % 4.0, instancePqr, modelBackgroundPqr]

      try:
        outFile = open (self.modelLog, "w")
        subprocess.check_call (command, stderr = outFile, stdout = outFile)
        outFile.close ()
      except:
        raise ContinuumElectrostaticsError ("Failed running command: %s" % " ".join (command))

    reader = MEADOutputFileReader (self.modelLog)
    reader.Parse ()

    if not hasattr (reader, "born") or not hasattr (reader, "back"):
      raise ContinuumElectrostaticsError ("Output file %s empty or corrupted. Empty the scratch directory and start anew." % self.modelLog)

    self.Gborn_model = reader.born
    self.Gback_model = reader.back


  def CalculateSiteInProtein (self, log = logFile):
    """Calculate Gborn and Gback of a site in protein environment.

    Also, calculate Wij between the site and the other sites.

    Use fpt-file for other instances of sites."""
    if os.path.exists (self.siteLog):
      pass
    else:
      model = self.parent.parent
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

    if not hasattr (reader, "born") or not hasattr (reader, "back") or not hasattr (reader, "interactions"):
      raise ContinuumElectrostaticsError ("Output file %s empty or corrupted. Empty the scratch directory and start anew." % self.modelLog)

    self.Gborn_protein = reader.born
    self.Gback_protein = reader.back

    sites         = []
    instances     = []
    prevSiteIndex = -1

    for siteIndex, instanceIndex, energy in reader.interactions:
      if siteIndex == (self.parent.siteID - 1):
        Wij = 0.0
      else:
        Wij = energy

      if siteIndex is not prevSiteIndex and siteIndex > 0:
        sites.append (instances)
        instances = []
      instances.append (Wij)
      prevSiteIndex = siteIndex

    if instances:
      sites.append (instances)

# Return it to the calling method? The matrix of interactions should be inside the MEADModel?
    self.interactions = sites


  def CalculateGintr (self, log = logFile):
    """Calculate Gintr of an instance of a site in the protein."""
    self.Gintr = self.Gmodel + (self.Gborn_protein - self.Gborn_model) + (self.Gback_protein - self.Gback_model)


  def PrintInteractions (self, sort = False, log = logFile):
    """Print interactions of the instance of a site with other instances of other sites."""
    if LogFileActive (log):
      site  = self.parent
      model = site.parent

      if model.isCalculated:
        instances = []

        for site in model.meadSites:
          for instance in site.instances:
            s   = site.siteID     - 1
            i   = instance.instID - 1
            Wij = self.interactions[s][i]
            instances.append ([Wij, site.segName, site.resName, site.resNum, instance.label])
        if sort: 
          instances.sort ()

        table = log.GetTable (columns = [6, 6, 6, 6, 16])
        table.Start ()
        table.Heading ("Instance of a site", columnSpan = 4)
        table.Heading ("Wij")

        for Wij, segName, resName, resNum, label in instances:
          table.Entry (segName)
          table.Entry (resName)
          table.Entry (resNum)
          table.Entry (label)
          table.Entry ("%16.4f" % Wij)
        table.Stop ()


  def TableEntry (self, table = None):
    """Report calculated energies in a table."""
    if table:
      site = self.parent
      table.Entry (site.segName)
      table.Entry (site.resName)
      table.Entry (site.resNum)
      table.Entry (self.label)
      table.Entry ("%16.4f" % self.Gborn_model)
      table.Entry ("%16.4f" % self.Gback_model)
      table.Entry ("%16.4f" % self.Gborn_protein)
      table.Entry ("%16.4f" % self.Gback_protein)
      table.Entry ("%16.4f" % self.Gmodel)
      table.Entry ("%16.4f" % self.Gintr)


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
