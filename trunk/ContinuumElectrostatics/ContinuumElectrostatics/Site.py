#-------------------------------------------------------------------------------
# . File      : Site.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADSite is a class representing a titratable site."""

__lastchanged__ = "$Id$"

from pCore       import Vector3
from Error       import ContinuumElectrostaticsError
from Instance    import MEADInstance

import os


#-------------------------------------------------------------------------------
class MEADSite (object):
  """Titratable site.
  Each site has at least two instances (protonated and deprotonated)."""

  defaultAttributes = {
                  "parent"           : None , # <--This should point to the MEAD model
                  "siteIndex"        : None ,
                  "segName"          : None ,
                  "resName"          : None ,
                  "resSerial"        : None , # <--Keep it as integer
                  "instances"        : None ,
                  "center"           : None ,
                  "siteAtomIndices"  : None ,
                  "modelAtomIndices" : None ,
                      }

  @property
  def label (self):
    return "%s_%s%s" % (self.segName, self.resName, self.resSerial)


  #===============================================================================
  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


  #===============================================================================
  def CalculateCenterOfGeometry (self, system, centralAtom=None):
    """Calculate center of geometry of a site."""
    if centralAtom:
      found = False
      for index in self.siteAtomIndices:
        atom = system.atoms[index]
        if atom.label == centralAtom:
          found = True
          break
      if found:
        center = system.coordinates3[index]
      else:
        raise ContinuumElectrostaticsError ("Cannot find central atom %s in site %s %s %d" % (centralAtom, self.segName, self.resName, self.resSerial))

    else:
      center = Vector3 ()
      natoms = len (self.siteAtomIndices)
      for atomIndex in self.siteAtomIndices:
        center.AddScaledVector3 (1., system.coordinates3[atomIndex])
      center.Scale (1. / natoms)

    # Set the center of geometry
    self.center = center


  #===============================================================================
  def GetMostProbableInstance (self):
    """Return the index, label and probability of the most probable instance of a site."""
    model = self.parent
    if not model.isProbability:
      raise ContinuumElectrostaticsError ("First calculate probabilities.")

    mostProbValue = 0.
    for instance in self.instances:
      if instance.probability > mostProbValue:
        mostProbValue = instance.probability
        mostProbIndex = instance.instIndex
        mostProbLabel = instance.label
    return (mostProbValue, mostProbIndex, mostProbLabel)


  #===============================================================================
  def GetSortedIndices (self):
    """Get a list of indices of instances sorted by increasing probability."""
    model = self.parent
    if not model.isProbability:
      raise ContinuumElectrostaticsError ("First calculate probabilities.")

    instances = []
    for instance in self.instances:
      instances.append ([instance.probability, instance.instIndex])
    instances.sort ()
    indices = [index for probability, index in instances]
    return indices 


  #===============================================================================
  def _CreateFilename (self, prefix, label, postfix):
    model = self.parent
    if model.splitToDirectories:
        return os.path.join (model.scratch, self.segName, "%s%d" % (self.resName, self.resSerial), "%s_%s.%s" % (prefix, label, postfix))
    else:
        return os.path.join (model.scratch, "%s_%s_%s_%d_%s.%s" % (prefix, self.segName, self.resName, self.resSerial, label, postfix))


  #===============================================================================
  def _CreateInstances (self, templatesOfInstances, instIndexGlobal):
    """Create instances of a site."""
    protons    = []
    instances  = []
    meadModel  = self.parent

    for instIndex, instance in enumerate (templatesOfInstances):
      # Set the protons later because they come from the central array
      newInstance = MEADInstance (
          parent          = self                                                           ,
          instIndex       = instIndex                                                      ,
          instIndexGlobal = instIndexGlobal                                                ,
          label           = instance [ "label"   ]                                         ,
          charges         = instance [ "charges" ]                                         ,
          Gmodel          = instance [ "Gmodel"  ] * meadModel.temperature / 300.          ,
          modelPqr        = self._CreateFilename ("model", instance [ "label"   ], "pqr")  ,
          modelLog        = self._CreateFilename ("model", instance [ "label"   ], "out")  ,
          modelGrid       = self._CreateFilename ("model", instance [ "label"   ], "mgm")  ,
          sitePqr         = self._CreateFilename ("site",  instance [ "label"   ], "pqr")  ,
          siteLog         = self._CreateFilename ("site",  instance [ "label"   ], "out")  ,
          siteGrid        = self._CreateFilename ("site",  instance [ "label"   ], "ogm")  ,
                                 )
      instances.append (newInstance)
      instIndexGlobal = instIndexGlobal + 1

      # Remember the number of protons of the current instance
      nprotons  = instance [ "protons" ]
      protons.append (nprotons)

    self.instances = instances
    return (protons, instIndexGlobal)


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
