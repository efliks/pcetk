#-------------------------------------------------------------------------------
# . File      : Site.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADSite is a class representing a titratable site."""

__lastchanged__ = "$Id$"


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
                  "modelAtomIndices" : None ,
                  "siteAtomIndices"  : None ,
                      }

  @property
  def label (self): return "%s_%s%s" % (self.segName, self.resName, self.resSerial)


  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


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
# Testing
#===============================================================================
if __name__ == "__main__": pass
