#-------------------------------------------------------------------------------
# . File      : Toolbox.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Different tools and helper functions."""

__lastchanged__ = "$Id$"

from Error          import ContinuumElectrostaticsError
from StateVector    import StateVector


def FormatEntry (items, header = False):
  """Generating headers or format strings for entries in files W.dat and gintr.dat"""
  if header:
    string = "#"
    for label, (width, digits) in items:
      string = string + (("%%%ds" % width) % label)
  else:
    string = ""
    for label, (width, digits) in items:
      if   digits < 0:
        string = string + "%%%ds" % width
      elif digits < 1:
        string = string + "%%%dd" % width
      else:
        string = string + "%%%d.%df" % (width, digits)
  return "%s\n" % string


def ConvertAttribute (attr):
  """Converting attributes to strings in Summary methods."""
  if isinstance (attr, bool):
    if attr:
      attrstring = "True"
    else:
      attrstring = "False"
  elif isinstance (attr, float):
    attrstring = "%g" % (attr)
  elif isinstance (attr, basestring):
    attrstring =  attr
  else:
    attrstring = str (attr)
  return attrstring


def StateVector_FromProbabilities (meadModel):
  """Create a state vector from the previously calculated probabilities."""
  if meadModel.isProbability:
    vector = StateVector (meadModel)

    for siteIndex, site in enumerate (meadModel.meadSites):
      pairs = []
      for instanceIndex, instance in enumerate (site.instances):
        pair = (instance.probability, instanceIndex)
        pairs.append (pair)
      maxProbPair = max (pairs)
      probability, instanceIndex = maxProbPair
      vector[siteIndex] = instanceIndex
    return vector
  else:
    raise ContinuumElectrostaticsError ("First calculate probabilities at a chosen pH.")


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
