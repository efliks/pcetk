#-------------------------------------------------------------------------------
# . File      : Site.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADSite is a class representing a titratable site."""

__lastchanged__ = "$Id :$"


class MEADSite (object):
  """Titratable site.

  Each site has at least two instances (protonated and deprotonated)."""

  defaultAttributes = {
                  "parent"           : None , # <--This should point to the MEAD model
                  "siteID"           : None ,
                  "segName"          : None ,
                  "resName"          : None ,
                  "resSerial"        : None ,
                  "instances"        : None ,
                  "center"           : None ,
                  "modelAtomIndices" : None ,
                  "siteAtomIndices"  : None ,
                  "label"            : None ,
                      }

  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)

    self.label = "%s_%s%s" % (self.segName, self.resName, self.resSerial)


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
