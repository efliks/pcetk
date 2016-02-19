#-------------------------------------------------------------------------------
# . File      : CEModelDefault.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore     import logFile, LogFileActive
from CEModel   import CEModel


class CEModelDefault (CEModel):
    """A class to represent a built-in continuum electrostatic model."""
    defaultAttributes = {
        }
    defaultAttributes.update (CEModel.defaultAttributes)

    defaultAttributeNames = {
        }
    defaultAttributeNames.update (CEModel.defaultAttributeNames)


    def __init__ (self, system, customFiles=None, log=logFile, **keywordArguments):
        """Constructor."""
        super (CEModelDefault, self).__init__ (system, customFiles=customFiles, log=log, **keywordArguments)

    @property
    def label (self):
        return "Default"


    #-------------------------------------------------------------------------------
    def _CreateSite (self, **keywordArguments):
        """Create a site and its instances specific to the built-in CE model."""
        pass


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
