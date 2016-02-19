#-------------------------------------------------------------------------------
# . File      : Model.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore         import logFile, LogFileActive
from CEModelMEAD   import CEModelMEAD

import os, glob


class MEADModel (CEModelMEAD):
    """A class to maintain backward compatibility."""

    def __init__ (self, system, log=logFile, **keywordArguments):
        """Constructor."""
        workDir     = os.getcwd ()
        customFiles = []
        for extension in ("*.yaml", "*.YAML", "*.est", "*.EST"):
            customFiles.extend (glob.glob (os.path.join (workDir, extension)))

        # . Call the constructor from the parent class
        super (MEADModel, self).__init__ (system, customFiles=customFiles, log=log, **keywordArguments)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
