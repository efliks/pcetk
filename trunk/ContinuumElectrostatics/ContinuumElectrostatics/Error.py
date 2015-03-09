#-------------------------------------------------------------------------------
# . File      : Error.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Error handling in the ContinuumElectrostatics module."""

__lastchanged__ = "$Id$"

import exceptions


class ContinuumElectrostaticsError (exceptions.StandardError):
    """A class for handling errors in the ContinuumElectrostatics module."""
    pass


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
