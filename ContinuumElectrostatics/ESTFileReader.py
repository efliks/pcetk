#-------------------------------------------------------------------------------
# . File      : ESTFileReader.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Reading EST parameter files."""

from pCore      import logFile, LogFileActive, TextFileReader
from Constants  import *


class ESTFileReader (TextFileReader):
  """A class for reading EST files."""

  def __init__ (self, name):
    """Constructor."""
    self.siteLabel = "XXX"
    TextFileReader.__init__ (self, name)


  def Parse (self, temperature = 300.0, log = logFile):
    """Parse the data on the file."""
    if not self.QPARSED:
      if LogFileActive (log): self.log = log

      self.Open ()

      line    = None
      convert = -1.0 / (CONSTANT_MOLAR_GAS_KCAL_MOL * CONSTANT_LN10 * temperature)

      try:
        while True:
          if not line:
            line = self.GetLine (QWARNING = False)

          if line.startswith ("label"):
            pass
          if line.startswith ("Gmodel"):
            pass
          if line.startswith ("proton"):
            pass
          if line.startswith ("center"):
            pass
          if line.startswith (self.siteLabel):
            pass
      except EOFError:
        pass

      self.WarningStop ()
      self.Close ()
      self.log     = None
      self.QPARSED = True


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
