#-------------------------------------------------------------------------------
# . File      : InputFileWriter.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""A simple writer of text files from lists of strings."""

__lastchanged__ = "$Id$"

from pCore      import TextFileWriter


class InputFileWriter (TextFileWriter):
  """A simple class for writing files from lists of strings."""

  def __init__ (self, name):
    """Constructor."""
    TextFileWriter.__init__ (self, name)


  def Write (self, listOfStrings):
    """Write a list of strings."""
    self.Open ()

    self.file.writelines (listOfStrings)
    self.Close ()


#===============================================================================
# Helper functions
#===============================================================================
def WriteInputFile (filename, listOfLines, addLineBreaks=False):
  """Writing input files from previously generated lists of strings."""
  out = InputFileWriter (filename)
  if addLineBreaks:
    out.Write (map (lambda line: "%s\n" % line, listOfLines))
  else:
    out.Write (listOfLines)


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
