#-------------------------------------------------------------------------------
# . File      : InputFileWriter.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore  import TextFileWriter


class InputFileWriter (TextFileWriter):
    """A class for writing text files."""

    def __init__ (self, name):
        """Constructor."""
        super (InputFileWriter, self).__init__ (name)


    def Write (self, listOfLines):
        """Write a list of lines."""
        self.Open  ()
        self.file.writelines (listOfLines)
        self.Close ()


#===============================================================================
# . Helper functions
#===============================================================================
def WriteInputFile (filename, listOfLines, addLineBreaks=False):
    """Write a file."""
    output = InputFileWriter (filename)
    if addLineBreaks:
        templ = []
        for line in listOfLines:        
            templ.append ("%s\n" % line)
        listOfLines = templ
    output.Write (listOfLines)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
