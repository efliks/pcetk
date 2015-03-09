#-------------------------------------------------------------------------------
# . File      : GMCTOutputFileReader.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Reading output files from GMCT."""

__lastchanged__ = "$Id$"

from pCore      import logFile, LogFileActive, TextFileReader
from Constants  import *


class GMCTOutputFileReader (TextFileReader):
    """A class for reading output files from GMCT."""

    def __init__ (self, name):
        """Constructor."""
        TextFileReader.__init__ (self, name)


    def Parse (self, temperature = 300.0, log = logFile):
        """Parse the data on the file."""
        if not self.QPARSED:
            if LogFileActive (log):
                self.log = log

            self.Open ()
            convert       = -1.0 / (CONSTANT_MOLAR_GAS_KCAL_MOL * CONSTANT_LN10 * temperature)
            line          = None
            pHtable       = []
            probabilities = {}

            try:
                while True:
                    if not line:
                        line = self.GetLine (QWARNING = False)

                    if line.startswith ("chemical potential"):
                        tokens = line.split ()
                        mu     = float (tokens[2])
                        pHtable.append (mu * convert)

                        self.GetLine ()
                        self.GetLine ()

                        while True:
                            line = self.GetLine (QWARNING = False)
                            if line.startswith ("chemical potential"):
                                break

                            tokens = line.split ()
                            label, probability, mu, protons, vlabel, vmemb = tokens

                            if probabilities.has_key (label):
                                entries = probabilities[label]
                            else:
                                entries = []
                            entries.append (float (probability))
                            probabilities[label] = entries
            except EOFError:
                pass

            self.pHtable       = pHtable
            self.probabilities = probabilities
            self.WarningStop ()
            self.Close ()
            self.log     = None
            self.QPARSED = True


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
