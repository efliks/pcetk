#-------------------------------------------------------------------------------
# . File      : ESTFileReader.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Reading EST parameter files."""

__lastchanged__ = "$Id$"

from pCore import logFile, LogFileActive, TextFileReader

import os


class ESTFileReader (TextFileReader):
  """A class for reading EST files."""

  def __init__ (self, name):
    """Constructor."""
    self.siteLabel = os.path.splitext (os.path.basename (name))[0]
    TextFileReader.__init__ (self, name)


  def Parse (self, log = logFile):
    """Parse the data on the file."""
    if not self.QPARSED:
      if LogFileActive (log): self.log = log

      self.Open ()
      atoms  = []
      line   = None
      # In some files center atom is absent?
      center = None

      try:
        while True:
          line   = self.GetLine (QWARNING = False)
          tokens = line.split ()
          if len (tokens) > 0:
            if tokens[0] == "Gmodel" : Gmodels = map (float, tokens[1:])
            if tokens[0] == "proton" : protons = map (int, tokens[1:])
            if tokens[0] == "label"  : labels  = tokens[1:]
            if tokens[0] == "center" : center  = tokens[1]
  
            if tokens[0] == self.siteLabel:
              label   = tokens[1]
              charges = map (float, tokens[2:])
              atoms.append ((label, charges))
      except EOFError:
        pass
      self.WarningStop ()
      self.Close ()
      self.log      = None
      self.QPARSED  = True

      instances = []
      for instanceIndex, (instanceLabel, instanceGmodel, instanceProtons) in enumerate (zip (labels, Gmodels, protons)):
        instanceAtoms   = []
        instanceCharges = []
        for atomLabel, atomCharges in atoms:
          instanceAtoms.append (atomLabel)
          instanceCharges.append (atomCharges[instanceIndex])

        instances.append ({"label" : instanceLabel, "Gmodel" : instanceGmodel, "protons" : instanceProtons, "charges" : instanceCharges})
      self.siteAtoms     = instanceAtoms
      self.siteInstances = instances
      self.siteCenter    = center


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
