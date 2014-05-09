#-------------------------------------------------------------------------------
# . File      : PQRFileWriter.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""
A module that handles writing PQR files. These files are needed as an input to extended MEAD.

This module is based on the PDBFileWriter module.
"""

__revision__ = "$Revision$"

# import itertools, time
# from pMolecule    import PDYNAMO_VERSION, PeriodicTable, 
# from ExportImport import _Exporter


from pCore        import Coordinates3, TextFileWriter
from pMolecule    import Sequence, System

_ATOMLINEFORMAT = "%-6s%5i %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%10.5f%8.3f%7d%7d%7d%7d\n"


#-------------------------------------------------------------------------------
class PQRFileWriter (TextFileWriter):
  """PQRFileWriter is the class for writing PQR files."""

  def WriteSystem (self, system, data = None, selection = None, charges = None, radii = None):
    """Write out a system.

    |system|    is the system to be written.
    |data|      is the data to be written. The coordinates of system are written if |data| is absent.
    |selection| gives the atoms to write. The default is to write all atoms.
    |charges|   is the table containing atomic charges. The charges are set to zero if |charges| is absent.
    |radii|     is the table containing atomic radii. The radii are set to zero if |radii| is absent.
    """

    # Code taken from PDBFileWriter
    if not isinstance (system, System): 
      raise TypeError ("Invalid |system| argument.")

    sequence = getattr (system, "sequence", None)
    if sequence is None: 
      sequence = Sequence.FromAtomContainer (system.atoms, componentLabel = "UNK.1")

    if isinstance (data, Coordinates3): 
      xyz = data
    else:                                 
      xyz = system.coordinates3

    if xyz is None: 
      raise TypeError ("Unable to obtain coordinate data from |system| or |data| arguments.")


    natoms = len (system.atoms)

    if natoms != xyz.rows: 
      raise TypeError ("The PQR model and coordinate data are of different lengths.")

    if selection is None: 
      towrite = xrange (natoms)
    else:                 
      towrite = selection

    if charges is None:
      charges = [0.] * natoms
    else:
      if natoms != len (charges):
        raise TypeError ("The PQR model and charge data are of different lengths.")

    if radii is None:
      radii = [0.] * natoms
    else:
      if natoms != len (radii):
        raise TypeError ("The PQR model and radii data are of different lengths.")
   
 
    self.Open ()

    for natm, iatom in enumerate (towrite):
      atom   = system.atoms [iatom]
      charge = charges      [iatom]
      radius = radii        [iatom]
      
      resName, resSeq, iCode = sequence.ParseLabel (atom.parent.label, fields = 3)
      segID     = ""
      chainID   = ""
      atomIndex = 0

#       element = PeriodicTable.Symbol (atom.atomicNumber).upper ()
#       if element == "*": 
#         element = ""
#       entityLabel = atom.parent.parent.label
#       if useSegmentEntityLabels:
#         chainID = ""
#         segID   = entityLabel[0:4]
#       else:
#         chainID = entityLabel[0:1]
#         segID   = ""

      label = atom.label
      if   len (label) >= 4: 
        outputlabel = label[0:4]
      elif label[0:1].isdigit (): 
        outputlabel = label
      else:                        
        outputlabel = " " + label

      x = xyz[iatom, 0]
      y = xyz[iatom, 1] 
      z = xyz[iatom, 2]

      self.file.write (_ATOMLINEFORMAT % ("ATOM", atomIndex, outputlabel, " ", resName[0:3], chainID, resSeq[:4], iCode, x, y, z, charge, radius, 0, 0, 0, 0))

    self.Close ()


#-------------------------------------------------------------------------------
def PQRFile_FromSystem (filename, system, selection = None, charges = None, radii = None):
  """Helper function that writes a system to a PQR file."""
  outfile = PQRFileWriter (filename)
  outfile.WriteSystem (system, selection = selection, charges = charges, radii = radii)

#  _Exporter.AddHandler ( { System : PQRFile_FromSystem } , [ "PQR" ], "Extended MEAD" )


#===============================================================================
if __name__ == "__main__" : pass
