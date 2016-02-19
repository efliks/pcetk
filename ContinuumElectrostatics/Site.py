#-------------------------------------------------------------------------------
# . File      : Site.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore  import Vector3
from Error  import ContinuumElectrostaticsError

import os


class Site (object):
    """Base class for a titratable site.

    This class should not be used directly."""
    defaultAttributes = {
        }
    # parent  siteIndex  segName  resName  resSerial  instances  center  siteAtomIndices  modelAtomIndices 

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, val) in keywordArguments.iteritems ():
            setattr (self, key, val)


    @property
    def label (self):
        checks = (
            hasattr (self, "segName"  ),
            hasattr (self, "resName"  ),
            hasattr (self, "resSerial"), )
        if all (checks):
            return "%s_%s%d" % (self.segName, self.resName, self.resSerial)
        return ""

    @property
    def ninstances (self):
        if hasattr (self, "instances"):
            return len (self.instances)
        return 0

    @property
    def charge (self):
        """Get the current charge of a site."""
        probability, index, label = self.GetMostProbableInstance ()
        instance = self.instances[index]
        charge   = sum (instance.charges)
        return charge


    def _CalculateCenter (self, centralAtom=None):
        """Calculate center of geometry of a site."""
        ceModel = self.parent
        system  = ceModel.owner
        if centralAtom:
            found = False
            for index in self.siteAtomIndices:
                atom = system.atoms[index]
                if atom.label == centralAtom:
                    found = True
                    break
            if found:
                center = system.coordinates3[index]
            else:
                raise ContinuumElectrostaticsError ("Cannot find central atom %s in site %s %s %d" % (centralAtom, self.segName, self.resName, self.resSerial))
        else:
            center = Vector3 ()
            natoms = len (self.siteAtomIndices)
            for atomIndex in self.siteAtomIndices:
                center.AddScaledVector3 (1., system.coordinates3[atomIndex])
            center.Scale (1. / natoms)
        self.center = center


    def GetMostProbableInstance (self):
        """Return the index, label and probability of the most probable instance of a site."""
        model = self.parent
        if not model.isProbability:
            raise ContinuumElectrostaticsError ("First calculate probabilities.")

        mostProbValue = 0.
        for instance in self.instances:
            if instance.probability > mostProbValue:
                mostProbValue = instance.probability
                mostProbIndex = instance.instIndex
                mostProbLabel = instance.label
        return (mostProbValue, mostProbIndex, mostProbLabel)


    def GetSortedIndices (self):
        """Get a list of indices of instances sorted by increasing probability."""
        model     = self.parent
        if not model.isProbability:
            raise ContinuumElectrostaticsError ("First calculate probabilities.")
        instances = []
        for instance in self.instances:
            instances.append ([instance.probability, instance.instIndex])
        instances.sort ()

        indices = [index for probability, index in instances]
        return indices


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
