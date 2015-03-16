#-------------------------------------------------------------------------------
# . File      : Site.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADSite is a class representing a titratable site."""

__lastchanged__ = "$Id$"

from pCore            import Vector3, Selection, Clone
from Error            import ContinuumElectrostaticsError
from Instance         import MEADInstance
from PQRFileWriter    import PQRFile_FromSystem
from InputFileWriter  import WriteInputFile

import os


#-------------------------------------------------------------------------------
class MEADSite (object):
    """Titratable site.
    Each site has at least two instances (protonated and deprotonated)."""
    defaultAttributes = { "parent"           : None ,
                          "siteIndex"        : None ,
                          "segName"          : None ,
                          "resName"          : None ,
                          "resSerial"        : None ,
                          "instances"        : None ,
                          "center"           : None ,
                          "siteAtomIndices"  : None ,
                          "modelAtomIndices" : None }

    @property
    def label (self):
        return "%s_%s%s" % (self.segName, self.resName, self.resSerial)

    @property
    def ninstances (self):
        if self.instances:
            return len (self.instances)
        else:
            return 0

    @property
    def charge (self):
        """Get the current charge of a site."""
        probability, index, label = self.GetMostProbableInstance ()
        instance = self.instances[index]
        charge   = sum (instance.charges)
        return charge


    #===============================================================================
    def __init__ (self, *arguments, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


    #===============================================================================
    def CalculateCenterOfGeometry (self, system, centralAtom=None):
        """Calculate center of geometry of a site."""
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

        # Set the center of geometry
        self.center = center


    #===============================================================================
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


    #===============================================================================
    def GetSortedIndices (self):
        """Get a list of indices of instances sorted by increasing probability."""
        model = self.parent
        if not model.isProbability:
            raise ContinuumElectrostaticsError ("First calculate probabilities.")

        instances = []
        for instance in self.instances:
            instances.append ([instance.probability, instance.instIndex])
        instances.sort ()
        indices = [index for probability, index in instances]
        return indices


    #===============================================================================
    def _WriteMEADFiles (self, system, systemCharges, systemRadii):
        """For each instance of each site, write:
           - PQR file for the site    - PQR file for the model compound
           - OGM file for the site    - MGM file for the model compound"""
        grids = []
        model = self.parent
        for stepIndex, (nodes, resolution) in enumerate (model.focussingSteps):
            if stepIndex < 1:
                grids.append ("ON_GEOM_CENT %d %f\n" % (nodes, resolution))
            else:
                x, y, z = self.center
                grids.append ("(%f %f %f) %d %f\n"% (x, y, z, nodes, resolution))

        selectSite  = Selection (self.siteAtomIndices)
        selectModel = Selection (self.modelAtomIndices)


        # In the PQR file of the model compound, charges of the site atoms must be set to zero (requirement of the my_2diel_solver program)
        chargesZeroSite = Clone (systemCharges)
        for atomIndex in self.siteAtomIndices:
            chargesZeroSite[atomIndex] = 0.

        for instance in self.instances:
            PQRFile_FromSystem (instance.modelPqr, system, selection=selectModel, charges=chargesZeroSite, radii=systemRadii)

            # Update system charges with instance charges
            chargesInstance = Clone (systemCharges)
            for chargeIndex, atomIndex in enumerate (self.siteAtomIndices):
                chargesInstance[atomIndex] = instance.charges[chargeIndex]

            PQRFile_FromSystem (instance.sitePqr, system, selection=selectSite, charges=chargesInstance, radii=systemRadii)
            del chargesInstance

            # Write OGM and MGM files (they have the same content)
            for fileGrid in (instance.modelGrid, instance.siteGrid):
                WriteInputFile (fileGrid, grids)

        del chargesZeroSite


    #===============================================================================
    def _CreateFilename (self, prefix, label, postfix):
        model = self.parent
        if model.splitToDirectories:
            return os.path.join (model.pathScratch, self.segName, "%s%d" % (self.resName, self.resSerial), "%s_%s.%s" % (prefix, label, postfix))
        else:
            return os.path.join (model.pathScratch, "%s_%s_%s_%d_%s.%s" % (prefix, self.segName, self.resName, self.resSerial, label, postfix))


    #===============================================================================
    def _CreateInstances (self, templatesOfInstances, instIndexGlobal):
        """Create instances of a site."""
        instances  = []
        meadModel  = self.parent

        for instIndex, instance in enumerate (templatesOfInstances):
            newInstance = MEADInstance (
                parent           = self                                                           ,
                instIndex        = instIndex                                                      ,
                _instIndexGlobal = instIndexGlobal                                                ,
                label            = instance [ "label"   ]                                         ,
                charges          = instance [ "charges" ]                                         ,
                modelPqr         = self._CreateFilename ("model", instance [ "label"   ], "pqr")  ,
                modelLog         = self._CreateFilename ("model", instance [ "label"   ], "out")  ,
                modelGrid        = self._CreateFilename ("model", instance [ "label"   ], "mgm")  ,
                sitePqr          = self._CreateFilename ("site",  instance [ "label"   ], "pqr")  ,
                siteLog          = self._CreateFilename ("site",  instance [ "label"   ], "out")  ,
                siteGrid         = self._CreateFilename ("site",  instance [ "label"   ], "ogm")  ,
                                       )
            Gmodel   = instance [ "Gmodel"  ] * meadModel.temperature / 300.
            meadModel.energyModel.SetGmodel  (instIndexGlobal, Gmodel)
            nprotons = instance [ "protons" ]
            meadModel.energyModel.SetProtons (instIndexGlobal, nprotons)

            instances.append (newInstance)
            instIndexGlobal = instIndexGlobal + 1

        self.instances = instances
        return instIndexGlobal


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
