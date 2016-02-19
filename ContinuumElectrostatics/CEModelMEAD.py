#-------------------------------------------------------------------------------
# . File      : CEModelMEAD.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore            import logFile, LogFileActive, YAMLUnpickle, Selection

from Error            import ContinuumElectrostaticsError
from Constants        import YAMLPATHIN
from CEModel          import CEModel
from SiteMEAD         import SiteMEAD
from InstanceThread   import InstanceThread
from PQRFileWriter    import PQRFile_FromSystem
from InputFileWriter  import WriteInputFile

import os, time


_DEFAULT_THREADS        =  1
_DEFAULT_PATH_MEAD      =  os.path.join ("usr", "local", "bin")
_DEFAULT_PATH_SCRATCH   =  os.getenv ("PDYNAMO_SCRATCH")


class CEModelMEAD (CEModel):
    """A class to represent a continuum electrostatic model based on MEAD."""
    defaultAttributes = {
        "nthreads"             :   _DEFAULT_THREADS       ,
        "pathMEAD"             :   _DEFAULT_PATH_MEAD     ,
        "pathScratch"          :   _DEFAULT_PATH_SCRATCH  ,
        "deleteJobFiles"       :   False                  ,
        "splitToDirectories"   :   True                   ,
        }
    defaultAttributes.update (CEModel.defaultAttributes)


    defaultAttributeNames = {
        "Threads"              :  "nthreads"              ,
        "Split Directories"    :  "splitToDirectories"    ,
        "Delete Job Files"     :  "deleteJobFiles"        ,
        }
    defaultAttributeNames.update (CEModel.defaultAttributeNames)

    @property
    def label (self):
        return "MEAD"


    #-------------------------------------------------------------------------------
    def __init__ (self, system, customFiles=None, log=logFile, **keywordArguments):
        """Constructor."""
        super (CEModelMEAD, self).__init__ (system, customFiles=customFiles, log=log, **keywordArguments)

        # . Prepare filenames (do not write actual files)
        generate = (
                ("pathPqrProtein" ,  "protein.pqr"),
                ("pathPqrBack"    ,  "back.pqr"   ),
                ("pathFptSites"   ,  "site.fpt"   ), )
        for attribute, filename in generate:
            setattr (self, attribute, os.path.join (self.pathScratch, filename))


    #-------------------------------------------------------------------------------
    def _CreateSite (self, **keywordArguments):
        """Create a site and its instances specific to the MEAD-based CE model."""
        newSite = SiteMEAD (
            parent            =  self                                       ,
            siteIndex         =  keywordArguments  [ "siteIndex"        ]   ,
            segName           =  keywordArguments  [ "segName"          ]   ,
            resName           =  keywordArguments  [ "resName"          ]   ,
            resSerial         =  keywordArguments  [ "resSerial"        ]   ,
            siteAtomIndices   =  keywordArguments  [ "siteAtomIndices"  ]   ,
            modelAtomIndices  =  keywordArguments  [ "modelAtomIndices" ]   ,)
        # . Initialize instances
        libSite            = keywordArguments [ "libSite"         ]
        instIndexGlobal    = keywordArguments [ "instIndexGlobal" ]
        updatedIndexGlobal = newSite._CreateInstances (libSite.instances, instIndexGlobal)

        # . Calculate center of geometry
        newSite._CalculateCenter (centralAtom=libSite.center)

        # . Add the site to the list of sites
        self.sites.append (newSite)

        # . Finalize
        return updatedIndexGlobal


    #-------------------------------------------------------------------------------
    def CalculateElectrostaticEnergies (self, calculateETA=False, asymmetricTolerance=0.05, asymmetricSummary=False, log=logFile):
        if self.isFilesWritten:
            ninstances = self.ninstances
            times      = []
            tab        = None

            if LogFileActive (log):
                if self.nthreads < 2:
                    log.Text ("\nStarting serial run.\n")
                else:
                    log.Text ("\nStarting parallel run on %d CPUs.\n" % self.nthreads)

                heads = [("Instance of a site" , 4),
                         ("Gborn_model"        , 0),
                         ("Gback_model"        , 0),
                         ("Gborn_protein"      , 0),
                         ("Gback_protein"      , 0),
                         ("Gmodel"             , 0),
                         ("Gintr"              , 0),]
                columns = [6, 6, 6, 6, 16, 16, 16, 16, 16, 16]
                if calculateETA:
                    heads.append (("ETA", 0))
                    columns.append (16)
                tab = log.GetTable (columns = columns)
                tab.Start ()
                for head, span in heads:
                    if span > 0:
                        tab.Heading (head, columnSpan = span)
                    else:
                        tab.Heading (head)


            if self.nthreads < 2:
                for site in self.sites:
                    for instance in site.instances:
                        time0 = time.time ()
                        instance.CalculateModelCompound (log)
                        instance.CalculateProtein (log)
                        instance.CalculateGintr (log)

                        if calculateETA:
                            times.append (time.time () - time0)
                            averageTimePerInstance = sum (times) / len (times)
                            ninstances = ninstances - 1
                            instance._TableEntry (tab, secondsToCompletion = averageTimePerInstance * ninstances)
                        else:
                            instance._TableEntry (tab)
            else:
                batches = []
                threads = []
                limit   = self.nthreads - 1

                for site in self.sites:
                    for instance in site.instances:
                        if len (threads) > limit:
                            batches.append (threads)
                            threads = []
                        thread = InstanceThread (instance, log)
                        threads.append (thread)

                if threads:
                    batches.append (threads)

                for batch in batches:
                    for thread in batch: thread.start ()
                    for thread in batch: thread.join ()

                    secondsToCompletion = None
                    if calculateETA:
                        # . Collect times of execution
                        nthreads = len (batch)
                        for thread in batch:
                            times.append (thread.time)

                        averageTimePerInstance = sum (times) / len (times) / nthreads
                        ninstances = ninstances - nthreads
                        secondsToCompletion = averageTimePerInstance * ninstances

                    # . Print the results at the end of each batch, otherwise they come in random order
                    for thread in batch:
                        instance = thread.instance
                        instance._TableEntry (tab, secondsToCompletion = secondsToCompletion)
            if tab:
                tab.Stop ()
                log.Text ("\nCalculating electrostatic energies complete.\n")

            # . Check for symmetricity of the matrix of interactions
            self._CheckIfSymmetric (tolerance=asymmetricTolerance, printSummary=asymmetricSummary, log=log)

            # . Symmetrize interaction energies inside the matrix of interactions
            self.energyModel.SymmetrizeInteractions (log=log)

            # . Finalize
            self.isCalculated = True


    #-------------------------------------------------------------------------------
    def _CheckIfSymmetric (self, tolerance=0.05, printSummary=False, log=logFile):
        """This method is a wrapper for the EnergyModel's CheckIfSymmetric method.

        The wrapper is able to print summaries."""
        isSymmetric, maxDeviation = self.energyModel.CheckIfSymmetric (tolerance=tolerance)

        if LogFileActive (log):
            if isSymmetric:
                log.Text ("\nInteractions are symmetric within the given tolerance (%0.4f kcal/mol).\n" % tolerance)
            else:
                if not printSummary:
                    log.Text ("\nWARNING: Maximum deviation of interactions is %0.4f kcal/mol.\n" % maxDeviation)
                else:
                    heads = [("Instance of a site A" , 4),
                             ("Instance of a site B" , 4),
                             ("Deviation"            , 0),]
                    columns = (7, 7, 7, 7, 7, 7, 7, 7, 12)
                    gaps = ("%7s", "%7s", "%7d", "%7s")

                    tab = log.GetTable (columns=columns)
                    tab.Start ()
                    tab.Title ("Deviations of interactions")
                    for head, span in heads:
                        if span > 0:
                            tab.Heading (head, columnSpan=span)
                        else:
                            tab.Heading (head)

                    # . This fragment should be rewritten to work faster
                    report = []
                    for rowSite in self.sites:
                        for rowInstance in rowSite.instances:
                            for columnSite in self.sites:
                                for columnInstance in columnSite.instances:
                                    deviation = self.energyModel.GetDeviation (rowInstance._instIndexGlobal, columnInstance._instIndexGlobal)
                                    if abs (deviation) > tolerance:
                                        report.append ([rowInstance, columnInstance, deviation])

                    for ainstance, binstance, deviation in report:
                        asite = ainstance.parent
                        for gap, content in zip (gaps, (asite.segName, asite.resName, asite.resSerial, ainstance.label)):
                            tab.Entry (gap % content)
                        bsite = binstance.parent
                        for gap, content in zip (gaps, (bsite.segName, bsite.resName, bsite.resSerial, binstance.label)):
                            tab.Entry (gap % content)

                        tab.Entry ("%0.4f" % deviation)
                    tab.Stop ()
        return isSymmetric


    #-------------------------------------------------------------------------------
    def WriteJobFiles (self, log=logFile):
        """Write files: PQR, FPT, OGM and MGM."""
        if self.isInitialized:
            # . Get atomic charges and radii for the system
            system = self.owner

            systemCharges = system.AtomicCharges ()
            systemRadii   = []
            systemTypes   = system.energyModel.mmAtoms.AtomTypes ()
            radii         = YAMLUnpickle ("%s/%s" % (YAMLPATHIN, "radii.yaml"))

            for atomType in systemTypes:
                if radii.has_key (atomType):
                    radius = radii[atomType]
                else:
                    generalAtomType = "%s*" % atomType[0]
                    if radii.has_key (generalAtomType):
                        radius = radii[generalAtomType]
                    else:
                        raise ContinuumElectrostaticsError ("Cannot find atomic radius for atom type %s" % atomType)
                systemRadii.append (radius)

            # . Prepare scratch space
            if not os.path.exists (self.pathScratch):
                try:
                    os.makedirs (self.pathScratch)
                except:
                    raise ContinuumElectrostaticsError ("Cannot create scratch directory %s" % self.pathScratch)

            # . Create subdirectories, if necessary
            if self.splitToDirectories:
                for site in self.sites:
                    sitePqr   = site.instances[0].sitePqr
                    directory = os.path.dirname (sitePqr)
                    if not os.path.exists (directory):
                        try:
                            os.makedirs (directory)
                        except:
                            raise ContinuumElectrostaticsError ("Cannot create directory %s" % directory)

            # . Write PQR, OGM and MGM files of all instances of all sites
            for site in self.sites:
                site._WriteMEADFiles (system, systemCharges, systemRadii)

            # . Write background PQR file
            PQRFile_FromSystem (self.pathPqrBack, system, selection=Selection (self.backAtomIndices), charges=systemCharges, radii=systemRadii)

            # . Write full-protein PQR file (to be used as eps2set_region)
            PQRFile_FromSystem (self.pathPqrProtein, system, selection=Selection (self.proteinAtomIndices), charges=systemCharges, radii=systemRadii)

            # . Write FPT-file
            lines = []
            for siteIndex, site in enumerate (self.sites):
                for instanceIndex, instance in enumerate (site.instances):
                    for atomIndex, charge in zip (site.siteAtomIndices, instance.charges):
                        x, y, z = system.coordinates3[atomIndex]
                        line    = "%d %d %f %f %f %f\n" % (siteIndex, instanceIndex, x, y, z, charge)
                        lines.append (line)
            WriteInputFile (self.pathFptSites, lines)

            self.isFilesWritten = True


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
