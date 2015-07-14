#-------------------------------------------------------------------------------
# . File      : Model.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADModel is a class representing the continuum electrostatic model."""

__lastchanged__ = "$Id$"


import os, glob, time

from pCore           import logFile, LogFileActive, Selection, YAMLUnpickle
from Constants       import *
from Error           import ContinuumElectrostaticsError
from Site            import MEADSite
from Instance        import InstanceThread
from Utils           import FormatEntry, ConvertAttribute
from EnergyModel     import EnergyModel
from MCModelGMCT     import MCModelGMCT
from MCModelDefault  import MCModelDefault

# File handling
from ESTFileReader   import ESTFileReader
from InputFileWriter import WriteInputFile
from PQRFileWriter   import PQRFile_FromSystem


_DefaultTemperature     =  300.
_DefaultIonicStrength   =     .1
_DefaultPathScratch     =  os.getenv ("PDYNAMO_SCRATCH")
_DefaultPathMEAD        =  "/usr/local/bin"
_DefaultThreads         =  1
_DefaultFocussingSteps  =  ((121, 2.), (101, 1.), (101, .5), (101, .25))


class MEADModel (object):
    """Continuum electrostatic model."""

    defaultAttributes = {
        "nthreads"             :   _DefaultThreads          ,
        "temperature"          :   _DefaultTemperature      ,
        "ionicStrength"        :   _DefaultIonicStrength    ,
        "splitToDirectories"   :   True                     ,
        "deleteJobFiles"       :   False                    ,
        "isInitialized"        :   False                    ,
        "isFilesWritten"       :   False                    ,
        "isCalculated"         :   False                    ,
        "isProbability"        :   False                    ,
        "proteinAtomIndices"   :   None                     ,
        "backAtomIndices"      :   None                     ,
        "pathMEAD"             :   _DefaultPathMEAD         ,
        "pathScratch"          :   _DefaultPathScratch      ,
        "pathFptSites"         :   None                     ,
        "pathPqrProtein"       :   None                     ,
        "pathPqrBack"          :   None                     ,
        "focussingSteps"       :   _DefaultFocussingSteps   ,
        "librarySites"         :   None                     ,
        "meadSites"            :   None                     ,
        "energyModel"          :   None                     ,
        "sampler"              :   None                     ,
        "owner"                :   None                     ,
            }

    defaultAttributeNames = {
        "Temperature"          :  "temperature"             ,
        "Ionic Strength"       :  "ionicStrength"           ,
        "Threads"              :  "nthreads"                ,
        "Split Directories"    :  "splitToDirectories"      ,
        "Delete Job Files"     :  "deleteJobFiles"          ,
        "Initialized"          :  "isInitialized"           ,
        "Files Written"        :  "isFilesWritten"          ,
        "Calculated"           :  "isCalculated"            ,
        "Calculated Prob."     :  "isProbability"           ,
            }

    @property
    def ninstances (self):
        if self.meadSites:
            return sum (map (lambda site: site.ninstances, self.meadSites))
        else:
            return 0

    @property
    def nsites (self):
        if self.meadSites:
            return len (self.meadSites)
        else:
            return 0


    #===============================================================================
    def __del__ (self):
        """Deallocation."""
        if self.deleteJobFiles: self.DeleteJobFiles ()


    #===============================================================================
    def __init__ (self, system, log=logFile, *arguments, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in keywordArguments.iteritems (): setattr (self, key, value)

        self.owner       = system
        self.pathMEAD    = os.path.abspath (self.pathMEAD)
        self.pathScratch = os.path.abspath (self.pathScratch)
        self.LoadLibraryOfSites (log=log)


    #===============================================================================
    def LoadLibraryOfSites (self, log=logFile):
        """Load a set of YAML or EST files with parameters for titratable sites.
        If there are YAML or EST files in the current directory, they are loaded as well.
        If these additional files have names coinciding with the names from the library, the library parameters will be overwritten.
        Notice! EST files have priority over YAML files."""
        directory    = os.getcwd ()
        filesLibrary = glob.glob (os.path.join (YAMLPATHIN, "sites", "*.yaml"))
        filesExtra   = glob.glob (os.path.join (directory, "*.yaml"))
        filesEST     = glob.glob (os.path.join (directory, "*.est"))

        if LogFileActive (log):
            for fileExtra in (filesExtra + filesEST):
                log.Text ("\nIncluded custom file: %s\n" % os.path.basename (fileExtra))

        self.librarySites = {}
        for fileSite in (filesLibrary + filesExtra):
            site      = YAMLUnpickle (fileSite)
            name      = site [ "site"      ]
            atoms     = site [ "atoms"     ]
            instances = site [ "instances" ]
            self.librarySites[name] = {"atoms" : atoms, "instances" : instances, "center" : None}

        for fileSite in filesEST:
            reader    = ESTFileReader (fileSite)
            reader.Parse ()
            name      = reader.siteLabel
            atoms     = reader.siteAtoms
            instances = reader.siteInstances
            center    = reader.siteCenter
            self.librarySites[name] = {"atoms" : atoms, "instances" : instances, "center" : center}


    #===============================================================================
    def WriteW (self, filename="W.dat", precision=3, log=logFile):
        """Write an interaction matrix compatible with GMCT."""
        if precision < 3 or precision > 8:
            raise ContinuumElectrostaticsError ("Wrong value for precision (%d)." % precision)

        if self.isCalculated:
            spacing = precision * 2
            if spacing < 10:
                spacing = 10

            items = (
                ( "idSite1"  , (   8   ,      0    ) ),
                ( "idInst1"  , (   8   ,      0    ) ),
                ( "labSite1" , (  14   ,     -1    ) ),
                ( "labInst1" , (  14   ,     -1    ) ),
                ( "idSite2"  , (   8   ,      0    ) ),
                ( "idInst2"  , (   8   ,      0    ) ),
                ( "labSite2" , (  14   ,     -1    ) ),
                ( "labInst2" , (  14   ,     -1    ) ),
                ( "Wij_symm" , (spacing,  precision) ),
                ( "Wij"      , (spacing,  precision) ),
                ( "Wij_err"  , (spacing,  precision) ),
                    )
            header = FormatEntry (items, header = True)
            entry  = FormatEntry (items)
            lines  = [header]

            for asite in self.meadSites:
                for ainstance in asite.instances:
                    for bsite in self.meadSites:
                        for binstance in bsite.instances:
                            interSymmetric = self.energyModel.GetInteractionSymmetric (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            interaction    = self.energyModel.GetInteraction          (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            deviation      = self.energyModel.GetDeviation            (ainstance._instIndexGlobal, binstance._instIndexGlobal)

                            lines.append (entry % (asite.siteIndex + 1, ainstance.instIndex + 1, asite.label, ainstance.label, bsite.siteIndex + 1, binstance.instIndex + 1, bsite.label, binstance.label, interSymmetric, interaction, deviation))
            WriteInputFile (filename, lines)


    #===============================================================================
    def WriteGintr (self, filename="gintr.dat", precision=3, log=logFile):
        """Iterate over instances and write a gintr.dat file compatible with GMCT.

        This file contains Gintr of each instance of each site."""
        if precision < 3 or precision > 8:
            raise ContinuumElectrostaticsError ("Wrong value for precision (%d)." % precision)

        if self.isCalculated:
            spacing = precision * 2
            if spacing < 10:
                spacing = 10
            items = (( "siteID"    , (  12   ,      0    ) ),
                     ( "instID"    , (  12   ,      0    ) ),
                     ( "siteLabel" , (  16   ,     -1    ) ),
                     ( "instLabel" , (  16   ,     -1    ) ),
                     ( "Gintr"     , (spacing,  precision) ),
                     ( "protons"   , (  12   ,      0    ) ),)
            header = FormatEntry (items, header = True)
            entry  = FormatEntry (items)
            lines  = [header]

            for site in self.meadSites:
                for instance in site.instances:
                    lines.append (entry % (site.siteIndex + 1, instance.instIndex + 1, site.label, instance.label, instance.Gintr, instance.protons))
            WriteInputFile (filename, lines)


    #===============================================================================
    def CalculateMicrostateEnergy (self, stateVector, pH=7.0):
        """Wrapper function to calculate microstate energy."""
        return self.energyModel.CalculateMicrostateEnergy (stateVector, pH=pH)


    #===============================================================================
    def DefineMCModel (self, sampler, log=logFile):
        """Define Monte Carlo model."""
        if isinstance (sampler, MCModelDefault) or isinstance (sampler, MCModelGMCT):
            self.sampler = sampler
            self.sampler.Initialize (self)
            self.sampler.PrintPairs (log=log)
            self.sampler.Summary    (log=log)
        elif sampler is None:
            self.sampler = None
        else:
            raise ContinuumElectrostaticsError ("Cannot define MC model.")


    #===============================================================================
    def CalculateProbabilities (self, pH=7.0, unfolded=False, isCalculateCurves=False, logFrequency=-1, log=logFile):
        """Calculate probabilities."""
        nstates = -1
        sites   = None
        if       self.sampler and     unfolded:
            raise ContinuumElectrostaticsError ("Monte Carlo sampling of unfolded proteins unsupported.")
        elif     self.sampler and not unfolded:
            self.sampler.CalculateOwnerProbabilities (pH=pH, logFrequency=logFrequency, log=log)
        elif not self.sampler and     unfolded:
            nstates = self.energyModel.CalculateProbabilitiesAnalyticallyUnfolded (pH=pH)
        elif not self.sampler and not unfolded:
            nstates = self.energyModel.CalculateProbabilitiesAnalytically (pH=pH)

        if isCalculateCurves:
            sites = []
            for site in self.meadSites:
                instances = []
                for instance in site.instances:
                    instances.append (instance.probability)
                sites.append (instances)

        if nstates > 0:
            if LogFileActive (log):
                log.Text ("\nCalculated %d protonation states.\n" % nstates)

        self.isProbability = True
        return sites


    #===============================================================================
    def CalculateElectrostaticEnergies (self, calculateETA=False, asymmetricTolerance=0.05, asymmetricSummary=False, log=logFile):
        """
        Calculate for each instance of each site:
        - self (Born) energy in the model compound
        - interaction energy   between the site and the background charge set of the model compound

        - self (Born) energy in the protein
        - interaction energy   between the site and the background charge set of the protein
        - interaction energies between the site and the other sites in their different protonation forms

        Finally, use the calculated heterotransfer energies to calculate Gintr from Gmodel.
        """
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
                for meadSite in self.meadSites:
                    for instance in meadSite.instances:
                        time0 = time.time ()
                        instance.CalculateSiteInModelCompound (log)
                        instance.CalculateSiteInProtein (log)
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

                for meadSite in self.meadSites:
                    for instance in meadSite.instances:
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
                        # Collect times of execution
                        nthreads = len (batch)
                        for thread in batch:
                            times.append (thread.time)

                        averageTimePerInstance = sum (times) / len (times) / nthreads
                        ninstances = ninstances - nthreads
                        secondsToCompletion = averageTimePerInstance * ninstances

                    # Print the results at the end of each batch, otherwise they come in random order
                    for thread in batch:
                        instance = thread.instance
                        instance._TableEntry (tab, secondsToCompletion = secondsToCompletion)
            if tab:
                tab.Stop ()
                log.Text ("\nCalculating electrostatic energies complete.\n")


            # Check for symmetricity of the matrix of interactions
            self._CheckIfSymmetric (tolerance=asymmetricTolerance, printSummary=asymmetricSummary, log=log)

            # Symmetrize interaction energies inside the matrix of interactions
            self.energyModel.SymmetrizeInteractions (log=log)

            self.isCalculated = True


    #===============================================================================
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

                    # This fragment should be rewritten to work faster
                    report = []
                    for rowSite in self.meadSites:
                        for rowInstance in rowSite.instances:
                            for columnSite in self.meadSites:
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


    #===============================================================================
    def _CheckResidue (self, excludeResidues, segmentName, residueName, residueSerial, log=logFile):
        """Check if the residue should be included."""
        includeResidue = True
        if excludeResidues:
            for exclSegmentName, exclResidueName, exclResidueSerial in excludeResidues:
                if   (    exclSegmentName) and (    exclResidueName) and (    exclResidueSerial):
                    if exclSegmentName == segmentName and exclResidueName == residueName and exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (    exclSegmentName) and (    exclResidueName) and (not exclResidueSerial):
                    if exclSegmentName == segmentName and exclResidueName == residueName:
                        includeResidue = False
                        break

                elif (    exclSegmentName) and (not exclResidueName) and (    exclResidueSerial):
                    if exclSegmentName == segmentName and exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (    exclSegmentName) and (not exclResidueName) and (not exclResidueSerial):
                    if exclSegmentName == segmentName:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (    exclResidueName) and (    exclResidueSerial):
                    if exclResidueName == residueName and exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (    exclResidueName) and (not exclResidueSerial):
                    if exclResidueName == residueName:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (not exclResidueName) and (    exclResidueSerial):
                    if exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (not exclResidueName) and (not exclResidueSerial):
                    includeResidue = False
                    break

            if not includeResidue:
                if LogFileActive (log):
                    log.Text ("\nExcluding residue: %s %s %d\n" % (segmentName, residueName, residueSerial))
        return includeResidue


    #===============================================================================
    def Initialize (self, excludeSegments=None, excludeResidues=None, includeTermini=False, log=logFile):
        """Decompose the system into model compounds, sites and a background charge set.

        |excludeSegments| is a sequence of segment names to exclude from the model, usually segments of water molecules.
        |excludeResidues| is a sequence of three-element sequences (segmentName, residueName, residueSerial).

        It is possible to leave some of the elements blank, for example ("PRTA", "CYS", "") means exclude all cysteines in segment PRTA.
        """
        if not self.isInitialized:
            system = self.owner

            # Check for the CHARMM energy model
            if system.energyModel.mmModel.label is not "CHARMM":
                raise ContinuumElectrostaticsError ("The energy model of the system is different from CHARMM.")

            # Perform a dry run to calculate the numbers of sites and instances
            totalSites, totalInstances = self._SplitModel (excludeSegments=excludeSegments, excludeResidues=excludeResidues, includeTermini=includeTermini, dryRun=True, log=log)
    
            # Allocate arrays of Gmodels, protons, intrinsic energies, interaction energies and probabilities
            self.energyModel = EnergyModel (self, totalSites, totalInstances)
    
            # Perform the actual initialization
            self._SplitModel (excludeSegments=excludeSegments, excludeResidues=excludeResidues, includeTermini=includeTermini, dryRun=False, log=None)
    
            # Complete the initialization of the energy model
            self.energyModel.Initialize ()
   
            # Construct the background set of charges and the protein (to be used as eps2set_region)
            self._SetupBackground ()
 
            # Prepare the filenames (does not write the actual files)
            self.pathPqrProtein     = os.path.join (self.pathScratch, "protein.pqr")
            self.pathPqrBack        = os.path.join (self.pathScratch, "back.pqr")
            self.pathFptSites       = os.path.join (self.pathScratch, "site.fpt")
    
            # Finish up
            self.isInitialized = True


    #===============================================================================
    def _GetResidueInfo (self, residue):
        system        = self.owner
        ParseLabel    = system.sequence.ParseLabel
        residueName, residueSerial = ParseLabel (residue.label, fields=2)
        residueSerial = int (residueSerial)

        segment       = residue.parent
        segmentName   = segment.label
        return (segmentName, residueName, residueSerial)


    #===============================================================================
    def _GetIndices (self, residue, atomLabels, check=True):
        missingLabels = []
        atomIndices   = []
        atoms         = residue.children
        for label in atomLabels:
            index = -1
            for atom in atoms:
                if label == atom.label:
                    index = atom.index
                    break
            if index >= 0:
                atomIndices.append (index)
            else:
                missingLabels.append (label)

            if check and missingLabels:
                segmentName, residueName, residueSerial = self._GetResidueInfo (residue)
                raise ContinuumElectrostaticsError ("Cannot include residue %s %s %d because of missing atoms: %s" % (segmentName, residueName, residueSerial, " ".join (missingLabels)))
        return atomIndices


    #===============================================================================
    def _SetupBackground (self):
        allSiteAtomIndices  = []
        for site in self.meadSites:
            allSiteAtomIndices.extend (site.siteAtomIndices)
        system              = self.owner
        segments            = system.sequence.children
        backAtomIndices     = []
        proteinAtomIndices  = []

        #============ Iterate segments ============
        for segment in segments:
            residues = segment.children
        
            #============ Iterate residues ============
            for residue in residues:
                foo, residueName, residueSerial = self._GetResidueInfo (residue)
        
                # Remove residues not defined in PROTEIN_RESIDUES, usually waters and ions
                if residueName not in REMOVE_RESIDUES:
                    atoms = residue.children
        
                    #============ Iterate atoms ============
                    for atom in atoms:
                        proteinAtomIndices.append (atom.index)
                        if atom.index not in allSiteAtomIndices:
                            backAtomIndices.append (atom.index)
        self.proteinAtomIndices = proteinAtomIndices
        self.backAtomIndices    = backAtomIndices


    #===============================================================================
    def _SetupSites (self, residue, prevResidue=None, nextResidue=None, terminal=None, log=logFile):
        segmentName, residueName, residueSerial = self._GetResidueInfo (residue)
        setupSites    = []

        if residueName in PROTEIN_RESIDUES:
            if terminal:
                if   terminal == "N":
                    if   residueName == "GLY":
                        libTerm = self.librarySites["NGL"]
                    elif residueName == "PRO":
                        libTerm = self.librarySites["NPR"]
                    else:
                        libTerm = self.librarySites["NTR"]
                elif terminal == "C":
                    libTerm = self.librarySites["CTR"]
                termAtomIndices = self._GetIndices (residue, libTerm["atoms"])

        # Check if titratable residue
        if residueName in self.librarySites:
            libSite      = self.librarySites[residueName]
            siteIndices  = self._GetIndices (residue, libSite["atoms"])
            modelIndices = []
            for atom in residue.children:
                modelIndices.append (atom.index)

            prevIndices  = []
            if prevResidue:
                foo, prevResidueName, prevResidueSerial = self._GetResidueInfo (prevResidue)
                if prevResidueName in PROTEIN_RESIDUES:
                    prevIndices = self._GetIndices (prevResidue, PREV_RESIDUE, check=False) 

            nextIndices  = []
            if nextResidue:
                foo, nextResidueName, nextResidueSerial = self._GetResidueInfo (nextResidue)
                if nextResidueName in PROTEIN_RESIDUES:
                    if   nextResidueName == "GLY":
                        nextLabels = NEXT_RESIDUE_GLY
                    elif nextResidueName == "PRO":
                        nextLabels = NEXT_RESIDUE_PRO
                    else:
                        nextLabels = NEXT_RESIDUE
                    nextIndices = self._GetIndices (nextResidue, nextLabels, check=False)

            setupSites.append ([libSite, siteIndices, prevIndices + modelIndices + nextIndices])
        return setupSites


    #===============================================================================
    def _SplitModel (self, excludeSegments=None, excludeResidues=None, includeTermini=False, dryRun=True, log=logFile):
        totalSites       =  0
        totalInstances   =  0
        siteIndex        =  0
        instIndexGlobal  =  0
        if not dryRun:
            self.meadSites = []

        if excludeSegments is None:
            excludeSegments = ["WATA", "WATB", "WATC", "WATD", ]
        system     = self.owner
        segments   = system.sequence.children


        #============ Iterate segments ============
        for segment in segments:
            segmentName = segment.label

            # Include segment?
            if segmentName not in excludeSegments:
                residues  = segment.children
                nresidues = len (residues)

                #============ Iterate residues ============
                for residueIndex, residue in enumerate (residues):
                    foo, residueName, residueSerial = self._GetResidueInfo (residue)
                    includeResidue = self._CheckResidue (excludeResidues, segmentName, residueName, residueSerial, log=log)
                    if not includeResidue:
                        continue
                    if residueIndex < 1:
                        prevResidue = None
                        terminal    = "N"
                    else:
                        prevResidue = residues[residueIndex - 1]
                        terminal    = None
                    if residueIndex > (nresidues - 2):
                        nextResidue = None
                        terminal    = "C"
                    else:
                        nextResidue = residues[residueIndex + 1]
                        terminal    = None

                    setupSites = self._SetupSites (residue, prevResidue, nextResidue, terminal=(terminal if includeTermini else None), log=log)
                    for libSite, siteAtomIndices, modelAtomIndices in setupSites:
                        if not dryRun:
                            newSite = MEADSite (
                                parent           = self              ,
                                siteIndex        = siteIndex         ,
                                segName          = segmentName       ,
                                resName          = residueName       ,
                                resSerial        = residueSerial     ,
                                siteAtomIndices  = siteAtomIndices   ,
                                modelAtomIndices = modelAtomIndices  ,
                                               )
                            # Calculate center of geometry
                            newSite.CalculateCenterOfGeometry (system, libSite["center"])

                            # Add instances to the newly created site
                            instIndexGlobal = newSite._CreateInstances (libSite["instances"], instIndexGlobal)

                            # Add the newly created site to the list of sites
                            self.meadSites.append (newSite)
                        siteIndex      = siteIndex      + 1
                        totalSites     = totalSites     + 1
                        totalInstances = totalInstances + len (libSite["instances"])

        # Return the calculated numbers of sites and instances (useful for dryRun)
        return (totalSites, totalInstances)


    #===============================================================================
    def Summary (self, log=logFile):
        """Summary."""
        if LogFileActive (log):
            summary = log.GetSummary ()
            summary.Start ("Continuum Electrostatic Model MEAD")

            keys = self.__class__.defaultAttributeNames.keys ()
            keys.sort ()
            for key in keys:
                attr     = getattr (self, self.__class__.defaultAttributeNames[key])
                attrConv = ConvertAttribute (attr)
                summary.Entry (key, attrConv)

            summary.Entry ("Number Of Sites",     "%s" % self.nsites)
            summary.Entry ("Number Of Instances", "%s" % self.ninstances)
            summary.Stop ()


    #===============================================================================
    def SummarySites (self, log=logFile):
        """List titratable residues."""
        if LogFileActive (log):
            if self.isInitialized:
                tab = log.GetTable (columns = [8, 8, 8, 8, 10, 10, 10, 10])
                tab.Start ()
                tab.Heading ("SiteID"),
                tab.Heading ("Site", columnSpan = 3),
                tab.Heading ("Instances"),
                tab.Heading ("Center", columnSpan = 3),

                for site in self.meadSites:
                    entries = (( "%d"      %  (site.siteIndex + 1) ),
                               ( "%s"      %  site.segName         ),
                               ( "%s"      %  site.resName         ),
                               ( "%d"      %  site.resSerial       ),
                               ( "%d"      %  len (site.instances) ),
                               ( "%10.3f"  %  site.center[0]       ),
                               ( "%10.3f"  %  site.center[1]       ),
                               ( "%10.3f"  %  site.center[2]       ),)
                    for entry in entries:
                        tab.Entry (entry)
                tab.Stop ()


    #===============================================================================
    def SummaryProbabilities (self, reportOnlyUnusual=False, maxProbThreshold=0.75, log=logFile):
        """List probabilities of occurance of instances."""
        unusualProtonations = {
            "HIS" : ("HSE", "HSD", "fd"),
            "ARG" : ("d", ),
            "ASP" : ("p", ),
            "CYS" : ("d", ),
            "GLU" : ("p", ),
            "LYS" : ("d", ),
            "TYR" : ("d", ), }

        if LogFileActive (log):
            if self.isProbability:
                maxinstances = 0
                for site in self.meadSites:
                    ninstances = len (site.instances)
                    if ninstances > maxinstances: maxinstances = ninstances

                tab = log.GetTable (columns=[6, 6, 6] + [8, 8] * maxinstances)
                tab.Start ()
                tab.Heading ("Site", columnSpan=3)
                tab.Heading ("Probabilities of instances", columnSpan=maxinstances * 2)

                for site in self.meadSites:
                    maxProb = 0.
                    for instance in site.instances:
                        if instance.probability > maxProb:
                            maxLabel = instance.label
                            maxIndex = instance.instIndex
                            maxProb  = instance.probability

                    skipSite = False
                    if reportOnlyUnusual:
                        if site.resName in unusualProtonations:
                            labels = unusualProtonations[site.resName]
                            if maxLabel not in labels:
                                skipSite = True

                    # If the maximum probability is lower than a certain threshold, do not skip the site in the report
                    if maxProb < maxProbThreshold:
                        skipSite = False

                    if not skipSite:
                        tab.Entry ("%6s" % site.segName)
                        tab.Entry ("%6s" % site.resName)
                        tab.Entry ("%6d" % site.resSerial)

                        for instance in site.instances:
                            if instance.instIndex == maxIndex:
                                label = "*%s" % instance.label
                            else:
                                label = instance.label
                            tab.Entry ("%8s"   % label)
                            tab.Entry ("%8.4f" % instance.probability)

                        for filler in range (0, maxinstances - len (site.instances)):
                            tab.Entry ("")
                            tab.Entry ("")
                tab.Stop ()


    #===============================================================================
    def SedScript_FromProbabilities (self, filename="his_repl.sed", overwrite=False, putPath=False, log=logFile):
        """Generate a sed-script for substituting histidines in the source PDB file based on the calculated probabilities."""
        if not self.isProbability:
            raise ContinuumElectrostaticsError ("First calculate probabilities.")

        if not overwrite:
            if os.path.exists (filename):
                raise ContinuumElectrostaticsError ("File %s already exists." % filename)

        lines = []
        if putPath:
            lines.append ("# %s\n" % os.path.abspath (filename))

        # First take care of histidines
        for site in self.meadSites:
            if site.resName in ("HIS", "HSP"):
                mostProbValue, mostProbIndex, mostProbLabel = site.GetMostProbableInstance ()
                lines.append ("/H.. .%4d/  s/H../%3s/  # %.4f\n" % (site.resSerial, mostProbLabel, mostProbValue))

        # Then everything else
        unusualProtonations = {
            "ARG" : ("d", ),
            "ASP" : ("p", ),
            "CYS" : ("d", ),
            "GLU" : ("p", ),
            "LYS" : ("d", ),
            "TYR" : ("d", ), }
        translateLabels = {
            "p" :   "protonated",
            "d" : "deprotonated", }
        warnings = []

        for site in self.meadSites:
            if site.resName not in ("HIS", "HSP"):
                mostProbValue, mostProbIndex, mostProbLabel = site.GetMostProbableInstance ()

                if site.resName in unusualProtonations:
                    unusualLabels = unusualProtonations[site.resName]
                    if mostProbLabel in unusualLabels:
                        if   site.resName == "ASP":
                            lines.append ("# patch ASPP %4s %4d setup  ! %.4f\n" % (site.segName, site.resSerial, mostProbValue))
                        elif site.resName == "GLU":
                            lines.append ("# patch GLUP %4s %4d setup  ! %.4f\n" % (site.segName, site.resSerial, mostProbValue))
                        else:
                            warnings.append ("# Warning: %4s %3s %4d is %s with probability of %.4f\n" % (site.segName, site.resName, site.resSerial, translateLabels[mostProbLabel], mostProbValue))
                else:
                    warnings.append ("# Unknown residue: %4s %3s %4d (most probable instance is \"%s\" with probability of %.4f)\n" % (site.segName, site.resName, site.resSerial, mostProbLabel, mostProbValue))
        lines.extend (warnings)

        # Write the sed script to disk
        WriteInputFile (filename, lines)

        if LogFileActive (log):
            log.Text ("\nWrote file: %s\n" % filename)


    #===============================================================================
    def DeleteJobFiles (self):
        """Delete job files."""
        files = []
        files.append (self.pathPqrProtein)
        files.append (self.pathPqrBack)
        files.append (self.pathFptSites)

        for site in self.meadSites:
            for instance in site.instances:
                files.append (instance.sitePqr)
                files.append (instance.siteLog)
                files.append (instance.siteGrid)
                files.append (instance.modelPqr)
                files.append (instance.modelLog)
                files.append (instance.modelGrid)

        for f in files:
            if os.path.exists (f): os.remove (f)

        if self.splitToDirectories:
            # Remove directories
            pass


    #===============================================================================
    def WriteJobFiles (self, log=logFile):
        """Write files: PQR, FPT, OGM and MGM."""
        if self.isInitialized:
            # Get atomic charges and radii for the system
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

            # Prepare scratch space
            if not os.path.exists (self.pathScratch):
                try:
                    os.mkdir (self.pathScratch)
                except:
                    raise ContinuumElectrostaticsError ("Cannot create scratch directory %s" % self.pathScratch)

            # Create subdirectories, if necessary
            if self.splitToDirectories:
                for meadSite in self.meadSites:
                    sitePqr   = meadSite.instances[0].sitePqr
                    directory = os.path.dirname (sitePqr)
                    if not os.path.exists (directory):
                        try:
                            os.makedirs (directory)
                        except:
                            raise ContinuumElectrostaticsError ("Cannot create directory %s" % directory)

            # Write PQR, OGM and MGM files of all instances of all sites
            for meadSite in self.meadSites:
                meadSite._WriteMEADFiles (system, systemCharges, systemRadii)

            # Write background PQR file
            PQRFile_FromSystem (self.pathPqrBack, system, selection=Selection (self.backAtomIndices), charges=systemCharges, radii=systemRadii)

            # Write full-protein PQR file (to be used as eps2set_region)
            PQRFile_FromSystem (self.pathPqrProtein, system, selection=Selection (self.proteinAtomIndices), charges=systemCharges, radii=systemRadii)

            # Write FPT-file
            lines = []
            for siteIndex, meadSite in enumerate (self.meadSites):
                for instanceIndex, instance in enumerate (meadSite.instances):
                    for atomIndex, charge in zip (meadSite.siteAtomIndices, instance.charges):
                        x, y, z = system.coordinates3[atomIndex]
                        line    = "%d %d %f %f %f %f\n" % (siteIndex, instanceIndex, x, y, z, charge)
                        lines.append (line)
            WriteInputFile (self.pathFptSites, lines)

            self.isFilesWritten = True


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
