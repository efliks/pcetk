#-------------------------------------------------------------------------------
# . File      : CEModel.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore             import logFile, LogFileActive, Selection, YAMLUnpickle
from pMolecule         import System

from Error             import ContinuumElectrostaticsError
from Constants         import TERM_REMOVE, PROTEIN_RESIDUES, NEXT_RESIDUE_GLY, NEXT_RESIDUE_PRO, NEXT_RESIDUE, PREV_RESIDUE, REMOVE_RESIDUES
from EnergyModel       import EnergyModel
from InputFileWriter   import WriteInputFile
from TemplatesLibrary  import TemplatesLibrary
from MCModelGMCT       import MCModelGMCT
from MCModelDefault    import MCModelDefault

import os


_DEFAULT_FOCUSING_STEPS    =  ((121, 2.),  (101, 1.),  (101, .5),  (101, .25))
_DEFAULT_EXCLUDE_SEGMENTS  =  ("WATA", "WATB", "WATC", "WATD", )

_DEFAULT_TEMPERATURE       =  300.      # 300 K
_DEFAULT_IONIC_STRENGTH    =     .1     # 100 mM = 0.1 M
_DEFAULT_EPSILON_WATER     =   80.
_DEFAULT_EPSILON_PROTEIN   =    4.

# CEAtom = collections.namedtuple ("CEAtom", "label  x  y  z  charge  radii")


class CEModel (object):
    """Base class for continuum electrostatic models.

    This class should not be used directly."""
    defaultAttributes = {
        "temperature"      :   _DEFAULT_TEMPERATURE      ,
        "ionicStrength"    :   _DEFAULT_IONIC_STRENGTH   ,
        "focusingSteps"    :   _DEFAULT_FOCUSING_STEPS   ,
        "epsilonWater"     :   _DEFAULT_EPSILON_WATER    ,
        "epsilonProtein"   :   _DEFAULT_EPSILON_PROTEIN  ,
        }

    defaultAttributeNames = {
        "Temperature"           :   "temperature"     ,
        "Ionic Strength"        :   "ionicStrength"   ,
        "Initialized"           :   "isInitialized"   ,
        "Files Written"         :   "isFilesWritten"  ,
        "Calculated"            :   "isCalculated"    ,
        "Calculated Prob."      :   "isProbability"   ,
        "Focusing Steps"        :   "focusingSteps"   ,
        "Water   Diel. Const."  :   "epsilonWater"    ,
        "Protein Diel. Const."  :   "epsilonProtein"  ,
        }

    @property
    def nsites (self):
        if hasattr (self, "sites"):
            return len (self.sites)
        return 0

    @property
    def ninstances (self):
        if hasattr (self, "sites"):
            counter = 0
            for site in self.sites:
                counter += site.ninstances
            return counter
        return 0

    @property
    def label (self):
        return "Base model"


    #-------------------------------------------------------------------------------
    def __init__ (self, system, customFiles=None, log=logFile, **keywordArguments):
        """Constructor."""
        # . Perform initial checks
        if not isinstance (system, System):
            raise ContinuumElectrostaticsError ("Cannot assign owner.")

        if system.energyModel.mmModel.label != "CHARMM":
            raise ContinuumElectrostaticsError ("Cannot use a non-CHARMM energy model.")


        # . Set owner and library
        self.owner   = system
        self.library = TemplatesLibrary (customFiles=customFiles, log=log)

        # . Set attributes
        attributes = self.__class__.defaultAttributes
        for (key, val) in attributes.iteritems ():
            setattr (self, key, val)

        for (key, val) in keywordArguments.iteritems ():
            if not attributes.has_key (key):
                raise ContinuumElectrostaticsError ("Unknown attribute: %s" % key)
            setattr (self, key, val)

        # . Set status attributes
        for attribute in ("isInitialized", "isFilesWritten", "isCalculated", "isProbability"):
            setattr (self, attribute, False)


    #-------------------------------------------------------------------------------
    def __del__ (self):
        """Deallocation."""
        pass


    #-------------------------------------------------------------------------------
    def Initialize (self, excludeSegments=_DEFAULT_EXCLUDE_SEGMENTS, excludeResidues=None, includeTermini=False, log=logFile):
        """Decompose the system into model compounds, sites and a background charge set.

        |excludeSegments| is a sequence of segment names to exclude from the model, usually segments of water molecules.
        |excludeResidues| is a sequence of three-element sequences (segmentName, residueName, residueSerial).

        It is possible to leave some of the elements blank, for example ("PRTA", "CYS", "") means exclude all cysteines in segment PRTA.
        """
        if not self.isInitialized:
            # . Perform a dry run to calculate the numbers of sites and instances
            totalSites, totalInstances = self._SplitModel (excludeSegments=excludeSegments, excludeResidues=excludeResidues, includeTermini=includeTermini, dryRun=True, log=log)

            # . Allocate arrays of Gmodels, protons, intrinsic energies, interaction energies and probabilities
            self.energyModel = EnergyModel (self, totalSites, totalInstances)

            # . Perform the actual initialization
            self._SplitModel (excludeSegments=excludeSegments, excludeResidues=excludeResidues, includeTermini=includeTermini, dryRun=False, log=None)
    
            # . Complete the initialization of the energy model
            self.energyModel.Initialize ()
   
            # . Construct the background set of charges and the protein (to be used as eps2set_region)
            self._SetupBackground ()
 
            # . Finish up
            self.isInitialized = True


    #-------------------------------------------------------------------------------
    def CalculateElectrostaticEnergies (self, **keywordArguments):
        """
        Calculate for each instance of each site:
        - self (Born) energy in the model compound
        - interaction energy   between the site and the background charge set of the model compound

        - self (Born) energy in the protein
        - interaction energy   between the site and the background charge set of the protein
        - interaction energies between the site and the other sites in their different protonation forms

        Finally, use the calculated heterotransfer energies to calculate Gintr from Gmodel.
        """
        pass


    #-------------------------------------------------------------------------------
    def Summary (self, log=logFile):
        """Summary."""
        if LogFileActive (log):
            summary = log.GetSummary ()
            # . Write header
            if hasattr (self, "label"):
                label = self.label
            else:
                label = "Unknown"
            summary.Start ("Continuum Electrostatic Model (%s)" % label)
            # . Write attributes
            names = self.__class__.defaultAttributeNames.keys ()
            names.sort ()
            for name in names:
                attribute = getattr (self, self.__class__.defaultAttributeNames[name])
                if   isinstance (attribute, bool):
                    convert = "True" if attribute else "False"
                elif isinstance (attribute, float):
                    convert = "%g" % attribute
                elif isinstance (attribute, basestring):
                    convert = attribute
                elif isinstance (attribute, tuple):
                    convert = "%d" % len (attribute)
                else:
                    convert = str (attribute)
                summary.Entry (name, convert)

            summary.Entry ("Number Of Sites",     "%s" % self.nsites)
            summary.Entry ("Number Of Instances", "%s" % self.ninstances)
            summary.Stop ()


    #-------------------------------------------------------------------------------
    def SummarySites (self, log=logFile):
        """List titratable residues."""
        if LogFileActive (log):
            if self.isInitialized:
                tab = log.GetTable (columns=[8, 8, 8, 8, 10, 10, 10, 10])
                # . Write header
                tab.Start ()
                tab.Heading ("SiteID"),
                tab.Heading ("Site", columnSpan=3),
                tab.Heading ("Instances"),
                tab.Heading ("Center", columnSpan=3),
                # . Write sites
                for site in self.sites:
                    entries = (
                        ( "%d"      %  (site.siteIndex + 1) ),
                        ( "%s"      %  site.segName         ),
                        ( "%s"      %  site.resName         ),
                        ( "%d"      %  site.resSerial       ),
                        ( "%d"      %  site.ninstances      ),
                        ( "%10.3f"  %  site.center[0]       ),
                        ( "%10.3f"  %  site.center[1]       ),
                        ( "%10.3f"  %  site.center[2]       ), )
                    for entry in entries:
                        tab.Entry (entry)
                tab.Stop ()


    #-------------------------------------------------------------------------------
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
                for site in self.sites:
                    ninstances = len (site.instances)
                    if ninstances > maxinstances: maxinstances = ninstances

                tab = log.GetTable (columns=[6, 6, 6] + [8, 8] * maxinstances)
                tab.Start ()
                tab.Heading ("Site", columnSpan=3)
                tab.Heading ("Probabilities of instances", columnSpan=maxinstances * 2)
                for site in self.sites:
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

                    # . If the maximum probability is lower than a certain threshold, do not skip the site in the report
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
                        for filler in range (maxinstances - site.ninstances):
                            tab.Entry ("")
                            tab.Entry ("")
                tab.Stop ()


    #-------------------------------------------------------------------------------
    def SedScript_FromProbabilities (self, filename="his_repl.sed", overwrite=False, putPath=False, log=logFile):
        """Generate a sed script for substituting histidines in the source PDB file based on the calculated probabilities."""
        if not self.isProbability:
            raise ContinuumElectrostaticsError ("First calculate probabilities.")

        if not overwrite:
            if os.path.exists (filename):
                if LogFileActive (log):
                    log.Text ("\nFile %s already exists, skipping.\n" % filename)
                    return
        lines = []
        if putPath:
            lines.append ("# %s\n" % os.path.abspath (filename))

        # . First write histidines
        for site in self.sites:
            if site.resName in ("HIS", "HSP"):
                mostProbValue, mostProbIndex, mostProbLabel = site.GetMostProbableInstance ()
                lines.append ("/H.. .%4d/  s/H../%3s/  # %.4f\n" % (site.resSerial, mostProbLabel, mostProbValue))
        # . Then everything else
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

        for site in self.sites:
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

        # . Write a sed script
        WriteInputFile (filename, lines)
        # . Summarize
        if LogFileActive (log):
            log.Text ("\nWrote file: %s\n" % filename)


    #-------------------------------------------------------------------------------
    def PrintInteractions (self, log=logFile):
        """Print the matrix of interactions in a readable form.

        It only remains readable as long as there are a few titratable sites with two instances each."""
        if self.isCalculated:
            if LogFileActive (log):
                header = "%16s" % ""
                for asite in self.sites:
                    header = "%s%16s" % (header, asite.label.center (8))
                logFile.Text (header + "\n")
                header = "%16s" % ""
                for asite in self.sites:
                    for ainstance in asite.instances:
                        header = "%s%8s" % (header, ainstance.label.center (8))
                logFile.Text (header + "\n")
                # . Outer loop
                for asite in self.sites:
                    for ainstance in asite.instances:
                        line = "%16s" % (asite.label + " " + ainstance.label)
                        # . Inner loop
                        for bsite in self.sites:
                            for binstance in bsite.instances:
                                Wij  = self.energyModel.GetInteractionSymmetric (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                                line = "%s%8.2f" % (line, Wij)
                        logFile.Text (line + "\n")


    #-------------------------------------------------------------------------------
    def WriteW (self, filename="W.dat", precision=3, log=logFile):
        """Write a GMCT-compatible matrix of interactions."""
        if self.isCalculated:
            items = (
                ( "idSiteA"  ,   8  ,   0 ),
                ( "idInstA"  ,   8  ,   0 ),
                ( "labSiteA" ,  14  ,  -1 ),
                ( "labInstA" ,  10  ,  -1 ),
                ( "idSiteB"  ,   8  ,   0 ),
                ( "idInstB"  ,   8  ,   0 ),
                ( "labSiteB" ,  14  ,  -1 ),
                ( "labInstB" ,  10  ,  -1 ),
                ( "Wij_symm" ,  10  ,   4 ),
                ( "Wij"      ,  10  ,   4 ),
                ( "Wij_err"  ,  10  ,   4 ), )
            # . Prepare header
            header = "# "
            for label, width, digits in items:
                header = "%s%s" % (header, label.center (width))
            lines  = ["%s\n" % header, ]
            # . Prepare formating string
            form = ""
            for label, width, digits in items:
                if   digits < 0:
                    form = "%s%%%ds"    % (form, width)
                elif digits < 1:
                    form = "%s%%%dd"    % (form, width)
                else:
                    form = "%s%%%d.%df" % (form, width, digits)
            form = "%s\n" % form
            # . Outer loop
            for asite in self.sites:
                for ainstance in asite.instances:
                    # . Inner loop
                    for bsite in self.sites:
                        for binstance in bsite.instances:
                            interSymmetric = self.energyModel.GetInteractionSymmetric (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            interaction    = self.energyModel.GetInteraction          (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            deviation      = self.energyModel.GetDeviation            (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            lines.append (form % (asite.siteIndex + 1, ainstance.instIndex + 1, asite.label, ainstance.label, bsite.siteIndex + 1, binstance.instIndex + 1, bsite.label, binstance.label, interSymmetric, interaction, deviation))
            # . Write to a file
            WriteInputFile (filename, lines)


    #-------------------------------------------------------------------------------
    def WriteGintr (self, filename="gintr.dat", precision=3, log=logFile):
        """Write a GMCT-compatible file containing intrinsic energies of each instance of each site."""
        if self.isCalculated:
            items = (
                ( "siteID"    ,   8  ,   0 ),
                ( "instID"    ,   8  ,   0 ),
                ( "siteLabel" ,  14  ,  -1 ),
                ( "instLabel" ,  10  ,  -1 ),
                ( "Gintr"     ,  12  ,   4 ),
                ( "protons"   ,   8  ,   0 ), )
            # . Prepare header
            header = "# "
            for label, width, digits in items:
                header = "%s%s" % (header, label.center (width))
            lines  = ["%s\n" % header, ]
            # . Prepare formating string
            form = ""
            for label, width, digits in items:
                if   digits < 0:
                    form = "%s%%%ds"    % (form, width)
                elif digits < 1:
                    form = "%s%%%dd"    % (form, width)
                else:
                    form = "%s%%%d.%df" % (form, width, digits)
            form = "%s\n" % form
            # . Prepare instances
            for site in self.sites:
                for instance in site.instances:
                    lines.append (form % (site.siteIndex + 1, instance.instIndex + 1, site.label, instance.label, instance.Gintr, instance.protons))
            # . Write to a file
            WriteInputFile (filename, lines)


    #-------------------------------------------------------------------------------
    def DefineMCModel (self, sampler, log=logFile):
        """Assign a Monte Carlo model to the continuum electrostatic model."""
        if sampler is not None:
            # . Perform initial checks
            checks = (
                isinstance (sampler, MCModelDefault) ,
                isinstance (sampler, MCModelGMCT)    ,)
            if not any (checks):
                raise ContinuumElectrostaticsError ("Cannot define MC model.")

            # . Assign and initialize MC model
            self.sampler = sampler
            self.sampler.Initialize (self)
            self.sampler.PrintPairs (log=log)
            self.sampler.Summary    (log=log)


    #-------------------------------------------------------------------------------
    def CalculateMicrostateEnergy (self, stateVector, pH=7.0):
        """Calculate microstate energy."""
        return self.energyModel.CalculateMicrostateEnergy (stateVector, pH=pH)


    #-------------------------------------------------------------------------------
    def CalculateProbabilities (self, pH=7.0, unfolded=False, isCalculateCurves=False, logFrequency=-1, trajectoryFilename="", log=logFile):
        """Calculate probabilities.

        Setting |trajectoryFilename| will cause writing energies of sampled states to a file."""
        nstates = -1
        sites   = None

        if       hasattr (self, "sampler")  and     unfolded:
            raise ContinuumElectrostaticsError ("Monte Carlo sampling of unfolded proteins unsupported.")
        elif     hasattr (self, "sampler")  and not unfolded:
            self.sampler.CalculateOwnerProbabilities (pH=pH, logFrequency=logFrequency, trajectoryFilename=trajectoryFilename, log=log)
        elif not hasattr (self, "sampler")  and     unfolded:
            if (trajectoryFilename != ""):
                raise ContinuumElectrostaticsError ("Writing trajectories of unfolded proteins unsupported.")
            nstates = self.energyModel.CalculateProbabilitiesAnalyticallyUnfolded (pH=pH)
        elif not hasattr (self, "sampler")  and not unfolded:
            # TODO !!!
            if (trajectoryFilename != ""):
                raise ContinuumElectrostaticsError ("Writing trajectories unsupported.")
            nstates = self.energyModel.CalculateProbabilitiesAnalytically (pH=pH)

        if isCalculateCurves:
            sites = []
            for site in self.sites:
                instances = []
                for instance in site.instances:
                    instances.append (instance.probability)
                sites.append (instances)
        if nstates > 0:
            if LogFileActive (log):
                log.Text ("\nCalculated %d protonation states.\n" % nstates)

        self.isProbability = True
        return sites


    #-------------------------------------------------------------------------------
    def _GetResidueInfo (self, residue):
        system        = self.owner
        ParseLabel    = system.sequence.ParseLabel
        # . Get segment label
        segment       = residue.parent
        segmentName   = segment.label
        # . Get residue label and serial
        residueName, residueSerial = ParseLabel (residue.label, fields=2)
        residueSerial = int (residueSerial)

        return (segmentName, residueName, residueSerial)


    #-------------------------------------------------------------------------------
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


    #-------------------------------------------------------------------------------
    def _SetupBackground (self):
        allSiteAtomIndices  = []
        for site in self.sites:
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
        
                # . Remove residues not defined in PROTEIN_RESIDUES, usually waters and ions
                if residueName not in REMOVE_RESIDUES:
                    atoms = residue.children
        
                    #============ Iterate atoms ============
                    for atom in atoms:
                        proteinAtomIndices.append (atom.index)
                        if atom.index not in allSiteAtomIndices:
                            backAtomIndices.append (atom.index)
        self.proteinAtomIndices = proteinAtomIndices
        self.backAtomIndices    = backAtomIndices


    #-------------------------------------------------------------------------------
    def _SetupSites (self, residue, prevResidue=None, nextResidue=None, terminal=None, log=logFile):
        segmentName, residueName, residueSerial = self._GetResidueInfo (residue)
        setupSites    = []
        if terminal:
            if   terminal == "N":
                if   residueName == "GLY":
                    libTerm = self.library["NGL"]
                elif residueName == "PRO":
                    libTerm = self.library["NPR"]
                else:
                    libTerm = self.library["NTR"]
            elif terminal == "C":
                libTerm = self.library["CTR"]
            termIndices = self._GetIndices (residue, libTerm.atomLabels)
            setupSites.append ([terminal, libTerm, termIndices, termIndices])

        # . Check for a titratable residue
        if residueName in self.library:
            libSite      = self.library[residueName]
            siteIndices  = self._GetIndices (residue, libSite.atomLabels)
            modelIndices = []
            for atom in residue.children:
                if terminal:
                    if (atom.label in TERM_REMOVE) or (atom.index in termIndices):
                        continue
                modelIndices.append (atom.index)

            if not terminal:
                prevIndices  = []
                if prevResidue:
                    foo, prevResidueName, prevResidueSerial = self._GetResidueInfo (prevResidue)
                    if prevResidueName in PROTEIN_RESIDUES:
                        prevIndices = self._GetIndices (prevResidue, PREV_RESIDUE, check=False)
                    modelIndices = prevIndices + modelIndices
    
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
                    modelIndices = modelIndices + nextIndices
            setupSites.append (["SITE", libSite, siteIndices, modelIndices])

            if terminal == "C":
                setupSites.reverse ()
        return setupSites


    #-------------------------------------------------------------------------------
    def _CreateSite (self, **keywordArguments):
        """Create a site and its instances specific to the CE model."""
        pass


    #-------------------------------------------------------------------------------
    def _SplitModel (self, excludeSegments=_DEFAULT_EXCLUDE_SEGMENTS, excludeResidues=None, includeTermini=False, dryRun=True, log=logFile):
        siteIndex        =  0
        instIndexGlobal  =  0
        totalInstances   =  0

        if not dryRun:
            self.sites = []
        system    = self.owner
        segments  = system.sequence.children


        #============ Iterate segments ============
        for segment in segments:
            segmentName = segment.label

            # . Include segment?
            if segmentName not in excludeSegments:
                residues  = segment.children
                nresidues = len (residues)

                #============ Iterate residues ============
                for residueIndex, residue in enumerate (residues):
                    foo, residueName, residueSerial = self._GetResidueInfo (residue)
                    includeResidue = self._CheckResidue (excludeResidues, segmentName, residueName, residueSerial, log=log)
                    if not includeResidue:
                        continue
                    if residueIndex > 0:
                        prevResidue = residues[residueIndex - 1]
                    else:
                        prevResidue = None
                    if residueIndex < (nresidues - 1):
                        nextResidue = residues[residueIndex + 1]
                    else:
                        nextResidue = None
                    if   residueIndex < 1:
                        terminal = "N"
                    elif residueIndex > (nresidues - 2):
                        terminal = "C"
                    else:
                        terminal = None

                    setupSites = self._SetupSites (residue, prevResidue, nextResidue, terminal=(terminal if (includeTermini and (residueName in PROTEIN_RESIDUES)) else None), log=log)
                    for siteType, libSite, siteAtomIndices, modelAtomIndices in setupSites:
                        if not dryRun:
                            if   siteType == "N":
                                updatedSerial = 998
                                updatedName   = libSite.label
                            elif siteType == "C":
                                updatedSerial = 999
                                updatedName   = libSite.label
                            else:
                                updatedSerial = residueSerial
                                updatedName   = residueName

                            # . Create a new site in the CE model
                            instIndexGlobal = self._CreateSite (
                                siteIndex        = siteIndex         ,
                                segName          = segmentName       ,
                                resName          = updatedName       ,
                                resSerial        = updatedSerial     ,
                                siteAtomIndices  = siteAtomIndices   ,
                                modelAtomIndices = modelAtomIndices  ,

                                libSite          = libSite           ,
                                instIndexGlobal  = instIndexGlobal   ,
                                )
                        # . Update the numbers of sites and instances
                        siteIndex      += 1
                        totalInstances += len (libSite.instances)

        # . Return the numbers of sites and instances
        totalSites = siteIndex
        return (totalSites, totalInstances)


    #-------------------------------------------------------------------------------
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
# . Main program
#===============================================================================
if __name__ == "__main__": pass
