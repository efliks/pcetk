#-------------------------------------------------------------------------------
# . File      : Model.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADModel is a class representing the continuum electrostatic model."""

__lastchanged__ = "$Id$"


import os, glob, subprocess, time

from pCore                 import logFile, LogFileActive, Selection, YAMLUnpickle
from Constants             import *
from Error                 import ContinuumElectrostaticsError
from Site                  import MEADSite
from Instance              import InstanceThread
from Utils                 import FormatEntry, ConvertAttribute
from EnergyModel           import EnergyModel

# File handling
from ESTFileReader         import ESTFileReader
from GMCTOutputFileReader  import GMCTOutputFileReader
from InputFileWriter       import WriteInputFile
from PQRFileWriter         import PQRFile_FromSystem


_DefaultTemperature        = 300.0

_DefaultIonicStrength      = 0.1

_DefaultPathScratch        = "/tmp"

_DefaultPathMEAD           = "/usr/bin"

_DefaultPathGMCT           = "/usr/bin"

_DefaultThreads            = 1

_DefaultDoubleFlip         = 2

_DefaultTripleFlip         = 3

_DefaultEquilibrationScans = 500

_DefaultProductionScans    = 20000

_DefaultFocussingSteps     = ((121, 2.00), (101, 1.00), (101, 0.50), (101, 0.25))



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
    "pathGMCT"             :   _DefaultPathGMCT         ,
    "pathScratch"          :   _DefaultPathScratch      ,
    "pathFptSites"         :   None                     ,
    "pathPqrProtein"       :   None                     ,
    "pathPqrBack"          :   None                     ,
    "focussingSteps"       :   _DefaultFocussingSteps   ,
    "librarySites"         :   None                     ,
    "meadSites"            :   None                     ,
    "energyModel"          :   None                     ,
    "owner"                :   None                     ,
        }

  defaultAttributeNames = {
    "Temperature"       : "temperature"        ,    "Initialized"       : "isInitialized"      ,
    "Ionic Strength"    : "ionicStrength"      ,    "Files Written"     : "isFilesWritten"     ,
    "Threads"           : "nthreads"           ,    "Calculated"        : "isCalculated"       ,
    "Split Directories" : "splitToDirectories" ,    "Calculated Prob."  : "isProbability"      ,
    "Delete Job Files"  : "deleteJobFiles"     ,
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
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)

    self.owner       = system
    self.pathMEAD    = os.path.abspath (self.pathMEAD)
    self.pathGMCT    = os.path.abspath (self.pathGMCT)
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
      items = (
         ( "siteID"    , (  12   ,      0    ) ),
         ( "instID"    , (  12   ,      0    ) ),
         ( "siteLabel" , (  16   ,     -1    ) ),
         ( "instLabel" , (  16   ,     -1    ) ),
         ( "Gintr"     , (spacing,  precision) ),
         ( "protons"   , (  12   ,      0    ) ),
              )
      header = FormatEntry (items, header = True)
      entry  = FormatEntry (items)
      lines  = [header]

      for site in self.meadSites:
        for instance in site.instances:
          lines.append (entry % (site.siteIndex + 1, instance.instIndex + 1, site.label, instance.label, instance.Gintr, instance.protons))
      WriteInputFile (filename, lines)


  #===============================================================================
  def CalculateProbabilitiesGMCT (self, pH=7.0, dryRun=False, doubleFlip=_DefaultDoubleFlip, tripleFlip=_DefaultTripleFlip, equilibrationScans=_DefaultEquilibrationScans, productionScans=_DefaultProductionScans, log=logFile):
    """Use GMCT to estimate probabilities of protonation states.

    With |dryRun=True|, GMCT is not called and only the directories and files are created. This is necessary in the parallel mode because the function mkdir does not work with multiple threads."""

    gmctSetupTemplate = """
     blab        1
     nconfflip   10
     limit       %d
     tlimit      %d
     itraj       0
     nmcequi     %d
     nmcfull     %d
     temp        %f
     icorr       0
     nmu         1
     mu          %f  %f  0.0  0  0
    """

    if self.isCalculated:
      potential = -CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10 * pH
      project   = "job"
      sites     = None

      # Prepare input files and directories for GMCT
      dirConf = os.path.join (self.pathScratch, "gmct", "conf")
      if not os.path.exists (dirConf): os.makedirs (dirConf)

      dirCalc = os.path.join (self.pathScratch, "gmct", "%s" % pH)
      if not os.path.exists (dirCalc): os.makedirs (dirCalc)

      fileGint = os.path.join (dirConf, "%s.gint" % project)
      if not os.path.exists (fileGint): self.WriteGintr (fileGint, precision = 8)

      fileInter = os.path.join (dirConf, "%s.inter" % project)
      if not os.path.exists (fileInter): self.WriteW (fileInter, precision = 8)

      fileConf = os.path.join (dirCalc, "%s.conf" % project)
      if not os.path.exists (fileConf): WriteInputFile (fileConf, ["conf  0.0  0.0  0.0\n"])

      fileSetup = os.path.join (dirCalc, "%s.setup" % project)
      if not os.path.exists (fileSetup): WriteInputFile (fileSetup, gmctSetupTemplate % (doubleFlip, tripleFlip, equilibrationScans, productionScans, self.temperature, potential, potential))

      linkname = os.path.join (dirCalc, "conf")
      if not os.path.exists (linkname): os.symlink ("../conf", linkname)


      if not dryRun:
        output = os.path.join (dirCalc, "%s.gmct-out" % project)
        error  = os.path.join (dirCalc, "%s.gmct-err" % project)

        if os.path.exists (os.path.join (dirCalc, output)):
          pass
        else:
          command = [os.path.join (self.pathGMCT, "gmct"), project]
          try:
            out = open (output, "w")
            err = open (error,  "w")
            subprocess.check_call (command, stderr = err, stdout = out, cwd = dirCalc)
            out.close ()
            err.close ()
          except:
            raise ContinuumElectrostaticsError ("Failed running command: %s" % " ".join (command))

        # Read probabilities from the output file
        reader = GMCTOutputFileReader (output)
        reader.Parse (temperature = self.temperature)

        # Construct a two-dimensional list of M-sites, each site N-instances, initiated with zeros
        sites = []
        for site in self.meadSites:
          instances = []
          for instance in site.instances:
            instances.append (0.)
          sites.append (instances)

        for siteIndex, site in enumerate (self.meadSites):
          for instanceIndex, instance in enumerate (site.instances):
            key                             = "conf_%s_%s%d_%s" % (site.segName, site.resName, site.resSerial, instance.label)
            probability                     = reader.probabilities[key][0]
            sites[siteIndex][instanceIndex] = probability
            instance.probability            = probability

      # The instances now contain calculated probabilities
      self.isProbability = True

      # Return the two-dimensional list (useful for calculating titration curves)
      return sites


  #===============================================================================
  # Wrapper methods for the calculation of energies and protonation states
  def CalculateMicrostateEnergy (self, stateVector, pH=7.0):
    """Calculate energy of a microstate defined by the state vector."""
    return self.energyModel.CalculateMicrostateEnergy (stateVector, pH=pH)


  def CalculateProbabilitiesMonteCarlo (self, pH=7.0, nequi=_DefaultEquilibrationScans, nprod=_DefaultProductionScans, log=logFile):
    """Calculate the probability of occurence of each instance of each site, using the in-house implementation of Metropolis Monte Carlo (experimental)."""
    self.energyModel.CalculateProbabilitiesMonteCarlo (pH=pH, nequi=nequi, nprod=nprod, log=log)
    return self._FinalizeProbabilities ()


  def CalculateProbabilitiesAnalytically (self, pH=7.0, log=logFile):
    """Calculate the probability of occurence of each instance of each site, using statistical mechanics."""
    nstates = self.energyModel.CalculateProbabilitiesAnalytically (pH=pH)

    if LogFileActive (log):
      log.Text ("\nCalculated %d protonation states.\n" % nstates)
    return self._FinalizeProbabilities ()


  def _FinalizeProbabilities (self):
    # The instances now contain calculated probabilities
    self.isProbability = True

    # It is necessary in parallel mode to return a list of probabilities
    sites = []
    for site in self.meadSites:
      instances = []
      for instance in site.instances:
        instances.append (instance.probability)
      sites.append (instances)
    return sites


  #===============================================================================
  def CalculateElectrostaticEnergies (self, calculateETA=True, asymmetricTolerance=0.05, asymmetricSummary=False, log=logFile):
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

        heads = [ ("Instance of a site" , 4),
                  ("Gborn_model"        , 0),
                  ("Gback_model"        , 0),
                  ("Gborn_protein"      , 0),
                  ("Gback_protein"      , 0),
                  ("Gmodel"             , 0),
                  ("Gintr"              , 0), ]
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
      self.energyModel.SymmetrizeInteractions ()

      if LogFileActive (log):
        log.Text ("\nSymmetrizing interactions complete.\n")

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
          heads = [ ("Instance of a site A" , 4),
                    ("Instance of a site B" , 4),
                    ("Deviation"            , 0), ]
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
  def Initialize (self, excludeSegments=None, excludeResidues=None, includeTermini=False, log=logFile):
    """Decompose the system into model compounds, sites and a background charge set.
    |excludeSegments| is a sequence of segment names to exclude from the model, usually segments of water molecules.
    |excludeResidues| is a sequence of three-element sequences (segmentName, residueName, residueSerial).
    It is possible to leave some of the elements blank, for example ("PRTA", "CYS", "") means exclude all cysteines in segment PRTA.
    Only C-terminus is supported and the support is experimental."""

    # Check for the CHARMM energy model
    system = self.owner

    if system.energyModel.mmModel.label is not "CHARMM":
      raise ContinuumElectrostaticsError ("The energy model of the system is different from CHARMM.")

    if not self.isInitialized:
      ParseLabel = system.sequence.ParseLabel
      segments   = system.sequence.children

      instIndexGlobal = 0
      siteIndex       = 0
      self.meadSites  = []

      # A temporary list of protons is needed because the central array of protons is initialized only after the initialization of all instances
      protons         = []

      if excludeSegments is None:
        excludeSegments = ["WATA", "WATB", "WATC", "WATD", ]


      #============ Go over segments ============
      for segment in segments:
        segmentName = segment.label

        # Include segment?
        if segmentName not in excludeSegments:
          residues  = segment.children
          nresidues = len (residues)


          #============ Go over residues ============
          for residueIndex, residue in enumerate (residues):
            residueName, residueSerial = ParseLabel (residue.label, fields = 2)
            residueSerial = int (residueSerial)

            # Include residue?
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


            if includeResidue:

              #============= Start experimental code =============
              terminus = None
              excludeTerminusNames = []

              if includeTermini:
                if residueName in PROTEIN_RESIDUES:

                  # Check for C-terminus
                  if residueIndex > (nresidues - 2):
                    terminus       = "CTER"
                    libTerminus    = self.librarySites[terminus]
                    terminusSerial = residueSerial + 1

                  # Check for N-terminus
                  #elif residueIndex < 1:
                  #  if   residueName == "GLY":
                  #    terminus = "GLYP"
                  #  elif residueName == "PRO":
                  #    terminus = "PROP"
                  #  else:
                  #    terminus = "NTER"
                  #  libTerminus    = self.librarySites[terminus]
                  #  terminusSerial = residueSerial - 1

                  if terminus:
                    excludeTerminusNames = libTerminus["atoms"]
              #============== End experimental code ==============


              # Titratable residue?
              if residueName in self.librarySites:
                prevIndices = []
                nextIndices = []

                # Include atoms from the previous residue to the model compound?
                if residueIndex > 1:
                  prevResidue = residues[residueIndex - 1]
                  prevResidueName, prevResidueSerial = ParseLabel (prevResidue.label, fields = 2)
                  prevResidueSerial = int (prevResidueSerial)

                  if prevResidueName in PROTEIN_RESIDUES:
                    prevNames = PREV_RESIDUE

                    for atom in prevResidue.children:
                      if atom.label in prevNames:
                        prevIndices.append (atom.index)

                # Include atoms from the next residue to the model compound?
                if residueIndex < (nresidues - 1):
                  nextResidue = residues[residueIndex + 1]
                  nextResidueName, nextResidueSerial = ParseLabel (nextResidue.label, fields = 2)
                  nextResidueSerial = int (nextResidueSerial)

                  if nextResidueName in PROTEIN_RESIDUES:
                    if   nextResidueName == "PRO":
                      nextNames = NEXT_RESIDUE_PRO
                    elif nextResidueName == "GLY":
                      nextNames = NEXT_RESIDUE_GLY
                    else:
                      nextNames = NEXT_RESIDUE

                    for atom in nextResidue.children:
                      if atom.label in nextNames:
                        nextIndices.append (atom.index)


                # Collect atom indices
                libSite          = self.librarySites[residueName]
                libSiteAtoms     = libSite["atoms"]
                atoms            = residue.children
                modelAtomIndices = prevIndices
                siteAtomIndices  = []
                missingNames     = []

                for libAtomName in libSiteAtoms:
                  index = -1
                  for atom in atoms:
                    if libAtomName == atom.label:
                      index = atom.index
                      break
                  if index >= 0:
                    siteAtomIndices.append (index)
                  else:
                    missingNames.append (libAtomName)
                if missingNames:
                  raise ContinuumElectrostaticsError ("Cannot include site %s %s %d because of missing atoms: %s" % (segmentName, residueName, residueSerial, " ".join (missingNames)))

                for atom in atoms:
                  if atom.label not in excludeTerminusNames:
                    modelAtomIndices.append (atom.index)
                modelAtomIndices.extend (nextIndices)


                # Create a site (add instances later)
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
                protonsOfInstances, updatedIndexGlobal = newSite._CreateInstances (libSite["instances"], instIndexGlobal)

                protons.extend (protonsOfInstances)
                instIndexGlobal = updatedIndexGlobal

                # Add the newly created site to the list of sites
                self.meadSites.append (newSite)
                siteIndex = siteIndex + 1


              #============= Start experimental code =============
              if terminus:
                libSiteAtoms     = libTerminus["atoms"]
                atoms            = residue.children
                siteAtomIndices  = []
                modelAtomIndices = []
                missingNames     = []

                for libAtomName in libSiteAtoms:
                  index = -1
                  for atom in atoms:
                    if libAtomName == atom.label:
                      index = atom.index
                      break
                  if index >= 0:
                    siteAtomIndices.append (index)
                  else:
                    missingNames.append (libAtomName)
                if missingNames:
                  raise ContinuumElectrostaticsError ("Cannot include terminus %s %s %d because of missing atoms: %s" % (segmentName, residueName, residueSerial, " ".join (missingNames)))

                # Construct a model compound
                if residueIndex > 1:
                  prevResidue = residues[residueIndex - 1]
                  prevNames = PREV_RESIDUE
                  for atom in prevResidue.children:
                    if atom.label in prevNames:
                      modelAtomIndices.append (atom.index)

                # Do not include atoms in the model compound which are already part of another site
                if self.librarySites.has_key (residueName):
                  excludeSiteNames = self.librarySites[residueName]["atoms"]
                else:
                  excludeSiteNames = []

                for atom in atoms:
                  if atom.label not in excludeSiteNames:
                    modelAtomIndices.append (atom.index)

                newSite = MEADSite (
                      parent           = self              ,
                      siteIndex        = siteIndex         ,
                      segName          = segmentName       ,
                      resName          = terminus          ,
                      resSerial        = terminusSerial    ,
                      siteAtomIndices  = siteAtomIndices   ,
                      modelAtomIndices = modelAtomIndices  ,
                                   )
                newSite.CalculateCenterOfGeometry (system, libSite["center"])
                protonsOfInstances, updatedIndexGlobal = newSite._CreateInstances (libSite["instances"], instIndexGlobal)

                protons.extend (protonsOfInstances)
                instIndexGlobal = updatedIndexGlobal

                self.meadSites.append (newSite)
                siteIndex = siteIndex + 1
              #============== End experimental code ==============



      # Construct the background set of charges and the protein (to be used as eps2set_region)
      allSiteAtomIndices  = []
      for site in self.meadSites:
        allSiteAtomIndices.extend (site.siteAtomIndices)

      backAtomIndices     = []
      proteinAtomIndices  = []

      #============ Go over segments ============
      for segment in segments:
        residues = segment.children

        #============ Go over residues ============
        for residue in residues:
          residueName, residueSerial = ParseLabel (residue.label, fields = 2)
          residueSerial = int (residueSerial)

          # Remove residues not defined in PROTEIN_RESIDUES, usually waters and ions
          if residueName not in REMOVE_RESIDUES:
            atoms = residue.children

            #============ Go over atoms ============
            for atom in atoms:
              proteinAtomIndices.append (atom.index)

              if atom.index not in allSiteAtomIndices:
                backAtomIndices.append (atom.index)


      self.proteinAtomIndices = proteinAtomIndices
      self.backAtomIndices    = backAtomIndices

      self.pathPqrProtein     = os.path.join (self.pathScratch, "protein.pqr")
      self.pathPqrBack        = os.path.join (self.pathScratch, "back.pqr")
      self.pathFptSites       = os.path.join (self.pathScratch, "site.fpt")

      # Allocate arrays of protons, intrinsic energies, interaction energies and probabilities
      self.energyModel = EnergyModel (self)

      # Initialize the array of protons
      for site in self.meadSites:
        for instance in site.instances:
          instance.protons = protons[instance._instIndexGlobal]
          # Or: self.energyModel.SetProtons (instance._instIndexGlobal, protons[instance._instIndexGlobal])

      # Finish up
      self.isInitialized = True


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
          entries = (
              ( "%d"      %  (site.siteIndex + 1) ),
              ( "%s"      %  site.segName         ),
              ( "%s"      %  site.resName         ),
              ( "%d"      %  site.resSerial       ),
              ( "%d"      %  len (site.instances) ),
              ( "%10.3f"  %  site.center[0]       ),
              ( "%10.3f"  %  site.center[1]       ),
              ( "%10.3f"  %  site.center[2]       ),
                    )
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
                        "TYR" : ("d", ),
                          }

    if LogFileActive (log):
      if self.isProbability:
        maxinstances = 0
        for site in self.meadSites:
          ninstances = len (site.instances)
          if ninstances > maxinstances: maxinstances = ninstances

        tab = log.GetTable (columns = [6, 6, 6] + [8, 8] * maxinstances)
        tab.Start ()
        tab.Heading ("Site", columnSpan = 3)
        tab.Heading ("Probabilities of instances", columnSpan = maxinstances * 2)

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
                        "TYR" : ("d", ),
                          }
    translateLabels = {
                    "p" :   "protonated",
                    "d" : "deprotonated",
                      }
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
