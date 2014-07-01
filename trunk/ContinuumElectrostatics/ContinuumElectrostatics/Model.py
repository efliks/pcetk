#-------------------------------------------------------------------------------
# . File      : Model.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADModel is a class representing the continuum electrostatic model.

CurveThread is a class for running parallel calculations of titration curves."""

__lastchanged__ = "$Id$"

import os, glob, math, threading, subprocess, time

from pCore                 import logFile, LogFileActive, Selection, Vector3, YAMLUnpickle, Clone
from Constants             import *
from Error                 import ContinuumElectrostaticsError
from Site                  import MEADSite
from Instance              import MEADInstance, InstanceThread 
from Toolbox               import FormatEntry, ConvertAttribute
from StateVector           import StateVector

# File handling
from ESTFileReader         import ESTFileReader
from GMCTOutputFileReader  import GMCTOutputFileReader
from InputFileWriter       import WriteInputFile
from PQRFileWriter         import PQRFile_FromSystem


_DefaultTemperature     = 300.0

_DefaultIonicStrength   = 0.1

_DefaultMeadPath        = "/usr/bin"

_DefaultGmctPath        = "/usr/bin"

_DefaultScratch         = "/tmp"

_DefaultThreads         = 1

_DefaultCleanUp         = False

_DefaultFocussingSteps  = ((121, 2.00), (101, 1.00), (101, 0.50), (101, 0.25))

_DefaultGmctSetup       = """
blab        1
nconfflip   10
tlimit      3
itraj       0
nmcfull     20000
temp        %f
icorr       0
limit       2
nmcequi     500
nmu         1
mu          %f  %f  0.0  0  0
"""


#-------------------------------------------------------------------------------
class CurveThread (threading.Thread):
  """Calculate each pH-step in a separate thread."""

  def __init__ (self, meadModel, isAnalytic, pH):
    threading.Thread.__init__ (self)
    self.model = meadModel
    self.sites = None
    self.pH    = pH
    if isAnalytic:
      self.calculate = meadModel.CalculateProbabilitiesAnalytically
    else:
      self.calculate = meadModel.CalculateProbabilitiesGMCT


  def run (self):
    """The method that runs the calculations."""
    self.sites = self.calculate (self.pH, log = None)


#-------------------------------------------------------------------------------
class MEADModel (object):
  """Continuum electrostatic model."""

  defaultAttributes = {
    "temperature"        :  _DefaultTemperature     ,
    "ionicStrength"      :  _DefaultIonicStrength   ,
    "meadPath"           :  _DefaultMeadPath        ,
    "gmctPath"           :  _DefaultGmctPath        ,
    "nthreads"           :  _DefaultThreads         ,
    "deleteJobFiles"     :  _DefaultCleanUp         ,
    "scratch"            :  _DefaultScratch         ,
    "focussingSteps"     :  _DefaultFocussingSteps  ,
    "librarySites"       :  None                    ,
    "meadSites"          :  None                    ,
    "backAtomIndices"    :  None                    ,
    "proteinAtomIndices" :  None                    ,
    "backPqr"            :  None                    ,
    "proteinPqr"         :  None                    ,
    "sitesFpt"           :  None                    ,
    "splitToDirectories" :  True                    ,
    "isInitialized"      :  False                   ,
    "isFilesWritten"     :  False                   ,
    "isCalculated"       :  False                   ,
    "isProbability"      :  False                   ,
        }

  defaultAttributeNames = {
    "Temperature"       : "temperature"        ,    "Initialized"       : "isInitialized"      ,
    "Ionic Strength"    : "ionicStrength"      ,    "Files Written"     : "isFilesWritten"     ,
    "Threads"           : "nthreads"           ,    "Calculated"        : "isCalculated"       ,
    "Split Directories" : "splitToDirectories" ,    "Calculated prob."  : "isProbability"      ,
    "Delete Job Files"  : "deleteJobFiles"     ,
        }

  def __del__ (self):
    """Deallocation."""
    if self.deleteJobFiles: self.DeleteJobFiles ()


  def __init__ (self, log = logFile, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)

    self.meadPath = os.path.abspath (self.meadPath)
    self.gmctPath = os.path.abspath (self.gmctPath)
    self.scratch  = os.path.abspath (self.scratch)
    self.LoadLibraryOfSites (log = log)


  def LoadLibraryOfSites (self, log = logFile):
    """Load a set of YAML or EST files with parameters for titratable sites.

    If there are YAML or EST files in the current directory, they are loaded as well.

    If these additional files have names coinciding with the names from the library, the library parameters will be overwritten.

    Notice! EST files overwrite YAML files."""
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


  def WriteW (self, filename = "W.dat", log = logFile):
    """Write an interaction matrix compatible with GMCT."""
    if self.isCalculated:
      items = (
         ("idSite1"  , ( 8,   0) ),
         ("idInst1"  , ( 8,   0) ),
         ("labSite1" , (14,  -1) ),
         ("labInst1" , (14,  -1) ),
         ("idSite2"  , ( 8,   0) ),
         ("idInst2"  , ( 8,   0) ),
         ("labSite2" , (14,  -1) ),
         ("labInst2" , (14,  -1) ),
         ("Wij_symm" , (16,   8) ),
         ("Wij"      , (16,   8) ),
         ("Wij_err"  , (16,   8) ),
              )
      header = FormatEntry (items, header = True)
      entry  = FormatEntry (items)
      lines  = [header]

      for sitea in self.meadSites:
        for insta in sitea.instances:
          for siteb in self.meadSites:
            for instb in siteb.instances:
              Wij      = insta.interactions[siteb.siteID - 1][instb.instID - 1]
              Wji      = instb.interactions[sitea.siteID - 1][insta.instID - 1]
              Wij_symm = 0.5 * (Wij + Wji)
              Wij_err  = Wij_symm - Wij
              line = entry % (
                      sitea.siteID, insta.instID, sitea.label, insta.label,
                      siteb.siteID, instb.instID, siteb.label, instb.label,
                      Wij_symm, Wij, Wij_err
                             )
              lines.append (line)

      WriteInputFile (filename, lines)


  def WriteGintr (self, filename = "gintr.dat", log = logFile):
    """Iterate over instances and write a gintr.dat file compatible with GMCT.

    This file contains Gintr of each instance of each site."""
    if self.isCalculated:
      items = (
         ("siteID"    , (12,  0) ),
         ("instID"    , (12,  0) ),
         ("siteLabel" , (16, -1) ),
         ("instLabel" , (16, -1) ),
         ("Gintr"     , (16,  8) ),
         ("protons"   , (12,  0) ),
              )
      header = FormatEntry (items, header = True)
      entry  = FormatEntry (items)
      lines  = [header]

      for site in self.meadSites:
        for instance in site.instances:
          line          = entry % (site.siteID, instance.instID, site.label, instance.label, instance.Gintr, instance.protons)
          lines.append (line)
      WriteInputFile (filename, lines)


  def CalculateCurves (self, isAnalytic = False, curveSampling = 0.5, curveStart = 0.0, curveStop = 14.0, directory = "curves", forceSerial = False, log = logFile):
    """Calculate titration curves."""
    if self.isCalculated:
      nsteps = int ((curveStop - curveStart) / curveSampling + 1)
      steps  = []
      tab    = None

      if LogFileActive (log):
        if self.nthreads < 2 or forceSerial:
          log.Text ("\nStarting serial run.\n")
        else:
          log.Text ("\nStarting parallel run on %d CPUs.\n" % self.nthreads)

        tab = log.GetTable (columns = [10, 10])
        tab.Start ()
        tab.Heading ("Step")
        tab.Heading ("pH")


      if self.nthreads < 2 or forceSerial:
        for step in range (0, nsteps):
          pH = curveStart + step * curveSampling
          if isAnalytic:
            result = self.CalculateProbabilitiesAnalytically (pH = pH, log = None)
          else:
            result = self.CalculateProbabilitiesGMCT (pH = pH, log = None)
          steps.append (result)
     
          if tab:
            tab.Entry ("%10d"   % step)
            tab.Entry ("%10.2f" % pH)
      else:
        limit   = self.nthreads - 1
        batches = []
        batch   = []
        for step in range (0, nsteps):
          batch.append (CurveThread (self, isAnalytic, curveStart + step * curveSampling))
          if len (batch) > limit:
            batches.append (batch)
            batch = []
        if batch:
          batches.append (batch)
    
        # If GMCT is to be used, perform a dry run in serial mode to create directories and files
        if not isAnalytic:
          for step in range (0, nsteps):
            self.CalculateProbabilitiesGMCT (pH = curveStart + step * curveSampling, dryRun = True, log = None)

        step = 0 
        for batch in batches:
          for thread in batch:
            thread.start ()
          for thread in batch:
            thread.join ()
    
          for thread in batch:
            steps.append (thread.sites)
            if tab:
              tab.Entry ("%10d"   % step)
              tab.Entry ("%10.2f" % (curveStart + step * curveSampling))
            step = step + 1

      if LogFileActive (log):
        tab.Stop ()
        log.Text ("\nCalculating titration curves complete.\n")

      # Set the flag back to False because the probabilities of instances become senseless after the calculation of curves
      self.isProbability = False
  
      # Write results to files
      if not os.path.exists (directory):
        try:
          os.mkdir (directory)
        except:
          raise ContinuumElectrostaticsError ("Cannot create directory %s" % directory)
  
      # For each instance of each site, write a curve file
      for siteIndex, site in enumerate (self.meadSites):
        for instanceIndex, instance in enumerate (site.instances):
          lines = []
          for step in range (0, nsteps):
            lines.append ("%f %f\n" % (curveStart + step * curveSampling, steps[step][siteIndex][instanceIndex]))
          filename = os.path.join (directory, "%s_%s.dat" % (site.label, instance.label))
          WriteInputFile (filename, lines)
  
      if LogFileActive (log):
        log.Text ("\nWriting curve files complete.\n")


  def CalculateProbabilitiesGMCT (self, pH = 7.0, dryRun = False, log = logFile):
    """Use GMCT to estimate protonation probabilities.

    With |dryRun = True|, GMCT is not called and only the directories and files are created. This is necessary in the parallel mode because the function mkdir does not work with multiple threads."""
    if self.isCalculated:
      potential = -CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10 * pH
      project   = "job"
      sites     = None

      # Prepare input files and directories for GMCT
      dirConf = os.path.join (self.scratch, "gmct", "conf")
      if not os.path.exists (dirConf): os.makedirs (dirConf)

      dirCalc = os.path.join (self.scratch, "gmct", "%s" % pH)
      if not os.path.exists (dirCalc): os.makedirs (dirCalc)

      fileGint = os.path.join (dirConf, "%s.gint" % project)
      if not os.path.exists (fileGint): self.WriteGintr (fileGint)

      fileInter = os.path.join (dirConf, "%s.inter" % project)
      if not os.path.exists (fileInter): self.WriteW (fileInter)

      fileConf = os.path.join (dirCalc, "%s.conf" % project)
      if not os.path.exists (fileConf): WriteInputFile (fileConf, ["conf  0.0  0.0  0.0\n"])

      fileSetup = os.path.join (dirCalc, "%s.setup" % project)
      if not os.path.exists (fileSetup): WriteInputFile (fileSetup, _DefaultGmctSetup % (self.temperature, potential, potential))

      linkname = os.path.join (dirCalc, "conf")
      if not os.path.exists (linkname): os.symlink ("../conf", linkname)


      if not dryRun:
        output = os.path.join (dirCalc, "%s.gmct-out" % project)
        error  = os.path.join (dirCalc, "%s.gmct-err" % project)
  
        if os.path.exists (os.path.join (dirCalc, output)):
          pass
        else:
          command = [os.path.join (self.gmctPath, "gmct"), project]
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


  def CalculateProbabilitiesAnalytically (self, pH = 7.0, log = logFile):
    """For each site, calculate the probability of occurance of each instance, using the Boltzmann weighted sum."""
    if self.isCalculated:
      nsites = len (self.meadSites)
      if nsites > ANALYTIC_SITES:
        raise ContinuumElectrostaticsError ("Too many sites for analytic treatment (%d)\n" % nsites)
  
      # Calculate all state energies
      increment     = True
      stateEnergies = []
      stateVector   = StateVector (self)
      stateVector.Reset ()

      while increment:
        energy    = self.CalculateMicrostateEnergy (stateVector, pH = pH)
        increment = stateVector.Increment ()
        stateEnergies.append (energy)

      # Find the minimum energy
      energyZero = min (stateEnergies)
  
      if LogFileActive (log):
        log.Text ("\nSearching for the minimum energy complete.\n")
  
      # Go over all calculated state energies and calculate Boltzmann factors
      bfactors = []
      beta     = -1.0 / (CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature)
      nstates  = len (stateEnergies)
      for stateIndex in range (0, nstates):
        bfactor = math.exp (beta * (stateEnergies[stateIndex] - energyZero))
        bfactors.append (bfactor)
  
      if LogFileActive (log):
        log.Text ("\nCalculating Boltzmann factors complete.\n")
 
      # Construct a two-dimensional list of M-sites, each site N-instances, initiated with zeros
      sites = []
      for site in self.meadSites:
        instances = []
        for instance in site.instances:
          instances.append (0.)
        sites.append (instances)

      # Calculate the probabilities 
      stateVector.Reset ()
      increment  = True
      stateIndex = 0

      while increment:
        for siteIndex, site in enumerate (self.meadSites):
          instanceIndex                   = stateVector[siteIndex]
          probability                     = sites[siteIndex][instanceIndex]
          probability                     = probability + bfactors[stateIndex]
          sites[siteIndex][instanceIndex] = probability
        increment  = stateVector.Increment ()
        stateIndex = stateIndex + 1

      bsum = 1.0 / sum (bfactors)
      for siteIndex, site in enumerate (sites):
        for instanceIndex, instance in enumerate (site):
          probability                     = sites[siteIndex][instanceIndex] * bsum
          sites[siteIndex][instanceIndex] = probability
  
      if LogFileActive (log):
        log.Text ("\nCalculating protonation probabilities complete.\n")

      # Copy the calculated probabilities into the MEADModel
      for siteIndex, site in enumerate (self.meadSites):
        for instanceIndex, instance in enumerate (site.instances):
          instance.probability = sites[siteIndex][instanceIndex]

      # The instances now contain calculated probabilities
      self.isProbability = True

      # Return the two-dimensional list (useful for calculating titration curves)
      return sites


  def CalculateMicrostateEnergy (self, stateVector, pH = 7.0):
    """Calculate energy of a protonation state (=microstate).

    The protonation state is defined by a state vector. 

    The energy is calculated at a given pH."""
    Gmicro = None
    if self.isCalculated:
      totalGintr    = 0.
      totalInteract = 0.
      nprotons      = 0
      nsites        = len (self.meadSites)

      for siteIndex in range (0, nsites):
        site          = self.meadSites [siteIndex]
        instanceIndex = stateVector    [siteIndex]
        instance      = site.instances [instanceIndex]
        Gintr         = instance.Gintr
        cprotons      = instance.protons
        interactions  = instance.interactions
        totalGintr    = totalGintr + Gintr
        nprotons      = nprotons + cprotons

        for siteIndexInner in range (0, siteIndex):
          instanceIndexInner = stateVector  [siteIndexInner]
          interaction        = interactions [siteIndexInner]
          totalInteract      = totalInteract + interaction [instanceIndexInner]
      Gmicro = totalGintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10 * pH) + totalInteract

    return Gmicro
# Slower but more accurate?
#        for siteIndexInner in range (0, nsites):
#          instanceIndexInner = stateVector  [siteIndexInner]
#          interaction        = interactions [siteIndexInner]
#          totalInteract      = totalInteract + interaction [instanceIndexInner]
#
#      protonChemicalPotential = -CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10 * pH
#      Gmicro = totalGintr - nprotons * protonChemicalPotential + 0.5 * totalInteract


  def CalculateElectrostaticEnergies (self, log = logFile):
    """
    Calculate for each instance of each site:
    - self (Born) energy in the model compound
    - interaction energy   between the site and the background charge set of the model compound

    - self (Born) energy in the protein
    - interaction energy   between the site and the background charge set of the protein
    - interaction energies between the site and the other sites in their different protonation forms

    Finally, use the calculated homotransfer energies to calculate Gintr from Gmodel.
    """
    if self.isFilesWritten:
      tab = None

      if LogFileActive (log):
        if self.nthreads < 2:
          log.Text ("\nStarting serial run.\n")
        else:
          log.Text ("\nStarting parallel run on %d CPUs.\n" % self.nthreads)

        tab = log.GetTable (columns = [6, 6, 6, 6, 16, 16, 16, 16, 16, 16, 16])
        tab.Start ()
        tab.Heading ("Instance of a site", columnSpan = 4)
        tab.Heading ("Gborn_model"  )
        tab.Heading ("Gback_model"  )
        tab.Heading ("Gborn_protein")
        tab.Heading ("Gback_protein")
        tab.Heading ("Gmodel"       )
        tab.Heading ("Gintr"        )
        tab.Heading ("ETA"          )

      # Total number of instances is needed for calculating ETA
      times      = []
      ninstances = 0
      for site in self.meadSites:
        ninstances = ninstances + len (site.instances)


      if self.nthreads < 2:
        for meadSite in self.meadSites:
          for instance in meadSite.instances:
            time0 = time.time ()
            instance.CalculateSiteInModelCompound (log)
            instance.CalculateSiteInProtein (log)
            instance.CalculateGintr (log)

            times.append (time.time () - time0)
            averageTimePerInstance = sum (times) / len (times)
            ninstances = ninstances - 1
            instance.TableEntry (tab, secondsToCompletion = averageTimePerInstance * ninstances)
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

          # Collect times of execution
          nthreads = len (batch)
          for thread in batch: times.append (thread.time)
         
          averageTimePerInstance = sum (times) / len (times) / nthreads
          ninstances = ninstances - nthreads
          secondsToCompletion = averageTimePerInstance * ninstances

          # Print the results at the end of each batch, otherwise they come in random order
          for thread in batch:
            instance = thread.instance
            instance.TableEntry (tab, secondsToCompletion = secondsToCompletion)

      if tab:
        tab.Stop ()
        log.Text ("\nCalculating electrostatic energies complete.\n")

      self.isCalculated = True


  def Initialize (self, system, excludeSegments = None, excludeResidues = None, log = logFile):
    """Decompose the system into model compounds, sites and a background charge set.

    |excludeSegments| is a sequence of segment names to exclude.

    |excludeResidues| is a sequence of three-element sequences (segmentName, residueName, residueSerial).

    It is possible to leave some of the elements blank, for example ("PRTA", "CYS", "") means exclude all cysteines in segment PRTA."""

    # Check for the CHARMM energy model
    if system.energyModel.mmModel.label is not "CHARMM":
      raise ContinuumElectrostaticsError ("The energy model of the system is different from CHARMM.")

    if not self.isInitialized:
      ParseLabel = system.sequence.ParseLabel
      segments   = system.sequence.children

      # instIndexGlobal becomes useful after the interactions are centralized in a common matrix
      instIndexGlobal = 0
      siteIndex       = 0
      self.meadSites  = []

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

                for libAtomName in libSiteAtoms:
                  for atom in atoms:
                    if libAtomName == atom.label:
                      siteAtomIndices.append (atom.index)

                for atom in atoms:
                  modelAtomIndices.append (atom.index)
                modelAtomIndices.extend (nextIndices)


                # Create instances
                libSiteInstances = libSite["instances"]
                instances        = []
  
                for instIndex, instance in enumerate (libSiteInstances):
                  label    = instance [ "label"   ]
                  protons  = instance [ "protons" ]
                  charges  = instance [ "charges" ]
                  Gmodel   = instance [ "Gmodel"  ] * 300.0 / self.temperature

                  if self.splitToDirectories:
                    modelPqr  = os.path.join (self.scratch, segmentName, "%s%d" % (residueName, residueSerial), "%s_%s.%s" % ("model", label, "pqr"))
                    modelLog  = os.path.join (self.scratch, segmentName, "%s%d" % (residueName, residueSerial), "%s_%s.%s" % ("model", label, "out"))
                    modelGrid = os.path.join (self.scratch, segmentName, "%s%d" % (residueName, residueSerial), "%s_%s.%s" % ("model", label, "mgm"))
                    sitePqr   = os.path.join (self.scratch, segmentName, "%s%d" % (residueName, residueSerial), "%s_%s.%s" % ("site",  label, "pqr"))
                    siteLog   = os.path.join (self.scratch, segmentName, "%s%d" % (residueName, residueSerial), "%s_%s.%s" % ("site",  label, "out"))
                    siteGrid  = os.path.join (self.scratch, segmentName, "%s%d" % (residueName, residueSerial), "%s_%s.%s" % ("site",  label, "ogm"))
                  else:
                    modelPqr  = os.path.join (self.scratch, "%s_%s_%s_%d_%s.%s" % ("model", segmentName, residueName, residueSerial, label, "pqr"))
                    modelLog  = os.path.join (self.scratch, "%s_%s_%s_%d_%s.%s" % ("model", segmentName, residueName, residueSerial, label, "out"))
                    modelGrid = os.path.join (self.scratch, "%s_%s_%s_%d_%s.%s" % ("model", segmentName, residueName, residueSerial, label, "mgm"))
                    sitePqr   = os.path.join (self.scratch, "%s_%s_%s_%d_%s.%s" % ("site",  segmentName, residueName, residueSerial, label, "pqr"))
                    siteLog   = os.path.join (self.scratch, "%s_%s_%s_%d_%s.%s" % ("site",  segmentName, residueName, residueSerial, label, "out"))
                    siteGrid  = os.path.join (self.scratch, "%s_%s_%s_%d_%s.%s" % ("site",  segmentName, residueName, residueSerial, label, "ogm"))

                  # Set the parent later
                  newInstance = MEADInstance (
                                 instID          = instIndex + 1   ,
                                 instIndexGlobal = instIndexGlobal ,
                                 label           = label           ,
                                 protons         = protons         ,
                                 charges         = charges         ,
                                 Gmodel          = Gmodel          ,
                                 modelPqr        = modelPqr        ,
                                 modelLog        = modelLog        ,
                                 modelGrid       = modelGrid       ,
                                 sitePqr         = sitePqr         ,
                                 siteLog         = siteLog         ,
                                 siteGrid        = siteGrid        ,
                                             )
                  instances.append (newInstance)
                  instIndexGlobal = instIndexGlobal + 1


                # Calculate center of geometry
                centralAtom = libSite["center"]
                if centralAtom:
                  atoms = residue.children

                  for atom in atoms:
                    if atom.label == centralAtom:
                      centralIndex = atom.index
                      break
                  center = system.coordinates3[centralIndex]

                else:
                  center = Vector3 ()
                  natoms = len (siteAtomIndices)

                  for atomIndex in siteAtomIndices:
                    center.AddScaledVector3 (1.0, system.coordinates3[atomIndex])
                  center.Scale (1.0 / natoms)


                # Create a site
                newSite = MEADSite (
                               parent           = self             ,
                               siteID           = siteIndex + 1    ,
                               segName          = segmentName      ,
                               resName          = residueName      ,
                               resSerial        = residueSerial    ,
                               modelAtomIndices = modelAtomIndices ,
                               siteAtomIndices  = siteAtomIndices  ,
                               instances        = instances        ,
                               center           = center           ,
                                   )
                for instance in newSite.instances:
                  instance.parent = newSite    # <--Setting parents

                self.meadSites.append (newSite)
                siteIndex = siteIndex + 1


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

      self.proteinPqr         = os.path.join (self.scratch, "protein.pqr")
      self.backPqr            = os.path.join (self.scratch, "back.pqr")

      # Define FPT-file
      self.sitesFpt = os.path.join (self.scratch, "site.fpt")

      self.isInitialized = True


  def Summary (self, log = logFile):
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

      nsites     = len (self.meadSites)
      ninstances = 0
      for site in self.meadSites:
        ninstances = ninstances + len (site.instances)

      summary.Entry ("Number Of Sites",     "%s" % nsites)
      summary.Entry ("Number Of Instances", "%s" % ninstances)
      summary.Stop ()


  def SummarySites (self, log = logFile):
    """List titratable residues."""
    if LogFileActive (log):
      if self.isInitialized:
        tab = log.GetTable (columns = [8, 8, 8, 8, 10, 10, 10, 10])
        tab.Start ()
        tab.Heading ("SiteID")
        tab.Heading ("Site", columnSpan = 3)
        tab.Heading ("Instances")
        tab.Heading ("Center", columnSpan = 3)

        for site in self.meadSites:
          tab.Entry ("%d" % site.siteID)
          tab.Entry (site.segName)
          tab.Entry (site.resName)
          tab.Entry ("%d" % site.resSerial)
          tab.Entry ("%d" % len (site.instances))
          tab.Entry ("%10.3f" % site.center[0])
          tab.Entry ("%10.3f" % site.center[1])
          tab.Entry ("%10.3f" % site.center[2])
        tab.Stop ()

 
  def SummaryProbabilities (self, reportOnlyUnusual = False, maxProbThreshold = 0.75, log = logFile):
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
      if self.isCalculated:
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
              maxProb  = instance.probability
              maxID    = instance.instID

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
              if instance.instID == maxID:
                label = "*%s" % instance.label
              else:
                label = instance.label
              tab.Entry ("%8s"   % label)
              tab.Entry ("%8.4f" % instance.probability)
  
            for filler in range (0, maxinstances - len (site.instances)):
              tab.Entry ("")
              tab.Entry ("")
        tab.Stop ()


  def DeleteJobFiles (self):
    """Delete job files."""
    files = []
    files.append (self.backPqr)
    files.append (self.proteinPqr)
    files.append (self.sitesFpt)

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


  def WriteJobFiles (self, system, log = logFile):
    """Write files: PQR, FPT, OGM and MGM."""
    if self.isInitialized:

      # Get atomic charges and radii for the system
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
      if not os.path.exists (self.scratch):
        try:
          os.mkdir (self.scratch)
        except:
          raise ContinuumElectrostaticsError ("Cannot create scratch directory %s" % self.scratch)


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


      # Write two PQR files for each instance of every site, first for the model compund and second for the site itself
      for meadSite in self.meadSites:
        model = Selection (meadSite.modelAtomIndices)
        site  = Selection (meadSite.siteAtomIndices)

        for instance in meadSite.instances:

          # In the PQR file of the model compound, charges of the site atoms must be set to zero (requirement of the my_2diel_solver program)
          chargesUpdated = Clone (systemCharges)
          for atomIndex in meadSite.siteAtomIndices:
            chargesUpdated[atomIndex] = 0.0

          PQRFile_FromSystem (instance.modelPqr, system, selection = model, charges = chargesUpdated, radii = systemRadii)


          chargesUpdated = Clone (systemCharges)
          for chargeIndex, atomIndex in enumerate (meadSite.siteAtomIndices):
            pickCharge                = instance.charges[chargeIndex]
            chargesUpdated[atomIndex] = pickCharge

          PQRFile_FromSystem (instance.sitePqr, system, selection = site,  charges = chargesUpdated, radii = systemRadii)


      # Write background PQR file
      PQRFile_FromSystem (self.backPqr, system, selection = Selection (self.backAtomIndices), charges = systemCharges, radii = systemRadii)

      # Write full-protein PQR file (to be used as eps2set_region)
      PQRFile_FromSystem (self.proteinPqr, system, selection = Selection (self.proteinAtomIndices), charges = systemCharges, radii = systemRadii)

      # Write FPT-file
      lines = []

      for siteIndex, meadSite in enumerate (self.meadSites):
        for instanceIndex, instance in enumerate (meadSite.instances):
          for atomIndex, charge in zip (meadSite.siteAtomIndices, instance.charges):
            x, y, z = system.coordinates3[atomIndex]
            line    = "%d %d %f %f %f %f\n" % (siteIndex, instanceIndex, x, y, z, charge)
            lines.append (line)

      WriteInputFile (self.sitesFpt, lines)


      # Write MGM and OGM files
      filesWritten = []

      for meadSite in self.meadSites:
        lines = []

        for stepIndex, (nodes, resolution) in enumerate (self.focussingSteps):
          if stepIndex < 1:
            lines.append ("ON_GEOM_CENT %d %f\n" % (nodes, resolution))
          else:
            x, y, z = meadSite.center
            lines.append ("(%f %f %f) %d %f\n"% (x, y, z, nodes, resolution))

        for instance in meadSite.instances:
          for fileGrid in (instance.modelGrid, instance.siteGrid):
            WriteInputFile (fileGrid, lines)
            filesWritten.append (fileGrid)

      self.isFilesWritten = True


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
