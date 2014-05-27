#-------------------------------------------------------------------------------
# . File      : Model.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""MEADModel is a class representing the continuum electrostatic model."""

import os, glob, math


from pCore             import logFile, LogFileActive, Selection, Vector3, YAMLUnpickle, Clone, UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE
                       
from Constants         import *
from Error             import ContinuumElectrostaticsError
from InputFileWriter   import WriteInputFile
from PQRFileWriter     import PQRFile_FromSystem
from Site              import MEADSite
from Instance          import MEADInstance, InstanceThread 
from Toolbox           import FormatEntry, ConvertAttribute
from StateVector       import StateVector


_DefaultTemperature     = 300.0

_DefaultIonicStrength   = 0.1

_DefaultMeadPath        = "/usr/bin"

_DefaultScratch         = "/tmp"

_DefaultThreads         = 1

_DefaultCleanUp         = False

_DefaultFocussingSteps  = [(121, 2.00), (101, 1.00), (101, 0.50), (101, 0.25)]


#-------------------------------------------------------------------------------
class MEADModel (object):
  """Continuum electrostatic model."""

  defaultAttributes = {
                 "temperature"        :  _DefaultTemperature    ,
                 "ionicStrength"      :  _DefaultIonicStrength  ,
                 "meadPath"           :  _DefaultMeadPath       ,
                 "nthreads"           :  _DefaultThreads        ,
                 "deleteJobFiles"     :  _DefaultCleanUp        ,
                 "scratch"            :  _DefaultScratch        ,
                 "focussingSteps"     :  Clone (_DefaultFocussingSteps) ,
                 "librarySites"       :  None     ,
                 "meadSites"          :  None     ,
                 "backAtomIndices"    :  None     ,
                 "proteinAtomIndices" :  None     ,
                 "backPqr"            :  None     ,
                 "proteinPqr"         :  None     ,
                 "sitesFpt"           :  None     ,
                 "isInitialized"      :  False    ,
                 "isFilesWritten"     :  False    ,
                 "isCalculated"       :  False    ,
                 "splitToDirectories" :  True     ,
                      }

  defaultAttributeNames = {
                  "Temperature"       : "temperature"        ,
                  "Ionic Strength"    : "ionicStrength"      ,
                  "Threads"           : "nthreads"           ,
                  "Delete Job Files"  : "deleteJobFiles"     ,
                  "Split Directories" : "splitToDirectories" ,
                          }
#                   "Scratch Directory" : "scratch"            ,
#                   "Path to MEAD"      : "meadPath"           ,
#                   "Initialized"       : "isInitialized"      ,
#                   "Files Written"     : "isFilesWritten"     ,
#                   "Calculated"        : "isCalculated"       ,


  def __del__ (self):
    """Deallocation."""
    if self.deleteJobFiles: self.DeleteJobFiles ()


  def __init__ (self, log = logFile, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)

    self.LoadLibraryOfSites (log = log)


  def LoadLibraryOfSites (self, log = logFile):
    """Load a set of YAML files with parameters for titratable sites.

    If there are YAML files in the current directory, they are loaded as well.

    If these additional files have names coinciding with the names from the library, the library parameters will be overwritten."""
    filesLibrary = glob.glob (os.path.join (YAMLPATHIN, "sites/", "*.yaml"))
    filesExtra   = glob.glob (os.path.join (os.getcwd (), "*.yaml"))
    filesLibrary.extend (filesExtra)

    if LogFileActive (log):
      for fileExtra in filesExtra:
        log.Text ("\nIncluded custom file: %s\n" % os.path.basename (fileExtra))

    self.librarySites = {}

    for fileSite in filesLibrary:
      site      = YAMLUnpickle (fileSite)
      name      = site [ "site"      ]
      atoms     = site [ "atoms"     ]
      instances = site [ "instances" ]
      self.librarySites[name] = {"atoms" : atoms, "instances" : instances}


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


  def CalculateCurvesGMCT (self, resolution = 0.5, log = logFile):
    """Calculate titration curves using GMCT."""
    pass


  def CalculateCurvesAnalytically (self, resolution = 0.5, log = logFile):
    """Calculate titration curves analytically."""
    if self.isCalculated:
      nsteps = int (14.0 / resolution + 1)
      sites  = []
      for site in self.meadSites:
        instances = []
        for instance in site.instances:
          instances.append ([None] * nsteps)
        sites.append (instances)
  
      # Go over pH values 0...14
      pH = 0.
      for step in range (0, nsteps):
        self.CalculateProbabilitiesAnalytically (pH = pH, log = None)
  
        for siteIndex, site in enumerate (self.meadSites):
          for instIndex, instance in enumerate (site.instances):
             sites[siteIndex][instIndex][step] = instance.probability
  
        if LogFileActive (log):
          log.Text ("Calculated pH = %.2f\n" % pH)
        pH = pH + resolution
  
      # Write results to files
      directory = os.path.join (self.scratch, "curves")
      if not os.path.exists (directory): os.mkdir (directory)
  
      for siteIndex, site in enumerate (self.meadSites):
        for instIndex, instance in enumerate (site.instances):
          probabilities = sites[siteIndex][instIndex]
          lines         = []
          for step in range (0, nsteps):
            lines.append ("%f %f\n" % (step * resolution, probabilities[step]))

          filename = os.path.join (directory, "%s_%s.dat" % (site.label, instance.label))
          WriteInputFile (filename, lines)


  def CalculateProbabilitiesGMCT (self, pH = 7.0, log = logFile):
    """Use GMCT to estimate probabilities."""
    pass


  def CalculateProbabilitiesAnalytically (self, pH = 7.0, log = logFile):
    """For each site, calculate the probability of occurance of each instance, using the Boltzmann weighted sum."""
    if self.isCalculated:
      if len (self.meadSites) > MAX_SITES:
        raise ContinuumElectrostaticsError ("Too many sites for analytic treatment (%d)\n" % nsites)
  
      # Calculate all state energies
      stateVector   = StateVector (self)
      stateEnergies = []
      increment     = True
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
 
      # Reset the probabilities 
      for site in self.meadSites:
        for instance in site.instances:
          instance.probability = 0.
 
      # Calculate the probabilities 
      increment  = True
      stateIndex = 0
      stateVector.Reset ()
  
      while increment:
        for siteIndex, site in enumerate (self.meadSites):
          instanceIndex        = stateVector[siteIndex]
          instance             = site.instances[instanceIndex]
          instance.probability = instance.probability + bfactors[stateIndex]
        increment  = stateVector.Increment ()
        stateIndex = stateIndex + 1
  
      bsum = 1.0 / sum (bfactors)
      for site in self.meadSites:
        for instance in site.instances:
          instance.probability = instance.probability * bsum
  
      if LogFileActive (log):
        log.Text ("\nCalculating protonation probabilities complete.\n")


  def CalculateMicrostateEnergy (self, stateVector, pH = 7.0):
    """Calculate energy of a protonation state (=microstate).

    The protonation state is defined by a state vector. 

    The energy is calculated at a given pH."""
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

      protonChemicalPotential = -CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10 * pH
      Gmicro = totalGintr - nprotons * protonChemicalPotential + totalInteract

# Slower but more accurate?
#        for siteIndexInner in range (0, nsites):
#          instanceIndexInner = stateVector  [siteIndexInner]
#          interaction        = interactions [siteIndexInner]
#          totalInteract      = totalInteract + interaction [instanceIndexInner]
#
#      protonChemicalPotential = -CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10 * pH
#      Gmicro = totalGintr - nprotons * protonChemicalPotential + 0.5 * totalInteract
    else:
      Gmicro = None

    return Gmicro


  def CalculateEnergies (self, log = logFile):
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
      table = None

      if LogFileActive (log):
        if self.nthreads < 2:
          log.Text ("\nStarting serial run.\n")
        else:
          log.Text ("\nStarting parallel run on %d CPUs.\n" % self.nthreads)

        table = log.GetTable (columns = [6, 6, 6, 6, 16, 16, 16, 16, 16, 16])
        table.Start ()
        table.Heading ("Instance of a site", columnSpan = 4)
        table.Heading ("Gborn_model"  )
        table.Heading ("Gback_model"  )
        table.Heading ("Gborn_protein")
        table.Heading ("Gback_protein")
        table.Heading ("Gmodel"       )
        table.Heading ("Gintr"        )

      if self.nthreads < 2:
        for meadSite in self.meadSites:
          for instance in meadSite.instances:
            instance.CalculateSiteInModelCompound (log)
            instance.CalculateSiteInProtein (log)
            instance.CalculateGintr (log)

            instance.TableEntry (table)
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
          for thread in batch:
            thread.start ()
          for thread in batch:
            thread.join ()

          # Print the results at the end of each batch, otherwise they come in random order
          for thread in batch:
            instance = thread.instance
            instance.TableEntry (table)

      if table:
        table.Stop ()
        log.Text ("\nCalculating electrostatic energies complete.\n")

      self.isCalculated = True


  def Initialize (self, system, excludeSegments = None, excludeResidues = None, log = logFile):
    """Decompose the system into model compounds, sites and a background charge set."""

    # Check for the CHARMM energy model
    if system.energyModel.mmModel.label is not "CHARMM":
      raise ContinuumElectrostaticsError ("The energy model of the system is different from CHARMM.")


    if not self.isInitialized:
      ParseLabel = system.sequence.ParseLabel
      segments   = system.sequence.children

      instIndexGlobal = 0
      siteIndex       = 0
      self.meadSites  = []

      if excludeSegments is None:
        excludeSegments = ["WATA", ]

      if excludeResidues is None:
        excludeResidues = []  # For example ["ARG", "LYS", ]


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

            # Include residue?
            if residueName not in excludeResidues:
 
              # Titratable residue? 
              if residueName in self.librarySites:
  
                prevIndices = []
                nextIndices = []
  
                # Include atoms from the previous residue to the model compound?
                if residueIndex > 1:
                  prevResidue = residues[residueIndex - 1]
                  prevResidueName, prevResidueSerial = ParseLabel (prevResidue.label, fields = 2)
  
                  if prevResidueName in PROTEIN_RESIDUES:
                    prevNames = PREV_RESIDUE
  
                    for atom in prevResidue.children:
                      if atom.label in prevNames:
                        prevIndices.append (atom.index)
  
                # Include atoms from the next residue to the model compound?
                if residueIndex < (nresidues - 1):
                  nextResidue = residues[residueIndex + 1]
                  nextResidueName, nextResidueSerial = ParseLabel (nextResidue.label, fields = 2)
  
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
                    modelPqr  = os.path.join (self.scratch, segmentName, "%s%s" % (residueName, residueSerial), "%s_%s.%s" % ("model", label, "pqr"))
                    modelLog  = os.path.join (self.scratch, segmentName, "%s%s" % (residueName, residueSerial), "%s_%s.%s" % ("model", label, "out"))
                    modelGrid = os.path.join (self.scratch, segmentName, "%s%s" % (residueName, residueSerial), "%s_%s.%s" % ("model", label, "mgm"))
                    sitePqr   = os.path.join (self.scratch, segmentName, "%s%s" % (residueName, residueSerial), "%s_%s.%s" % ("site",  label, "pqr"))
                    siteLog   = os.path.join (self.scratch, segmentName, "%s%s" % (residueName, residueSerial), "%s_%s.%s" % ("site",  label, "out"))
                    siteGrid  = os.path.join (self.scratch, segmentName, "%s%s" % (residueName, residueSerial), "%s_%s.%s" % ("site",  label, "ogm"))
                  else:
                    modelPqr  = os.path.join (self.scratch, "%s_%s_%s_%s_%s.%s" % ("model", segmentName, residueName, residueSerial, label, "pqr"))
                    modelLog  = os.path.join (self.scratch, "%s_%s_%s_%s_%s.%s" % ("model", segmentName, residueName, residueSerial, label, "out"))
                    modelGrid = os.path.join (self.scratch, "%s_%s_%s_%s_%s.%s" % ("model", segmentName, residueName, residueSerial, label, "mgm"))
                    sitePqr   = os.path.join (self.scratch, "%s_%s_%s_%s_%s.%s" % ("site",  segmentName, residueName, residueSerial, label, "pqr"))
                    siteLog   = os.path.join (self.scratch, "%s_%s_%s_%s_%s.%s" % ("site",  segmentName, residueName, residueSerial, label, "out"))
                    siteGrid  = os.path.join (self.scratch, "%s_%s_%s_%s_%s.%s" % ("site",  segmentName, residueName, residueSerial, label, "ogm"))

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


                # Calculate the center of geometry
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
        table = log.GetTable (columns = [8, 8, 8, 8, 10, 10, 10, 10])
        table.Start ()
        table.Heading ("SiteID")
        table.Heading ("Site", columnSpan = 3)
        table.Heading ("Instances")
        table.Heading ("Center", columnSpan = 3)

        for site in self.meadSites:
          table.Entry ("%d" % site.siteID)
          table.Entry (site.segName)
          table.Entry (site.resName)
          table.Entry (site.resSerial)
          table.Entry ("%d" % len (site.instances))
          table.Entry ("%10.3f" % site.center[0])
          table.Entry ("%10.3f" % site.center[1])
          table.Entry ("%10.3f" % site.center[2])
        table.Stop ()

 
  def SummaryProbabilities (self, log = logFile):
    """List probabilities of occurance of instances."""
    if LogFileActive (log):
      if self.isCalculated:
        
        maxinstances = 0
        for site in self.meadSites:
          ninstances = len (site.instances)
          if ninstances > maxinstances: maxinstances = ninstances

        table = log.GetTable (columns = [6, 6, 6] + [8, 8] * maxinstances)
        table.Start ()
        table.Heading ("Site", columnSpan = 3)
        table.Heading ("Probabilities of instances", columnSpan = maxinstances * 2)

        for site in self.meadSites:
          table.Entry ("%6s" % site.segName)
          table.Entry ("%6s" % site.resName)
          table.Entry ("%6s" % site.resSerial)

          maxProb = 0.
          for instance in site.instances:
            if instance.probability > maxProb: 
              maxProb = instance.probability
              maxID   = instance.instID

          for instance in site.instances:
            if instance.instID == maxID:
              label = "*%s" % instance.label
            else:
              label = instance.label
            table.Entry ("%8s"   % label)
            table.Entry ("%8.4f" % instance.probability)

          for filler in range (0, maxinstances - len (site.instances)):
            table.Entry ("")
            table.Entry ("")
        table.Stop ()


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
