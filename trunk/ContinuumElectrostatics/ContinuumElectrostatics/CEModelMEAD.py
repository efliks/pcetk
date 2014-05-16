#-------------------------------------------------------------------------------
# . File      : CEModelMEAD.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------

# FIXME: Move parts of WriteJobFiles to the instance class
# 
# FIXME: Rename variables containing filenames so that they start from "file"
# 
# FIXME: Coordinates from the FPT file should have their own data structure?
# 
# FIXME: Efficiency improvements during partitioning of the system and writing job files?
# 
# FIXME: Drop _Residue class, use components of the system instead (children)
# 
# FIXME: Use arrays instead of lists for interactions (Real1DArray or SymmetricMatrix?)
# 
# FIXME: Have a column of ETA (Estimated Time for Accomplishment) in MEAD calculations
# 
# FIXME: Optionally convert kcal/mol (MEAD units) to kJ/mol (pDynamo units)


import exceptions, os, subprocess, threading, glob


from pCore                  import logFile, LogFileActive, Selection, Vector3, YAMLUnpickle, Clone, UNITS_ENERGY_KILOCALORIES_PER_MOLE_TO_KILOJOULES_PER_MOLE
from pMolecule              import System
                            
from Constants              import *
from Tools                  import _Residue, _FormatEntry, _ConvertAttribute
from StateVector            import StateVector
from InputFileWriter        import WriteInputFile
from PQRFileWriter          import PQRFile_FromSystem
from MEADOutputFileReader   import MEADOutputFileReader


_DefaultTemperature     = 300.0

_DefaultIonicStrength   = 0.1

_DefaultMeadPath        = "/usr/bin"

_DefaultScratch         = "/tmp"

_DefaultThreads         = 1

_DefaultCleanUp         = False

_DefaultFocussingSteps  = [(121, 2.00), (101, 1.00), (101, 0.50), (101, 0.25)]


#-------------------------------------------------------------------------------
class CEModelMEADError (exceptions.StandardError):
  """A class for handling errors in CEModelMEAD."""
  pass


#-------------------------------------------------------------------------------
class _MEADThread (threading.Thread):
  """A class for running parallel calculations for sites both in model compounds and in the protein.

  It also calculates Gintr."""

  def __init__ (self, instance = None, log = logFile):
    """Constructor."""
    threading.Thread.__init__ (self)
    self.instance = instance
    self.log      = log


  def run (self):
    """The method that runs the calculations."""
    self.instance.CalculateSiteInModelCompound (log = self.log)
    self.instance.CalculateSiteInProtein       (log = self.log)
    self.instance.CalculateGintr               (log = self.log)


#   def run (self):
#     """The method that runs the calculations."""
#     if not self.complete:
#       tstart   = time.time ()
#       self.instance.CalculateSiteInModelCompound (log = self.log)
#
#       tmodel   = time.time ()
#       self.instance.CalculateSiteInProtein (log = self.log)
#
#       tprotein = time.time ()
#       self.instance.CalculateGintr (log = self.log)
#
#       tstop    = time.time ()
#       self.tmodel   = tmodel   - tstart
#       self.tprotein = tprotein - tmodel
#       self.ttotal   = tstop    - tstart
#       self.complete = True


#-------------------------------------------------------------------------------
class MEADInstance (object):
  """Instance of a titratable site.

  For the time being, the instances only differ in charges. There are no rotameric instances."""

  defaultAttributes = {
                  "parent"        : None , # <--This should point to the instance's site
                  "instID"        : None ,
                  "label"         : None ,
                  "protons"       : 0    ,
                  "charges"       : None ,
                  "modelPqr"      : None ,
                  "modelLog"      : None ,
                  "modelGrid"     : None ,
                  "sitePqr"       : None ,
                  "siteLog"       : None ,
                  "siteGrid"      : None ,
                  "Gmodel"        : None ,
                  "Gintr"         : None ,
                  "Gborn_model"   : None ,
                  "Gback_model"   : None ,
                  "Gborn_protein" : None ,
                  "Gback_protein" : None ,
                  "interactions"  : None ,
                      }

  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


  def CalculateSiteInModelCompound (self, log = logFile):
    """Calculate Gborn and Gback of a site in the model compound."""
    if os.path.exists (self.modelLog):
      pass
    else:
      instancePqr        = self.sitePqr  [:-4]
      modelBackgroundPqr = self.modelPqr [:-4]

      model   = self.parent.parent
      command = ["%s/my_2diel_solver" % model.meadPath, "-T", "%f" % model.temperature, "-ionicstr", "%f" % model.ionicStrength, "-epsin", "%f" % 4.0, instancePqr, modelBackgroundPqr]

      try:
        outFile = open (self.modelLog, "w")
        subprocess.check_call (command, stderr = outFile, stdout = outFile)

        outFile.close ()
      except:
        report = " ".join (command)
        raise CEModelMEADError ("Failed running command: %s" % report)

    reader = MEADOutputFileReader (self.modelLog)
    reader.Parse ()

    self.Gborn_model = reader.born
    self.Gback_model = reader.back


  def CalculateSiteInProtein (self, log = logFile):
    """Calculate Gborn and Gback of a site in protein environment.

    Also, calculate Wij between the site and the other sites.

    Use fpt-file for other instances of sites."""
    if os.path.exists (self.siteLog):
      pass
    else:
      model = self.parent.parent
      sitesFpt             = model.sitesFpt   [:-4]
      proteinPqr           = model.proteinPqr [:-4]
      proteinBackgroundPqr = model.backPqr    [:-4]
      instancePqr          = self.sitePqr     [:-4]

      # epsin1 is never used but must be given

      # eps2set defines the whole protein

      command = ["%s/my_3diel_solver" % model.meadPath, "-T", "%f" % model.temperature, "-ionicstr", "%f" % model.ionicStrength, "-epsin1", "%f" % 1.0, "-epsin2", "%f" % 4.0, "-eps2set", "%s" % proteinPqr, "-fpt", "%s" % sitesFpt, instancePqr, proteinBackgroundPqr]

      try:
        outFile = open (self.siteLog, "w")
        subprocess.check_call (command, stderr = outFile, stdout = outFile)

        outFile.close ()
      except:
        report = " ".join (command)
        raise CEModelMEADError ("Failed running command: %s" % report)

    reader = MEADOutputFileReader (self.siteLog)
    reader.Parse ()

    self.Gborn_protein = reader.born
    self.Gback_protein = reader.back

    sites         = []
    instances     = []
    prevSiteIndex = -1

    for siteIndex, instanceIndex, energy in reader.interactions:
      if siteIndex == (self.parent.siteID - 1):
        Wij = 0.0
      else:
        Wij = energy

      if siteIndex is not prevSiteIndex and siteIndex > 0:
        sites.append (instances)
        instances = []
      instances.append (Wij)
      prevSiteIndex = siteIndex

    if instances:
      sites.append (instances)

# Return it to the calling method? The matrix of interactions should be inside the CEModelMEAD?
    self.interactions = sites


  def CalculateGintr (self, log = logFile):
    """Calculate Gintr of an instance of a site in the protein."""
    self.Gintr = self.Gmodel + (self.Gborn_protein - self.Gborn_model) + (self.Gback_protein - self.Gback_model)


  def PrintInteractions (self, sort = False, log = logFile):
    """Print interactions of the instance of a site with other instances of other sites."""
    if LogFileActive (log):
      site  = self.parent
      model = site.parent

      if model.isCalculated:
        instances = []

        for site in model.meadSites:
          for instance in site.instances:
            s   = site.siteID     - 1
            i   = instance.instID - 1
            Wij = self.interactions[s][i]
            instances.append ([Wij, site.segName, site.resName, site.resNum, instance.label])
        if sort: instances.sort ()

        table = log.GetTable (columns = [6, 6, 6, 6, 16])
        table.Start ()
        table.Heading ("Instance of a site", columnSpan = 4)
        table.Heading ("Wij")

        for Wij, segName, resName, resNum, label in instances:
          table.Entry (segName)
          table.Entry (resName)
          table.Entry (resNum)
          table.Entry (label)
          table.Entry ("%16.4f" % Wij)
        table.Stop ()


  def TableEntry (self, table = None):
    """Report calculated energies in a table."""
    if table:
      site = self.parent
      table.Entry (site.segName)
      table.Entry (site.resName)
      table.Entry (site.resNum)
      table.Entry (self.label)
      table.Entry ("%16.4f" % self.Gborn_model)
      table.Entry ("%16.4f" % self.Gback_model)
      table.Entry ("%16.4f" % self.Gborn_protein)
      table.Entry ("%16.4f" % self.Gback_protein)
      table.Entry ("%16.4f" % self.Gmodel)
      table.Entry ("%16.4f" % self.Gintr)


#-------------------------------------------------------------------------------
class MEADSite (object):
  """Titratable site.

  Each site has at least two instances."""

# FIXME: Atom names are not necessary? Or also add coordinates and radii (modelAtomRadii, siteAtomRadii)
  defaultAttributes = {
                  "parent"           : None , # <--This should point to the MEAD model
                  "siteID"           : None ,
                  "segName"          : None ,
                  "resName"          : None ,
                  "resNum"           : None ,
                  "instances"        : None ,
                  "center"           : None ,
                  "modelAtomNames"   : None ,
                  "modelAtomIndices" : None ,
                  "siteAtomNames"    : None ,
                  "siteAtomIndices"  : None ,
                  "label"            : None ,
                      }

  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)

    self.label = "%s_%s%s" % (self.segName, self.resName, self.resNum)


#-------------------------------------------------------------------------------
class CEModelMEAD (object):
  """Continuum electrostatic model."""

  defaultAttributes = {
                 "temperature"        :  _DefaultTemperature      ,
                 "ionicStrength"      :  _DefaultIonicStrength    ,
                 "meadPath"           :  _DefaultMeadPath         ,
                 "nthreads"           :  _DefaultThreads          ,
                 "deleteJobFiles"     :  _DefaultCleanUp          ,
                 "scratch"            :  _DefaultScratch          ,
                 "focussingSteps"     :  _DefaultFocussingSteps   ,
                 "librarySites"       :  {}                       ,
                 "meadSites"          :  []                       ,
                 "backAtomIndices"    :  None                     ,
                 "proteinAtomIndices" :  None                     ,
                 "backPqr"            :  None                     ,
                 "proteinPqr"         :  None                     ,
                 "sitesFpt"           :  None                     ,
                 "system"             :  None                     ,
                 "isInitialized"      :  False                    ,
                 "isFilesWritten"     :  False                    ,
                 "isCalculated"       :  False                    ,
                 "splitToDirectories" :  True                     ,
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


  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


    if not isinstance (self.system, System):
      raise CEModelMEADError ("The system argument is mandatory.")

    if not self.system.energyModel.mmModel.label is "CHARMM":
      raise CEModelMEADError ("The energy model of the system is different from CHARMM.")

    if not os.path.exists (self.scratch):
      try:
        os.mkdir (self.scratch)
      except:
        raise CEModelMEADError ("Cannot create scratch directory %s" % self.scratch)

    # Load the library of sites
    search     = "%s/%s" % (YAMLPATHIN, "sites/")
    sitefiles  = glob.glob ("%s/*.yaml" % search)

    # Use custom files if they are present
    extrafiles = glob.glob ("*.yaml")
    if extrafiles:
      sitefiles.extend (extrafiles)

    for sitefile in sitefiles:
      site = YAMLUnpickle (sitefile)
      name = site["site"]
      self.librarySites[name] = {"atoms" : site["atoms"], "instances" : site["instances"]}


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
         ("Wij_symm" , (16,   10) ),
         ("Wij"      , (16,   10) ),
         ("Wij_err"  , (16,   10) ),
              )
      header = _FormatEntry (items, header = True)
      entry  = _FormatEntry (items)
      lines  = [header]

      for site1 in self.meadSites:
        for instance1 in site1.instances:
          for site2 in self.meadSites:
            for instance2 in site2.instances:
              s   = site2.siteID     - 1
              i   = instance2.instID - 1
              Wij = instance1.interactions[s][i]

              s   = site1.siteID     - 1
              i   = instance1.instID - 1
              Wji = instance2.interactions[s][i]

              Wij_symm = 0.5 * (Wij + Wji)
              Wij_err  = Wij_symm - Wij

              line = entry % (
                      site1.siteID, instance1.instID, site1.label, instance1.label,
                      site2.siteID, instance2.instID, site2.label, instance2.label,
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
         ("Gintr"     , (16,  10) ),
         ("protons"   , (12,  0) ),
              )
      header = _FormatEntry (items, header = True)
      entry  = _FormatEntry (items)
      lines  = [header]

      for site in self.meadSites:
        for instance in site.instances:
          line          = entry % (site.siteID, instance.instID, site.label, instance.label, instance.Gintr, instance.protons)
          lines.append (line)
      WriteInputFile (filename, lines)


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

        #for siteIndexInner in range (0, nsites):
        for siteIndexInner in range (0, siteIndex):
          instanceIndexInner = stateVector  [siteIndexInner]
          interaction        = interactions [siteIndexInner]
          totalInteract      = totalInteract + interaction [instanceIndexInner]

      protonChemicalPotential = -CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10 * pH

      #Gmicro = totalGintr - nprotons * protonChemicalPotential + 0.5 * totalInteract
      Gmicro = totalGintr - nprotons * protonChemicalPotential + totalInteract
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

    Finally, calculate Gintr from Gmodel and transfer energies.
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
            thread = _MEADThread (instance, log)
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


  def Initialize_Testing (self, excludeSegments = None, excludeResidues = None, log = logFile):
    """Decompose the system into model compounds, sites and background charge set."""

    if not self.isInitialized:
      system      = self.system
      ParseLabel  = system.sequence.ParseLabel
      segments    = system.sequence.children
      siteIndex   = 0

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
                atomIndicesModel = prevIndices
                atomIndicesSite  = []
  
                for atom in atoms:
                  atomIndex   = atom.index
                  atomName    = atom.label
  
                  if atomName in libSiteAtoms:
                    atomIndicesSite.append (atomIndex)

                  atomIndicesModel.append (atomIndex)
                atomIndicesModel.extend (nextIndices)


                # Create instances
                libSiteInstances = libSite["instances"]
                instances        = []
  
                for instIndex, instance in enumerate (libSiteInstances):
                  label    = instance [ "label"   ]
                  protons  = instance [ "protons" ]
                  charges  = instance [ "charges" ]
                  Gmodel   = instance [ "Gmodel"  ] * 300.0 / self.temperature
                  pattern  = "%s/%s_%s.%s"

                  if self.splitToDirectories:
                    foo = "%s/%s/%s%s" % (self.scratch, segmentName, residueName, residueSerial)
                    bar = label
                  else:
                    foo = self.scratch
                    bar = "%s_%s_%s_%s" % (segmentName, residueName, residueSerial, label)

                  modelPqr  = pattern % (foo, "model", bar, "pqr") 
                  modelLog  = pattern % (foo, "model", bar, "log") 
                  modelGrid = pattern % (foo, "model", bar, "mgm") 
                  sitePqr   = pattern % (foo, "site" , bar, "pqr") 
                  siteLog   = pattern % (foo, "site" , bar, "log") 
                  siteGrid  = pattern % (foo, "site" , bar, "ogm") 

                  # Note: set the instance's parent later
                  newInstance = MEADInstance (
                                 instID     =  instIndex + 1 ,
                                 label      =  label         ,
                                 protons    =  protons       ,
                                 charges    =  charges       ,
                                 Gmodel     =  Gmodel        ,
                                 modelPqr   =  modelPqr      ,
                                 modelLog   =  modelLog      ,
                                 modelGrid  =  modelGrid     ,
                                 sitePqr    =  sitePqr       ,
                                 siteLog    =  siteLog       ,
                                 siteGrid   =  siteGrid      ,
                                             )
                  instances.append (newInstance)


                # Calculate the center of geometry
                center = Vector3 ()
                for atomIndex in atomIndicesSite:
                  center.AddScaledVector3 (1.0, system.coordinates3[atomIndex])

                center.Scale (1.0 / len (atomIndicesSite))

                # Create a site
                newSite = MEADSite (
                               parent           = self             ,
                               siteID           = siteIndex + 1    ,
                               segName          = segmentName      ,
                               resName          = residueName      ,
                               resNum           = residueSerial    ,
                               modelAtomIndices = atomIndicesModel ,
                               siteAtomIndices  = atomIndicesSite  ,
                               instances        = instances        ,
                               center           = center           ,
                                   )
                for instance in newSite.instances:
                  instance.parent = newSite    # <--Setting parents
                
                self.meadSites.append (newSite)
                siteIndex = siteIndex + 1


      # Construct the background set of charges
      allSiteAtomIndices = []
      for site in self.meadSites:
        allSiteAtomIndices.extend (site.siteAtomIndices)

      backAtomIndices = []

      #============ Go over segments ============
      for segment in segments:
        residues = segment.children

        #============ Go over residues ============
        for residue in residues:
          atoms = residue.children

          #============ Go over atoms ============
#          backAtomIndices.extend (filter (lambda atom: atom.index not in allSiteAtomIndices, atoms))
          for atom in atoms:
            if atom.index not in allSiteAtomIndices:
              backAtomIndices.append (atom.index)

      self.backAtomIndices = backAtomIndices
      self.backPqr         = "%s/back.pqr" % self.scratch


      # Construct the protein (this means removing residues not defined in PROTEIN_RESIDUES, usually waters and ions)
      proteinAtomIndices = []

      #============ Go over segments ============
      for segment in segments:
        residues = segment.children

        #============ Go over residues ============
        for residue in residues:
          residueName, residueSerial = ParseLabel (residue.label, fields = 2)

          if residueName not in REMOVE_RESIDUES:
            atoms = residue.children
  
            #============ Go over atoms ============
#            proteinAtomIndices.extend (map (lambda atom: atom.index, atoms))
            for atom in atoms:
              proteinAtomIndices.append (atom.index)
      self.proteinAtomIndices = proteinAtomIndices


      # Define full-protein PQR file (to be used as eps2set_region)
      self.proteinPqr = "%s/protein.pqr" % self.scratch

      # Define FPT-file
      self.sitesFpt = "%s/site.fpt" % self.scratch

      self.isInitialized = True


  def Initialize (self, log = logFile):
    """Decompose the system into model compounds, sites and background charge set.

    Create the necessary data structures inside the model object.

    FIXME: Rewrite this method to use system.sequence.children and eliminate the _Residue class.
    """
    if not self.isInitialized:
      residues    = []
      atomNames   = []
      atomIndices = []

      residue     = ""
      prevResidue = ""

      for atomIndex, atom in enumerate (self.system.atoms):
        segName         = atom.parent.parent.label
        residue         = atom.parent.label
        resName, resNum = self.system.sequence.ParseLabel (residue, fields = 2)
        atomName        = atom.label

        if residue != prevResidue and prevResidue != "":
          segPrevName, resPrevName, resPrevNum = prevPath

          if resPrevName in PROTEIN_RESIDUES:
            r = _Residue (
                     segName     = segPrevName ,
                     resName     = resPrevName ,
                     resNum      = resPrevNum  ,
                     atomNames   = atomNames   ,
                     atomIndices = atomIndices
                         )
            residues.append (r)

          atomNames   = []
          atomIndices = []

        atomNames.append (atomName)
        atomIndices.append (atomIndex)

        prevResidue = residue
        prevPath    = [segName, resName, resNum]

      if resName in PROTEIN_RESIDUES:
        if atomNames:
          r = _Residue (
                   segName     = segPrevName ,
                   resName     = resPrevName ,
                   resNum      = resPrevNum  ,
                   atomNames   = atomNames   ,
                   atomIndices = atomIndices
                       )
          residues.append (r)


      # Out of all residues, select titratable residues and extend them by adding some atoms from neighbouring residues to create model compounds
      totalResidues = len (residues)
      siteID = 1

      for resIndex, cr in enumerate (residues):
        modelAtomIndices = []
        modelAtomNames   = []

        if cr.resName in self.librarySites:
          if resIndex > 1:
            pr = residues[resIndex - 1]

            for prAtomName, prAtomIndex in zip (pr.atomNames, pr.atomIndices):
              if prAtomName in PREV_RESIDUE:
                modelAtomIndices.append (prAtomIndex)
                modelAtomNames.append (prAtomName)


          modelAtomIndices.extend (cr.atomIndices)
          modelAtomNames.extend (cr.atomNames)

          if resIndex < (totalResidues - 2):
            nr = residues[resIndex + 1]

            if   nr.resName == "PRO":
              nextResidueAtomNames = NEXT_RESIDUE_PRO
            elif nr.resName == "GLY":
              nextResidueAtomNames = NEXT_RESIDUE_GLY
            else:
              nextResidueAtomNames = NEXT_RESIDUE

            for nrAtomName, nrAtomIndex in zip (nr.atomNames, nr.atomIndices):
              if nrAtomName in nextResidueAtomNames:
                modelAtomIndices.append (nrAtomIndex)
                modelAtomNames.append (nrAtomName)


          site            = self.librarySites[cr.resName]
          siteAtoms       = site["atoms"]
          siteInstances   = site["instances"]

          siteAtomIndices = []
          siteAtomNames   = []

          for atomName in siteAtoms:
            if atomName in cr.atomNames:
              siteAtomNames.append (atomName)

              atomIndex = cr.atomIndices[cr.atomNames.index (atomName)]
              siteAtomIndices.append (atomIndex)


          if self.splitToDirectories:
            direc = "%s/%s/%s%s/" % (self.scratch, cr.segName, cr.resName, cr.resNum)
            if not os.path.exists (direc):
              try:
                os.makedirs (direc)
              except:
                raise CEModelMEADError ("Cannot create directory: %s" % direc)

          instances = []

          for instID, instance in enumerate (siteInstances, 1):
            label      = instance["label"]
            protons    = instance["protons"]
            charges    = instance["charges"]

            Gmodel300K = instance["Gmodel"]
            Gmodel     = Gmodel300K * 300.0 / self.temperature

            if self.splitToDirectories:
              modelPqr  = "%s/model_%s.pqr" % (direc, label)
              modelLog  = "%s/model_%s.out" % (direc, label)
              modelGrid = "%s/model_%s.mgm" % (direc, label)
              sitePqr   = "%s/site_%s.pqr"  % (direc, label)
              siteLog   = "%s/site_%s.out"  % (direc, label)
              siteGrid  = "%s/site_%s.ogm"  % (direc, label)
            else:
              modelPqr  = "%s/model_%s_%s_%s_%s.pqr" % (self.scratch, cr.segName, cr.resName, cr.resNum, label)
              modelLog  = "%s/model_%s_%s_%s_%s.out" % (self.scratch, cr.segName, cr.resName, cr.resNum, label)
              modelGrid = "%s/model_%s_%s_%s_%s.mgm" % (self.scratch, cr.segName, cr.resName, cr.resNum, label)
              sitePqr   = "%s/site_%s_%s_%s_%s.pqr"  % (self.scratch, cr.segName, cr.resName, cr.resNum, label)
              siteLog   = "%s/site_%s_%s_%s_%s.out"  % (self.scratch, cr.segName, cr.resName, cr.resNum, label)
              siteGrid  = "%s/site_%s_%s_%s_%s.ogm"  % (self.scratch, cr.segName, cr.resName, cr.resNum, label)


            newInstance = MEADInstance (
                                    parent    = None       , # <--Do not forget to set this later
                                    instID    = instID     ,
                                    Gmodel    = Gmodel     ,
                                    label     = label      ,
                                    protons   = protons    ,
                                    charges   = charges    ,
                                    modelPqr  = modelPqr   ,
                                    modelLog  = modelLog   ,
                                    modelGrid = modelGrid  ,
                                    sitePqr   = sitePqr    ,
                                    siteLog   = siteLog    ,
                                    siteGrid  = siteGrid   ,
                                       )
            instances.append (newInstance)

          # For each site, calculate the center of geometry
          center = Vector3 ()
          for atomIndex in siteAtomIndices:
            center.AddScaledVector3 (1.0, self.system.coordinates3[atomIndex])

          center.Scale (1.0 / len (siteAtomIndices))


          newMeadSite = MEADSite (
                             parent           = self             ,
                             siteID           = siteID           ,
                             segName          = cr.segName       ,
                             resName          = cr.resName       ,
                             resNum           = cr.resNum        ,
                             modelAtomNames   = modelAtomNames   ,
                             modelAtomIndices = modelAtomIndices ,
                             siteAtomNames    = siteAtomNames    ,
                             siteAtomIndices  = siteAtomIndices  ,
                             instances        = instances        ,
                             center           = center           ,
                                 )
          for instance in newMeadSite.instances:
            instance.parent = newMeadSite    # <--Setting parents

          self.meadSites.append (newMeadSite)
          siteID += 1


      # Construct background set of charges
      allSiteAtomIndices = []

      for meadSite in self.meadSites:
        allSiteAtomIndices.extend (meadSite.siteAtomIndices)

      backAtomIndices = []

      for residue in residues:
        for atomIndex in residue.atomIndices:
          if not atomIndex in allSiteAtomIndices:
            backAtomIndices.append (atomIndex)

      self.backAtomIndices = backAtomIndices
      self.backPqr         = "%s/back.pqr" % self.scratch


      # Construct the protein (this means removing residues not defined in PROTEIN_RESIDUES, usually waters and ions)
      proteinAtomIndices = []

      for residue in residues:
        proteinAtomIndices.extend (residue.atomIndices)
      self.proteinAtomIndices = proteinAtomIndices


      # Define full-protein PQR file (to be used as eps2set_region)
      self.proteinPqr = "%s/protein.pqr" % self.scratch

      # Define FPT-file
      self.sitesFpt = "%s/site.fpt" % self.scratch

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
        attrConv = _ConvertAttribute (attr)
        summary.Entry (key, attrConv)

      nsites     = len (self.meadSites)
      ninstances = sum (map (lambda site: len (site.instances), self.meadSites))
#      ninstances = 0
#      for site in self.meadSites:
#        ninstances += len (site.instances)

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
          table.Entry (site.resNum)
          table.Entry ("%d" % len (site.instances))
          table.Entry ("%10.3f" % site.center[0])
          table.Entry ("%10.3f" % site.center[1])
          table.Entry ("%10.3f" % site.center[2])
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


  def WriteJobFiles (self, log = logFile):
    """Write files: PQR and FPT."""
    if self.isInitialized:

      # Get atomic charges and radii for the system
      systemCharges = self.system.AtomicCharges ()
      systemRadii   = []

      systemTypes   = self.system.energyModel.mmAtoms.AtomTypes ()
      radii         = YAMLUnpickle ("%s/%s" % (YAMLPATHIN, "radii.yaml"))

      for atomType in systemTypes:
        if radii.has_key (atomType):
          radius = radii[atomType]
        else:
          generalAtomType = "%s*" % atomType[0]

          if radii.has_key (generalAtomType):
            radius = radii[generalAtomType]
          else:
            raise CEModelMEADError ("Cannot find atomic radius for atom type %s" % atomType)
        systemRadii.append (radius)


      # Create subdirectories, if necessary
      if self.splitToDirectories:

        for meadSite in self.meadSites:
          sitePqr   = meadSite.instances[0].sitePqr
          directory = os.path.dirname (sitePqr)

          if not os.path.exists (directory):
            try:
              os.makedirs (directory)
            except:
              raise CEModelMEADError ("Cannot create directory %s" % directory)


      # Write two PQR files for each instance of every site, first for the model compund and second for the site itself
      for meadSite in self.meadSites:
        model = Selection (meadSite.modelAtomIndices)
        site  = Selection (meadSite.siteAtomIndices)

        for instance in meadSite.instances:

          # In the PQR file of the model compound, charges of the site atoms must be set to zero (requirement of the my_2diel_solver program)
          chargesUpdated = Clone (systemCharges)
          for atomIndex in meadSite.siteAtomIndices:
            chargesUpdated[atomIndex] = 0.0

          PQRFile_FromSystem (instance.modelPqr, self.system, selection = model, charges = chargesUpdated, radii = systemRadii)


          chargesUpdated = Clone (systemCharges)

# This fragment does not work when Initialize_Testing is used
#           for atomName, atomIndex in zip (meadSite.siteAtomNames, meadSite.siteAtomIndices):
#             pickCharge                = instance.charges[meadSite.siteAtomNames.index (atomName)]
#             chargesUpdated[atomIndex] = pickCharge

          for chargeIndex, atomIndex in enumerate (meadSite.siteAtomIndices):
            pickCharge                = instance.charges[chargeIndex]
            chargesUpdated[atomIndex] = pickCharge


          PQRFile_FromSystem (instance.sitePqr,  self.system, selection = site,  charges = chargesUpdated, radii = systemRadii)


      # Write background PQR file
      PQRFile_FromSystem (self.backPqr, self.system, selection = Selection (self.backAtomIndices), charges = systemCharges, radii = systemRadii)

      # Write full-protein PQR file (to be used as eps2set_region)
      PQRFile_FromSystem (self.proteinPqr, self.system, selection = Selection (self.proteinAtomIndices), charges = systemCharges, radii = systemRadii)

      # Write FPT-file
      lines = []

      for siteIndex, meadSite in enumerate (self.meadSites):
        for instanceIndex, instance in enumerate (meadSite.instances):
          for atomIndex, charge in zip (meadSite.siteAtomIndices, instance.charges):
            x, y, z = self.system.coordinates3[atomIndex]
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
