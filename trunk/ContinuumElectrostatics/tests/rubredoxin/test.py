# Test script for the new module

from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3, PDBFile_FromSystem

from pCore import Pickle, Unpickle, logFile

from ContinuumElectrostatics import CEModelMEAD, StateVector


logFile.Header ("Test case for MEAD calculations")


#===========================================
par_tab = ["toppar/ff/par_all27_prot_na.inp", "toppar/cluster/par_iron_sulfur.inp"]

mol  = CHARMMPSFFile_ToSystem ("rubredoxin_oxidized_exclude_iron/rubr_xplor.psf", isXPLOR = True, parameters = CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("rubredoxin_oxidized_exclude_iron/rubr.crd")


#===========================================
ce_model = CEModelMEAD (system = mol, meadPath = "/home/mikolaj/local/bin/", scratch = "scratch", nthreads = 8)

ce_model.Initialize ()

ce_model.Summary ()

ce_model.SummarySites ()

ce_model.WriteJobFiles ()

ce_model.CalculateEnergies ()

ce_model.WriteGintr ()

ce_model.WriteW ()


#===========================================
#vector    = StateVector (ce_model)
#increment = True
#
#while increment:
#  Gmicro = ce_model.CalculateMicrostateEnergy (vector, pH = 7.0)
#
#  increment = vector.Increment ()
#logFile.Footer ()


#===========================================
#site = ce_model.meadSites[0]
#inst = site.instances[0]
#inst.PrintInteractions (sort = False)


#vector.Increment ()
#vector.Print (ce_model)
#moreVectors = True
#iteration   = 0
#
#vector.Reset ()
# 
#while moreVectors:
#  moreVectors = vector.Increment ()
#
#  iteration = iteration + 1
#  if (iteration % 10000) == 0:
#    logFile.Text ("Iteration: %d\n" % iteration)
