# Test script for the new module

from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3, PDBFile_FromSystem

from pCore import logFile

from ContinuumElectrostatics import MEADModel


logFile.Header ("Test case for MEAD calculations")


#===========================================
par_tab = ["toppar/ff/par_all27_prot_na.inp", "toppar/cluster/par_iron_sulfur.inp"]

mol  = CHARMMPSFFile_ToSystem ("rubredoxin_oxidized/rubr_xplor.psf", isXPLOR = True, parameters = CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("rubredoxin_oxidized/rubr.crd")


#===========================================
ce_model = MEADModel (meadPath = "/home/mikolaj/local/bin/", scratch = "scratch", nthreads = 2)

ce_model.Initialize (mol)

ce_model.Summary ()

ce_model.SummarySites ()

ce_model.WriteJobFiles (mol)

#  ce_model.CalculateEnergies ()
#  
#  ce_model.WriteGintr ()
#  
#  ce_model.WriteW ()

logFile.Footer ()
