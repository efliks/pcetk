# Test script for the new module

from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3, PDBFile_FromSystem

from pCore import Pickle, Unpickle, logFile

from ContinuumElectrostatics import CEModelMEAD, StateVector


logFile.Header ("Calculate protonation states in lysozyme")


#===========================================
par_tab = ["../charmm/par_all27_prot_na.inp", ]

mol  = CHARMMPSFFile_ToSystem ("../charmm/lysozyme_xplor.psf", isXPLOR = True, parameters = CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("../charmm/lysozyme.crd")


#===========================================
ce_model = CEModelMEAD (system = mol, meadPath = "/home/mikolaj/local/bin/", scratch = "scratch", nthreads = 8)

del ce_model.standardSites["ARG"]
del ce_model.standardSites["CYS"]

ce_model.Initialize ()

ce_model.Summary ()

ce_model.SummarySites ()

ce_model.WriteJobFiles ()

ce_model.CalculateEnergies ()

ce_model.WriteGintr ()

ce_model.WriteW ()


#===========================================
v = StateVector (ce_model)

v.Print (ce_model)
Gmicro = ce_model.CalculateMicrostateEnergy (v, pH = 7.0)

logFile.Text ("Gmicro = %f\n" % Gmicro)


v.Increment ()

v.Print (ce_model)
Gmicro = ce_model.CalculateMicrostateEnergy (v, pH = 7.0)

logFile.Text ("Gmicro = %f\n" % Gmicro)





#===========================================
logFile.Footer ()
