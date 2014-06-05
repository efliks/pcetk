# Test script for the new module

from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3, PDBFile_FromSystem

from pCore import Pickle, Unpickle, logFile

from ContinuumElectrostatics import MEADModel


logFile.Header ("Testing a small protein (only chain A)")


#===========================================
par_tab = ["charmm/par_all27_prot_na.inp", ]

mol  = CHARMMPSFFile_ToSystem ("charmm/defensin_xplor.psf", isXPLOR = True, parameters = CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/defensin.crd")


#===========================================
ce_model = MEADModel (meadPath = "/home/mikolaj/local/bin/", gmctPath = "/home/mikolaj/local/bin/", scratch = "scratch", nthreads = 8)

ce_model.Initialize (mol)

ce_model.Summary ()

ce_model.SummarySites ()

ce_model.WriteJobFiles (mol)

ce_model.CalculateEnergies ()



logFile.Text ("\n*** Calculating curves analytically in serial mode ***\n")
ce_model.nthreads = 1
ce_model.CalculateCurves (isAnalytic = True,  directory = "curves_analytic_serial")


logFile.Text ("\n*** Calculating curves using GMCT in serial mode ***\n")
ce_model.nthreads = 1
ce_model.CalculateCurves (isAnalytic = False,  directory = "curves_gmct_serial")


logFile.Text ("\n*** Calculating curves analytically in parallel mode ***\n")
ce_model.nthreads = 4
ce_model.CalculateCurves (isAnalytic = True,  directory = "curves_analytic_parallel")


logFile.Text ("\n*** Calculating curves using GMCT in parallel mode ***\n")
ce_model.nthreads = 4
ce_model.CalculateCurves (isAnalytic = False,  directory = "curves_gmct_parallel")


#===========================================
logFile.Footer ()
