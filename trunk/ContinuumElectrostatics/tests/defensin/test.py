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


logFile.Text ("\n*** Calculating titration curves with GMCT ***\n")

ce_model.CalculateCurvesSerial (directory = "curves")


#logFile.Text ("\n*** Calculating titration curves analytically ***\n")

#ce_model.CalculateCurves (analytically = True, directory = "curves_analytic")



#  logFile.Text ("\n*** Calculating probabilities analytically ***\n")
#  
#  ce_model.CalculateProbabilitiesAnalytically ()
#  
#  ce_model.SummaryProbabilities ()
#  
#  
#  logFile.Text ("\n*** Calculating probabilities in GMCT ***\n")
#  
#  ce_model.CalculateProbabilitiesGMCT ()
#  
#  ce_model.SummaryProbabilities ()

#===========================================
logFile.Footer ()
