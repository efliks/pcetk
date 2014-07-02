# Example script for the ContinuumElectrostatics module

from pBabel   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore    import logFile

from ContinuumElectrostatics import MEADModel, StateVector_FromProbabilities, MEADSubstate


logFile.Header ("Calculate protonation states in chain A of defensin")


#===========================================
par_tab = ["charmm/toppar/par_all27_prot_na.inp", ]

mol  = CHARMMPSFFile_ToSystem ("charmm/defensin_xplor.psf", isXPLOR = True, parameters = CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/defensin.crd")


#===========================================
ce_model = MEADModel (meadPath = "/home/mikolaj/local/bin/", gmctPath = "/home/mikolaj/local/bin/", scratch = "scratch", nthreads = 4)

ce_model.Initialize (mol)

ce_model.Summary ()

ce_model.SummarySites ()

ce_model.WriteJobFiles (mol)

ce_model.CalculateElectrostaticEnergies ()


#===========================================
logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 analytically ***\n")

ce_model.CalculateProbabilitiesAnalytically ()

ce_model.SummaryProbabilities ()



logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 using GMCT ***\n")

ce_model.CalculateProbabilitiesGMCT ()

ce_model.SummaryProbabilities ()


#===========================================
sites = (
("PRTA", "ASP" , 2),
("PRTA", "GLU", 14),
("PRTA", "ARG", 15),
("PRTA", "CYS", 31),
)

substate = MEADSubstate (ce_model, sites, pH = 7.0)

substate.CalculateSubstateEnergies ()

substate.Summary ()



#===========================================
#  logFile.Text ("\n*** Calculating titration curves analytically ***\n")
#  
#  ce_model.CalculateCurves (isAnalytic = True, forceSerial = True, directory = "curves_analytic")
#  
#  
#  
#  logFile.Text ("\n*** Calculating titration curves using GMCT in serial mode ***\n")
#  
#  ce_model.CalculateCurves (forceSerial = True, directory = "curves_gmct")
#  
#  
#  
#  logFile.Text ("\n*** Calculating titration curves using GMCT in parallel mode ***\n")
#  
#  ce_model.CalculateCurves (directory = "curves_gmct_parallel")


#===========================================
logFile.Footer ()
