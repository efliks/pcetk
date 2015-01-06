# Example script for the ContinuumElectrostatics module

from pBabel   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore    import logFile

from ContinuumElectrostatics import MEADModel, MEADSubstate


logFile.Header ("Calculate protonation states in chain A of defensin")

#===========================================
par_tab = ["charmm/toppar/par_all27_prot_na.inp", ]

mol = CHARMMPSFFile_ToSystem ("charmm/defensin_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (par_tab))
mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/defensin.crd")


#===========================================
cem = MEADModel (system=mol, pathMEAD= "/home/mikolaj/local/bin/", pathGMCT="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)

cem.Initialize ()
cem.Summary ()
cem.SummarySites ()
cem.WriteJobFiles ()
cem.CalculateElectrostaticEnergies ()


#===========================================
logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 analytically ***\n")
cem.CalculateProbabilitiesAnalytically ()
cem.SummaryProbabilities ()


logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 using GMCT ***\n")
cem.CalculateProbabilitiesGMCT ()
cem.SummaryProbabilities ()


#===========================================
sites = (
("PRTA", "ASP" , 2),
("PRTA", "GLU", 14),
("PRTA", "ARG", 15),
)

substate = MEADSubstate (cem, sites, pH=7.0)
substate.CalculateSubstateEnergies ()
substate.Summary ()


#===========================================
#  logFile.Text ("\n*** Calculating titration curves analytically ***\n")
#  cem.CalculateCurves (isAnalytic = True, forceSerial = True, directory = "curves_analytic")
#  
#  
#  logFile.Text ("\n*** Calculating titration curves using GMCT in serial mode ***\n")
#  cem.CalculateCurves (forceSerial = True, directory = "curves_gmct")
#  
#  
#  logFile.Text ("\n*** Calculating titration curves using GMCT in parallel mode ***\n")
#  cem.CalculateCurves (directory = "curves_gmct_parallel")


#===========================================
logFile.Footer ()
