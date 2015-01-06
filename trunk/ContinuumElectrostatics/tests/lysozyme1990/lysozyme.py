# Example script for the ContinuumElectrostatics module
# http://www.poissonboltzmann.org/apbs/examples/pka-calculations/lysozyme-pka-example

from pBabel   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore    import logFile

from ContinuumElectrostatics import MEADModel, MEADSubstate, StateVector


logFile.Header ("Calculate protonation states in lysozyme")


#===========================================
par_tab = ["toppar/par_all27_prot_na.inp", ]

mol  = CHARMMPSFFile_ToSystem ("lysozyme1990_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (par_tab))
mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("lysozyme1990.crd")


#===========================================
cem = MEADModel (system=mol, pathMEAD="/home/mikolaj/local/bin/", pathGMCT="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)

# Exclude cysteines involved in disulfide bonds
# Also exclude all arginines (otherwise the system has too many sites for analytic treatment)
exclusions = (
("PRTA", "CYS",   6),
("PRTA", "CYS", 127),
("PRTA", "CYS",  30),
("PRTA", "CYS", 115),
("PRTA", "CYS",  64),
("PRTA", "CYS",  80),
("PRTA", "CYS",  76),
("PRTA", "CYS",  94),
("PRTA", "ARG",   0),
)

cem.Initialize (excludeResidues=exclusions)
cem.Summary ()
cem.SummarySites ()
cem.WriteJobFiles ()
cem.CalculateElectrostaticEnergies (calculateETA=False, asymmetricSummary=True, asymmetricThreshold=0.01)

cem.WriteW ()
cem.WriteGintr ()


#===========================================
logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 using GMCT ***\n")
cem.CalculateProbabilitiesGMCT ()
cem.SummaryProbabilities ()


logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 analytically ***\n")
cem.CalculateProbabilitiesAnalytically ()
cem.SummaryProbabilities ()


#===========================================
sites = (
("PRTA", "GLU" , 35),
("PRTA", "ASP" , 52),
#  ("PRTA", "ASP" , 66),
#  ("PRTA", "HIS" , 15),
)

substate = MEADSubstate (cem, sites)
substate.CalculateSubstateEnergies ()
substate.Summary ()


#===========================================
# vector = StateVector (cem)
# go = True
# 
# while go:
#     cem.CalculateMicrostateEnergy (vector)
#     go = vector.Increment ()


#===========================================
# logFile.Text ("\n*** Calculating titration curves using GMCT ***\n")
# cem.CalculateCurves (directory="curves_gmct")
#  
# logFile.Text ("\n*** Calculating titration curves analytically ***\n")
# cem.CalculateCurves (directory="curves_analytic", isAnalytic=True, forceSerial=True)


#===========================================
logFile.Footer ()
