# Example script for the ContinuumElectrostatics module
#
# http://www.poissonboltzmann.org/apbs/examples/pka-calculations/lysozyme-pka-example

from pBabel   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore    import logFile

from ContinuumElectrostatics import MEADModel, MEADSubstate


logFile.Header ("Calculate protonation states in lysozyme")


#===========================================
par_tab = ["toppar/par_all27_prot_na.inp", ]

mol  = CHARMMPSFFile_ToSystem ("lysozyme1990_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("lysozyme1990.crd")


#===========================================
ce_model = MEADModel (meadPath="/home/mikolaj/local/bin/", gmctPath="/home/mikolaj/local/bin/", scratch="mead", nthreads=6)

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

ce_model.Initialize (mol, excludeResidues=exclusions)
ce_model.Summary ()
ce_model.SummarySites ()
ce_model.WriteJobFiles (mol)
ce_model.CalculateElectrostaticEnergies (calculateETA=False)


#===========================================
logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 using GMCT ***\n")
ce_model.CalculateProbabilitiesGMCT ()
ce_model.SummaryProbabilities ()


logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 analytically ***\n")
ce_model.CalculateProbabilitiesAnalytically ()
ce_model.SummaryProbabilities ()


#===========================================
sites = (
("PRTA", "GLU" , 35),
("PRTA", "ASP" , 52),
#  ("PRTA", "ASP" , 66),
#  ("PRTA", "HIS" , 15),
)

substate = MEADSubstate (ce_model, sites)
substate.CalculateSubstateEnergies ()
substate.Summary ()


#===========================================
# logFile.Text ("\n*** Calculating titration curves using GMCT ***\n")
# ce_model.CalculateCurves (directory="curves_gmct")
#  
# logFile.Text ("\n*** Calculating titration curves analytically ***\n")
# ce_model.CalculateCurves (directory="curves_analytic", isAnalytic=True, forceSerial=True)


#===========================================
logFile.Footer ()
