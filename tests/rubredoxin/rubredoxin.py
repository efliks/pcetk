# Example script for the ContinuumElectrostatics module

from pBabel   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore    import logFile

from ContinuumElectrostatics import MEADModel


logFile.Header ("Calculate protonation states in rubredoxin")


#===========================================
par_tab = ["charmm/toppar/par_all27_prot_na_fixed.inp", "charmm/toppar/par_iron_sulfur.inp"]

mol  = CHARMMPSFFile_ToSystem ("charmm/rubredoxin_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/rubredoxin.crd")


#===========================================
ce_model = MEADModel (meadPath="/home/mikolaj/local/bin/", gmctPath="/home/mikolaj/local/bin/", scratch="mead", nthreads=4)

# Do not take cysteines
exclusions = (
  ("", "CYS", ""),
)

ce_model.Initialize (mol, excludeResidues=exclusions)
ce_model.Summary ()
ce_model.SummarySites ()
ce_model.WriteJobFiles (mol)
ce_model.CalculateElectrostaticEnergies (calculateETA=False)


logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 using GMCT ***\n")
ce_model.CalculateProbabilitiesGMCT ()
ce_model.SummaryProbabilities ()

logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 analytically ***\n")
ce_model.CalculateProbabilitiesAnalytically ()
ce_model.SummaryProbabilities ()


logFile.Text ("\n*** Calculating titration curves using GMCT ***\n")
ce_model.CalculateCurves (directory="curves_gmct")

# Analytic curves may be slow on some computers
logFile.Text ("\n*** Calculating titration curves analytically ***\n")
ce_model.CalculateCurves (directory="curves_analytic", isAnalytic=True, forceSerial=True)


#===========================================
logFile.Footer ()
