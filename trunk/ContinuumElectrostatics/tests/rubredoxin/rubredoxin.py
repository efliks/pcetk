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
cem = MEADModel (system=mol, pathMEAD="/home/mikolaj/local/bin/", pathGMCT="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)

# Do not take cysteines
exclusions = (
  ("", "CYS", ""),
)

cem.Initialize (excludeResidues=exclusions)
cem.Summary ()
cem.SummarySites ()
cem.WriteJobFiles ()
cem.CalculateElectrostaticEnergies (calculateETA=False)


logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 using GMCT ***\n")
cem.CalculateProbabilitiesGMCT ()
cem.SummaryProbabilities ()

logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 analytically ***\n")
cem.CalculateProbabilitiesAnalytically ()
cem.SummaryProbabilities ()


logFile.Text ("\n*** Calculating titration curves using GMCT ***\n")
cem.CalculateCurves (directory="curves_gmct")

# Analytic curves may be slow on some computers
logFile.Text ("\n*** Calculating titration curves analytically ***\n")
cem.CalculateCurves (directory="curves_analytic", isAnalytic=True, forceSerial=True)


#===========================================
logFile.Footer ()
