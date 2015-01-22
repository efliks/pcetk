# Example script for the ContinuumElectrostatics module

from pBabel   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore    import logFile

from ContinuumElectrostatics import MEADModel, StateVector


logFile.Header ("Calculate protonation states of the histidine in a hypothetical peptide.")


#===========================================
par_tab = ["charmm/toppar/par_all27_prot_na.inp", ]

mol = CHARMMPSFFile_ToSystem ("charmm/testpeptide_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (par_tab))
mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/testpeptide.crd")


#===========================================
cem = MEADModel (system=mol, pathMEAD="/home/mikolaj/local/bin/", pathGMCT="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)

cem.Initialize ()
cem.Summary ()
cem.SummarySites ()
cem.WriteJobFiles ()
cem.CalculateElectrostaticEnergies ()


logFile.Text ("\n*** Calculating protonation probabilities at pH=7 analytically ***\n")
cem.CalculateProbabilitiesAnalytically ()
cem.SummaryProbabilities ()

logFile.Text ("\n*** Calculating protonation probabilities at pH=7 using GMCT ***\n")
cem.CalculateProbabilitiesGMCT ()
cem.SummaryProbabilities ()


vector = StateVector (cem)
go = True

while go:
    Gmicro = cem.CalculateMicrostateEnergy (vector)
    vector.Print (title="Gmicro = %f" % Gmicro)
    go = vector.Increment ()


#===========================================
logFile.Footer ()
