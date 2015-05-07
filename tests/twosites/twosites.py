# Example script
from pCore                   import logFile
from pBabel                  import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from ContinuumElectrostatics import MEADModel, StateVector, MCModelDefault, TitrationCurves


logFile.Header ("Calculate protonation states of two titratable sites in a hypothetical peptide.")

parameters = ["charmm/toppar/par_all27_prot_na.inp", ]
mol = CHARMMPSFFile_ToSystem ("charmm/testpeptide_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (parameters))
mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/testpeptide.crd")


cem = MEADModel (system=mol, pathMEAD="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)
cem.Initialize ()
cem.Summary ()
cem.SummarySites ()
cem.WriteJobFiles ()
cem.CalculateElectrostaticEnergies ()


logFile.Text ("\n*** Calculating microstate energies of all states at pH=7 ***\n")

statevector = StateVector (cem)
increment   = True
while increment:
    Gmicro = cem.CalculateMicrostateEnergy (statevector, pH=7.0)
    statevector.Print (title="Gmicro = %f" % Gmicro)
    increment = statevector.Increment ()


logFile.Text ("\n*** Calculating protonation probabilities at pH=7 analytically ***\n")
cem.CalculateProbabilities (pH=7.0)
cem.SummaryProbabilities ()

logFile.Text ("\n*** Calculating protonation probabilities at pH=7 using in-house MC sampling ***\n")
mc = MCModelDefault (nprod=30000)
cem.DefineMCModel (mc)
cem.CalculateProbabilities ()
cem.SummaryProbabilities ()


#===========================================
logFile.Text ("\n*** Calculating titration curves analytically ***\n")
cmc = TitrationCurves (cem, curveSampling=0.5)
cmc.CalculateCurves ()
cmc.WriteCurves (directory="curves_analytic")

logFile.Text ("\n*** Calculating titration curves using in-house MC sampling ***\n")
cem.DefineMCModel (None)
ca = TitrationCurves (cem, curveSampling=0.5)
ca.CalculateCurves ()
ca.WriteCurves (directory="curves_mc")


#===========================================
logFile.Footer ()
