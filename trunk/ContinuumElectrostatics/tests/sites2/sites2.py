# Example script for the ContinuumElectrostatics module

from pBabel   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore    import logFile

from ContinuumElectrostatics import MEADModel, StateVector


logFile.Header ("Calculate protonation states of two titratable sites in a hypothetical peptide.")


#===========================================
par_tab = ["charmm/toppar/par_all27_prot_na.inp", ]

mol  = CHARMMPSFFile_ToSystem ("charmm/testpeptide_xplor.psf", isXPLOR = True, parameters = CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/testpeptide.crd")


#===========================================
ce_model = MEADModel (meadPath = "/home/mikolaj/local/bin/", gmctPath = "/home/mikolaj/local/bin/", scratch = "scratch", nthreads = 1)

ce_model.Initialize (mol)

ce_model.Summary ()

ce_model.SummarySites ()

ce_model.WriteJobFiles (mol)

ce_model.CalculateElectrostaticEnergies ()


#===========================================
logFile.Text ("\n*** Calculating microstate energies of all states at pH = 7 ***\n")

statevector = StateVector (ce_model)
increment   = True

while increment:
  Gmicro = ce_model.CalculateMicrostateEnergy (statevector, pH = 7.0)
  statevector.Print (ce_model, title = "Gmicro = %f" % Gmicro)
  increment = statevector.Increment ()


#===========================================
ce_model.WriteGintr ()

ce_model.WriteW ()


#===========================================
logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 analytically ***\n")

ce_model.CalculateProbabilitiesAnalytically ()

ce_model.SummaryProbabilities ()



logFile.Text ("\n*** Calculating protonation probabilities at pH = 7 using GMCT ***\n")

ce_model.CalculateProbabilitiesGMCT ()

ce_model.SummaryProbabilities ()


#===========================================
logFile.Text ("\n*** Calculating titration curves analytically ***\n")

ce_model.CalculateCurves (isAnalytic = True, forceSerial = True, directory = "curves_analytic")


logFile.Text ("\n*** Calculating titration curves using GMCT ***\n")

ce_model.CalculateCurves (directory = "curves_gmct")


#===========================================
logFile.Footer ()
