"""Example: Protonation states in lysozyme."""

from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from ContinuumElectrostatics import MEADModel, MCModelDefault, TitrationCurves

parameters = ["par_all27_prot_na.inp", ]

mol = CHARMMPSFFile_ToSystem ("lysozyme_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (parameters))
mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("lysozyme.crd")

cem = MEADModel (system=mol, pathMEAD="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)

exclusions = (
("PRTA", "CYS",   6),  ("PRTA", "CYS", 127),  ("PRTA", "CYS",  30),
("PRTA", "CYS", 115),  ("PRTA", "CYS",  64),  ("PRTA", "CYS",  80),
("PRTA", "CYS",  76),  ("PRTA", "CYS",  94),  ("PRTA", "ARG",   0), )

cem.Initialize (excludeResidues=exclusions)
cem.Summary ()
cem.SummarySites ()
cem.WriteJobFiles ()
cem.CalculateElectrostaticEnergies (calculateETA=False)

cem.CalculateProbabilities (pH=7.0)
cem.SummaryProbabilities ()

curves = TitrationCurves (cem)
curves.CalculateCurves ()
curves.WriteCurves (directory="curves_analytic")
curves.PrintHalfpKs ()

sampling = MCModelDefault ()
cem.DefineMCModel (sampling)

cem.CalculateProbabilities (pH=7.0)
cem.SummaryProbabilities ()

mcc = TitrationCurves (cem)
mcc.CalculateCurves ()
mcc.WriteCurves (directory="curves_mc")
mcc.PrintHalfpKs ()
