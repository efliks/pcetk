"""Example: Protonation states in IFP2.0."""

from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from ContinuumElectrostatics import MEADModel, MEADSubstate, MCModelDefault, TitrationCurves

parameters = ("toppar/par_all27_prot_na.inp", "toppar/blf.inp")

mol = CHARMMPSFFile_ToSystem ("setup/ifp20_xplor_separated.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (parameters))
mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("setup/ifp20.crd")
mol.Summary ()

mead = MEADModel (system=mol, pathMEAD="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)
mead.Initialize (excludeResidues=(("PRTA", "CYS", 24),))
mead.Summary ()

mead.WriteJobFiles ()
mead.CalculateElectrostaticEnergies (calculateETA=False)

sampling = MCModelDefault ()
mead.DefineMCModel (sampling)

tc = TitrationCurves (mead)
tc.CalculateCurves ()
tc.WriteCurves (directory="curves")

selection = (("PRTA", "HIS", 260), ("CHRO", "BLF", 1), ("CHRO", "ACB", 2), ("CHRO", "ACC", 3))
substate = MEADSubstate (mead, selection, pH=7.0)
substate.CalculateSubstateEnergies ()
substate.Summary ()
