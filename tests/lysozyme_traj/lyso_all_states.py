"""Example: Protonation states in lysozyme."""

import os, sys

# . Initialize pDynamo
pDynamo_version = "1.9.0"
pDynamo_root = os.path.join (os.environ["HOME"], "local", "opt", "pDynamo-%s" % pDynamo_version)

pDynamo_pbabel = os.path.join (pDynamo_root, "pBabel-%s" % pDynamo_version)
pDynamo_pcore = os.path.join (pDynamo_root, "pCore-%s" % pDynamo_version)
pDynamo_pmolecule = os.path.join (pDynamo_root, "pMolecule-%s" % pDynamo_version)
pDynamo_pmoleculescripts = os.path.join (pDynamo_root, "pMoleculeScripts-%s" % pDynamo_version)

for path in (pDynamo_pbabel, pDynamo_pcore, pDynamo_pmolecule, pDynamo_pmoleculescripts):
    sys.path.append (path)

pDynamo_parameters = os.path.join (pDynamo_root, "parameters")
pDynamo_style = os.path.join (pDynamo_parameters, "ccsStyleSheets", "defaultStyle.css")
pDynamo_scratch = "/tmp"

os.environ["PDYNAMO_PARAMETERS"] = pDynamo_parameters
os.environ["PDYNAMO_STYLE"] = pDynamo_style
os.environ["PDYNAMO_SCRATCH"] = pDynamo_scratch

# . Initialize Pcetk
pcetk_root = os.path.join (os.environ["HOME"], "devel", "pcetk_testing")
sys.path.append (pcetk_root)

os.environ["PDYNAMO_PCETK"] = pcetk_root


#-------------------------------------------------------------------------------
from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3

from ContinuumElectrostatics import MEADModel, MCModelDefault, TitrationCurves, StateVector


parameters = ["toppar/par_all27_prot_na.inp", ]

mol = CHARMMPSFFile_ToSystem ("setup/lysozyme1977_xplor.psf", isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (parameters))
mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("setup/lysozyme1977.crd")

cem = MEADModel (system=mol, pathMEAD="/home/mikolaj/local/bin/", pathScratch="mead", nthreads=1)

exclusions = (
("PRTA", "CYS",   6),  ("PRTA", "CYS", 127),  ("PRTA", "CYS",  30),
("PRTA", "CYS", 115),  ("PRTA", "CYS",  64),  ("PRTA", "CYS",  80),
("PRTA", "CYS",  76),  ("PRTA", "CYS",  94),  ("PRTA", "ARG",   0), )

cem.Initialize (excludeResidues=exclusions, includeTermini=True)
cem.Summary ()
cem.SummarySites ()
cem.WriteJobFiles ()
cem.CalculateElectrostaticEnergies (calculateETA=False)

statevector = StateVector (cem)
increment   = True
collect     = []

while increment:
    Gmicro = cem.CalculateMicrostateEnergy (statevector, pH=7.0)
    collect.append (Gmicro)
    increment = statevector.Increment ()

print ("\nCalculated %d states." % len (collect))

output = open ("stat_anl.dat", "w")
for gmicro in collect:
    output.write ("%f\n" % gmicro)
output.close ()
