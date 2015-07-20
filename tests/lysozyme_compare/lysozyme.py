# Lysozyme example with two alternative structures
from pCore                    import logFile
from pBabel                   import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from ContinuumElectrostatics  import MEADModel, MEADSubstate, StateVector, TitrationCurves, MCModelDefault, MCModelGMCT


# Exclude cysteines involved in disulfide bonds and arginines
exclusions = (
    ( "PRTA"  ,  "CYS"  ,    6 ),
    ( "PRTA"  ,  "CYS"  ,  127 ),
    ( "PRTA"  ,  "CYS"  ,   30 ),
    ( "PRTA"  ,  "CYS"  ,  115 ),
    ( "PRTA"  ,  "CYS"  ,   64 ),
    ( "PRTA"  ,  "CYS"  ,   80 ),
    ( "PRTA"  ,  "CYS"  ,   76 ),
    ( "PRTA"  ,  "CYS"  ,   94 ),
    ( "PRTA"  ,  "ARG"  ,    0 ), )
parameters = ["toppar/par_all27_prot_na.inp", ]
models     = {}

for label, proteinPsf, proteinCrd, meadDir, curveDir, message in (
    ( "old" , "setup/7LYZ/lysozyme1977_xplor.psf" , "setup/7LYZ/lysozyme1977.crd" , "mead/7LYZ", "curves/7LYZ" , "Old lysozyme" ),
    ( "new" , "setup/2LZT/lysozyme1990_xplor.psf" , "setup/2LZT/lysozyme1990.crd" , "mead/2LZT", "curves/2LZT" , "New lysozyme" ), ):

    logFile.Text ("\n*** Now calculating: %s ***\n" % message)

    protein = CHARMMPSFFile_ToSystem (proteinPsf, isXPLOR=True, parameters=CHARMMParameterFiles_ToParameters (parameters))
    protein.coordinates3 = CHARMMCRDFile_ToCoordinates3 (proteinCrd)

    model = MEADModel (system=protein, pathMEAD="/home/mikolaj/local/bin/", pathScratch=meadDir, nthreads=1)
    model.Initialize (excludeResidues=exclusions, includeTermini=True)
    model.Summary ()
    model.SummarySites ()
    model.WriteJobFiles ()
    model.CalculateElectrostaticEnergies ()

    mcModel = MCModelDefault ()
    model.DefineMCModel (mcModel)

    curves  = TitrationCurves (model, curveSampling=.5)
    curves.CalculateCurves ()
    curves.WriteCurves (directory=curveDir)
    curves.CalculateHalfpKs ()
    models[label] = [model, curves]


logFile.Text ("\n*** pK1/2 values for the old model ***\n")
model, curves = models["old"]
curves.PrintHalfpKs ()

logFile.Text ("\n*** pK1/2 values for the new model ***\n")
model, curves = models["new"]
curves.PrintHalfpKs ()

# End of script
logFile.Footer ()


# Table from the paper Miteva et al. Nucleic Acids Research, 2005, 33, p. W372-375
#        residue          pKintr    pKcalc    pKbashford        pKexper
# earlierResults = (
#   ( "NTR"   ,   998   ,    5.6   ,    5.1   ,    6.4   ,    7.8   ,    8.0  ), 
#   ( "HIS"   ,    15   ,    3.5   ,    2.4   ,    4.0   ,    5.8   ,    5.8  ), 
#   ( "GLU"   ,     7   ,    5.5   ,    3.2   ,    2.1   ,    2.6   ,    2.6  ), 
#   ( "GLU"   ,    35   ,    6.5   ,    5.7   ,    6.3   ,    6.1   ,    6.1  ), 
#   ( "ASP"   ,    18   ,    3.7   ,    1.6   ,    3.1   ,    2.8   ,    3.0  ), 
#   ( "ASP"   ,    48   ,    5.3   ,    2.5   ,    1.0   ,    4.3   ,    4.3  ), 
#   ( "ASP"   ,    52   ,    6.9   ,    7.4   ,    7.0   ,    3.5   ,    3.7  ), 
#   ( "ASP"   ,    66   ,    5.9   ,    1.5   ,    1.7   ,    1.5   ,    2.5  ), 
#   ( "ASP"   ,    87   ,    3.6   ,    1.9   ,    1.2   ,    3.5   ,    3.75 ), 
#   ( "ASP"   ,   101   ,    5.3   ,    4.3   ,    7.9   ,    4.0   ,    4.25 ), 
#   ( "ASP"   ,   119   ,    4.6   ,    3.6   ,    3.2   ,    2.2   ,    2.8  ), 
#   ( "TYR"   ,    20   ,   12.5   ,   12.7   ,   14.0   ,   10.3   ,   10.3  ), 
#   ( "TYR"   ,    23   ,   10.2   ,    9.5   ,   11.7   ,    9.8   ,    9.8  ), 
#   ( "TYR"   ,    53   ,   12.9   ,   16.0   ,   20.8   ,   12.1   ,   12.1  ), 
#   ( "LYS"   ,     1   ,    9.7   ,   11.2   ,    9.6   ,   10.7   ,   10.9  ), 
#   ( "LYS"   ,    13   ,    9.6   ,   12.9   ,   11.6   ,   10.4   ,   10.6  ), 
#   ( "LYS"   ,    33   ,   10.3   ,   10.0   ,    9.6   ,   10.5   ,   10.7  ), 
#   ( "LYS"   ,    96   ,   10.4   ,   10.7   ,   10.4   ,   10.7   ,   10.9  ), 
#   ( "LYS"   ,    97   ,   10.6   ,   10.9   ,   10.6   ,   10.2   ,   10.4  ), 
#   ( "LYS"   ,   116   ,   10.4   ,   10.3   ,    9.9   ,   10.3   ,   10.5  ), 
#   ( "CTR"   ,   999   ,    5.0   ,    2.7   ,    2.3   ,    2.7   ,    2.8  ),)
