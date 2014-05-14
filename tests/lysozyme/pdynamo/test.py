# Test script for the new module

from pBabel import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3, PDBFile_FromSystem

from pCore import Pickle, Unpickle, logFile

from ContinuumElectrostatics import CEModelMEAD, StateVector


logFile.Header ("Calculate protonation states in lysozyme")


#===========================================
par_tab = ["../charmm/par_all27_prot_na.inp", ]

mol  = CHARMMPSFFile_ToSystem ("../charmm/lysozyme_xplor.psf", isXPLOR = True, parameters = CHARMMParameterFiles_ToParameters (par_tab))

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("../charmm/lysozyme.crd")


#===========================================
ce_model = CEModelMEAD (system = mol, meadPath = "/home/mikolaj/local/bin/", scratch = "scratch", nthreads = 8)


del ce_model.librarySites["ARG"]
del ce_model.librarySites["CYS"]

ce_model.Initialize ()

ce_model.Summary ()

ce_model.SummarySites ()

ce_model.WriteJobFiles ()

ce_model.CalculateEnergies ()

ce_model.WriteGintr ()

ce_model.WriteW ()

Pickle ("ce_model.pkl", ce_model)


#===========================================
vector    = StateVector (ce_model)

#increment = True
#nstates   = 0
#
#while increment:
#  increment = vector.Increment ()
#  Gmicro    = ce_model.CalculateMicrostateEnergy (vector)
#
#  logFile.Text ("Gmicro = %g\n" % Gmicro)
#
#  nstates = nstates + 1
#  if nstates > 5000:
#    break


vector.Reset ()

vector.Print (ce_model)

energy = ce_model.CalculateMicrostateEnergy (vector)

logFile.Text ("Gmicro = %.6f\n" % energy)


vector.Increment ()

vector.Print (ce_model)

energy = ce_model.CalculateMicrostateEnergy (vector)

logFile.Text ("Gmicro = %.6f\n" % energy)


vector.Increment ()

vector.Print (ce_model)

energy = ce_model.CalculateMicrostateEnergy (vector)

logFile.Text ("Gmicro = %.6f\n" % energy)


vector.Increment ()

vector.Print (ce_model)

energy = ce_model.CalculateMicrostateEnergy (vector)

logFile.Text ("Gmicro = %.6f\n" % energy)

# vector[0]  = 0    #   PRTA_LYS1     p   1
# vector[1]  = 1    #   PRTA_GLU7     d   0
# vector[2]  = 1    #  PRTA_LYS13     d   0
# vector[3]  = 2    #  PRTA_HIS15   HSE   1
# vector[4]  = 1    #  PRTA_ASP18     d   0
# vector[5]  = 0    #  PRTA_TYR20     p   1
# vector[6]  = 0    #  PRTA_TYR23     p   1
# vector[7]  = 0    #  PRTA_LYS33     p   1
# vector[8]  = 0    #  PRTA_GLU35     p   1
# vector[9]  = 1    #  PRTA_ASP48     d   0
# vector[10] = 0    #  PRTA_ASP52     p   1
# vector[11] = 0    #  PRTA_TYR53     p   1
# vector[12] = 1    #  PRTA_ASP66     d   0
# vector[13] = 1    #  PRTA_ASP87     d   0
# vector[14] = 0    #  PRTA_LYS96     p   1
# vector[15] = 0    #  PRTA_LYS97     p   1
# vector[16] = 1    # PRTA_ASP101     d   0
# vector[17] = 0    # PRTA_LYS116     p   1
# vector[18] = 1    # PRTA_ASP119     d   0
# 
# vector.Print (ce_model)
# energy = ce_model.CalculateMicrostateEnergy (vector)
# logFile.Text ("Gmicro = %.6f\n" % energy)


#===========================================
logFile.Footer ()
