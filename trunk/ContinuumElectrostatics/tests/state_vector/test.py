# Test script for the new module

from pCore import Pickle, Unpickle, logFile

from ContinuumElectrostatics import CEModelMEAD, StateVector


logFile.Header ("Testing of the state vector")


ce_model = Unpickle ("ce_model.pkl")

ce_model.Summary ()

ce_model.SummarySites ()


#===========================================
v = StateVector (ce_model)

v.Print (ce_model)

r = True
iteration = 0

while r:
  r = v.Increment ()

  iteration = iteration + 1
  if (iteration % 10000 == 0):
    logFile.Text ("Iteration: %d\n" % iteration)

v.Print (ce_model)


v.ResetToMaximum ()

v.Print (ce_model)


#===========================================
logFile.Footer ()
