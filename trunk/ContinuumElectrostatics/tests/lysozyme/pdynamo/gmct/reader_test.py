#===== test =====

from pCore                       import logFile
from ContinuumElectrostatics     import GMCTOutputFileReader


logFile.Header ("test")


temperature = 300.0

reader = GMCTOutputFileReader ("lysozyme.gmct-out")
reader.Parse (temperature)

probs = reader.probabilities['conf_PRTA_ASP87_p']

for pH, prob in zip (reader.pHtable, probs):
  logFile.Text ("%8.3f    %8.3f\n" % (pH, prob))

logFile.Footer ()
