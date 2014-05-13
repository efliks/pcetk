#!/usr/bin/python

data  = open ("state_energy.out").readlines ()
table = []

for line in data:
  if line.startswith ("State"):
    tokens = line.split ()
    state  = int (tokens[1])
    energy = float (tokens[5])
    table.append ([energy, state])

table.sort ()

for energy, state in table:
  print ("%8d     %10.4f" % (state, energy))
