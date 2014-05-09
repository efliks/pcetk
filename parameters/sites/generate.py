#=========
_DefaultTitratableSites = {
    "ASP" : (( "CB", "HB1", "HB2", "CG", "OD1", "OD2", ),
  ( "p", -5.487135, 1, -0.21, 0.09, 0.09, 0.75, -0.36, -0.36, ),
  ( "d",  0.000000, 0, -0.28, 0.09, 0.09, 0.62, -0.76, -0.76, ),
  ),

  "TYR" : (( "CE1", "HE1", "CE2", "HE2", "CZ", "OH", "HH", ),
  ( "p",-13.169124, 1, -0.12, 0.12, -0.12, 0.12,  0.11, -0.54, 0.43, ),
  ( "d",  0.000000, 0, -0.17, 0.07, -0.17, 0.07, -0.11, -0.69, 0.00, ),
  ),

  "GLU" : (( "CG", "HG1", "HG2", "CD", "OE1", "OE2", ),
  ( "p", -6.035848, 1, -0.21, 0.09, 0.09, 0.75, -0.36, -0.36, ),
  ( "d",  0.000000, 0, -0.28, 0.09, 0.09, 0.62, -0.76, -0.76, ),
  ),

  "CYS" : (( "CB", "HB1", "HB2", "SG", "HG1", ),
  ( "p",-12.483232, 1, -0.11,  0.09,  0.09, -0.23,  0.16, ),
  ( "d",  0.000000, 0, -0.20, -0.05, -0.05, -0.70, -0.00, ),
  ),

  "LYS" : (( "CE", "HE1", "HE2", "NZ", "HZ1", "HZ2", "HZ3", ),
  ( "p",-14.266551, 1, 0.21, 0.05, 0.05, -0.30, 0.33, 0.33, 0.33, ),
  ( "d",  0.000000, 0, 0.10, 0.00, 0.00, -0.25, 0.05, 0.05, 0.05, ),
  ),

  "ARG" : (( "CZ", "CD", "HD1", "HD2", "NE", "HE", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22", ),
  ( "p",-16.461405, 1, 0.64, 0.20, 0.09, 0.09, -0.70, 0.44, -0.80, 0.46, 0.46, -0.80, 0.46, 0.46, ),
  ( "d",  0.000000, 0, 0.48, 0.06, 0.00, 0.00, -0.33, 0.17, -0.39, 0.10, 0.10, -0.39, 0.10, 0.10, ),
  ),

  "HIS" : (("NE2", "CG", "ND1", "HD1", "CE1", "HE2", "CD2", "CB", "HB1", "HB2", "HD2", "HE1", ),
  ( "HSP",  0.000000, 2, -0.51,  0.19, -0.51, 0.44, 0.32, 0.44,  0.19, -0.05, 0.09, 0.09, 0.13, 0.18, ),
  ( "HSD",  9.602486, 1, -0.70, -0.05, -0.36, 0.32, 0.25, 0.00,  0.22, -0.09, 0.09, 0.09, 0.10, 0.13, ),
  ( "HSE",  9.053770, 1, -0.36,  0.22, -0.70, 0.00, 0.25, 0.32, -0.05, -0.08, 0.09, 0.09, 0.09, 0.13, ),
  ( "fd" , 27.434000, 0, -0.70,  0.12, -0.70, 0.00, 0.12, 0.00, -0.05, -0.08, 0.09, 0.09, 0.09, 0.02, ),
  ),
}

site = """
# The Gmodel values were calculated for T = 300K
# They are appropriately recalculated for a given temperature of the model
# 
# Acids have their contribution to Gmicro negative, bases positive (is that true?)
---
site      : %s
atoms     : [%s]
instances :"""

inst ="""
  - label   : %s
    Gmodel  : %s
    protons : %s
    charges : [%s]
"""

items = _DefaultTitratableSites.iteritems ()

for (name, data) in items:
  atoms     = data[0]
  instances = data[1:] 
  output    = site % (name, ", ".join (atoms))
  
  for instance in instances:
    label   = instance[0]
    Gmodel  = "%.6f" % instance[1]
    protons = instance[2]
    charges = ", ".join (map (lambda c: "%5.2f" % c, instance[3:]))

    output  = "%s%s" % (output, inst % (label, Gmodel, protons, charges))

  output = "%s...\n" % output

  filename = "%s.yaml" % name
  f = open (filename, "w")
  f.writelines (output)
  f.close ()

