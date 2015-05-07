#!/usr/bin/python
# Detach acetate groups from the chromophore and treat them as separate titratable sites

"""
    5766 CHRO 1    BLF  CAC  CPM   -0.820000E-01   12.0110           0
"""

import sys


acetate_b = ("CBD" ,"HO3" ,"HO4" ,"CGD" ,"O2D" ,"O1D" ,"HO2D","CAD" ,"HO7" ,"HO8")
acetate_c = ("CBA" ,"HO1" ,"HO2" ,"CGA" ,"O2A" ,"O1A" ,"HO2A","CAA" ,"HO5" ,"HO6")

if len (sys.argv) < 2:
    psf = "parent_xplor.psf"
else:
    psf = sys.argv[1]


lines = open (psf).readlines ()
lines2 = []

format_line = "%8s%5s%2d%7s%5s%5s%16s%10s%12s\n"


for line in lines:
    if line.count ("CHRO"):
        tokens = line.split ()
        atom_serial, seg_label, res_serial, res_label, atom_label, atom_type, atom_charge, atom_radius, filler = tokens

        if atom_label in acetate_b:
            line = format_line % (atom_serial, seg_label, int (res_serial) + 1, "ACB", atom_label, atom_type, atom_charge, atom_radius, filler)

        if atom_label in acetate_c:
            line = format_line % (atom_serial, seg_label, int (res_serial) + 2, "ACC", atom_label, atom_type, atom_charge, atom_radius, filler)
    lines2.append (line)


w = open (psf.replace (".psf", "_separated.psf"), "w")
w.writelines (lines2)
w.close ()
