# Fit proteins
set selstr "protein and backbone and not hydrogen and resid 1 to 128"

# Template
set       mol      "../7LYZ/7LYZ"
mol       new      $mol.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set       ref      [atomselect 0 $selstr]

set       mol      "2LZT"
set       molid    1
mol       new      $mol.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set       fit      [atomselect $molid $selstr]
set       fit_all  [atomselect $molid all]

$fit_all  move     [measure fit $fit $ref]
$fit_all  writepdb "2LZT_temp.pdb"

# End of script
quit
