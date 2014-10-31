#!/usr/bin/python

import Cython.Compiler.Main as main
import os, glob, sys


# If getenv does not work, set the pCore directory manually, for example:
# pdynamo_pcore = "/home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0"
if len (sys.argv) < 2:
    pdynamo_pcore = os.getenv ("PDYNAMO_PCORE")
else:
    pdynamo_pcore = sys.argv[1]

current_directory  = os.getcwd ()
pxd_directories    = [current_directory, os.path.join (pdynamo_pcore, "extensions/pyrex")]

# Compile all files in the cython directory
sources = glob.glob (os.path.join (current_directory, "*.pyx"))

for source in sources:
    options = main.CompilationOptions (main.default_options)
    options.include_path.extend (pxd_directories)
    main.compile (source, options)
