#!/usr/bin/python

import Cython.Compiler.Main as main
import os, glob, sys


nargs = len (sys.argv)

# If getenv does not work, set the pCore directory manually, for example:
# pdynamo_pcore = "/home/mikolaj/local/opt/pDynamo-1.8.0/pCore-1.8.0"
if nargs < 2:
    pdynamo_pcore = os.getenv ("PDYNAMO_PCORE")
else:
    pdynamo_pcore = sys.argv[1]

current_directory  = os.getcwd ()
pxd_directories    = [current_directory, os.path.join (pdynamo_pcore, "extensions/pyrex")]

# Decide between taking all files in the current directory or the files from the command line
if nargs > 2:
    sources = [os.path.abspath (filename) for filename in sys.argv[2:]]
else:
    sources = glob.glob (os.path.join (current_directory, "*.pyx"))

# Now compile
for source in sources:
    options = main.CompilationOptions (main.default_options)
    options.include_path.extend (pxd_directories)
    main.compile (source, options)
