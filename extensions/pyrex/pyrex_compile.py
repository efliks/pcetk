#!/usr/bin/python

import Cython.Compiler.Main as main

import os, glob


pdynamo_pcore      =  os.getenv ("PDYNAMO_PCORE")
current_directory  =  os.getcwd ()

pxd_directories = [current_directory, os.path.join (pdynamo_pcore, "extensions/pyrex")]

# Compile all files in pyrex directory
sources = glob.glob (os.path.join (current_directory, "*.pyx"))

for source in sources:
    options = main.CompilationOptions (main.default_options)
    options.include_path.extend (pxd_directories)

    main.compile (source, options)
