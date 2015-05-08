#!/usr/bin/python

from   pkcalc import *


blva = ProtonationForm (label     = "a"                        ,
                        fileHess  = "data/a/a.out"             ,
                        fileCbs   = "data/a/cbs/a_cbs.out"     ,
                        fileCosmo = "data/a/cosmo/a_cosmo.out" ,)

blvb = ProtonationForm (label     = "b"                        ,
                        fileHess  = "data/b/b.out"             ,
                        fileCbs   = "data/b/cbs/b_cbs.out"     ,
                        fileCosmo = "data/b/cosmo/b_cosmo.out" ,)

blvc = ProtonationForm (label     = "c"                        ,
                        fileHess  = "data/c/c.out"             ,
                        fileCbs   = "data/c/cbs/c_cbs.out"     ,
                        fileCosmo = "data/c/cosmo/c_cosmo.out" ,)

blvd = ProtonationForm (label     = "d"                        ,
                        fileHess  = "data/d/d.out"             ,
                        fileCbs   = "data/d/cbs/d_cbs.out"     ,
                        fileCosmo = "data/d/cosmo/d_cosmo.out" ,)

blvf = ProtonationForm (label     = "full"                           ,
                        fileHess  = "data/full/full.out"             ,
                        fileCbs   = "data/full/cbs/full_cbs.out"     ,
                        fileCosmo = "data/full/cosmo/full_cosmo.out" ,)


for form in blva, blvb, blvc, blvd, blvf:
    form.ParseFiles     ()
    form.CalculateTerms ()

for form in blva, blvb, blvc, blvd:
    form.CalculateDeltas (blvf)
    form.CorrectDeltas   (a=0.623742, b=6.990237)
    form.Summary         ()
# pKa_lfit = 6.990237 + 0.623742 * pKa_calc
