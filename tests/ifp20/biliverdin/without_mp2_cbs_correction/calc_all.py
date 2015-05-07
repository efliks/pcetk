#!/usr/bin/python

# Fitting taken from
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html

from   pkcalc import *
import numpy  as np


acids = ( "4_bromo_anilinium"      ,
          "anilinium"              ,
          "246_collidinium"        ,
          "benzylammonium"         ,
          "triethylammonium"       ,
          "pyrrolidinium"          ,
          "guanidinium"            ,
          # Own compounds          
          "4_cyano_anilinium"      ,
          "lysine"                 ,
          # Previously excluded
          "25_dichloro_anilinium"  ,
          "p-anisidinium"          ,
          "26_dimethyl_pyridinium" ,
          "DMAP"                   , )
source  = "../cationic_acids"
results = {}

for method in ("b3lyp", "accopt"):
    for acid in acids:
        prot   = ProtonationForm (
            label     = "%s:%s:prot"                            % (        acid, method        ),
            fileHess  = "%s/%s/prot/%s/prot_%s.out"             % (source, acid, method, method),
            fileCbs   = "%s/%s/prot/%s/cbs/prot_%s_cbs.out"     % (source, acid, method, method),
            fileCosmo = "%s/%s/prot/%s/cosmo/prot_%s_cosmo.out" % (source, acid, method, method),)
        prot.ParseFiles      ()
        prot.CalculateTerms  (useLowLevelEel=True)

        deprot = ProtonationForm (
            label     = "%s:%s:deprot"                              % (        acid, method        ),
            fileHess  = "%s/%s/deprot/%s/deprot_%s.out"             % (source, acid, method, method),
            fileCbs   = "%s/%s/deprot/%s/cbs/deprot_%s_cbs.out"     % (source, acid, method, method),
            fileCosmo = "%s/%s/deprot/%s/cosmo/deprot_%s_cosmo.out" % (source, acid, method, method),)
        deprot.ParseFiles            ()
        deprot.CalculateTerms        (useLowLevelEel=True)
        deprot.CalculateDeltas       (prot)
        deprot.CalculateExperimental (pKa_exp=pKasExperimental[acid])
        deprot.Summary               ()

        if results.has_key (method):
            acids = results[method]
        else:
            acids = {}
        acids[acid]     = deprot
        results[method] = acids


# Generate a list of acids sorted depending on their experimental pKas
pairs = []

for acid in acids:
    deprot  = results["accopt"][acid]
    pKa_exp = deprot.pKa_exp
    pairs.append ([acid, pKa_exp])
pairs.sort (key=lambda k: k[1])

acids_sorted = []
for acid, pKa_exp in pairs:
    acids_sorted.append (acid)


#======================================
for method in ("b3lyp", "accopt"):
    exp  = []
    corr = []
    for acid in acids:
        deprot = results[method][acid]
        exp.append  (deprot.pKa_exp)
        corr.append (deprot.pKa_corr)
    x = np.array (exp)
    y = np.array (corr)

    A = np.vstack ([x, np.ones (len (x))]).T
    c, d = np.linalg.lstsq (A, y)[0]
    a = 1. + c / (1. - c)
    b = d / (1. - c)

    sigma = 0.
    lines = []
    for acid in acids_sorted:
        deprot    = results[method][acid]
        pKa_lfit  = a * deprot.pKa + b
        error     = pKa_lfit - deprot.pKa_exp
        sigma    += error ** 2
        lines.append ([deprot.pKa_exp, deprot.pKa, deprot.DeltaG_corr, pKa_lfit, error, acid])

    nacids = len (acids)
    sigma  = (1. / nacids * sigma) ** .5


    fo = open ("%s_fit_without_corr.dat" % method, "w")
    fo.write ("# pKa_lfit = %f + %f * pKa_calc\n" % (b, a))
    fo.write ("#%13s%14s%14s%14s%14s\n" % ("pKa_exp", "pKa_calc", "DeltaG_corr", "pKa_lfit", "error"))

    for pKa_exp, pKa_calc, DeltaG_corr, pKa_lfit, error, acid in lines:
        fo.write ("%14.1f%14.1f%14.1f%14.1f%14.1f     \"%s\"\n" % (pKa_exp, pKa_calc, DeltaG_corr, pKa_lfit, error, acid.replace ("_", "-")))
    fo.write ("# sigma=%f\n" % sigma)
    fo.close ()
