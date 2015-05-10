#!/usr/bin/python
# Fitting taken from
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html

from   pkcalc2 import *
import numpy


class TrainingSet (object):
    """Training set of acids."""

    def __init__ (self, rootDir, acidDirs, method, useLowLevelEel=False):
        """Constructor."""
        acids = []
        for acidDir in acidDirs:
            prot   = ProtonationForm (
                label     = "%s" % acidDir.replace ("_", "-"),
                fileHess  = "%s/%s/prot/%s/prot_%s.out"             % (rootDir, acidDir, method, method),
                fileCbs   = "%s/%s/prot/%s/cbs/prot_%s_cbs.out"     % (rootDir, acidDir, method, method),
                fileCosmo = "%s/%s/prot/%s/cosmo/prot_%s_cosmo.out" % (rootDir, acidDir, method, method),)
            prot.ParseFiles      ()
            prot.CalculateTerms  (useLowLevelEel=useLowLevelEel)
    
            deprot = ProtonationForm (
                label     = "%s" % acidDir.replace ("_", "-"),
                fileHess  = "%s/%s/deprot/%s/deprot_%s.out"             % (rootDir, acidDir, method, method),
                fileCbs   = "%s/%s/deprot/%s/cbs/deprot_%s_cbs.out"     % (rootDir, acidDir, method, method),
                fileCosmo = "%s/%s/deprot/%s/cosmo/deprot_%s_cosmo.out" % (rootDir, acidDir, method, method),)
            deprot.ParseFiles            ()
            deprot.CalculateTerms        (useLowLevelEel=useLowLevelEel)
            deprot.CalculateDeltas       (prot)
            deprot.CalculateExperimental (pKa_exp=pKasExperimental[acidDir])

            acids.append (deprot)
        self.acids  = acids
        self.nacids = len (acids)


    def CalculateCorrectionFunction (self):
        """Find parameters of correction function."""
        nacids = len (self.acids)
        calc   = numpy.arange (nacids, dtype=numpy.float)
        exp    = numpy.arange (nacids, dtype=numpy.float)
        for i, acid in enumerate (self.acids):
            calc[i] = acid.pKa_corr
            exp [i] = acid.pKa_exp
        A      = numpy.vstack ([exp, numpy.ones (len (exp))]).T
        c, d   = numpy.linalg.lstsq (A, calc)[0]
        self.a = c / (1. - c) + 1.
        self.b = d / (1. - c)


    def CorrectAcids (self):
        """Calculate corrected pKas of acids and standard deviation."""
        sigma  = 0.
        for acid in self.acids:
            acid.CalculateCorrections (self.a, self.b)
            sigma += acid.pKa_error ** 2
        nacids     = len (self.acids)
        self.sigma = (1. / nacids * sigma) ** .5


    def WriteTable (self, filename="fit.dat", sortByExperiment=False):
        """Write table with corrected pKas."""
        if sortByExperiment:
            pairs = []
            for iacid, acid in enumerate (self.acids):
                pairs.append ([iacid, acid.pKa_exp])
            pairs.sort (key=lambda k: k[1])

            sequence = []
            for iacid, foo in pairs:
                sequence.append (iacid)
        else:
            sequence = range (0, self.nacids)

        fo = open (filename, "w")
        fo.write ("# pKa_lfit = %f + %f * pKa_calc\n" % (self.b, self.a))
        fo.write ("#%13s%14s%14s%14s%14s\n" % ("pKa_exp", "pKa_calc", "DeltaG_corr", "pKa_lfit", "error"))

        for i in sequence:
            acid = self.acids[i]
            fo.write ("%14.1f%14.1f%14.1f%14.1f%14.1f     \"%s\"\n" % (acid.pKa_exp, acid.pKa, acid.DeltaG_corr, acid.pKa_lfit, acid.pKa_error, acid.label))
        fo.write ("# sigma=%f\n" % self.sigma)
        fo.close ()


#======================================
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

trainAccopt = TrainingSet ("../cationic_acids", acids, "accopt", useLowLevelEel=True)
trainAccopt.CalculateCorrectionFunction ()
trainAccopt.CorrectAcids ()
trainAccopt.WriteTable (filename="accopt_fit_without_corr.dat", sortByExperiment=True)

trainB3LYP = TrainingSet ("../cationic_acids", acids, "b3lyp", useLowLevelEel=True)
trainB3LYP.CalculateCorrectionFunction ()
trainB3LYP.CorrectAcids ()
trainB3LYP.WriteTable (filename="b3lyp_fit_without_corr.dat", sortByExperiment=True)
