#-------------------------------------------------------------------------------
# . File      : ProtonationForm.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2015)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from Constants import *


class ProtonationForm (object):
    """Protonation form of a model compound."""

    defaultAttributes = {
        "label"          :   "ProtonationForm" ,
        "temperature"    :   TEMPERATURE       ,
        "fileCbs"        :   "cbs.out"         ,
        "fileHess"       :   "hess.out"        ,
        "fileCosmo"      :   "cosmo.out"       ,
        "isParsed"       :   None              ,
        "isCalcTerms"    :   None              ,
        "isCalcDeltas"   :   None              ,
        "isCorrected"    :   None              ,
        "isExperimental" :   None              ,
                        }

    def __init__ (self, *arguments, **keywordArguments):
        """Constructor."""
        for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
        for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


    def _Parser (self, f, items):
        """Parse ORCA output file."""
        lines = open (f).readlines ()
        ret   = {}
        for line in lines:
            for item, item_token_index in items:
                if line.count (item):
                    tokens = line.split ()
                    ret[item] = float (tokens[item_token_index])
        return ret


    def ParseFiles (self):
        """Get information from ORCA files."""
        if not self.isParsed:
            # Level: B3LYP/6-311G+(d,p) or ACCOPT
            pattern        = "FINAL SINGLE POINT ENERGY"
            self.Eel       = self._Parser (self.fileHess  , ((pattern, -1), ))[pattern] * UNITS_ENERGY_HARTREES_TO_KILOCALORIES_PER_MOLE
            pattern        = "G-E(el)"
            self.Gcorr     = self._Parser (self.fileHess  , ((pattern,  2), ))[pattern] * UNITS_ENERGY_HARTREES_TO_KILOCALORIES_PER_MOLE
            pattern        = "Total Energy after outlying charge correction"
            self.Eel_cosmo = self._Parser (self.fileCosmo , ((pattern, -2), ))[pattern] * UNITS_ENERGY_HARTREES_TO_KILOCALORIES_PER_MOLE
            # Level: MP2/CBS
            pattern        = "FINAL SINGLE POINT ENERGY"
            self.Eel_cbs   = self._Parser (self.fileCbs   , ((pattern, -1), ))[pattern] * UNITS_ENERGY_HARTREES_TO_KILOCALORIES_PER_MOLE
            self.isParsed  = True


    def CalculateTerms (self, useLowLevelEel=False):
        """Calculate Esolv and Gaq."""
        if self.isParsed:
            self.Esol = self.Eel_cosmo - self.Eel
            if useLowLevelEel : self.Gaq = self.Eel     + self.Esol + self.Gcorr
            else              : self.Gaq = self.Eel_cbs + self.Esol + self.Gcorr
            self.isCalcTerms = True


    def CalculateDeltas (self, other):
        """Calculate DeltaG(aq) and pKa(aq)."""
        if self.isCalcTerms:
            if not other.isCalcTerms:
                if not other.isParsed:
                    other.ParseFiles ()
                other.CalculateTerms ()
            self.DeltaGaq     = self.Gaq + G_PROTON_AQ - other.Gaq
            self.pKa          = self.DeltaGaq / (CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10)
            self.isCalcDeltas = True


    def CalculateExperimental (self, pKa_exp):
        """Add experimental pKa and calculate difference pKa_exp-pKa."""
        if self.isCalcDeltas:
            self.pKa_exp        = pKa_exp
            self.pKa_corr       = pKa_exp - self.pKa
            self.DeltaG_corr    = self.pKa_corr * CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10
            self.isExperimental = True


    def CalculateCorrections (self, a=DEFAULT_A, b=DEFAULT_B):
        """Correct pKa by linear fitting."""
        if hasattr (self, "pKa"):
            self.pKa_lfit = b + a * self.pKa
            if hasattr (self, "pKa_exp"):
                self.pKa_error = self.pKa_lfit - self.pKa_exp
            self.Gmodel      = self.pKa_lfit * CONSTANT_MOLAR_GAS_KCAL_MOL * self.temperature * CONSTANT_LN10
            self.isCorrected = True


    def DoAll (self, other, a=DEFAULT_A, b=DEFAULT_B):
        """Run all calculations."""
        self.ParseFiles      ()
        self.CalculateTerms  ()
        self.CalculateDeltas (other)
        self.CorrectDeltas   (a=a, b=b)


    def Summary (self):
        """Print a summary."""
        print ("*** %s ***"       % self.label   )
        print ("pKa      = %6.2f" % self.pKa     )
        if self.isCorrected:
            print ("pKa_lfit = %6.2f" % self.pKa_lfit)
            print ("Gmodel   = %6.2f" % self.Gmodel  )
            print ("=================")


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass
