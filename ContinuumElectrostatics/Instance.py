#-------------------------------------------------------------------------------
# . File      : Instance.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014-2016)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
from pCore  import logFile, LogFileActive
from Error  import ContinuumElectrostaticsError


class Instance (object):
    """Base class for an instance of a titratable site.

    This class should not be used directly."""
    defaultAttributes = {
        }
    # parent  instIndex  _instIndexGlobal  label  charges  Gborn_model  Gback_model  Gborn_protein  Gback_protein

    def __init__ (self, **keywordArguments):
        """Constructor."""
        for (key, val) in keywordArguments.iteritems ():
            setattr (self, key, val)


    def _GetEnergyModel (self):
        if hasattr (self, "parent"):
            site        = self.parent
            ceModel     = site.parent
            energyModel = ceModel.energyModel
        else:
            raise ContinuumElectrostaticsError ("Energy model is undefined.")
        return energyModel


    @property
    def Gmodel (self):
        em = self._GetEnergyModel ()
        return em.GetGmodel (self._instIndexGlobal)

    @Gmodel.setter
    def Gmodel (self, value):
        em = self._GetEnergyModel ()
        em.SetGmodel (self._instIndexGlobal, value)

    @property
    def Gintr (self):
        em = self._GetEnergyModel ()
        return em.GetGintr (self._instIndexGlobal)

    @Gintr.setter
    def Gintr (self, value):
        em = self._GetEnergyModel ()
        em.SetGintr (self._instIndexGlobal, value)

    @property
    def protons (self):
        em = self._GetEnergyModel ()
        return em.GetProtons (self._instIndexGlobal)

    @protons.setter
    def protons (self, value):
        em = self._GetEnergyModel ()
        em.SetProtons (self._instIndexGlobal, value)

    @property
    def probability (self):
        em = self._GetEnergyModel ()
        return em.GetProbability (self._instIndexGlobal)

    @probability.setter
    def probability (self, value):
        em = self._GetEnergyModel ()
        em.SetProbability (self._instIndexGlobal, value)

    @property
    def interactions (self):
        # . Return a list of interactions of the current instance with other instances
        if not hasattr (self, "parent"):
            raise ContinuumElectrostaticsError ("Energy model is undefined.")
        site        = self.parent
        ceModel     = site.parent
        energyModel = ceModel.energyModel
        energies    = []
        for indexOther in range (ceModel.ninstances):
            energies.append (energyModel.GetInteractionSymmetric (self._instIndexGlobal, indexOther))
        return energies


    def CalculateModelCompound (self, log=logFile):
        """Calculate Gborn and Gback of a site in a model compound."""
        pass


    def CalculateProtein (self, log=logFile):
        """Calculate Gborn, Gback and Wij of a site in protein environment."""
        pass


    def CalculateGintr (self, log=logFile):
        """Calculate Gintr of an instance of a site in a protein."""
        checks = (
            hasattr (self , "Gborn_protein") ,
            hasattr (self , "Gback_protein") ,
            hasattr (self , "Gborn_model"  ) ,
            hasattr (self , "Gback_model"  ) ,)
        if all (checks):
            self.Gintr = self.Gmodel + (self.Gborn_protein - self.Gborn_model) + (self.Gback_protein - self.Gback_model)


    def PrintInteractions (self, sort=False, log=logFile):
        """Print interactions of an instance of a site with other instances of other sites."""
        if LogFileActive (log):
            site         = self.parent
            model        = site.parent
            interactions = self.interactions

            if model.isCalculated:
                instances = []
                for site in model.sites:
                    for instance in site.instances:
                        wij = interactions[instance._instIndexGlobal]
                        instances.append ([wij, site.segName, site.resName, site.resSerial, instance.label])
                if sort:
                    instances.sort ()

                tab = log.GetTable (columns=[6, 6, 6, 6, 16])
                tab.Start ()
                tab.Heading ("Instance of a site", columnSpan=4)
                tab.Heading ("Wij")
                for wij, segName, resName, resSerial, label in instances:
                    entries = ( ( "%s"     % segName   ),
                                ( "%s"     % resName   ),
                                ( "%d"     % resSerial ),
                                ( "%s"     % label     ),
                                ( "%16.4f" % wij       ), )
                    for entry in entries:
                        tab.Entry (entry)
                tab.Stop ()


    def _TableEntry (self, tab=None, secondsToCompletion=None):
        """Report calculated energies in a table.

        Optionally, include Estimated Time for Accomplishment (ETA).

        ETA has to be calculated outside of this method."""
        if tab:
            site = self.parent
            entries = (( "%s"     % site.segName       ),
                       ( "%s"     % site.resName       ),
                       ( "%d"     % site.resSerial     ),
                       ( "%s"     % self.label         ),
                       ( "%16.4f" % self.Gborn_model   ),
                       ( "%16.4f" % self.Gback_model   ),
                       ( "%16.4f" % self.Gborn_protein ),
                       ( "%16.4f" % self.Gback_protein ),
                       ( "%16.4f" % self.Gmodel        ),
                       ( "%16.4f" % self.Gintr         ), )
            for entry in entries:
                tab.Entry (entry)

            if isinstance (secondsToCompletion, float):
                minutes, seconds = divmod (secondsToCompletion, 60)
                hours, minutes   = divmod (minutes, 60)
                tab.Entry ("%16s" % ("%d:%02d:%02d" % (hours, minutes, seconds)))


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass
