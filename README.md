# Pcetk
Pcetk (a pDynamo-based continuum electrostatic toolkit) is a Python module
extending the pDynamo library with a Poisson-Boltzmann continuum electrostatic
model that allows for protonation state calculations in proteins.

The module links pDynamo to the external solver of the Poisson-Boltzmann
equation, extended-MEAD, which is used for the calculation of electrostatic
energy terms. The calculation of protonation states and titration curves is
done by using the module's own analytic or Monte Carlo routines or through
an interface to the external sampling program, GMCT.

The key features of the module include:
  * Improved treatment of multiprotic sites, such as histidines or multiprotic 
        ligands
  * Energy model that allows for independent treatment of sites
  * Fast calculations thanks to the parallelization
  * Analytic calculation of protonation states (for small proteins) or Monte 
        Carlo sampling (in-house method or by using GMCT)
  * Automatic generation of titration curves
  * Calculations of substate energies
  * Easily extensible scripting environment based on Python
  * Time-consuming parts of code written in C and Cython
  * Close integration with pDynamo

More information about the theoretical methods and algorithms used by Pcetk can be found
in our article:
Feliks, M.; Field, M. "*Pcetk*: A pDynamo-based Toolkit for Protonation State Calculations in Proteins."
*J. Chem. Inf. Model.* **2015**, *55*, pp. 2288-2296 (http://dx.doi.org/10.1021/acs.jcim.5b00262)


Auxiliary software:
 * pDynamo (Martin Field)
    http://www.pdynamo.org

 * Extended-MEAD (Donald Bashford, Timm Essigke, Thomas R. Ullmann)
    http://www.bisb.uni-bayreuth.de/People/ullmannt/index.php?name=extended-mead

 * GMCT (Matthias G. Ullmann, Thomas R. Ullmann)
    http://www.bisb.uni-bayreuth.de/People/ullmannt/index.php?name=gmct-gcem
