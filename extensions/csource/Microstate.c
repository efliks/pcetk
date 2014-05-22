/*------------------------------------------------------------------------------
! . File      : MicrostateEnergy.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "Microstate.h"

Real Microstate_Energy (Real *gintr, Integer *protons, Real **interactions, StateVector *vector, Real protonChemicalPotential) {
  Integer site;
  Integer *v = vector->vector;
  Real    W = 0.0, Gintr = 0.0;

  for (site = 0; site < vector->lenght; site++, v++) {
    instance = *v;

    for (siteInnter = 0; siteInner < vector->lenght; siteInner++) {

    }
  }

  return Gintr - nprotons * protonChemicalPotential + 0.5 * W;
}
