/*------------------------------------------------------------------------------
! . File      : EnergyModel.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "EnergyModel.h"


EnergyModel *EnergyModel_Allocate (const Integer nsites, const Integer ninstances, Status *status) {
  EnergyModel *self = NULL;

  MEMORY_ALLOCATE (self, EnergyModel);
  if (self == NULL) {
    Status_Set (status, Status_MemoryAllocationFailure);
  }
  else {
    self->vector           =  NULL  ;
    self->protons          =  NULL  ;
    self->intrinsic        =  NULL  ;
    self->interactions     =  NULL  ;
    self->probabilities    =  NULL  ;
    self->symmetricmatrix  =  NULL  ;

    if (ninstances > 0) {
      self->vector          = StateVector_Allocate     (nsites,     status)             ;
      self->protons         = Integer1DArray_Allocate  (ninstances, status)             ;
      self->intrinsic       = Real1DArray_Allocate     (ninstances, status)             ;
      self->interactions    = Real2DArray_Allocate     (ninstances, ninstances, status) ;
      self->probabilities   = Real1DArray_Allocate     (ninstances, status)             ;
      self->symmetricmatrix = SymmetricMatrix_Allocate (ninstances)                     ;
      if (self->symmetricmatrix == NULL)
        Status_Set (status, Status_MemoryAllocationFailure);

      if (*status != Status_Continue)
        EnergyModel_Deallocate (self);
    }
  }
  return self;
}

void EnergyModel_Deallocate (EnergyModel *self) {
  if ( self->symmetricmatrix != NULL)  SymmetricMatrix_Deallocate ( &self->symmetricmatrix ) ;
  if ( self->probabilities   != NULL)  Real1DArray_Deallocate     ( &self->probabilities   ) ;
  if ( self->interactions    != NULL)  Real2DArray_Deallocate     ( &self->interactions    ) ;
  if ( self->intrinsic       != NULL)  Real1DArray_Deallocate     ( &self->intrinsic       ) ;
  if ( self->protons         != NULL)  Integer1DArray_Deallocate  ( &self->protons         ) ;
  if ( self->vector          != NULL)  StateVector_Deallocate     (  self->vector          ) ;
  if (self != NULL) MEMORY_DEALLOCATE (self);
}

Boolean EnergyModel_CheckInteractionsSymmetric (const EnergyModel *self, Real tolerance, Real *maxDeviation) {
  return Real2DArray_IsSymmetric (self->interactions, &tolerance, maxDeviation);
}

void EnergyModel_SymmetrizeInteractions (const EnergyModel *self, Status *status) {
  SymmetricMatrix_CopyFromReal2DArray (self->symmetricmatrix, self->interactions, status);
}

/*=============================================================================
  Getters
=============================================================================*/
Real EnergyModel_GetGintr (const EnergyModel *self, const Integer instIndexGlobal) {
  return Real1DArray_Item (self->intrinsic, instIndexGlobal);
}

Integer EnergyModel_GetProtons (const EnergyModel *self, const Integer instIndexGlobal) {
  return Integer1DArray_Item (self->protons, instIndexGlobal);
}

Real EnergyModel_GetProbability (const EnergyModel *self, const Integer instIndexGlobal) {
  return Real1DArray_Item (self->probabilities, instIndexGlobal);
}

Real EnergyModel_GetInteraction (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
  return Real2DArray_Item (self->interactions, instIndexGlobalA, instIndexGlobalB);
}

Real EnergyModel_GetInteractionSymmetric (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
  Real value;

  if (instIndexGlobalA >= instIndexGlobalB) {
    value = SymmetricMatrix_Item (self->symmetricmatrix, instIndexGlobalA, instIndexGlobalB);
  }
  else {
    value = SymmetricMatrix_Item (self->symmetricmatrix, instIndexGlobalB, instIndexGlobalA);
  }
  return value;
}

Real EnergyModel_GetDeviation (const EnergyModel *self, const Integer i, const Integer j) {
  Real wij, wji, deviation;

  wij = Real2DArray_Item (self->interactions, i, j);
  wji = Real2DArray_Item (self->interactions, j, i);
  deviation = (wij + wji) * .5 - wij;
  return deviation;
}

/*=============================================================================
  Setters
=============================================================================*/
void EnergyModel_SetGintr (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
  Real1DArray_Item (self->intrinsic, instIndexGlobal) = value;
}

void EnergyModel_SetProtons (const EnergyModel *self, const Integer instIndexGlobal, const Integer value) {
  Integer1DArray_Item (self->protons, instIndexGlobal) = value;
}

void EnergyModel_SetProbability (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
  Real1DArray_Item (self->probabilities, instIndexGlobal) = value;
}

void EnergyModel_SetInteraction (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB, const Real value) {
  Real2DArray_Item (self->interactions, instIndexGlobalA, instIndexGlobalB) = value;
}

/*=============================================================================
  Calculating microstate energy
=============================================================================*/
Real EnergyModel_CalculateMicrostateEnergy (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature) {
  Real      Gintr, W, *data;
  Integer   nprotons, i, j;
  TitrSite *site, *siteInner;

  W        = 0. ;
  Gintr    = 0. ;
  nprotons = 0  ;

  site = vector->sites;
  for (i = 0; i < vector->nsites; i++, site++) {
    Gintr     += Real1DArray_Item    (self->intrinsic , site->indexActive);
    nprotons  += Integer1DArray_Item (self->protons   , site->indexActive);

    data       = &self->symmetricmatrix->data[(site->indexActive * (site->indexActive + 1)) >> 1];
    siteInner  = vector->sites;
    for (j = 0; j <= i; j++, siteInner++) {
      W += *(data + (siteInner->indexActive));
    }
  }
  return (Gintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * temperature * CONSTANT_LN10 * pH) + W);
}

/*=============================================================================
  Calculating probabilities analytically
=============================================================================*/
void EnergyModel_CalculateProbabilitiesAnalytically (const EnergyModel *self, const Real pH, const Real temperature, Status *status) {
  Real         *bfactor, G, Gzero, bsum;
  Real1DArray  *bfactors;
  Integer       i, j;
  TitrSite     *ts;

  bfactors = Real1DArray_Allocate (self->nstates, status);
  if (*status != Status_Continue) {
    return;
  }

  /* Calculate state energies of ALL states */
  StateVector_Reset (self->vector);
  i       = self->nstates;
  bfactor = Real1DArray_Data (bfactors);
  Gzero   = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH, temperature);

  for (; i > 0; i--, bfactor++) {
    G = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH, temperature);
    if (G < Gzero) {
      Gzero = G;
    }
    *bfactor = G;
    StateVector_Increment (self->vector);
  }

  /* Convert energies to Boltzmann factors */
  Real1DArray_AddScalar (bfactors, -Gzero);
  Real1DArray_Scale (bfactors, -1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * temperature));
  Real1DArray_Exp (bfactors);

  /* Calculate probabilities */
  Real1DArray_Set (self->probabilities, 0.);
  StateVector_Reset (self->vector);

  i = self->nstates;
  bfactor = Real1DArray_Data (bfactors);
  for (; i > 0; i--, bfactor++) {
    j  = self->vector->nsites ;
    ts = self->vector->sites  ;
    for (; j > 0; j--, ts++) {
      Real1DArray_Item (self->probabilities, ts->indexActive) += *bfactor;
    }
    StateVector_Increment (self->vector);
  }

  bsum = Real1DArray_Sum (bfactors);
  Real1DArray_Scale (self->probabilities, 1. / bsum);

  Real1DArray_Deallocate (&bfactors);
}

/*=============================================================================
  Calculating probabilities using Metropolis Monte Carlo
=============================================================================*/
Boolean EnergyModel_Metropolis (const Real GdeltaRT) {
  Boolean metropolis;
  Real ran;

  /* Prepare */
  ran = rand () / (Real) RAND_MAX;

  /* Apply the Metropolis criterion; based on GMCT */
  if (GdeltaRT < 0.)
    metropolis = True;
  else {
    if (-GdeltaRT < TOO_SMALL)
      metropolis = False;
    else if (ran < exp (-GdeltaRT))
      metropolis = True;
    else
      metropolis = False;
  }
  return metropolis;
}

/* The purpose of a scan is to generate a state vector representing a low-energy, statistically relevant protonation state */
Real EnergyModel_MCScan (const EnergyModel *self, const Real pH, const Real temperature, Integer nmoves) {
  Real      G, Gnew, GdeltaRT, beta;
  Integer   site, oldIndexActive;
  Boolean   accept;
  TitrSite *ts;

  G = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH, temperature);
  beta = 1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * temperature);

  for (; nmoves > 0; nmoves--) {
    /* Perform the move */
    StateVector_Move (self->vector, &site, &oldIndexActive);
    Gnew     = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH, temperature);
    GdeltaRT = (Gnew - G) * beta;
    accept   = EnergyModel_Metropolis (GdeltaRT);

    if (accept == True) {
      /* Accept the move */
      G = Gnew;
    }
    else {
      /* Reject the move and revert to the previous state */
      ts = &self->vector->sites[site];
      ts->indexActive = oldIndexActive;
    }
  }
  /* Return the last accepted energy (only for info) */
  return G;
}

void EnergyModel_CalculateProbabilitiesMonteCarlo (const EnergyModel *self, const Real pH, const Real temperature, const Boolean equil, Integer nscans, Status *status) {
  Real    Gfinal, scale;
  Integer nmoves;

  /* The number of moves is proportional to the number of sites */
  nmoves = self->vector->nsites;
  scale  = 1. / nscans;

  /* Equilibration phase? */
  if (equil) {
    StateVector_Randomize (self->vector);

    /* Run the scans */
    for (; nscans > 0; nscans--)
      Gfinal = EnergyModel_MCScan (self, pH, temperature, nmoves);
  }

  /* Production phase? */
  else {
    /* Reset the probabilities */
    Real1DArray_Set (self->probabilities, 0.);

    /* Run the scans */
    for (; nscans > 0; nscans--) {
      Gfinal = EnergyModel_MCScan (self, pH, temperature, nmoves);

      /* Update the counts */
      EnergyModel_UpdateProbabilities (self);
    }
    /* Average the probabilities */
    Real1DArray_Scale (self->probabilities, scale);
  }
}

void EnergyModel_UpdateProbabilities (const EnergyModel *self) {
  Real     *counter;
  TitrSite *ts;
  Integer   i;
  ts = self->vector->sites  ;
  i  = self->vector->nsites ;

  for (; i > 0; i--, ts++) {
    counter = Real1DArray_ItemPointer (self->probabilities, ts->indexActive);
    (*counter)++;
  }
}
