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
  if (self == NULL)
    goto failSet;

  self->vector           =  NULL  ;
  self->generator        =  NULL  ;

  self->protons          =  NULL  ;
  self->intrinsic        =  NULL  ;
  self->interactions     =  NULL  ;
  self->probabilities    =  NULL  ;
  self->symmetricmatrix  =  NULL  ;

  /* Initialize RNG */
  self->generator = RandomNumberGenerator_Allocate (RandomNumberGeneratorType_MersenneTwister);
  if (self->generator == NULL)
    goto failSetDealloc;
  RandomNumberGenerator_SetSeed (self->generator, (Cardinal) time (NULL));

  if (nsites > 0) {
    self->vector = StateVector_Allocate (nsites, status);
    if (*status != Status_Continue)
      goto failDealloc;
  }

  if (ninstances > 0) {
    self->protons         = Integer1DArray_Allocate  (ninstances, status)             ;
    self->intrinsic       = Real1DArray_Allocate     (ninstances, status)             ;
    self->interactions    = Real2DArray_Allocate     (ninstances, ninstances, status) ;
    self->probabilities   = Real1DArray_Allocate     (ninstances, status)             ;
    if (*status != Status_Continue)
      goto failDealloc;

    self->symmetricmatrix = SymmetricMatrix_Allocate (ninstances);
    if (self->symmetricmatrix == NULL)
      goto failSetDealloc;
  }
  return self;


failSetDealloc:
  Status_Set (status, Status_MemoryAllocationFailure);

failDealloc:
  EnergyModel_Deallocate (self);
  return NULL;

failSet:
  Status_Set (status, Status_MemoryAllocationFailure);
  return NULL;
}

void EnergyModel_Deallocate (EnergyModel *self) {
  if ( self->generator       != NULL)  RandomNumberGenerator_Deallocate ( &self->generator       ) ;
  if ( self->symmetricmatrix != NULL)  SymmetricMatrix_Deallocate       ( &self->symmetricmatrix ) ;
  if ( self->probabilities   != NULL)  Real1DArray_Deallocate           ( &self->probabilities   ) ;
  if ( self->interactions    != NULL)  Real2DArray_Deallocate           ( &self->interactions    ) ;
  if ( self->intrinsic       != NULL)  Real1DArray_Deallocate           ( &self->intrinsic       ) ;
  if ( self->protons         != NULL)  Integer1DArray_Deallocate        ( &self->protons         ) ;
  if ( self->vector          != NULL)  StateVector_Deallocate           (  self->vector          ) ;
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

Real EnergyModel_GetInterSymmetric (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
  return EnergyModel_GetW (self, instIndexGlobalA, instIndexGlobalB);
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
  Real      Gintr, W, *interact;
  Integer   nprotons, i, j;
  TitrSite *site, *siteInner;

  W        = 0. ;
  Gintr    = 0. ;
  nprotons = 0  ;

  site = vector->sites;
  for (i = 0; i < vector->nsites; i++, site++) {
    Gintr     += Real1DArray_Item    (self->intrinsic , site->indexActive);
    nprotons  += Integer1DArray_Item (self->protons   , site->indexActive);

    interact   = EnergyModel_RowPointer (self, site->indexActive);
    siteInner  = vector->sites;
    for (j = 0; j <= i; j++, siteInner++) {
      W += *(interact + (siteInner->indexActive));
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
  Monte Carlo-related functions
=============================================================================*/
Boolean EnergyModel_Metropolis (const Real GdeltaRT, const RandomNumberGenerator *generator) {
  Boolean metropolis;

  /* Apply the Metropolis criterion; based on GMCT */
  if (GdeltaRT < 0.)
    metropolis = True;
  else {
    if (-GdeltaRT < TOO_SMALL)
      metropolis = False;
    else if (RandomNumberGenerator_NextReal (generator) < exp (-GdeltaRT))
      metropolis = True;
    else
      metropolis = False;
  }
  return metropolis;
}

/*
 * The purpose of a scan is to generate a state vector representing a low-energy, statistically relevant protonation state
 */
Real EnergyModel_MCScan (const EnergyModel *self, const Real pH, const Real temperature, Integer nmoves) {
  Integer    site, siteOther, oldActive, oldActiveOther, selection, select;
  Real       G, Gnew, GdeltaRT, beta;
  Boolean    accept, isDouble;
  TitrSite  *ts;

  G         = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH, temperature);
  beta      = 1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * temperature);
  selection = self->vector->nsites + self->vector->npairs;

  for (; nmoves > 0; nmoves--) {
    select = RandomNumberGenerator_NextCardinal (self->generator) % selection;

    if (select < self->vector->nsites) {
      /* Perform a single move */
      isDouble = False;
      StateVector_Move (self->vector, &site, &oldActive, self->generator);
    }
    else {
      /* Perform a double move */
      isDouble = True;
      StateVector_DoubleMove (self->vector, &site, &siteOther, &oldActive, &oldActiveOther, self->generator);
    }

    Gnew     = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH, temperature);
    GdeltaRT = (Gnew - G) * beta;
    accept   = EnergyModel_Metropolis (GdeltaRT, self->generator);

    if (accept) {
      /* Accept the move */
      G = Gnew;
    }
    else {
      /* Revert the single move */
      ts              = &self->vector->sites[site];
      ts->indexActive = oldActive;

      /* Revert the double move */
      if (isDouble) {
        ts              = &self->vector->sites[siteOther];
        ts->indexActive = oldActiveOther;
      }
    }
  }
  return G;
}

/*
 * Increase the counts of "active" instances.
 * These counts, after scaling, will give the probabilities of occurrence of instances.
 */
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

/*
 * Finds maximum absolute interaction energy between two sites.
 */
Real EnergyModel_FindMaxInteraction (const EnergyModel *self, const TitrSite *site, const TitrSite *other) {
  Integer index, indexOther;
  Real W, Wmax;

  Wmax  = 0.;
  index = site->indexFirst;
  for (; index <= site->indexLast; index++) {

    indexOther = other->indexFirst;
    for (; indexOther <= other->indexLast; indexOther++) {
      W = fabs (EnergyModel_GetW (self, index, indexOther));
      if (W > Wmax) Wmax = W;
    }
  }
  return Wmax;
}

/*
 * Finds pairs of sites whose interaction energy is greater than the given limit.
 * 
 * If npairs < 1, dry run is assumed and only nfound is returned.
 * The value of npairs is used in the second run to allocate and fill out the pairs.
 */
Integer EnergyModel_FindPairs (const EnergyModel *self, const Real limit, const Integer npairs, Status *status) {
  TitrSite *site, *siteInner;
  Integer i, j, nfound;
  Real Wmax;

  if (npairs > 0) {
    StateVector_AllocatePairs (self->vector, npairs, status);
    if (*status != Status_Continue)
      return -1;
  }

  nfound = 0;
  site   = self->vector->sites;
  for (i = 0; i < self->vector->nsites; i++, site++) {

    siteInner = self->vector->sites;
    for (j = 0; j <= i; j++, siteInner++) {

      Wmax = EnergyModel_FindMaxInteraction (self, site, siteInner);
      if (Wmax >= limit) {
        if (npairs > 0) {
          StateVector_SetPair (self->vector, nfound, site->indexSite, siteInner->indexSite, Wmax, status);
        }
        nfound++;
      }
    }
  }
  return nfound;
}

/*
 * Runs a Monte Carlo simulation.
 *
 * If equil is true, only the equilibration run is performed. Otherwise the production run is done.
 */
void EnergyModel_CalculateProbabilitiesMonteCarlo (const EnergyModel *self, const Real pH, const Real temperature, const Boolean equil, Integer nscans, Status *status) {
  Real    Gfinal, scale;
  Integer nmoves;

  /* The number of moves is proportional to the number of sites and pairs */
  nmoves = self->vector->nsites + self->vector->npairs;
  scale  = 1. / nscans;

  /* Equilibration phase? */
  if (equil) {
    StateVector_Randomize (self->vector, self->generator);
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
