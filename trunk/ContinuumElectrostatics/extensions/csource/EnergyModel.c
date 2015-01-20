/*------------------------------------------------------------------------------
! . File      : EnergyModel.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "EnergyModel.h"


EnergyModel *EnergyModel_Allocate (const Integer ninstances, Status *status) {
  EnergyModel *self = NULL;

  MEMORY_ALLOCATE (self, EnergyModel);
  if (self == NULL) {
    Status_Set (status, Status_MemoryAllocationFailure);
  }
  else {
    self->protons          =  NULL  ;
    self->intrinsic        =  NULL  ;
    self->deviations       =  NULL  ;
    self->interactions     =  NULL  ;
    self->probabilities    =  NULL  ;
    self->symmetricmatrix  =  NULL  ;
    self->nstates          =  0     ;
    self->ninstances       =  0     ;

    if (ninstances > 0) {
      self->protons         = Integer1DArray_Allocate  (ninstances, status)             ;
      self->intrinsic       = Real1DArray_Allocate     (ninstances, status)             ;
      self->deviations      = Real2DArray_Allocate     (ninstances, ninstances, status) ;
      self->interactions    = Real2DArray_Allocate     (ninstances, ninstances, status) ;
      self->probabilities   = Real1DArray_Allocate     (ninstances, status)             ;
      self->symmetricmatrix = SymmetricMatrix_Allocate (ninstances)                     ;
      if (self->symmetricmatrix == NULL) {
        Status_Set (status, Status_MemoryAllocationFailure);
      }

      if (*status != Status_Continue) {
        EnergyModel_Deallocate (self);
      }
    }
  }
  return self;
}

void EnergyModel_Deallocate (EnergyModel *self) {
  if ( self->symmetricmatrix != NULL)  SymmetricMatrix_Deallocate ( &self->symmetricmatrix ) ;
  if ( self->probabilities   != NULL)  Real1DArray_Deallocate     ( &self->probabilities   ) ;
  if ( self->interactions    != NULL)  Real2DArray_Deallocate     ( &self->interactions    ) ;
  if ( self->deviations      != NULL)  Real2DArray_Deallocate     ( &self->deviations      ) ;
  if ( self->intrinsic       != NULL)  Real1DArray_Deallocate     ( &self->intrinsic       ) ;
  if ( self->protons         != NULL)  Integer1DArray_Deallocate  ( &self->protons         ) ;
  if (self != NULL) MEMORY_DEALLOCATE (self);
}

void EnergyModel_CalculateDeviations (const EnergyModel *self) {
  Integer i, j;
  Real    wij, wji, deviation;

  for (i = 0; i < self->ninstances; i++)
    for (j = 0; j < self->ninstances; j++) {
      wij = Real2DArray_Item (self->interactions, i, j);
      wji = Real2DArray_Item (self->interactions, j, i);

      deviation = (wij + wji) * .5 - wij;
      Real2DArray_Item (self->deviations, i, j) = deviation;
  }
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

Real EnergyModel_GetDeviation (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
  return Real2DArray_Item (self->deviations, instIndexGlobalA, instIndexGlobalB);
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
  Integer   nprotons, siteIndex, siteIndexInner;
  Integer  *instanceIndex, *instanceIndexInner;

  W        = 0. ;
  Gintr    = 0. ;
  nprotons = 0  ;

  for (siteIndex = 0, instanceIndex = vector->vector; siteIndex < vector->nsites; siteIndex++, instanceIndex++) {
    nprotons += Integer1DArray_Item (self->protons, *instanceIndex);
    Gintr += Real1DArray_Item (self->intrinsic, *instanceIndex);

    /* # define SymmetricMatrix_Item( self, i, j ) ( (self)->data[( (i) * ( i + 1 ) ) / 2 + j] ) */
    data = &self->symmetricmatrix->data[(*instanceIndex * (*instanceIndex + 1)) >> 1];

    for (siteIndexInner = 0, instanceIndexInner = vector->vector; siteIndexInner <= siteIndex; siteIndexInner++, instanceIndexInner++) {
      /* W += SymmetricMatrix_Item (self->symmetricmatrix, *instanceIndex, *instanceIndexInner); */
      W += *(data + (*instanceIndexInner));
    }
  }
  return (Gintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * temperature * CONSTANT_LN10 * pH) + W);
}


/*=============================================================================
  Calculating probabilities analytically
=============================================================================*/
void EnergyModel_CalculateProbabilitiesAnalytically (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature, Status *status) {
  Real1DArray *bfactors;
  Real        *bfactor;
  Real         energy, energyZero, bsum;
  Integer     *activeInstanceGlobalIndex;
  Integer      stateIndex, siteIndex;

  bfactors = Real1DArray_Allocate (self->nstates, status);
  if (*status != Status_Continue) {
    return;
  }
  StateVector_Reset (vector);

  for (stateIndex = 0, bfactor = bfactors->data; stateIndex < self->nstates; stateIndex++, bfactor++) {
    energy = EnergyModel_CalculateMicrostateEnergy (self, vector, pH, temperature);

    if (stateIndex < 1) {
      energyZero = energy;
    }
    else {
      if (energy < energyZero) {
        energyZero = energy;
      }
    }

    *bfactor = energy;
    StateVector_Increment (vector);
  }
  Real1DArray_AddScalar (bfactors, -energyZero);
  Real1DArray_Scale (bfactors, -1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * temperature));
  Real1DArray_Exp (bfactors);
  Real1DArray_Set (self->probabilities, 0.);
  StateVector_Reset (vector);

  for (stateIndex = 0, bfactor = bfactors->data; stateIndex < self->nstates; stateIndex++, bfactor++) {
    for (siteIndex = 0, activeInstanceGlobalIndex = vector->vector; siteIndex < vector->nsites; siteIndex++, activeInstanceGlobalIndex++) {
      self->probabilities->data[*activeInstanceGlobalIndex] += *bfactor;
    }
    StateVector_Increment (vector);
  }
  bsum = Real1DArray_Sum (bfactors);
  Real1DArray_Scale (self->probabilities, 1. / bsum);

  Real1DArray_Deallocate (&bfactors);
}


/*=============================================================================
  Calculating probabilities using Metropolis Monte Carlo
=============================================================================*/
/* The purpose of a scan is to generate a state vector representing a low-energy, statistically relevant protonation state */
Real EnergyModel_MCScan (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature, Integer nmoves) {
  Real      G, Gnew, GdeltaRT, ran, beta;
  Integer   site, instanceBefore, instanceAfter;
  Boolean   accept;

  G = EnergyModel_CalculateMicrostateEnergy (self, vector, pH, temperature);
  beta = 1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * temperature);

  for (; nmoves >= 0; nmoves--) {
    /* Perform the move */
    StateVector_Move (vector, &site, &instanceBefore, &instanceAfter);
    Gnew = EnergyModel_CalculateMicrostateEnergy (self, vector, pH, temperature);
    GdeltaRT = (Gnew - G) * beta;

    /* Prepare */
    ran = rand () / (Real) RAND_MAX;

    /* Apply the Metropolis criterion; based on GMCT */
    if (GdeltaRT < 0.)
      accept = True;
    else {
      if (-GdeltaRT < TOO_SMALL)
        accept = False;
      else if (ran < exp (-GdeltaRT))
        accept = True;
      else
        accept = False; 
    }

    /* Accept or reject the move? */
    if (accept == True)
      G = Gnew;
    else
      vector->vector[site] = instanceBefore;
  }
  /* Return the last accepted energy (only for info) */
  return G;
}

void EnergyModel_CalculateProbabilitiesMonteCarlo (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature, const Boolean equil, Integer nscans, Status *status) {
  Real      Gfinal, scale;
  Real     *counter;
  Integer  *activeInstance;
  Integer   site, nmoves;

  /* The number of moves is proportional to the number of sites */
  nmoves  = vector->nsites;
  scale   = 1. / nscans;

  /* Equilibration phase? */
  if (equil) {
    StateVector_Randomize (vector);

    /* Run the scans */
    for (; nscans >= 0; nscans--) {
      Gfinal = EnergyModel_MCScan (self, vector, pH, temperature, nmoves);
    }
  }

  /* Production phase? */
  else {
    /* Reset the probabilities */
    Real1DArray_Set (self->probabilities, 0.);

    /* Run the scans */
    for (; nscans >= 0; nscans--) {
      Gfinal = EnergyModel_MCScan (self, vector, pH, temperature, nmoves);

      /* Update the counts */
      for (site = 0, activeInstance = vector->vector; site < vector->nsites; site++, activeInstance++) {
        counter = Real1DArray_ItemPointer (self->probabilities, *activeInstance);
        (*counter)++;
      }
    }
    /* Average the probabilities */
    Real1DArray_Scale (self->probabilities, scale);
  }
}

/* Function to be called from Cython */
void EnergyModel_UpdateProbabilities (const EnergyModel *self, const StateVector *vector) {
  Integer site, *activeInstance;
  Real *counter;

  for (site = 0, activeInstance = vector->vector; site < vector->nsites; site++, activeInstance++) {
    counter = Real1DArray_ItemPointer (self->probabilities, *activeInstance);
    (*counter)++;
  }
}
