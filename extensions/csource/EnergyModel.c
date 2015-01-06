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

/* Getters */
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

/* Setters */
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


Boolean EnergyModel_CheckInteractionsSymmetric (const EnergyModel *self, const Real threshold, Real *maxDeviate) {
  Integer row, column;
  Real symmetry, itemRow, itemColumn;
  Real deviation, absoluteDeviation, maxDeviation = 0.;
  Boolean isSymmetric = True;
  /* 
  for (row = 0; row < self->ninstances; row++) {
    for (column = 0; column < self->ninstances; column++) {
    }
  }
  */

  for (row = 0; row < Real2DArray_Rows (self->interactions); row++) {
    for (column = 0; column < Real2DArray_Columns (self->interactions); column++) {
      itemRow    = Real2DArray_Item (self->interactions, row,    column);
      itemColumn = Real2DArray_Item (self->interactions, column, row);

      symmetry  = .5 * (itemRow + itemColumn);
      deviation = symmetry - itemRow;
      if (deviation < 0.) {
        absoluteDeviation = -1. * deviation;
      }
      else {
        absoluteDeviation = deviation;
      }
      Real2DArray_Item (self->deviations, row, column) = absoluteDeviation;

      if (absoluteDeviation > threshold) {
        isSymmetric = False;
      }

      if (absoluteDeviation > maxDeviation) {
        maxDeviation = absoluteDeviation;
      }
    }
  }

  *maxDeviate = maxDeviation;
  return isSymmetric;

  /*
    # define SYMMETRIC_TOLERANCE 1.0e-10
    Boolean Real2DArray_IsSymmetric ( const Real2DArray *self, const Real *tolerance, Real *deviation )
  */
}

Boolean EnergyModel_SymmetrizeInteractions (const EnergyModel *self, Status *status) {
  Boolean success;
  SymmetricMatrix_CopyFromReal2DArray (self->symmetricmatrix, self->interactions, status);

  success = (*status != Status_Continue)? False : True;
  return success;
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

  for (siteIndex = 0, instanceIndex = vector->vector; siteIndex < vector->length; siteIndex++, instanceIndex++) {
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
Boolean EnergyModel_CalculateProbabilitiesAnalytically (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature, Status *status) {
  Real1DArray *bfactors;
  Real        *bfactor;
  Real         energy, energyZero, bsum;
  Integer     *activeInstanceGlobalIndex;
  Integer      stateIndex, siteIndex;

  bfactors = Real1DArray_Allocate (self->nstates, status);
  if (*status != Status_Continue) {
    return False;
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
    for (siteIndex = 0, activeInstanceGlobalIndex = vector->vector; siteIndex < vector->length; siteIndex++, activeInstanceGlobalIndex++) {
      self->probabilities->data[*activeInstanceGlobalIndex] += *bfactor;
    }
    StateVector_Increment (vector);
  }
  bsum = Real1DArray_Sum (bfactors);
  Real1DArray_Scale (self->probabilities, 1. / bsum);

  Real1DArray_Deallocate (&bfactors);

  return True;
}
