/*------------------------------------------------------------------------------
! . File      : StateVector.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "StateVector.h"


/*StateVector *StateVector_Allocate (const Integer length, Status *status) {*/
StateVector *StateVector_Allocate (const Integer length) {
  StateVector *self = NULL;

  MEMORY_ALLOCATE (self, StateVector);
  if (self != NULL) {
    self->vector    = NULL;
    self->maxvector = NULL;
    self->length    = length;

    if (length > 0) {
      MEMORY_ALLOCATEARRAY (self->vector, length, Integer);
      if (self->vector == NULL) {
        MEMORY_DEALLOCATE (self);
      }

      MEMORY_ALLOCATEARRAY (self->maxvector, length, Integer);
      if (self->maxvector == NULL) {
        MEMORY_DEALLOCATE (self->vector);
        MEMORY_DEALLOCATE (self);
      }
    }
  }
  return self;
}

void StateVector_Deallocate (StateVector *self) {
  if (self != NULL) {
    MEMORY_DEALLOCATE (self->maxvector);
    MEMORY_DEALLOCATE (self->vector);
    MEMORY_DEALLOCATE (self);
  }
}

void StateVector_Reset (const StateVector *self) {
  Integer   i;
  Integer   *v = self->vector;

  for (i = 0; i < self->length; i++, v++) {
    *v = 0;
  }
}

void StateVector_ResetToMaximum (const StateVector *self) {
  Integer   i;
  Integer   *v = self->vector, *m = self->maxvector;

  for (i = 0; i < self->length; i++, v++, m++) {
    *v = *m;
  }
}

Integer StateVector_GetItem (const StateVector *self, const Integer index) {
  if (index < 0 || index > (self->length - 1)) {
    return -1000;
  }
  else {
    return self->vector[index];
  }
}

Boolean StateVector_SetItem (const StateVector *self, const Integer index, const Integer value) {
  if (index < 0 || index > (self->length - 1)) {
    return False;
  }
  else if (value < 0 || value > self->maxvector[index]) {
    return False;
  }
  else {
    self->vector[index] = value;
    return True;
  }
}

/* Incrementation algorithm by Timm Essigke */
Boolean StateVector_Increment (const StateVector *self) {
  Integer i;
  Integer *v = self->vector, *m = self->maxvector;
/*
  This prevents zeroing the state vector after the last iteration
  incr = False;
  for (i = 0; i < self->length; i++, v++, m++) {
    if ((*v) < (*m)) {
      incr = True;
      break;
    }
  }
  if (incr == False) {
    return False;
  }
*/

  for (i = 0; i < self->length; i++, v++, m++) {
    if ((*v) < (*m)) {
      (*v)++;
      return True;
    }
    else {
      *v = 0;
    }
  }
  return False;
}

/*
Real Microstate_Energy (Real *gintr, Integer *protons, Real **interactions, StateVector *vector, Real pH, Real temperature) {
  Integer siteIndex, siteIndexInner;
  Integer *v = vector->vector;
  Real    W = 0.0, Gintr = 0.0;

  for (siteIndex = 0; siteIndex < vector->lenght; siteIndex++, v++) {
    instanceIndex = *v;

    for (siteIndexInner = 0; siteIndexInner < vector->lenght; siteIndexInner++) {

    }
  }

  protonChemicalPotential = -CONSTANT_MOLAR_GAS_KCAL_MOL * temperature * CONSTANT_LN10 * pH
  return Gintr - nprotons * protonChemicalPotential + 0.5 * W;
}
*/
