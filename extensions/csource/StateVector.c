/*------------------------------------------------------------------------------
! . File      : StateVector.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "StateVector.h"


/*=============================================================================
  Allocation, deallocation, copying, etc.
=============================================================================*/
StateVector *StateVector_Allocate (const Integer nsites, Status *status) {
  StateVector *self = NULL;

  MEMORY_ALLOCATE (self, StateVector);
  if (self != NULL) {
    self->vector    = NULL;
    self->minvector = NULL;
    self->maxvector = NULL;
    self->substate  = NULL;
    self->nsites    = nsites;
    self->nssites   = 0;

    if (nsites > 0) {
      MEMORY_ALLOCATEARRAY (self->vector, nsites, Integer);
      if (self->vector == NULL) {
        MEMORY_DEALLOCATE (self);
      }
      else {
        MEMORY_ALLOCATEARRAY (self->minvector, nsites, Integer);
        if (self->minvector == NULL) {
          MEMORY_DEALLOCATE (self->vector);
          MEMORY_DEALLOCATE (self);
        }
        else {
          MEMORY_ALLOCATEARRAY (self->maxvector, nsites, Integer);
          if (self->maxvector == NULL) {
            MEMORY_DEALLOCATE (self->minvector);
            MEMORY_DEALLOCATE (self->vector);
            MEMORY_DEALLOCATE (self);
          }
        }
      }
    }
  }

  if (self == NULL) {
    Status_Set (status, Status_MemoryAllocationFailure);
  }
  return self;
}

void StateVector_Deallocate (StateVector *self) {
  if (self != NULL) {
    MEMORY_DEALLOCATE (self->maxvector);
    MEMORY_DEALLOCATE (self->minvector);
    MEMORY_DEALLOCATE (self->vector);

    /* Deallocate substate, if exists */
    if (self->substate != NULL) {
      MEMORY_DEALLOCATE (self->substate);
    }

    MEMORY_DEALLOCATE (self);
  }
}

StateVector *StateVector_Clone (const StateVector *self, Status *status) {
  StateVector *clone = NULL;

  if (self != NULL) {
    clone = StateVector_Allocate (self->nsites, status);
    if (*status == Status_Continue) {
      StateVector_CopyTo (self, clone, status);
    }
  }

  return clone;
}

Boolean StateVector_CopyTo (const StateVector *self, StateVector *other, Status *status) {
  /* Check for different nsitess */
  if (self->nsites != other->nsites) {
    StateVector_Deallocate (other);
    other = StateVector_Allocate (self->nsites, status);

    if (*status != Status_Continue) {
      return False;
    }
  }

  /* Copy */
  memcpy (other->vector    , self->vector    , other->nsites * sizeof (Integer));
  memcpy (other->minvector , self->minvector , other->nsites * sizeof (Integer));
  memcpy (other->maxvector , self->maxvector , other->nsites * sizeof (Integer));

  /* Copy substate? */
  if (self->substate != NULL) {
    StateVector_AllocateSubstate (other, self->nssites, status);
    if (*status != Status_Continue) {
      return False;
    }
    memcpy (other->substate, self->substate, other->nssites * sizeof (Integer));
  }

  return True;
}

/*=============================================================================
  Functions relatated to items
=============================================================================*/
void StateVector_Reset (const StateVector *self) {
  Integer   i;
  Integer   *v = self->vector, *m = self->minvector;
  for (i = 0; i < self->nsites; i++, v++, m++) {
    *v = *m;
  }
}

void StateVector_ResetToMaximum (const StateVector *self) {
  Integer   i;
  Integer   *v = self->vector, *m = self->maxvector;
  for (i = 0; i < self->nsites; i++, v++, m++) {
    *v = *m;
  }
}

/*-----------------------------------------------------------------------------
 Get the local index of an instance of a site, usually 0 and 1 for 
 most sites or 0, 1, 2, 3 for histidines
-----------------------------------------------------------------------------*/
Integer StateVector_GetItem (const StateVector *self, const Integer index, Status *status) {
  if (index < 0 || index > (self->nsites - 1)) {
    Status_Set (status, Status_IndexOutOfRange);
    return -1;
  }
  else {
    return (self->vector[index] - self->minvector[index]);
  }
}

Boolean StateVector_SetItem (const StateVector *self, const Integer index, const Integer value, Status *status) {
  Integer valueActual;

  if (index < 0 || index > (self->nsites - 1)) {
    Status_Set (status, Status_IndexOutOfRange);
    return False;
  }
  else {
    valueActual = value + self->minvector[index];
    if (valueActual < self->minvector[index] || valueActual > self->maxvector[index]) {
      Status_Set (status, Status_ValueError);
      return False;
    }
    else {
      self->vector[index] = valueActual;
      return True;
    }
  }
}

/*-----------------------------------------------------------------------------
 Get the actual content of the state vector, i.e. global index of 
 an instance in the central arrays (_protons, _intrinsic, _interactions)
-----------------------------------------------------------------------------*/
Integer StateVector_GetActualItem (const StateVector *self, const Integer index, Status *status) {
  if (index < 0 || index > (self->nsites - 1)) {
    Status_Set (status, Status_IndexOutOfRange);
    return -1;
  }
  else {
    return self->vector[index];
  }
}

Boolean StateVector_SetActualItem (const StateVector *self, const Integer index, const Integer value, Status *status) {
  if (index < 0 || index > (self->nsites - 1)) {
    Status_Set (status, Status_IndexOutOfRange);
    return False;
  }
  else if (value < self->minvector[index] || value > self->maxvector[index]) {
    Status_Set (status, Status_ValueError);
    return False;
  }
  else {
    self->vector[index] = value;
    return True;
  }
}

/*-----------------------------------------------------------------------------
 Incrementation algorithm by Timm Essigke 

 One could write a code that prevents zeroing the state vector after 
 the last iteration
-----------------------------------------------------------------------------*/
Boolean StateVector_Increment (const StateVector *self) {
  Integer i;
  Integer *v = self->vector, *minv = self->minvector, *maxv = self->maxvector;

  for (i = 0; i < self->nsites; i++, v++, minv++, maxv++) {
    if ((*v) < (*maxv)) {
      (*v)++;
      return True;
    }
    else {
      *v = *minv;
    }
  }
  return False;
}

/*=============================================================================
  Functions related to substate
=============================================================================*/
Boolean StateVector_AllocateSubstate (StateVector *self, const Integer nsites, Status *status) {
  if (self->substate != NULL) {
    /* Substate already allocated */
    Status_Set (status, Status_MemoryAllocationFailure);
    return False;
  }
  else {
    MEMORY_ALLOCATEARRAY (self->substate, nsites, Integer);
    if (self->substate == NULL) {
      /* Substate allocation failed */
      Status_Set (status, Status_MemoryAllocationFailure);
      return False;
    }
    self->nssites = nsites;
    return True;
  }
}

/*-----------------------------------------------------------------------------
 |selectedSiteIndex| is an index of the site to increment within the substate

 |index| is an index in the substate's array of selectedSiteIndices
-----------------------------------------------------------------------------*/
Boolean StateVector_SetSubstateItem (const StateVector *self, const Integer selectedSiteIndex, const Integer index, Status *status) {
  if (index < 0 || index > (self->nssites - 1)) {
    Status_Set (status, Status_IndexOutOfRange);
    return False;
  }
  else if (selectedSiteIndex < 0 || (selectedSiteIndex > self->nsites - 1)) {
    Status_Set (status, Status_ValueError);
    return False;
  }
  else {
    self->substate[index] = selectedSiteIndex;
    return True;
  }
}

Integer StateVector_GetSubstateItem (const StateVector *self, const Integer index, Status *status) {
  if (index < 0 || index > (self->nssites - 1)) {
    Status_Set (status, Status_IndexOutOfRange);
    return -1;
  }
  else {
    return self->substate[index];
  }
}

void StateVector_ResetSubstate (const StateVector *self) {
  Integer i;
  Integer *siteIndex = self->substate;

  if (self->substate != NULL) {
    for (i = 0; i < self->nssites; i++, siteIndex++) {
      self->vector[*siteIndex] = self->minvector[*siteIndex];
    }
  }
}

Boolean StateVector_IncrementSubstate (const StateVector *self) {
  Integer i, site, maxsite;
  Integer *siteIndex = self->substate;

  if (self->substate != NULL) {
    for (i = 0; i < self->nssites; i++, siteIndex++) {
      site    = self->vector    [*siteIndex];
      maxsite = self->maxvector [*siteIndex];
  
      if (site < maxsite) {
        site++;
        self->vector[*siteIndex] = site;
        return True;
      }
      else {
        self->vector[*siteIndex] = self->minvector[*siteIndex];
      }
    }
    return False;
  }
  else {
    return False;
  }
}

/*=============================================================================
  Monte Carlo-related functions
=============================================================================*/
void StateVector_Move (const StateVector *self, Integer *site, Integer *before, Integer *after) {
  Integer component, value;
  /* Choose component and save its current value */
  component = rand () % self->nsites;
  *site     = component;
  *before   = self->vector[component];
  /* Choose a new value */
  value     = rand () % (self->maxvector[component] - self->minvector[component]) + self->minvector[component];

  /* Why is this part necessary? */
  if (value == self->vector[component]) {
    value++;
    if (value > self->maxvector[component])
      value = self->minvector[component];
  }

  /* Set the new value */
  self->vector[component] = value;
  *after = value;
}

/* Generate a random vector */
void StateVector_Randomize (const StateVector *self) {
  Integer i;
  static Boolean first = True;
  if (first) {
    srandom ((unsigned int) time (NULL));
    first = False;
  }

  for (i = 0; i < self->nsites; i++) {
    self->vector[i] = rand () % (self->maxvector[i] - self->minvector[i]) + self->minvector[i];
  }
}
