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

    self->substate  = NULL;
    self->slength   = 0;

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

Boolean StateVector_AllocateSubstate (StateVector *self, const Integer nsites) {
  if (self->substate != NULL) {
    /* Substate already allocated */
    return False;
  }
  else {
    MEMORY_ALLOCATEARRAY (self->substate, nsites, Integer);
    if (self->substate == NULL) {
      /* Substate allocation failed */
      return False;
    }
    self->slength = nsites;
    return True;
  }
}

void StateVector_Deallocate (StateVector *self) {
  if (self != NULL) {
    MEMORY_DEALLOCATE (self->maxvector);
    MEMORY_DEALLOCATE (self->vector);

    /* Deallocate substate, if exists */
    if (self->substate != NULL) {
      MEMORY_DEALLOCATE (self->substate);
    }

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
    return -1;
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
    |selectedSiteIndex| is an index of the site to increment within the substate

    |index| is an index in the substate's array of selectedSiteIndices
*/
Boolean StateVector_SetSubstateItem (const StateVector *self, const Integer selectedSiteIndex, const Integer index) {
  if (index < 0 || index > (self->slength - 1)) {
    return False;
  }
  else if (selectedSiteIndex < 0 || (selectedSiteIndex > self->length - 1)) {
    return False;
  }
  else {
    self->substate[index] = selectedSiteIndex;
    return True;
  }
}

Integer StateVector_GetSubstateItem (const StateVector *self, const Integer index) {
  if (index < 0 || index > (self->slength - 1)) {
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
    for (i = 0; i < self->slength; i++, siteIndex++) {
      self->vector[*siteIndex] = 0;
    }
  }
}

/* Maybe this is not the fastest solution */
Boolean StateVector_IncrementSubstate (const StateVector *self) {
  Integer i, site, maxsite;
  Integer *siteIndex = self->substate;

  if (self->substate != NULL) {
    for (i = 0; i < self->slength; i++, siteIndex++) {
      site    = self->vector    [*siteIndex];
      maxsite = self->maxvector [*siteIndex];
  
      if (site < maxsite) {
        site++;
        self->vector[*siteIndex] = site;
        return True;
      }
      else {
        self->vector[*siteIndex] = 0;
      }
    }
    return False;
  }
  else {
    return False;
  }
}
