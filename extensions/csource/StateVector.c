/*------------------------------------------------------------------------------
! . File      : StateVector.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/

# include "StateVector.h"


//StateVector *StateVector_Allocate (const Integer length, Status *status) {
StateVector *StateVector_Allocate (const Integer length) {
  StateVector *self = NULL;
  MEMORY_ALLOCATE (self, StateVector);

  //TO DO: Check if memory allocation was successfull, use Status_Set

  if (self != NULL) {
    self->vector    = NULL;
    self->maxvector = NULL;
    self->length    = length;

    if (length > 0) {
      MEMORY_ALLOCATEARRAY (self->vector, length, Integer);
      //if (self->vector == NULL) StateVector_Deallocate (&self);

      MEMORY_ALLOCATEARRAY (self->maxvector, length, Integer);
      //if (self->maxvector == NULL) StateVector_Deallocate (&self);
    }
  }

  return self;
}


// Use **self ???
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


Boolean StateVector_Increment (const StateVector *self) {
/*
  Integer   i;
  Integer   value, maxvalue;

  for (i = 0; i < self->length; i++) {
    value    = self->vector[i];
    maxvalue = self->maxvector[i];

    if (value < maxvalue) {
      value++;
      self->vector[i] = value;
      return True;
    }
    else {
      self->vector[i] = 0;
    }
  }
  return False;
*/

/*
for index, (value, maxvalue) in enumerate (zip (self.vector, self.maxvector)):
  if value < maxvalue:
    value = value + 1
    self.vector[index] = value
    return True
  else:
    value = 0
    self.vector[index] = value
return False
*/
  Integer i;
  Integer *v = self->vector, *m = self->maxvector;

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
