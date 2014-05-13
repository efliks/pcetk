/*------------------------------------------------------------------------------
! . File      : StateVector.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#ifndef _STATEVECTOR
#define _STATEVECTOR

#include "Boolean.h"
#include "Integer.h"
#include "Memory.h"
#include "Status.h"


typedef struct {
  Integer *vector;
  Integer *maxvector;
  Integer length;
} StateVector;


//extern StateVector *StateVector_Allocate (const Integer length, Status *status);
extern StateVector *StateVector_Allocate       (const Integer length);
extern void         StateVector_Deallocate     (      StateVector *self);
extern void         StateVector_Reset          (const StateVector *self);
extern void         StateVector_ResetToMaximum (const StateVector *self);
extern Integer      StateVector_GetItem        (const StateVector *self, const Integer index);
extern Boolean      StateVector_SetItem        (const StateVector *self, const Integer index, const Integer value);
extern Boolean      StateVector_Increment      (const StateVector *self);

#endif
