/*------------------------------------------------------------------------------
! . File      : StateVector.h
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#ifndef _STATEVECTOR
#define _STATEVECTOR

/* memcpy comes from here */
#include <stdlib.h>
/* Needed for memcpy, otherwise warnings */
#include <string.h>
/* Needed for random seed */
#include <time.h>
/* Needed for random */
#include <stdio.h>

#include "Boolean.h"
#include "Integer.h"
#include "Real.h"
#include "Memory.h"
#include "Status.h"


typedef struct {
  Integer  *vector     ;
  Integer  *minvector  ;
  Integer  *maxvector  ;
  Integer  *substate   ;
  Integer   length     ;
  Integer   slength    ;
} StateVector;

extern StateVector *StateVector_Allocate          (const Integer length, Status *status);
extern StateVector *StateVector_Clone             (const StateVector *self, Status *status);
extern Boolean      StateVector_CopyTo            (const StateVector *self, StateVector *other, Status *status);
extern void         StateVector_Deallocate        (      StateVector *self);
extern void         StateVector_Reset             (const StateVector *self);
extern void         StateVector_ResetToMaximum    (const StateVector *self);
extern Integer      StateVector_GetItem           (const StateVector *self, const Integer index, Status *status);
extern Boolean      StateVector_SetItem           (const StateVector *self, const Integer index, const Integer value, Status *status);
extern Integer      StateVector_GetActualItem     (const StateVector *self, const Integer index, Status *status);
extern Boolean      StateVector_SetActualItem     (const StateVector *self, const Integer index, const Integer value, Status *status);
extern Boolean      StateVector_Increment         (const StateVector *self);

/* Substate-related functions */
extern void         StateVector_ResetSubstate     (const StateVector *self);
extern Boolean      StateVector_AllocateSubstate  (      StateVector *self, const Integer nsites, Status *status);
extern Boolean      StateVector_IncrementSubstate (const StateVector *self);
extern Boolean      StateVector_SetSubstateItem   (const StateVector *self, const Integer selectedSiteIndex, const Integer index, Status *status);
extern Integer      StateVector_GetSubstateItem   (const StateVector *self, const Integer index, Status *status);

/* Monte Carlo functions */
extern void         StateVector_Move              (const StateVector *self);
extern void         StateVector_Randomize         (const StateVector *self);

#endif
