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

#define CONSTANT_MOLAR_GAS_KCAL_MOL  0.001987165392
#define CONSTANT_LN10                2.302585092994


typedef struct {
  Integer *vector, *maxvector, *substate;
  Integer length, slength;
} StateVector;


/*extern StateVector *StateVector_Allocate (const Integer length, Status *status);*/
extern StateVector *StateVector_Allocate          (const Integer length);
extern void         StateVector_Deallocate        (      StateVector *self);
extern void         StateVector_Reset             (const StateVector *self);
extern void         StateVector_ResetToMaximum    (const StateVector *self);
extern Integer      StateVector_GetItem           (const StateVector *self, const Integer index);
extern Boolean      StateVector_SetItem           (const StateVector *self, const Integer index, const Integer value);
extern Boolean      StateVector_Increment         (const StateVector *self);

/* Substate-related functions */
extern Boolean      StateVector_AllocateSubstate  (      StateVector *self, const Integer nsites);
extern Boolean      StateVector_SetSubstateItem   (const StateVector *self, const Integer index, const Integer siteIndex);
extern Boolean      StateVector_IncrementSubstate (const StateVector *self);
extern void         StateVector_ResetSubstate     (const StateVector *self);

#endif
