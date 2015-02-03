/*------------------------------------------------------------------------------
! . File      : StateVector.h
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#ifndef _STATEVECTOR
#define _STATEVECTOR

#include <stdlib.h>   /* Needed for calloc */
#include <stdio.h>    /* Needed for random */
#include <time.h>     /* Needed for random seed */

#include "Real.h"
#include "Boolean.h"
#include "Integer.h"
#include "Cardinal.h"
#include "Memory.h"
#include "Status.h"


#define MEMORY_ALLOCATEARRAY_POINTERS( object, length, type ) { object = (type**) calloc ((CSize) length, sizeof (type*)) ; }

typedef struct {
  /* Site belongs to a substate */
  Boolean isSubstate;
  /* Index of the site itself */
  Integer indexSite;
  /* Global indices of instances of the site */
  Integer indexActive, indexFirst, indexLast;
} TitrSite;

typedef struct {
  /* Pointers to sites that make up a pair */
  TitrSite *a, *b;
} PairSite;

typedef struct {
  TitrSite  *sites, **substateSites;
  Integer    nsites, nssites;
  /* Handled by the EnergyModel module */
  PairSite  *pairs;
  Integer    npairs;
} StateVector;


/* Allocation and deallocation */
extern StateVector *StateVector_Allocate          (const Integer nsites, Status *status);
extern void         StateVector_AllocateSubstate  (      StateVector *self, const Integer nssites, Status *status);
extern void         StateVector_AllocatePairs     (      StateVector *self, const Integer npairs, Status *status);
extern void         StateVector_Deallocate        (      StateVector *self);

/* Copying and cloning */
extern StateVector *StateVector_Clone             (const StateVector *self, Status *status);
extern void         StateVector_CopyTo            (const StateVector *self, StateVector *other, Status *status);

/* Functions for setting all items at once */
extern void         StateVector_Reset             (const StateVector *self);
extern void         StateVector_ResetSubstate     (const StateVector *self);
extern void         StateVector_ResetToMaximum    (const StateVector *self);
extern void         StateVector_Randomize         (const StateVector *self);

/* Functions for accessing items */
extern void         StateVector_SetSite           (const StateVector *self, const Integer indexSite, const Integer indexFirst, const Integer indexLast, Status *status);
extern void         StateVector_SetPair           (const StateVector *self, const Integer indexPair, const Integer indexFirstSite, const Integer indexSecondSite, Status *status);
extern void         StateVector_GetPair           (const StateVector *self, const Integer indexPair, Integer *indexFirstSite, Integer *indexSecondSite, Status *status);
extern Boolean      StateVector_IsSubstate        (const StateVector *self, const Integer siteIndex, Status *status);
extern Integer      StateVector_GetItem           (const StateVector *self, const Integer siteIndex, Status *status);
extern void         StateVector_SetItem           (const StateVector *self, const Integer siteIndex, const Integer instanceLocalIndex, Status *status);
extern Integer      StateVector_GetActualItem     (const StateVector *self, const Integer siteIndex, Status *status);
extern void         StateVector_SetActualItem     (const StateVector *self, const Integer siteIndex, const Integer instanceGlobalIndex, Status *status);
extern Integer      StateVector_GetSubstateItem   (const StateVector *self, const Integer index, Status *status);
extern void         StateVector_SetSubstateItem   (const StateVector *self, const Integer selectedSiteIndex, const Integer index, Status *status);

/* Incrementation */
extern Boolean      StateVector_Increment         (const StateVector *self);
extern Boolean      StateVector_IncrementSubstate (const StateVector *self);

/* Monte Carlo-related functions */
extern void         StateVector_Move              (const StateVector *self, Integer *siteIndex, Integer *oldIndexActive);

#endif
