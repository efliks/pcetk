/*------------------------------------------------------------------------------
! . File      : StateVector.h
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#ifndef _STATEVECTOR
#define _STATEVECTOR

#include "Boolean.h"
#include "Integer.h"
#include "Real.h"
#include "Memory.h"
#include "Status.h"
#include "Integer1DArray.h"
#include "Real1DArray.h"
#include "Real2DArray.h"


#define CONSTANT_MOLAR_GAS_KCAL_MOL  0.001987165392
#define CONSTANT_LN10                2.302585092994


typedef struct {
  Integer *vector, *minvector, *maxvector, *substate;
  Integer length, slength;
} StateVector;

/*
extern StateVector *StateVector_Clone             (const StateVector *self);
extern void         StateVector_CopyTo            (const StateVector *self, StateVector *other);
*/

extern StateVector *StateVector_Allocate          (const Integer length);
extern void         StateVector_Deallocate        (      StateVector *self);
extern void         StateVector_Reset             (const StateVector *self);
extern void         StateVector_ResetToMaximum    (const StateVector *self);
extern Integer      StateVector_GetItem           (const StateVector *self, const Integer index);
extern Boolean      StateVector_SetItem           (const StateVector *self, const Integer index, const Integer value);
extern Integer      StateVector_GetActualItem     (const StateVector *self, const Integer index);
extern Boolean      StateVector_SetActualItem     (const StateVector *self, const Integer index, const Integer value);
extern Boolean      StateVector_Increment         (const StateVector *self);

/* Substate-related functions */
extern void         StateVector_ResetSubstate     (const StateVector *self);
extern Boolean      StateVector_AllocateSubstate  (      StateVector *self, const Integer nsites);
extern Boolean      StateVector_IncrementSubstate (const StateVector *self);
extern Boolean      StateVector_SetSubstateItem   (const StateVector *self, const Integer selectedSiteIndex, const Integer index);
extern Integer      StateVector_GetSubstateItem   (const StateVector *self, const Integer index);

/* Calculating microstate energy */
extern Real StateVector_CalculateMicrostateEnergy (const StateVector *self, const Integer1DArray *protons, const Real1DArray *intrinsic, const Real2DArray *interactions, const Real pH, const Real temperature);

/* Calculating probabilities of protonation states analytically */
extern Boolean StateVector_CalculateProbabilitiesAnalytically (const StateVector *self, const Integer1DArray *protons, const Real1DArray *intrinsic, const Real2DArray *interactions, const Real pH, const Real temperature, const Integer nstates, Real1DArray *probabilities);

#endif
