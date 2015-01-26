/*------------------------------------------------------------------------------
! . File      : EnergyModel.h
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#ifndef _ENERGYMODEL
#define _ENERGYMODEL

/* Needed for memset and random */
#include <stdio.h>
/* Needed for exp */
#include <math.h>

#include "Boolean.h"
#include "Integer.h"
#include "Real.h"
#include "Memory.h"
#include "Status.h"
#include "Real1DArray.h"
#include "Real2DArray.h"
#include "Integer1DArray.h"
#include "SymmetricMatrix.h"
#include "StateVector.h"


#define CONSTANT_MOLAR_GAS_KCAL_MOL  0.001987165392
#define CONSTANT_LN10                2.302585092994
/* Taken from GMCT */
#define TOO_SMALL                   -500.0

typedef struct {
  Integer1DArray    *protons          ;
  Real1DArray       *intrinsic        ;
  Real2DArray       *interactions     ;
  SymmetricMatrix   *symmetricmatrix  ;
  Real1DArray       *probabilities    ;
  Integer            nstates          ;
  Integer            ninstances       ;
} EnergyModel;


/* Allocation and deallocation */
extern EnergyModel *EnergyModel_Allocate (const Integer ninstances, Status *status);
extern void         EnergyModel_Deallocate (EnergyModel *self);

/* Handling of the interaction matrix */
extern void         EnergyModel_SymmetrizeInteractions     (const EnergyModel *self, Status *status);
extern Boolean      EnergyModel_CheckInteractionsSymmetric (const EnergyModel *self, Real tolerance, Real *maxDeviation);

/* Calculation functions */
extern Real         EnergyModel_CalculateMicrostateEnergy          (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature);
extern void         EnergyModel_CalculateProbabilitiesAnalytically (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature, Status *status);

/* Monte Carlo functions */
extern void         EnergyModel_CalculateProbabilitiesMonteCarlo   (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature, const Boolean equil, Integer nscans, Status *status);
extern Real         EnergyModel_MCScan                             (const EnergyModel *self, const StateVector *vector, const Real pH, const Real temperature, Integer nmoves);
extern void         EnergyModel_UpdateProbabilities                (const EnergyModel *self, const StateVector *vector);
extern Boolean      EnergyModel_Metropolis                         (const Real GdeltaRT);

/* Functions for accessing items */
extern Real         EnergyModel_GetGintr                (const EnergyModel *self, const Integer instIndexGlobal);
extern Integer      EnergyModel_GetProtons              (const EnergyModel *self, const Integer instIndexGlobal);
extern Real         EnergyModel_GetProbability          (const EnergyModel *self, const Integer instIndexGlobal);
extern Real         EnergyModel_GetInteraction          (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);
extern Real         EnergyModel_GetInteractionSymmetric (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);
extern Real         EnergyModel_GetDeviation            (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);

extern void         EnergyModel_SetGintr                (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void         EnergyModel_SetProtons              (const EnergyModel *self, const Integer instIndexGlobal, const Integer value);
extern void         EnergyModel_SetProbability          (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void         EnergyModel_SetInteraction          (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB, const Real value);

#endif
