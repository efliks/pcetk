/*------------------------------------------------------------------------------
! . File      : EnergyModel.h
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#ifndef _ENERGYMODEL
#define _ENERGYMODEL

/* Data types */
#include "Real.h"
#include "Boolean.h"
#include "Integer.h"
#include "Cardinal.h"

/* Arrays */
#include "Real1DArray.h"
#include "Real2DArray.h"
#include "Integer1DArray.h"
#include "SymmetricMatrix.h"

/* Other */
#include "Memory.h"
#include "Status.h"

/* Own modules */
#include "StateVector.h"


#define CONSTANT_MOLAR_GAS_KCAL_MOL  0.001987165392
#define CONSTANT_LN10                2.302585092994

/* Macros */
#define EnergyModel_RowPointer(self, i) (&self->symmetricmatrix->data[(i * (i + 1) >> 1)])

#define EnergyModel_GetW(self, i, j) (i >= j ? self->symmetricmatrix->data[(i * (i + 1) >> 1) + j] : self->symmetricmatrix->data[(j * (j + 1) >> 1) + i])

typedef struct {
    /* Number of bound protons of each instance */
    Integer1DArray   *protons;
    /* Gmodel of each instance (needed for energies of unfolded proteins) */
    Real1DArray      *models;
    /* Gintr of each instance */
    Real1DArray      *intrinsic;
    /* Interactions between instances before symmetrization */
    Real2DArray      *interactions;
    /* Symmetrized interactions */
    SymmetricMatrix  *symmetricmatrix;
    /* Probability of occurrence of each instance */
    Real1DArray      *probabilities;
    /* Private state vector of the energy model */
    StateVector      *vector;
    /* Total number of possible protonation states, no greater than ANALYTIC_STATES */
    Integer           nstates;
    /* Total number of instances */
    Integer           ninstances;
    /* Temperature at which the MEAD part was done */
    Real              temperature;
} EnergyModel;


/* Allocation and deallocation */
extern EnergyModel *EnergyModel_Allocate (const Integer nsites, const Integer ninstances, Status *status);
extern void         EnergyModel_Deallocate (EnergyModel *self);

/* Miscellaneous functions */
extern void    EnergyModel_SymmetrizeInteractions       (const EnergyModel *self, Status *status);
extern Boolean EnergyModel_CheckInteractionsSymmetric   (const EnergyModel *self, Real tolerance, Real *maxDeviation);
extern void    EnergyModel_ResetInteractions            (const EnergyModel *self);
extern void    EnergyModel_StateVectorFromProbabilities (const EnergyModel *self, StateVector *vector, Status *status);

/* Calculation of energies */
extern Real EnergyModel_CalculateMicrostateEnergy         (const EnergyModel *self, const StateVector *vector, const Real pH);
extern Real EnergyModel_CalculateMicrostateEnergyUnfolded (const EnergyModel *self, const StateVector *vector, const Real pH);

/* Calculation of partition functions */
extern Real EnergyModel_CalculateZ (const EnergyModel *self, Real (*EnergyFunction)(const EnergyModel*, const StateVector*, const Real), const Real pH, const Real Gzero, Real1DArray *bfactors);
extern Real EnergyModel_CalculateZunfolded (const EnergyModel *self, const Real pH, const Real Gzero, Status *status);
extern Real EnergyModel_CalculateZfolded (const EnergyModel *self, const Real pH, const Real Gzero, Status *status);

/* Calculation of probabilities */
extern void EnergyModel_CalculateProbabilitiesFromZ (const EnergyModel *self, const Real Z, const Real1DArray *bfactors);
extern void EnergyModel_CalculateProbabilitiesAnalytically (const EnergyModel *self, const Real pH, Status *status);
extern void EnergyModel_CalculateProbabilitiesAnalyticallyUnfolded (const EnergyModel *self, const Real pH, Status *status);

/* Functions for accessing items */
extern Real    EnergyModel_GetGmodel         (const EnergyModel *self, const Integer instIndexGlobal);
extern Real    EnergyModel_GetGintr          (const EnergyModel *self, const Integer instIndexGlobal);
extern Integer EnergyModel_GetProtons        (const EnergyModel *self, const Integer instIndexGlobal);
extern Real    EnergyModel_GetProbability    (const EnergyModel *self, const Integer instIndexGlobal);
extern Real    EnergyModel_GetInteraction    (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);
extern Real    EnergyModel_GetInterSymmetric (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);
extern Real    EnergyModel_GetDeviation      (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);

extern void    EnergyModel_SetGmodel         (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void    EnergyModel_SetGintr          (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void    EnergyModel_SetProtons        (const EnergyModel *self, const Integer instIndexGlobal, const Integer value);
extern void    EnergyModel_SetProbability    (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void    EnergyModel_SetInteraction    (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB, const Real value);

#endif
