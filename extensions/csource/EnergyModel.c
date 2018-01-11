/*------------------------------------------------------------------------------
! . File      : EnergyModel.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "EnergyModel.h"

/*
 * Allocate the energy model.
 * Other attributes of the model (nstates, ninstances, temperature) are set from the Cython level.
 */
EnergyModel *EnergyModel_Allocate (const Integer nsites, const Integer ninstances, Status *status) {
    EnergyModel *self = NULL;

    MEMORY_ALLOCATE (self, EnergyModel);
    if (self == NULL)
        goto failSet;

    self->vector           =  NULL  ;
    self->models           =  NULL  ;
    self->protons          =  NULL  ;
    self->intrinsic        =  NULL  ;
    self->interactions     =  NULL  ;
    self->probabilities    =  NULL  ;
    self->symmetricmatrix  =  NULL  ;

    if (nsites > 0) {
        self->vector = StateVector_Allocate (nsites, status);
        if (*status != Status_Continue)
            goto failDealloc;
    }
    if (ninstances > 0) {
        self->models          = Real1DArray_Allocate     (ninstances, status)  ;
        self->protons         = Integer1DArray_Allocate  (ninstances, status)  ;
        self->intrinsic       = Real1DArray_Allocate     (ninstances, status)  ;
        self->interactions    = Real2DArray_Allocate     (ninstances, ninstances, status) ;
        self->probabilities   = Real1DArray_Allocate     (ninstances, status)  ;
        if (*status != Status_Continue)
            goto failDealloc;

        self->symmetricmatrix = SymmetricMatrix_Allocate (ninstances);
        if (self->symmetricmatrix == NULL)
            goto failSetDealloc;
    }
    return self;


failSetDealloc:
    Status_Set (status, Status_MemoryAllocationFailure);

failDealloc:
    EnergyModel_Deallocate (self);
    return NULL;

failSet:
    Status_Set (status, Status_MemoryAllocationFailure);
    return NULL;
}

/*
 * Deallocate the energy model.
 */
void EnergyModel_Deallocate (EnergyModel *self) {
    if ( self->symmetricmatrix != NULL)  SymmetricMatrix_Deallocate ( &self->symmetricmatrix ) ;
    if ( self->probabilities   != NULL)  Real1DArray_Deallocate     ( &self->probabilities   ) ;
    if ( self->interactions    != NULL)  Real2DArray_Deallocate     ( &self->interactions    ) ;
    if ( self->intrinsic       != NULL)  Real1DArray_Deallocate     ( &self->intrinsic       ) ;
    if ( self->protons         != NULL)  Integer1DArray_Deallocate  ( &self->protons         ) ;
    if ( self->models          != NULL)  Real1DArray_Deallocate     ( &self->models          ) ;
    if ( self->vector          != NULL)  StateVector_Deallocate     (  self->vector          ) ;
    if (self != NULL) MEMORY_DEALLOCATE (self);
}

/*
 * Check if the array of interactions is symmetric within the given tolerance (kcal/mol).
 */
Boolean EnergyModel_CheckInteractionsSymmetric (const EnergyModel *self, Real tolerance, Real *maxDeviation) {
    return Real2DArray_IsSymmetric (self->interactions, &tolerance, maxDeviation);
}

/*
 * Symmetrize the array of interactions into a symmetric matrix.
 */
void EnergyModel_SymmetrizeInteractions (const EnergyModel *self, Status *status) {
    SymmetricMatrix_CopyFromReal2DArray (self->symmetricmatrix, self->interactions, status);
}

/*
 * Set all interactions to zero.
 */
void EnergyModel_ResetInteractions (const EnergyModel *self) {
    SymmetricMatrix_Set (self->symmetricmatrix, 0.);
}

/*
 * Scale interactions.
 */
void EnergyModel_ScaleInteractions (const EnergyModel *self, Real scale) {
    SymmetricMatrix_Scale (self->symmetricmatrix, scale);
}

/*
 * Generate the lowest energy state vector.
 * If "vector" is NULL, use the EnergyModel's private vector.
 */
void EnergyModel_StateVectorFromProbabilities (const EnergyModel *self, StateVector *vector, Status *status) {
    TitrSite *ts;
    Integer   i, index, maxi;
    Real      probability, maxp; 

    if (vector != NULL) {
        ts = vector->sites  ;
        i  = vector->nsites ;
        if (i != self->vector->nsites)
            goto fail;
    }
    else {
        ts = self->vector->sites  ;
        i  = self->vector->nsites ;
    }
    for (; i >= 0; i--, ts++) {
        index = ts->indexFirst;
        maxi  = index;
        maxp  = -1.;
        do {
            probability = Real1DArray_Item (self->probabilities, index);
            if (probability > maxp) {
                maxi = index;
                maxp = probability;
            }
        } while (++index <= ts->indexLast);
        ts->indexActive = maxi;
    }
    return;

fail:
    Status_Set (status, Status_ArrayNonConformableSizes);
}

/*
 * Getters.
 */
Real EnergyModel_GetGmodel (const EnergyModel *self, const Integer instIndexGlobal) {
    return Real1DArray_Item (self->models, instIndexGlobal);
}

Real EnergyModel_GetGintr (const EnergyModel *self, const Integer instIndexGlobal) {
    return Real1DArray_Item (self->intrinsic, instIndexGlobal);
}

Integer EnergyModel_GetProtons (const EnergyModel *self, const Integer instIndexGlobal) {
    return Integer1DArray_Item (self->protons, instIndexGlobal);
}

Real EnergyModel_GetProbability (const EnergyModel *self, const Integer instIndexGlobal) {
    return Real1DArray_Item (self->probabilities, instIndexGlobal);
}

Real EnergyModel_GetInteraction (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
    return Real2DArray_Item (self->interactions, instIndexGlobalA, instIndexGlobalB);
}

Real EnergyModel_GetInterSymmetric (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
    return EnergyModel_GetW (self, instIndexGlobalA, instIndexGlobalB);
}

Real EnergyModel_GetDeviation (const EnergyModel *self, const Integer i, const Integer j) {
    Real wij, wji, deviation;

    wij = Real2DArray_Item (self->interactions, i, j);
    wji = Real2DArray_Item (self->interactions, j, i);
    deviation = (wij + wji) * .5 - wij;
    return deviation;
}

/*
 * Setters.
 */
void EnergyModel_SetGmodel (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
    Real1DArray_Item (self->models, instIndexGlobal) = value;
}

void EnergyModel_SetGintr (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
    Real1DArray_Item (self->intrinsic, instIndexGlobal) = value;
}

void EnergyModel_SetProtons (const EnergyModel *self, const Integer instIndexGlobal, const Integer value) {
    Integer1DArray_Item (self->protons, instIndexGlobal) = value;
}

void EnergyModel_SetProbability (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
    Real1DArray_Item (self->probabilities, instIndexGlobal) = value;
}

void EnergyModel_SetInteraction (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB, const Real value) {
    Real2DArray_Item (self->interactions, instIndexGlobalA, instIndexGlobalB) = value;
}

/*
 * Calculate the energy of a microstate defined by the state vector.
 */
Real EnergyModel_CalculateMicrostateEnergy (const EnergyModel *self, const StateVector *vector, const Real pH) {
    Real      Gintr, W, *interact;
    Integer   nprotons, i, j;
    TitrSite *site, *siteInner;

    W        = 0. ;
    Gintr    = 0. ;
    nprotons = 0  ;
    site     = vector->sites;
    for (i = 0; i < vector->nsites; i++, site++) {
        Gintr     +=    Real1DArray_Item (self->intrinsic , site->indexActive);
        nprotons  += Integer1DArray_Item (self->protons   , site->indexActive);

        interact   = EnergyModel_RowPointer (self, site->indexActive);
        siteInner  = vector->sites;
        for (j = 0; j < i; j++, siteInner++) {
            W += *(interact + (siteInner->indexActive));
        }
    }
    return (Gintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature * CONSTANT_LN10 * pH) + W);
}

/*
 * Calculate the energy of a microstate in an unfolded (=denaturated) protein.
 * In the unfolded state, Gintr become Gmodel and all interactions are set to zero.
 * 
 * Reference: Yang A.-S., Honig B., J. Mol. Biol. 1993, 231, 459-474
 */
Real EnergyModel_CalculateMicrostateEnergyUnfolded (const EnergyModel *self, const StateVector *vector, const Real pH) {
    Integer   nprotons, i;
    TitrSite *site;
    Real      Gmodel;

    Gmodel   = 0. ;
    nprotons = 0  ;
    site     = vector->sites;
    for (i = 0; i < vector->nsites; i++, site++) {
        Gmodel    +=    Real1DArray_Item (self->models  , site->indexActive);
        nprotons  += Integer1DArray_Item (self->protons , site->indexActive);
    }
    return (Gmodel - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature * CONSTANT_LN10 * pH));
}

/*
 * Calculate partition function and Boltzmann factors using a custom energy function.
 *
 * Note: bfactors should be allocated beforehand.
 */
Real EnergyModel_CalculateZ (const EnergyModel *self, Real (*EnergyFunction)(const EnergyModel*, const StateVector*, const Real), const Real pH, const Real Gzero, Real1DArray *bfactors) {
    Real     *bfactor, G, Gmin, Z;
    Integer   i;

    StateVector_Reset (self->vector);
    i       = self->nstates;
    bfactor = Real1DArray_Data (bfactors);
    Gmin    = EnergyFunction (self, self->vector, pH) - Gzero;
    for (; i > 0; i--, bfactor++) {
        /* TODO
         * Optimize: do not calculate Gmicro after every increment of the state 
         * vector, calculate deltas like in MC moves.
         */
        G = EnergyFunction (self, self->vector, pH) - Gzero;
        if (G < Gmin) {
            Gmin = G;
        }
        *bfactor = G;
        StateVector_Increment (self->vector);
    }
    Real1DArray_AddScalar (bfactors, -Gmin);
    Real1DArray_Scale (bfactors, -1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature));
    Real1DArray_Exp (bfactors);

    Z = Real1DArray_Sum (bfactors);
    return Z;
}

/*
 * Calculate the statistical mechanical partition function of an unfolded protein.
 */
Real EnergyModel_CalculateZunfolded (const EnergyModel *self, const Real pH, const Real Gzero, Status *status) {
    Real1DArray *bfactors;
    Real Z;

    bfactors = Real1DArray_Allocate (self->nstates, status);
    if (*status != Status_Continue) {
        return -1.;
    }
    Z = EnergyModel_CalculateZ (self, EnergyModel_CalculateMicrostateEnergyUnfolded, pH, Gzero, bfactors);
    Real1DArray_Deallocate (&bfactors);
    return Z;
}

/*
 * Calculate the statistical mechanical partition function of a folded protein.
 */
Real EnergyModel_CalculateZfolded (const EnergyModel *self, const Real pH, const Real Gzero, Status *status) {
    Real1DArray *bfactors;
    Real Z;

    bfactors = Real1DArray_Allocate (self->nstates, status);
    if (*status != Status_Continue) {
        return -1.;
    }
    Z = EnergyModel_CalculateZ (self, EnergyModel_CalculateMicrostateEnergy, pH, Gzero, bfactors);
    Real1DArray_Deallocate (&bfactors);
    return Z;
}

/*
 * Calculate protonation state probabilities from the statistical mechanical partition function.
 */
void EnergyModel_CalculateProbabilitiesFromZ (const EnergyModel *self, const Real Z, const Real1DArray *bfactors) {
    Real      *bfactor;
    Integer    i, j;
    TitrSite  *ts;

    Real1DArray_Set (self->probabilities, 0.);
    StateVector_Reset (self->vector);

    i = self->nstates;
    bfactor = Real1DArray_Data (bfactors);
    for (; i > 0; i--, bfactor++) {
        j  = self->vector->nsites ;
        ts = self->vector->sites  ;
        for (; j > 0; j--, ts++) {
            Real1DArray_Item (self->probabilities, ts->indexActive) += *bfactor;
        }
        StateVector_Increment (self->vector);
    }
    Real1DArray_Scale (self->probabilities, 1. / Z);
}

/*
 * Analytic evaluation of protonation state probabilities.
 */
void EnergyModel_CalculateProbabilitiesAnalytically (const EnergyModel *self, const Real pH, Status *status) {
    Real1DArray *bfactors;
    Real         Z;

    bfactors = Real1DArray_Allocate (self->nstates, status);
    if (*status != Status_Continue) {
        return;
    }
    Z = EnergyModel_CalculateZ (self, EnergyModel_CalculateMicrostateEnergy, pH, 0., bfactors);
    EnergyModel_CalculateProbabilitiesFromZ (self, Z, bfactors);

    Real1DArray_Deallocate (&bfactors);
}

/*
 * Analytic evaluation of protonation state probabilities (unfolded protein).
 */
void EnergyModel_CalculateProbabilitiesAnalyticallyUnfolded (const EnergyModel *self, const Real pH, Status *status) {
    Real1DArray *bfactors;
    Real         Z;

    bfactors = Real1DArray_Allocate (self->nstates, status);
    if (*status != Status_Continue) {
        return;
    }
    Z = EnergyModel_CalculateZ (self, EnergyModel_CalculateMicrostateEnergyUnfolded, pH, 0., bfactors);
    EnergyModel_CalculateProbabilitiesFromZ (self, Z, bfactors);

    Real1DArray_Deallocate (&bfactors);
}
