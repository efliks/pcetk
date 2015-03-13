/*------------------------------------------------------------------------------
! . File      : EnergyModel.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "EnergyModel.h"

/*
 * Allocate the energy model and initialize its random number generator.
 * Other attributes of the model (nstates, ninstances, temperature) are set from the Cython level.
 */
EnergyModel *EnergyModel_Allocate (const Integer nsites, const Integer ninstances, Status *status) {
    EnergyModel *self = NULL;

    MEMORY_ALLOCATE (self, EnergyModel);
    if (self == NULL)
        goto failSet;

    self->vector           =  NULL  ;
    self->generator        =  NULL  ;

    self->protons          =  NULL  ;
    self->models           =  NULL  ;
    self->intrinsic        =  NULL  ;
    self->interactions     =  NULL  ;
    self->probabilities    =  NULL  ;
    self->symmetricmatrix  =  NULL  ;

    self->generator = RandomNumberGenerator_Allocate (RandomNumberGeneratorType_MersenneTwister);
    if (self->generator == NULL)
        goto failSetDealloc;
    RandomNumberGenerator_SetSeed (self->generator, (Cardinal) time (NULL));

    if (nsites > 0) {
        self->vector = StateVector_Allocate (nsites, status);
        if (*status != Status_Continue)
            goto failDealloc;
    }

    if (ninstances > 0) {
        self->protons         = Integer1DArray_Allocate  (ninstances, status)             ;
        self->models          = Real1DArray_Allocate     (ninstances, status)             ;
        self->intrinsic       = Real1DArray_Allocate     (ninstances, status)             ;
        self->interactions    = Real2DArray_Allocate     (ninstances, ninstances, status) ;
        self->probabilities   = Real1DArray_Allocate     (ninstances, status)             ;
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
    if ( self->generator       != NULL)  RandomNumberGenerator_Deallocate ( &self->generator       ) ;
    if ( self->symmetricmatrix != NULL)  SymmetricMatrix_Deallocate       ( &self->symmetricmatrix ) ;
    if ( self->probabilities   != NULL)  Real1DArray_Deallocate           ( &self->probabilities   ) ;
    if ( self->interactions    != NULL)  Real2DArray_Deallocate           ( &self->interactions    ) ;
    if ( self->intrinsic       != NULL)  Real1DArray_Deallocate           ( &self->intrinsic       ) ;
    if ( self->models          != NULL)  Real1DArray_Deallocate           ( &self->models          ) ;
    if ( self->protons         != NULL)  Integer1DArray_Deallocate        ( &self->protons         ) ;
    if ( self->vector          != NULL)  StateVector_Deallocate           (  self->vector          ) ;
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
 * Evaluate the probabilities of occurrence of instances from the statistical mechanical partition function.
 */
void EnergyModel_CalculateProbabilitiesAnalytically (const EnergyModel *self, const Real pH, Status *status) {
    Real         *bfactor, G, Gzero, bsum;
    Real1DArray  *bfactors;
    Integer       i, j;
    TitrSite     *ts;

    bfactors = Real1DArray_Allocate (self->nstates, status);
    if (*status != Status_Continue) {
        return;
    }

    StateVector_Reset (self->vector);
    i       = self->nstates;
    bfactor = Real1DArray_Data (bfactors);
    Gzero   = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH);

    for (; i > 0; i--, bfactor++) {
        G = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH);
        if (G < Gzero) {
            Gzero = G;
        }
        *bfactor = G;
        StateVector_Increment (self->vector);
    }

    Real1DArray_AddScalar (bfactors, -Gzero);
    Real1DArray_Scale (bfactors, -1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature));
    Real1DArray_Exp (bfactors);

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

    bsum = Real1DArray_Sum (bfactors);
    Real1DArray_Scale (self->probabilities, 1. / bsum);

    Real1DArray_Deallocate (&bfactors);
}

/*
 * The Metropolis criterion.
 *
 * Function adopted from GMCT. RT-units of energy are used, instead of the usual kcal/mol.
 */
Boolean EnergyModel_Metropolis (const Real GdeltaRT, const RandomNumberGenerator *generator) {
    Boolean metropolis;

    if (GdeltaRT < 0.)
        metropolis = True;
    else {
        if (-GdeltaRT < TOO_SMALL)
            metropolis = False;
        else if (RandomNumberGenerator_NextReal (generator) < exp (-GdeltaRT))
            metropolis = True;
        else
            metropolis = False;
    }
    return metropolis;
}

/*
 * Choose a random site and change its "active" instance.
 */
Boolean EnergyModel_Move (const EnergyModel *self, const Real pH, const Real G, Real *Gnew) {
    Integer    site, instance, nprotons, i;
    Real       Gintr, W, Gdelta, GdeltaRT, beta;
    TitrSite  *ts, *tsOther;
    Boolean    accept;

    site = RandomNumberGenerator_NextCardinal (self->generator) % self->vector->nsites;
    ts   = &self->vector->sites[site];
    do {
        instance = RandomNumberGenerator_NextCardinal (self->generator) % (ts->indexLast - ts->indexFirst + 1) + ts->indexFirst;
    } while (instance == ts->indexActive);

    Gintr    =    Real1DArray_Item (self->intrinsic , instance) -    Real1DArray_Item (self->intrinsic , ts->indexActive);
    nprotons = Integer1DArray_Item (self->protons   , instance) - Integer1DArray_Item (self->protons   , ts->indexActive);

    W       = 0.;
    i       = self->vector->nsites ;
    tsOther = self->vector->sites  ;
    for (; i > 0; i--, tsOther++) {
        W += (EnergyModel_GetW (self, instance, tsOther->indexActive) - EnergyModel_GetW (self, ts->indexActive, tsOther->indexActive));
    }

    Gdelta   = Gintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature * CONSTANT_LN10 * pH) + W;
    beta     = 1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature);
    GdeltaRT = Gdelta * beta;
    accept   = EnergyModel_Metropolis (GdeltaRT, self->generator);
    if (accept) {
        ts->indexActive = instance;
        *Gnew = G + Gdelta;
    }
    return accept;
}

/*
 * Choose a random pair of sites and change their "active" instances.
 */
Boolean EnergyModel_DoubleMove (const EnergyModel *self, const Real pH, const Real G, Real *Gnew) {
    Integer    newPair, indexSa, indexSb, nprotons, i;
    Real       Gintr, W, Gdelta, GdeltaRT, beta;
    TitrSite  *sa, *sb, *tsOther;
    PairSite  *pair;
    Boolean    accept;

    newPair = RandomNumberGenerator_NextCardinal (self->generator) % self->vector->npairs;
    pair    = &self->vector->pairs[newPair];
    sa      = pair->a;
    sb      = pair->b;

    do {
        indexSa = RandomNumberGenerator_NextCardinal (self->generator) % (sa->indexLast - sa->indexFirst + 1) + sa->indexFirst;
    } while (indexSa == sa->indexActive);
    do {
        indexSb = RandomNumberGenerator_NextCardinal (self->generator) % (sb->indexLast - sb->indexFirst + 1) + sb->indexFirst;
    } while (indexSb == sb->indexActive);

    Gintr =        Real1DArray_Item (self->intrinsic , indexSa) -    Real1DArray_Item (self->intrinsic , sa->indexActive)
                 + Real1DArray_Item (self->intrinsic , indexSb) -    Real1DArray_Item (self->intrinsic , sb->indexActive);
    nprotons =  Integer1DArray_Item (self->protons   , indexSa) - Integer1DArray_Item (self->protons   , sa->indexActive)
              + Integer1DArray_Item (self->protons   , indexSb) - Integer1DArray_Item (self->protons   , sb->indexActive);

    W       = EnergyModel_GetW (self, indexSa, indexSb) - EnergyModel_GetW (self, sa->indexActive, sb->indexActive);
    tsOther = self->vector->sites;
    for (i = 0; i < self->vector->nsites; i++, tsOther++) {
        if (i != sa->indexSite && i != sb->indexSite) {
            W +=  EnergyModel_GetW (self, indexSa, tsOther->indexActive) - EnergyModel_GetW (self, sa->indexActive, tsOther->indexActive)
                + EnergyModel_GetW (self, indexSb, tsOther->indexActive) - EnergyModel_GetW (self, sb->indexActive, tsOther->indexActive);
        }
    }

    Gdelta   = Gintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature * CONSTANT_LN10 * pH) + W;
    beta     = 1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature);
    GdeltaRT = Gdelta * beta;
    accept   = EnergyModel_Metropolis (GdeltaRT, self->generator);

    if (accept) {
        sa->indexActive = indexSa;
        sb->indexActive = indexSb;
        *Gnew = G + Gdelta;
    }
    return accept;
}

/*
 * Generate a state vector representing a low-energy, statistically relevant protonation state.
 * The energy of the vector is returned only for information.
 */
Real EnergyModel_MCScan (const EnergyModel *self, const Real pH, Integer nmoves) {
    Integer  selection, select;
    Real     G, Gnew;
    Boolean  accept;

    G = EnergyModel_CalculateMicrostateEnergy (self, self->vector, pH);
    selection = self->vector->nsites + self->vector->npairs;

    for (; nmoves > 0; nmoves--) {
        select = RandomNumberGenerator_NextCardinal (self->generator) % selection;
        if (select < self->vector->nsites)
            accept = EnergyModel_Move (self, pH, G, &Gnew);
        else
            accept = EnergyModel_DoubleMove (self, pH, G, &Gnew);

        if (accept)
            G = Gnew;
    }
    return G;
}

/*
 * Increase the counts of "active" instances.
 * These counts, after scaling, will give the probabilities of occurrence of instances.
 */
void EnergyModel_UpdateProbabilities (const EnergyModel *self) {
    Real     *counter;
    TitrSite *ts;
    Integer   i;
    ts = self->vector->sites  ;
    i  = self->vector->nsites ;

    for (; i > 0; i--, ts++) {
        counter = Real1DArray_ItemPointer (self->probabilities, ts->indexActive);
        (*counter)++;
    }
}

/*
 * Find maximum absolute interaction energy between two sites.
 */
Real EnergyModel_FindMaxInteraction (const EnergyModel *self, const TitrSite *site, const TitrSite *other) {
    Integer index, indexOther;
    Real W, Wmax;

    Wmax  = 0.;
    index = site->indexFirst;
    for (; index <= site->indexLast; index++) {

        indexOther = other->indexFirst;
        for (; indexOther <= other->indexLast; indexOther++) {
            W = fabs (EnergyModel_GetW (self, index, indexOther));
            if (W > Wmax) Wmax = W;
        }
    }
    return Wmax;
}

/*
 * Find pairs of sites whose interaction energy is greater than the given limit.
 *
 * If npairs < 1, dry run is assumed and only nfound is returned.
 * The value of npairs is used in the second run to allocate and fill out the pairs.
 */
Integer EnergyModel_FindPairs (const EnergyModel *self, const Real limit, const Integer npairs, Status *status) {
    TitrSite *site, *siteInner;
    Integer i, j, nfound;
    Real Wmax;

    if (npairs > 0) {
        if (self->vector->npairs > 0)
            StateVector_ReallocatePairs (self->vector, npairs, status);
        else
            StateVector_AllocatePairs (self->vector, npairs, status);
        if (*status != Status_Continue)
            return -1;
    }

    nfound = 0;
    site   = self->vector->sites;
    for (i = 0; i < self->vector->nsites; i++, site++) {

        siteInner = self->vector->sites;
        for (j = 0; j < i; j++, siteInner++) {

            Wmax = EnergyModel_FindMaxInteraction (self, site, siteInner);
            if (Wmax >= limit) {
                if (npairs > 0) {
                    StateVector_SetPair (self->vector, nfound, site->indexSite, siteInner->indexSite, Wmax, status);
                }
                nfound++;
            }
        }
    }
    return nfound;
}

/*
 * Run a Monte Carlo simulation.
 *
 * If equil is true, only the equilibration run is performed.
 *
 * Otherwise, the production run is done. The resulting state vectors are not accumulated.
 * Instead, they are immediately used to update the probabilities.
 *
 * The number of moves is proportional to the number of sites and pairs.
 */
void EnergyModel_CalculateProbabilitiesMonteCarlo (const EnergyModel *self, const Real pH, const Boolean equil, Integer nscans, Status *status) {
    Real    Gfinal, scale;
    Integer nmoves;

    nmoves = self->vector->nsites + self->vector->npairs;
    scale  = 1. / nscans;
    if (equil) {
        StateVector_Randomize (self->vector, self->generator);
        for (; nscans > 0; nscans--)
            Gfinal = EnergyModel_MCScan (self, pH, nmoves);
    }
    else {
        Real1DArray_Set (self->probabilities, 0.);

        for (; nscans > 0; nscans--) {
            Gfinal = EnergyModel_MCScan (self, pH, nmoves);
            EnergyModel_UpdateProbabilities (self);
        }
        Real1DArray_Scale (self->probabilities, scale);
    }
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
 * Calculate the statistical mechanical partition function using a custom energy function.
 */
static Real EnergyModel_PartitionFunction (const EnergyModel *self, Real (*EnergyFunction)(const EnergyModel*, const StateVector*, const Real), const Real pH, const Real Gneutral, Status *status) {
    Real1DArray  *bfactors;
    Real         *bf, Z, G, Gzero;
    Integer       nstates;

    nstates = 0;
    StateVector_Reset (self->vector);
    do {
        nstates++;
    } while (StateVector_Increment (self->vector));

    bfactors = Real1DArray_Allocate (nstates, status);
    if (*status != Status_Continue) {
        return -1.;
    }

    StateVector_Reset (self->vector);
    Gzero = EnergyFunction (self, self->vector, pH);
    bf    = Real1DArray_Data (bfactors);
    do {
        G = EnergyFunction (self, self->vector, pH) - Gneutral;
        if (G < Gzero)
            Gzero = G;
        *(bf++) = G;
    } while (StateVector_Increment (self->vector));

    Real1DArray_AddScalar (bfactors, -Gzero);
    Real1DArray_Scale (bfactors, -1. / (CONSTANT_MOLAR_GAS_KCAL_MOL * self->temperature));
    Real1DArray_Exp (bfactors);
    Z = Real1DArray_Sum (bfactors);

    Real1DArray_Deallocate (&bfactors);
    return Z;
}

/*
 * Calculate the statistical mechanical partition function of an unfolded protein.
 */
Real EnergyModel_PartitionFunctionUnfolded (const EnergyModel *self, const Real pH, const Real Gneutral, Status *status) {
    return EnergyModel_PartitionFunction (self, EnergyModel_CalculateMicrostateEnergyUnfolded, pH, Gneutral, status);
}

/*
 * Calculate the statistical mechanical partition function of a folded protein.
 */
Real EnergyModel_PartitionFunctionFolded (const EnergyModel *self, const Real pH, const Real Gneutral, Status *status) {
    return EnergyModel_PartitionFunction (self, EnergyModel_CalculateMicrostateEnergy, pH, Gneutral, status);
}
