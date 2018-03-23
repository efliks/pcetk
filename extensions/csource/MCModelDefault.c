/*------------------------------------------------------------------------------
! . File      : MCModelDefault.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "MCModelDefault.h"

/*
 * Allocate the default Monte Carlo model and initialize its random number generator.
 */
MCModelDefault *MCModelDefault_Allocate (const Real limit, const Integer nequil, 
                                         const Integer nprod, const Integer randomSeed, 
                                         Status *status) {
    MCModelDefault *self = NULL;

    MEMORY_ALLOCATE (self, MCModelDefault);
    if (self == NULL) {
        goto failSet;
    }
    self->limit       = limit  ;
    self->nprod       = nprod  ;
    self->nequil      = nequil ;
    self->vector      = NULL   ;
    self->energyModel = NULL   ;

    self->generator   = RandomNumberGenerator_Allocate (RandomNumberGeneratorType_MersenneTwister);
    if (self->generator == NULL) {
        goto failSetDealloc;
    }
    if (randomSeed <= 0) {
        RandomNumberGenerator_SetSeed (self->generator, (Cardinal) time (NULL));
    }
    else {
        RandomNumberGenerator_SetSeed (self->generator, (Cardinal) randomSeed);
    }
    return self;

failSetDealloc:
    MCModelDefault_Deallocate (self);

failSet:
    Status_Set (status, Status_MemoryAllocationFailure);
    return NULL;
}

/*
 * Deallocate Monte Carlo model.
 */
void MCModelDefault_Deallocate (MCModelDefault *self) {
    if ( self->generator       != NULL)  RandomNumberGenerator_Deallocate ( &self->generator ) ;
    if ( self->vector          != NULL)  StateVector_Deallocate           (  self->vector    ) ;
    if (self != NULL) MEMORY_DEALLOCATE (self);
}

/*
 * Link Monte Carlo and energy models.
 */
void MCModelDefault_LinkToEnergyModel (MCModelDefault *self, EnergyModel *energyModel, 
                                       Status *status) {
    self->energyModel = energyModel;
    self->vector      = StateVector_Clone (energyModel->vector, status);
}

/*
 * The Metropolis criterion.
 *
 * Function adopted from GMCT. RT-units of energy are used, instead of the usual kcal/mol.
 */
Boolean MCModelDefault_Metropolis (const Real GdeltaRT, 
                                   const RandomNumberGenerator *generator) {
    if (GdeltaRT < 0.0f) {
        return True;
    }
    else {
        if (-GdeltaRT < TOO_SMALL) {
            return False;
        }
        else if (RandomNumberGenerator_NextReal (generator) < exp (-GdeltaRT)) {
            return True;
        }
    }
    return False;
}

/*
 * Choose a random site and change its "active" instance.
 */
Boolean MCModelDefault_Move (const MCModelDefault *self, const Real pH, 
                             const Real G, Real *Gnew) {
    Integer    site, instance, nprotons, i;
    Real       Gintr, W, Gdelta, GdeltaRT, beta;
    TitrSite  *ts, *tsOther;
    Boolean    accept;

    site = RandomNumberGenerator_NextCardinal (self->generator) % self->vector->nsites;
    ts   = &self->vector->sites[site];
    do {
        instance = RandomNumberGenerator_NextCardinal (self->generator) % (ts->indexLast - ts->indexFirst + 1) + ts->indexFirst;
    } while (instance == ts->indexActive);

    Gintr    =    Real1DArray_Item (self->energyModel->intrinsic , instance) -    Real1DArray_Item (self->energyModel->intrinsic , ts->indexActive);
    nprotons = Integer1DArray_Item (self->energyModel->protons   , instance) - Integer1DArray_Item (self->energyModel->protons   , ts->indexActive);

    W       = 0.0f;
    i       = self->vector->nsites ;
    tsOther = self->vector->sites  ;
    for (; i > 0; i--, tsOther++) {
        W += (EnergyModel_GetW (self->energyModel, instance, tsOther->indexActive) - EnergyModel_GetW (self->energyModel, ts->indexActive, tsOther->indexActive));
    }

    Gdelta   = Gintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * self->energyModel->temperature * CONSTANT_LN10 * pH) + W;
    beta     = 1.0f / (CONSTANT_MOLAR_GAS_KCAL_MOL * self->energyModel->temperature);
    GdeltaRT = Gdelta * beta;
    accept   = MCModelDefault_Metropolis (GdeltaRT, self->generator);
    if (accept) {
        ts->indexActive = instance;
        *Gnew = G + Gdelta;
    }
    return accept;
}

/*
 * Choose a random pair of sites and change their "active" instances.
 */
Boolean MCModelDefault_DoubleMove (const MCModelDefault *self, const Real pH, 
                                   const Real G, Real *Gnew) {
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

    Gintr =        Real1DArray_Item (self->energyModel->intrinsic , indexSa) -    Real1DArray_Item (self->energyModel->intrinsic , sa->indexActive)
                 + Real1DArray_Item (self->energyModel->intrinsic , indexSb) -    Real1DArray_Item (self->energyModel->intrinsic , sb->indexActive);
    nprotons =  Integer1DArray_Item (self->energyModel->protons   , indexSa) - Integer1DArray_Item (self->energyModel->protons   , sa->indexActive)
              + Integer1DArray_Item (self->energyModel->protons   , indexSb) - Integer1DArray_Item (self->energyModel->protons   , sb->indexActive);

    W       = EnergyModel_GetW (self->energyModel, indexSa, indexSb) - EnergyModel_GetW (self->energyModel, sa->indexActive, sb->indexActive);
    tsOther = self->vector->sites;
    for (i = 0; i < self->vector->nsites; i++, tsOther++) {
        if (i != sa->indexSite && i != sb->indexSite) {
            W +=  EnergyModel_GetW (self->energyModel, indexSa, tsOther->indexActive) - EnergyModel_GetW (self->energyModel, sa->indexActive, tsOther->indexActive)
                + EnergyModel_GetW (self->energyModel, indexSb, tsOther->indexActive) - EnergyModel_GetW (self->energyModel, sb->indexActive, tsOther->indexActive);
        }
    }

    Gdelta   = Gintr - nprotons * (-CONSTANT_MOLAR_GAS_KCAL_MOL * self->energyModel->temperature * CONSTANT_LN10 * pH) + W;
    beta     = 1.0f / (CONSTANT_MOLAR_GAS_KCAL_MOL * self->energyModel->temperature);
    GdeltaRT = Gdelta * beta;
    accept   = MCModelDefault_Metropolis (GdeltaRT, self->generator);

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
Real MCModelDefault_MCScan (const MCModelDefault *self, const Real pH, Integer nmoves, 
                            Integer *movesDone, Integer *movesAccepted, Integer *flipsDone, 
                            Integer *flipsAccepted) {
    Integer  selection, select;
    Real     G, Gnew;
    Boolean  accept, move;

    G = EnergyModel_CalculateMicrostateEnergy (self->energyModel, self->vector, pH);
    selection = self->vector->nsites + self->vector->npairs;

    for (; nmoves > 0; nmoves--) {
        select = RandomNumberGenerator_NextCardinal (self->generator) % selection;
        move   = select < self->vector->nsites;
        if (move) {
            (*movesDone)++;
            accept = MCModelDefault_Move (self, pH, G, &Gnew);
            if (accept) {
                (*movesAccepted)++;
                G = Gnew;
            }
        }
        else {
            (*flipsDone)++;
            accept = MCModelDefault_DoubleMove (self, pH, G, &Gnew);
            if (accept) {
                (*flipsAccepted)++;
                G = Gnew;
            }
        }
    }
    return G;
}

/*
 * Increase the counts of "active" instances.
 * These counts, after scaling, will give the probabilities of occurrence of instances.
 */
void MCModelDefault_UpdateProbabilities (const MCModelDefault *self) {
    Real     *counter;
    TitrSite *ts;
    Integer   i;
    ts = self->vector->sites  ;
    i  = self->vector->nsites ;

    for (; i > 0; i--, ts++) {
        counter = Real1DArray_ItemPointer (self->energyModel->probabilities, ts->indexActive);
        (*counter)++;
    }
}

/*
 * Find maximum absolute interaction energy between two sites.
 */
Real MCModelDefault_FindMaxInteraction (const MCModelDefault *self, 
                                        const TitrSite *site, const TitrSite *other) {
    Integer index, indexOther;
    Real W, Wmax;

    Wmax  = 0.0f;
    index = site->indexFirst;
    for (; index <= site->indexLast; index++) {
        indexOther = other->indexFirst;
        for (; indexOther <= other->indexLast; indexOther++) {
            W = fabs (EnergyModel_GetW (self->energyModel, index, indexOther));
            if (W > Wmax) {
                Wmax = W;
            }
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
Integer MCModelDefault_FindPairs (const MCModelDefault *self, const Integer npairs, 
                                  Status *status) {
    TitrSite *site, *siteInner;
    Integer i, j, nfound;
    Real Wmax;

    if (npairs > 0) {
        if (self->vector->npairs > 0) {
            StateVector_ReallocatePairs (self->vector, npairs, status);
        }
        else {
            StateVector_AllocatePairs (self->vector, npairs, status);
        }
        if (*status != Status_Continue) {
            return -1;
        }
    }
    nfound = 0;
    site   = self->vector->sites;
    for (i = 0; i < self->vector->nsites; i++, site++) {
        siteInner = self->vector->sites;
        for (j = 0; j < i; j++, siteInner++) {
            Wmax = MCModelDefault_FindMaxInteraction (self, site, siteInner);
            if (Wmax >= self->limit) {
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
 * Run a Monte Carlo production.
 *
 * The number of moves during each scan is proportional to the number of sites and pairs.
 *
 * The resulting state vectors are not accumulated. Instead, they are immediately 
 * used to update the probabilities.
 */
void MCModelDefault_Production (const MCModelDefault *self, const Real pH) {
    Integer  nmoves, nscans, moves, movesAcc, flips, flipsAcc;
    Real     scale, Gfinal;

    Real1DArray_Set (self->energyModel->probabilities, 0.0f);
    nmoves = self->vector->nsites + self->vector->npairs;
    nscans = self->nprod;

    for (; nscans > 0; nscans--) {
        Gfinal = MCModelDefault_MCScan (self, pH, nmoves, &moves, &movesAcc, &flips, &flipsAcc);
        MCModelDefault_UpdateProbabilities (self);
    }
    scale = 1.0f / self->nprod;
    Real1DArray_Scale (self->energyModel->probabilities, scale);
}

/*
 * Run a Monte Carlo equilibration.
 */
void MCModelDefault_Equilibration (const MCModelDefault *self, const Real pH) {
    Integer  nmoves, nscans, moves, movesAcc, flips, flipsAcc;

    StateVector_Randomize (self->vector, self->generator);
    nmoves = self->vector->nsites + self->vector->npairs;
    nscans = self->nequil;

    for (; nscans > 0; nscans--) {
        MCModelDefault_MCScan (self, pH, nmoves, &moves, &movesAcc, &flips, &flipsAcc);
    }
}
