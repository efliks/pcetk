/*------------------------------------------------------------------------------
! . File      : StateVector.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "StateVector.h"


/*
 * Allocate the state vector.
 */
StateVector *StateVector_Allocate (const Integer nsites, Status *status) {
    StateVector *self;

    MEMORY_ALLOCATE (self, StateVector);
    if (self == NULL) {
        goto fail;
    }
    self->sites         = NULL  ;
    self->substateSites = NULL  ;
    self->pairs         = NULL  ;
    self->nsites        = 0     ;
    self->nssites       = 0     ;
    self->npairs        = 0     ;
    if (nsites > 0) {
        MEMORY_ALLOCATEARRAY (self->sites, nsites, TitrSite);
        if (self->sites == NULL) {
            MEMORY_DEALLOCATE (self);
            goto fail;
        }
    }
    self->nsites = nsites;
    return self;

fail:
    Status_Set (status, Status_MemoryAllocationFailure);
    return NULL;
}

/*
 * Allocate a substate within the state vector.
 */
void StateVector_AllocateSubstate (StateVector *self, const Integer nssites, Status *status) {
    if (self->substateSites != NULL) {
        Status_Set (status, Status_MemoryAllocationFailure);
    }
    else {
        MEMORY_ALLOCATEARRAY_POINTERS (self->substateSites, nssites, TitrSite);
        if (self->substateSites != NULL) {
            self->nssites = nssites;
        }
        else {
            Status_Set (status, Status_MemoryAllocationFailure);
        }
    }
}

/*
 * Allocate an array of pairs within the state vector.
 * The number of pairs and their contents are decided by the MCModelDefault module.
 */
void StateVector_AllocatePairs (StateVector *self, const Integer npairs, Status *status) {
    if (self->pairs != NULL) {
        Status_Set (status, Status_MemoryAllocationFailure);
    }
    else {
        MEMORY_ALLOCATEARRAY (self->pairs, npairs, PairSite);
        if (self->pairs != NULL) {
            self->npairs = npairs;
        }
        else {
            Status_Set (status, Status_MemoryAllocationFailure);
        }
    }
}

/*
 * Deallocate the old pairs and allocate the new ones.
 */
void StateVector_ReallocatePairs (StateVector *self, const Integer npairs, Status *status) {
    MEMORY_DEALLOCATE (self->pairs);
    self->npairs = 0;

    StateVector_AllocatePairs (self, npairs, status);
}

/*
 * Deallocate the state vector, including the optional arrays of pairs and substate.
 */
void StateVector_Deallocate (StateVector *self) {
    if (self != NULL) {
        if (self->pairs         != NULL) MEMORY_DEALLOCATE (self->pairs);
        if (self->substateSites != NULL) MEMORY_DEALLOCATE (self->substateSites);
        MEMORY_DEALLOCATE (self->sites);
        MEMORY_DEALLOCATE (self);
    }
}

/*
 * Clone a state vector.
 */
StateVector *StateVector_Clone (const StateVector *self, Status *status) {
    StateVector *clone = NULL;

    if (self != NULL) {
        clone = StateVector_Allocate (self->nsites, status);
        if (*status == Status_Continue) {
            StateVector_CopyTo (self, clone, status);
        }
    }
    return clone;
}

/*
 * Copy a state vector to another vector.
 */
void StateVector_CopyTo (const StateVector *self, StateVector *other, Status *status) {
    memcpy (other->sites, self->sites, self->nsites * sizeof (TitrSite));
}

/*
 * Set all sites of the vector to their initial instances.
 */
void StateVector_Reset (const StateVector *self) {
    TitrSite *site = self->sites;
    Integer   i = self->nsites;

    for (; i > 0; i--, site++) {
        site->indexActive = site->indexFirst;
    }
}

/*
 * Set all sites of the substate to their initial instances.
 */
void StateVector_ResetSubstate (const StateVector *self) {
    TitrSite *site, **pointToSite = self->substateSites;
    Integer   i = self->nssites;

    for (; i > 0; i--, pointToSite++) {
        site = *pointToSite;
        site->indexActive = site->indexFirst;
    }
}

/*
 * Set all sites of the vector to their final instances.
 */
void StateVector_ResetToMaximum (const StateVector *self) {
    TitrSite *site = self->sites;
    Integer   i = self->nsites;

    for (; i > 0; i--, site++) {
        site->indexActive = site->indexLast;
    }
}

/*
 * Set all sites of the vector to randomized instances.
 */
void StateVector_Randomize (const StateVector *self, const RandomNumberGenerator *generator) {
    TitrSite *site = self->sites;
    Integer   i = self->nsites;

    for (; i > 0; i--, site++) {
        site->indexActive = RandomNumberGenerator_NextCardinal (generator) % (site->indexLast - site->indexFirst + 1) + site->indexFirst;
    }
}

/*
 * Set a state vector site.
 */
void StateVector_SetSite (const StateVector *self, const Integer indexSite, 
                          const Integer indexFirst, const Integer indexLast, Status *status) {
    TitrSite *site;

    if (indexSite < 0 || indexSite >= self->nsites) {
        Status_Set (status, Status_IndexOutOfRange);
    }
    else {
        site = &self->sites[indexSite];
        site->isSubstate   =  False       ;
        site->indexSite    =  indexSite   ;
        site->indexLast    =  indexLast   ;
        site->indexFirst   =  indexFirst  ;
        site->indexActive  =  indexFirst  ;
    }
}

/*
 * Set a pair of strongly interacting sites.
 */
void StateVector_SetPair (const StateVector *self, const Integer indexPair, 
                          const Integer indexFirstSite, const Integer indexSecondSite, 
                          const Real Wmax, Status *status) {
    PairSite *pair;

    if (indexPair < 0 || indexPair >= self->npairs) {
        Status_Set (status, Status_IndexOutOfRange);
    }
    else {
        pair       = &self->pairs[indexPair];
        pair->a    = &self->sites[indexFirstSite];
        pair->b    = &self->sites[indexSecondSite];
        pair->Wmax = Wmax;
    }
}

/*
 * Get indices and maximum interaction energy of a pair of strongly interacting sites.
 */
void StateVector_GetPair (const StateVector *self, const Integer indexPair, 
                          Integer *indexFirstSite, Integer *indexSecondSite, Real *Wmax, 
                          Status *status) {
    PairSite *pair;
    TitrSite *site;

    if (indexPair < 0 || indexPair >= self->npairs) {
        Status_Set (status, Status_IndexOutOfRange);
    }
    else {
        pair  = &self->pairs[indexPair];
        site  = pair->a;
        *indexFirstSite  = site->indexSite;
        site  = pair->b;
        *indexSecondSite = site->indexSite;
        *Wmax = pair->Wmax;
    }
}

/*
 * Return true if the site belongs to a substate.
 */
Boolean StateVector_IsSubstate (const StateVector *self, const Integer siteIndex, Status *status) {
    TitrSite *site;

    if (siteIndex < 0 || siteIndex >= self->nsites) {
        Status_Set (status, Status_IndexOutOfRange);
        return False;
    }
    site = &self->sites[siteIndex];
    return site->isSubstate;
}

/*
 * Get the current protonation of a site, i.e. the local index of its currently "active" instance.
 */
Integer StateVector_GetItem (const StateVector *self, const Integer siteIndex, Status *status) {
    TitrSite *site;
    Integer instanceLocalIndex;

    if (siteIndex < 0 || siteIndex >= self->nsites) {
        Status_Set (status, Status_IndexOutOfRange);
        return -1;
    }
    site = &self->sites[siteIndex];
    instanceLocalIndex = site->indexActive - site->indexFirst;

    return instanceLocalIndex;
}

/*
 * Set the protonation of a site by defining a local index of its "active" instance.
 */
void StateVector_SetItem (const StateVector *self, const Integer siteIndex, const Integer instanceLocalIndex, 
                          Status *status) {
    TitrSite *site;
    Integer instanceGlobalIndex;

    if (siteIndex < 0 || siteIndex >= self->nsites) {
        Status_Set (status, Status_IndexOutOfRange);
    }
    else {
        site = &self->sites[siteIndex];
        /* Translate local index to global index */
        instanceGlobalIndex = instanceLocalIndex + site->indexFirst;

        if (instanceGlobalIndex < site->indexFirst || instanceGlobalIndex > site->indexLast) {
            Status_Set (status, Status_ValueError);
        }
        else {
            site->indexActive = instanceGlobalIndex;
        }
    }
}

/*
 * Get the current protonation of a site, i.e. the global index of its currently "active" instance.
 */
Integer StateVector_GetActualItem (const StateVector *self, const Integer siteIndex, Status *status) {
    TitrSite *site;
    Integer instanceGlobalIndex;

    if (siteIndex < 0 || siteIndex >= self->nsites) {
        Status_Set (status, Status_IndexOutOfRange);
        return -1;
    }
    site = &self->sites[siteIndex];
    instanceGlobalIndex = site->indexActive;

    return instanceGlobalIndex;
}

/*
 * Set the protonation of a site by defining a global index of its "active" instance.
 */
void StateVector_SetActualItem (const StateVector *self, const Integer siteIndex, const Integer instanceGlobalIndex, 
                                Status *status) {
    TitrSite *site;

    if (siteIndex < 0 || siteIndex >= self->nsites) {
        Status_Set (status, Status_IndexOutOfRange);
    }
    else {
        site = &self->sites[siteIndex];
        if (instanceGlobalIndex < site->indexFirst || instanceGlobalIndex > site->indexLast) {
            Status_Set (status, Status_ValueError);
        }
        else {
            site->indexActive = instanceGlobalIndex;
        }
    }
}

/*
 * Get the index of a site belonging to a substate.
 */
Integer StateVector_GetSubstateItem (const StateVector *self, const Integer index, Status *status) {
    TitrSite *site;

    if (index < 0 || index >= self->nssites) {
        Status_Set (status, Status_IndexOutOfRange);
        return -1;
    }
    site = self->substateSites[index];
    return site->indexSite;
}

/*
 * Attach the selected site to a substate by passing its index.
 */
void StateVector_SetSubstateItem (const StateVector *self, const Integer selectedSiteIndex, const Integer index, 
                                  Status *status) {
    TitrSite *site;

    if (index < 0 || index >= self->nssites) {
        Status_Set (status, Status_IndexOutOfRange);
    }
    else {
        if (selectedSiteIndex < 0 || selectedSiteIndex >= self->nsites) {
            Status_Set (status, Status_ValueError);
        }
        else {
            site = &self->sites[selectedSiteIndex];
            site->isSubstate = True;
            self->substateSites[index] = site;
        }
    }
}

/*
 * Increment the state vector.
 * After reaching the last vector, false is returned and the vector is back in its initial state.
 * True is returned as long as there are more vectors ahead.
 *
 * Incrementation algorithm by Timm Essigke.
 */
Boolean StateVector_Increment (const StateVector *self) {
    TitrSite *site = self->sites;
    Integer   i    = self->nsites;

    for (; i > 0; i--, site++) {
        if (site->indexActive < site->indexLast) {
            site->indexActive++;
            return True;
        }
        else {
            site->indexActive = site->indexFirst;
        }
    }
    return False;
}

/*
 * Increment only within the substate of sites of the vector.
 */
Boolean StateVector_IncrementSubstate (const StateVector *self) {
    TitrSite *site, **pointToSite = self->substateSites;
    Integer   i = self->nssites;

    for (; i > 0; i--, pointToSite++) {
        site = *pointToSite;
        if (site->indexActive < site->indexLast) {
            site->indexActive++;
            return True;
        }
        else {
            site->indexActive = site->indexFirst;
        }
    }
    return False;
}
