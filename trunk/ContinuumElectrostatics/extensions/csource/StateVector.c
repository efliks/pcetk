/*------------------------------------------------------------------------------
! . File      : StateVector.c
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#include "StateVector.h"


/*=============================================================================
  Allocation and deallocation
=============================================================================*/
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

void StateVector_AllocateSubstate (StateVector *self, const Integer nssites, Status *status) {
  if (self->substateSites != NULL) {
    /* Substate already allocated */
    Status_Set (status, Status_MemoryAllocationFailure);
  }
  else {
    /* Allocate an array of pointers */
    MEMORY_ALLOCATEARRAY_POINTERS (self->substateSites, nssites, TitrSite);

    if (self->substateSites != NULL) 
      self->nssites = nssites;
    else
      Status_Set (status, Status_MemoryAllocationFailure);
  }
}

void StateVector_AllocatePairs (StateVector *self, const Integer npairs, Status *status) {
  if (self->substateSites != NULL) {
    /* Pairs already allocated */
    Status_Set (status, Status_MemoryAllocationFailure);
  }
  else {
    MEMORY_ALLOCATEARRAY (self->pairs, npairs, PairSite);

    if (self->pairs != NULL) 
      self->npairs = npairs;
    else
      Status_Set (status, Status_MemoryAllocationFailure);
  }
}

void StateVector_Deallocate (StateVector *self) {
  if (self != NULL) {
    /* First deallocate optional tables */
    if (self->pairs         != NULL) MEMORY_DEALLOCATE (self->pairs);
    if (self->substateSites != NULL) MEMORY_DEALLOCATE (self->substateSites);

    MEMORY_DEALLOCATE (self->sites);
    MEMORY_DEALLOCATE (self);
  }
}

/*=============================================================================
  Copying and cloning
=============================================================================*/
/*StateVector *StateVector_Clone (const StateVector *self, Status *status) {
  StateVector *clone = NULL;

  if (self != NULL) {
    clone = StateVector_Allocate (self->nsites, status);
    if (*status == Status_Continue) {
      StateVector_CopyTo (self, clone, status);
    }
  }

  return clone;
}

void StateVector_CopyTo (const StateVector *self, StateVector *other, Status *status) {
  if (self->nsites != other->nsites) {
    StateVector_Deallocate (other);
    other = StateVector_Allocate (self->nsites, status);

    if (*status != Status_Continue) {
      return;
    }
  }
  memcpy (other->sites, self->sites, other->nsites * sizeof (TitrSite*));

  if (self->substateSites != NULL) {
    StateVector_AllocateSubstate (other, self->nssites, status);
    if (*status != Status_Continue) {
      return;
    }
    memcpy (other->substateSites, self->substateSites, other->nssites * sizeof (TitrSite*));
  }
}*/

/*=============================================================================
  Functions for setting all items at once
=============================================================================*/
void StateVector_Reset (const StateVector *self) {
  TitrSite *site = self->sites;
  Integer   i = self->nsites;

  for (; i > 0; i--, site++) {
    site->indexActive = site->indexFirst;
  }
}

void StateVector_ResetSubstate (const StateVector *self) {
  TitrSite *site, **pointToSite = self->substateSites;
  Integer   i = self->nssites;

  for (; i > 0; i--, pointToSite++) {
    site = *pointToSite;
    site->indexActive = site->indexFirst;
  }
}

void StateVector_ResetToMaximum (const StateVector *self) {
  TitrSite *site = self->sites;
  Integer   i = self->nsites;

  for (; i > 0; i--, site++) {
    site->indexActive = site->indexLast;
  }
}

void StateVector_Randomize (const StateVector *self, const RandomNumberGenerator *generator) {
  TitrSite *site = self->sites;
  Integer   i = self->nsites;

  for (; i > 0; i--, site++) {
    site->indexActive = RandomNumberGenerator_NextCardinal (generator) % (site->indexLast - site->indexFirst + 1) + site->indexFirst;
  }
}

/*=============================================================================
  Functions for accessing items
=============================================================================*/
void StateVector_SetSite (const StateVector *self, const Integer indexSite, const Integer indexFirst, const Integer indexLast, Status *status) {
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

void StateVector_SetPair (const StateVector *self, const Integer indexPair, const Integer indexFirstSite, const Integer indexSecondSite, const Real Wmax, Status *status) {
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

void StateVector_GetPair (const StateVector *self, const Integer indexPair, Integer *indexFirstSite, Integer *indexSecondSite, Real *Wmax, Status *status) {
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

Boolean StateVector_IsSubstate (const StateVector *self, const Integer siteIndex, Status *status) {
  TitrSite *site;

  if (siteIndex < 0 || siteIndex >= self->nsites) {
    Status_Set (status, Status_IndexOutOfRange);
    return False;
  }
  site = &self->sites[siteIndex];
  return site->isSubstate;
}

Integer StateVector_GetItem (const StateVector *self, const Integer siteIndex, Status *status) {
  TitrSite *site;
  Integer instanceLocalIndex;

  if (siteIndex < 0 || siteIndex >= self->nsites) {
    Status_Set (status, Status_IndexOutOfRange);
    return -1;
  }
  site = &self->sites[siteIndex];
  /* Translate global index to local index */
  instanceLocalIndex = site->indexActive - site->indexFirst;

  return instanceLocalIndex;
}

void StateVector_SetItem (const StateVector *self, const Integer siteIndex, const Integer instanceLocalIndex, Status *status) {
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

void StateVector_SetActualItem (const StateVector *self, const Integer siteIndex, const Integer instanceGlobalIndex, Status *status) {
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

Integer StateVector_GetSubstateItem (const StateVector *self, const Integer index, Status *status) {
  TitrSite *site;

  if (index < 0 || index >= self->nssites) {
    Status_Set (status, Status_IndexOutOfRange);
    return -1;
  }
  site = self->substateSites[index];
  return site->indexSite;
}

void StateVector_SetSubstateItem (const StateVector *self, const Integer selectedSiteIndex, const Integer index, Status *status) {
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

/*=============================================================================
 Incrementation algorithm by Timm Essigke 
=============================================================================*/
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
  /* Return false after reaching the last vector.
    The vector goes back to its initial state.  */
  return False;
}

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

/*=============================================================================
  Monte Carlo-related functions
=============================================================================*/
void StateVector_Move (const StateVector *self, Integer *site, Integer *oldActive, const RandomNumberGenerator *generator) {
  TitrSite *ts;
  Integer randomSite, randomInstance;

  randomSite = RandomNumberGenerator_NextCardinal (generator) % self->nsites;
  ts = &self->sites[randomSite];
  do {
    randomInstance = RandomNumberGenerator_NextCardinal (generator) % (ts->indexLast - ts->indexFirst + 1) + ts->indexFirst;
  } while (randomInstance == ts->indexActive);

  *site           = randomSite;
  *oldActive      = ts->indexActive;
  ts->indexActive = randomInstance;
}

void StateVector_DoubleMove (const StateVector *self, Integer *site, Integer *siteOther, Integer *oldActive, Integer *oldActiveOther, const RandomNumberGenerator *generator) {
  Integer    randomPair, randomInstance;
  TitrSite  *sa, *sb;
  PairSite  *pair;

  randomPair = RandomNumberGenerator_NextCardinal (generator) % self->npairs;
  pair = &self->pairs[randomPair];
  sa   = pair->a;
  sb   = pair->b;

  do {
    randomInstance = RandomNumberGenerator_NextCardinal (generator) % (sa->indexLast - sa->indexFirst + 1) + sa->indexFirst;
  } while (randomInstance == sa->indexActive);
  *site            =  sa->indexSite    ;
  *oldActive       =  sa->indexActive  ;
  sa->indexActive  =  randomInstance   ;

  do {
    randomInstance = RandomNumberGenerator_NextCardinal (generator) % (sb->indexLast - sb->indexFirst + 1) + sb->indexFirst;
  } while (randomInstance == sb->indexActive);
  *siteOther       =  sb->indexSite    ;
  *oldActiveOther  =  sb->indexActive  ;
  sb->indexActive  =  randomInstance   ;
}
