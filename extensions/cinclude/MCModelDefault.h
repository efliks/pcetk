/*------------------------------------------------------------------------------
! . File      : MCModelDefault.h
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
!                          Mikolaj J. Feliks (2014-2015)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
#ifndef _MCMODELDEFAULT
#define _MCMODELDEFAULT

/* Needed for random seed */
#include <time.h>
/* Needed for exp */
#include <math.h>

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
#include "RandomNumberGenerator.h"

/* Own modules */
#include "StateVector.h"
#include "EnergyModel.h"


/* Taken from GMCT */
#define TOO_SMALL -500.0

typedef struct {
    /* Energy limit for double moves */
    Real                   limit;
    /* Number of equlibration scans */
    Integer                nequil;
    /* Number of production scans */
    Integer                nprod;
    /* Pointer to the energy model */
    EnergyModel           *energyModel;
    /* Private state vector of the Monte Carlo model */
    StateVector           *vector;
    /* Mersenne Twister generator */
    RandomNumberGenerator *generator;
} MCModelDefault;


/* Allocation and deallocation */
extern MCModelDefault *MCModelDefault_Allocate          (const Real limit, const Integer nequil, const Integer nprod, const Integer randomSeed, Status *status);
extern void            MCModelDefault_Deallocate        (MCModelDefault *self);
extern void            MCModelDefault_LinkToEnergyModel (MCModelDefault *self, EnergyModel *energyModel, Status *status);

/* Monte Carlo-related functions */
extern Boolean MCModelDefault_Metropolis             (const Real GdeltaRT, const RandomNumberGenerator *generator);
extern Boolean MCModelDefault_Move                   (const MCModelDefault *self, const Real pH, const Real G, Real *Gnew);
extern Boolean MCModelDefault_DoubleMove             (const MCModelDefault *self, const Real pH, const Real G, Real *Gnew);
extern Real    MCModelDefault_MCScan                 (const MCModelDefault *self, const Real pH, Integer nmoves);
extern void    MCModelDefault_UpdateProbabilities    (const MCModelDefault *self);
extern Real    MCModelDefault_FindMaxInteraction     (const MCModelDefault *self, const TitrSite *site, const TitrSite *other);
extern Integer MCModelDefault_FindPairs              (const MCModelDefault *self, const Integer npairs, Status *status);
extern void    MCModelDefault_CalculateProbabilities (const MCModelDefault *self, const Real pH, const Boolean equil);

#endif
