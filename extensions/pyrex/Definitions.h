/*------------------------------------------------------------------------------
! . File      : Definitions.h
! . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
! . Copyright : CEA, CNRS, Martin J. Field (2007-2013)
! . License   : CeCILL French Free Software License     (http://www.cecill.info)
!-----------------------------------------------------------------------------*/
/*==================================================================================================================================
! . This file contains general definitions used throughout the cDynamo program.
!=================================================================================================================================*/
# ifndef _DEFINITIONS
# define _DEFINITIONS

# include <limits.h>
# include <stddef.h>

/* . The standard cardinal and integer need to be at least 32 bits. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Primitive value types.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Default types. */
# include "Boolean.h"
# include "Cardinal.h"
# include "Integer.h"
# include "Real.h"

/* . Specific types. */
typedef short unsigned int Cardinal16 ;
typedef       unsigned int Cardinal32 ;
typedef short          int Integer16  ;
typedef                int Integer32  ;
typedef float              Real32     ;
typedef double             Real64     ;

typedef char               Character  ;
typedef size_t             CSize      ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# include "Macros.h"

# endif
