/************************************************************************/
/**

   \file       sstruc.h
   
   \version    V1.0
   \date       10.07.15
   \brief      Header for secondary structure calculation
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1988-2015
   \author     Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College London,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is NOT IN THE PUBLIC DOMAIN, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <math.h>
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"

/************************************************************************/
/* Defines and macros
*/

/* Character symbols for secondary structure types                      */
#define SECSTR_COIL            '-'

#define SECSTR_BEND_START      '>'
#define SECSTR_BEND_END        '<'
#define SECSTR_BEND_BOTH       'X'

#define SECSTR_BEND            'S'
#define SECSTR_BRIDGE_FWD      'b'
#define SECSTR_BRIDGE_BACKWD   'B'
#define SECSTR_SHEET           'E'
#define SECSTR_SHEET_SMALL     'e'
#define SECSTR_TURN            'T'
#define SECSTR_BREAK           '!'
#define SECSTR_BULGE           '*'
#define SECSTR_HELIX           'H'
#define SECSTR_HELIX_SMALL     'h'
#define SECSTR_PI              'I'
#define SECSTR_PI_SMALL        'i'
#define SECSTR_3_10            'G'
#define SECSTR_3_10_SMALL      'g'

#define SECSTR_ERR_NOERR       0
#define SECSTR_ERR_NOMEM       (-1)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int blCalcSecStrucPDB(PDB *pdbStart, PDB *pdbStop, BOOL quiet);
