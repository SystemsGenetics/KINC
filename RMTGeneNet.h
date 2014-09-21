#ifndef _RMTGENENET_
#define _RMTGENENET_

/**
 * Description:
 * ------------
 * Calculates a Pearson correlation or Mutual Information matrix from an n x m
 * expression matrix where the columns are the samples and rows are
 * the genes or probesets. Cells represent expression levels. The first
 * column of the expression matrix should be the name of the gene or
 * probeset.
 *
 *
 * Change Log / History:
 * ----------
 * 6/2014 (by Stephen Ficklin)
 * 1) Added support for calculation of mutual information matrix
 * 2) Use getopt_long for --xxxx and -x style input arguments
 * 3) Re-organized the code
 *
 * 10/2014 (by Stephen Ficklin)
 * 1) Input files no longer require an implied .txt extesion but the full filename
 *    should be specified
 * 2) Removed the spaces in the output file name and renamed file.
 * 3) Incorporated GSL Pearson calculation function because the pre-calculating
 *    of sum of row and sum of squared of the row couldn't account for missing values
 *
 * 10/2011
 * Created by Scott Gibson at Clemson University under direction of
 * Dr. Melissa Smith in Collaboration with Alex Feltus, Feng Luo and Stephen
 * Ficklin
 */

#include <stdio.h>
#include <string.h>
#include <setjmp.h>

#include "threshold.h"
#include "similarity.h"
#include "extract.h"

/**
 * Definitions for mimicing a try, catch block.
 */
#define TRY do{ jmp_buf ex_buf__; if( !setjmp(ex_buf__) ){
#define CATCH } else {
#define ETRY } }while(0)
#define THROW longjmp(ex_buf__, 1)


/**
 * Function prototypes
 */
void print_usage();

#endif
