#ifndef _RMTGENENET_
#define _RMTGENENET_

#include <stdio.h>
#include <string.h>
#include <setjmp.h>

#include "threshold.h"
#include "similarity.h"

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
