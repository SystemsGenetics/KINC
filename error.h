#ifndef _ERROR_
#define _ERROR_

#include <stdio.h>
#include <setjmp.h>
#include <stdlib.h>


void handle_error(char * message);
void handle_warning(char * message);

#endif
