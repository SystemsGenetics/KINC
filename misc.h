#ifndef _MISC_
#define _MISC_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    unsigned long size, resident, share, text, lib, data, dt;
} statm_t;


statm_t * memory_get_usage();

int is_numeric(char * string);

#endif
