#ifndef _MISC_
#define _MISC_

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    unsigned long size, resident, share, text, lib, data, dt;
} statm_t;


statm_t * memory_get_usage();

#endif
