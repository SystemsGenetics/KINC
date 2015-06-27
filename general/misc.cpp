#include "misc.h"

/**
 * Returns a statm_t object containing information about memory usage.
 */
statm_t * memory_get_usage() {

  const char* statm_path = "/proc/self/statm";

  statm_t * result = (statm_t *) malloc(sizeof(statm_t));
  result->size = 0;
  result->resident = 0;
  result->share = 0;
  result->text = 0;
  result->lib = 0;
  result->data = 0;
  result->dt = 0;

  FILE *f = fopen(statm_path,"r");
  if(f){
    fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
      &result->size,
      &result->resident,
      &result->share,
      &result->text,
      &result->lib,
      &result->data,
      &result->dt
    );
    fclose(f);
  }
  return result;
}

/**
 *
 */
int is_numeric(char * string) {
   // ASCII value of
   // '-': 45
   // '.': 46
   // '0': 48
   // '9': 57
   // 'e': 101
   int num_period = 0;
   unsigned int i;
   int good = 0;
   // the position of the e in scientific notation
   unsigned int epos = 0;

   // all remaining characters can be numeric but with only one period
   for (i = 0; i < strlen(string); i++) {
     good = 0;
     if (string[i] >= 48 && string[i] <= 57) {
       good = 1;
     }
     // periods are acceptable, but only one.
     if (string[i] == 46) {
       num_period++;
       good = 1;
       // we can only have one period
       if (num_period > 1) {
         good = 0;
       }
     }
     // the minus sign is acceptable if it appears first in the string or
     // after an e in scientific notation
     if (string[i] == 45) {
       if (i == 0) {
         good = 1;
       }
       if (i == epos + 1) {
         good = 1;
       }
     }
     if (string[i] == 101) {
       good = 1;
       epos = i;
     }
     if (!good) {
       return 0;
     }
   }
   return 1;
}
