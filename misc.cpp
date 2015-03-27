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
    int status = fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
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
