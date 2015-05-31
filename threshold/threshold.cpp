#include "threshold.h"



/**
 * The function to call when running the 'similarity' command.
 */
int do_threshold(int argc, char *argv[]) {

  RMTThreshold * rmt = new RMTThreshold(argc, argv);

  // Find the RMT threshold.
  rmt->findThreshold();

  printf("Done.\n");
  return 1;
}


