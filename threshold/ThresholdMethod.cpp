#include "ThresholdMethod.h"



/**
 * DRArgs constructor.
 */
ThresholdMethod::ThresholdMethod(EMatrix *ematrix, char * method) {
  this->ematrix = ematrix;
  this->method = method;

  bin_dir = (char *) malloc(sizeof(char) * strlen(ematrix->getInfileName()));
  if (strcmp(method, "mi") == 0) {
    strcpy(bin_dir, "MI");
  }
  else if (strcmp(method, "pc") == 0) {
    strcpy(bin_dir, "Pearson");
  }
  else if (strcmp(method, "sc") == 0) {
    strcpy(bin_dir, "Spearman");
  }
}

/**
 * DRArgs destructor.
 */
ThresholdMethod::~ThresholdMethod() {
  free(bin_dir);
}
