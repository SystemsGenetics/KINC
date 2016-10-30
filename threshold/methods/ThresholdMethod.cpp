#include "ThresholdMethod.h"



/**
 * DRArgs constructor.
 */
ThresholdMethod::ThresholdMethod(EMatrix *ematrix, char ** method, int num_methods,
    char * th_method, char * clustering, int min_cluster_size, int max_missing,
    int max_modes, float min_range) {

  this->ematrix = ematrix;
  this->method = method;
  this->clustering = clustering;
  this->min_cluster_size = min_cluster_size;
  this->max_missing = max_missing;
  this->max_modes = max_modes;
  this->min_range = min_range;
  this->num_methods = num_methods;
  this->th_method = th_method;

  // Find the index of the th_method in the methods array
  for (int i = 0; i < this->num_methods; i++) {
    if (strcmp(method[i], th_method) == 0) {
      this->th_method_index = i;
    }
  }

  // For the binary file format:
  bin_dir = (char *) malloc(sizeof(char) * strlen(ematrix->getInfileName()));
  if (strcmp(th_method, "mi") == 0) {
    strcpy(bin_dir, "MI");
  }
  else if (strcmp(th_method, "pc") == 0) {
    strcpy(bin_dir, "Pearson");
  }
  else if (strcmp(th_method, "sc") == 0) {
    strcpy(bin_dir, "Spearman");
  }
}

/**
 * DRArgs destructor.
 */
ThresholdMethod::~ThresholdMethod() {
  free(bin_dir);
}

/**
 *
 */
float ** ThresholdMethod::parseScores(char * scores_str) {
  // Split the method into as many parts
  char * tmp;
  int i = 0;
  tmp = strstr(scores_str, ",");
  float ** scores = (float **) malloc(sizeof(float *) * this->num_methods);
  char tmp_score[255];

  while (tmp) {
    strncpy((char *) &tmp_score, scores_str, (tmp - scores_str));
    scores[i] = (float *) malloc(sizeof(float) * 1);
    *(scores[i]) = atof(tmp_score);
    scores_str = tmp + 1;
    tmp = strstr(scores_str, ",");
    i++;
  }
  // Get the last element of the methods_str.
  strncpy((char *) &tmp_score, scores_str, strlen(scores_str));
  scores[i] = (float *) malloc(sizeof(float) * 1);
  *(scores[i]) = atof(tmp_score);

  return scores;
}
