#include "clusters.h"


/**
 * Initializes the head PairWiseCluster object.
 */
PairWiseCluster::PairWiseCluster() {

  // Initialize the class members.
  samples = NULL;
  num_samples = 0;
  cluster_size = 0;
  score = 0;
  neighbor = NULL;
}

PairWiseClusterSet::PairWiseClusterSet(int g1, int g2, const char * m) {
  gene1 = g1;
  gene2 = g2;
  num_clusters = 0;
  method = (char *) malloc(sizeof(char) * strlen(m) + 1);
  strcpy(method, m);
  head = NULL;
}

void PairWiseClusterSet::addCluster(PairWiseCluster * pwc) {

  // Check the list to see if it is empty. If so, then make this item the
  // new head.
  if (head == NULL) {
    head = pwc;
    return;
  }

  // Traverse the list to the end and then add the new item
  PairWiseCluster * curr = head;
  while (curr->neighbor != NULL) {
    curr = curr->neighbor;
  }
  curr->neighbor = pwc;
}

/**
 * Constructor
 */
PairWiseClusterWriter::PairWiseClusterWriter(const char * m, const char * fp, int i) {
  id = i;
  method = (char *) malloc(sizeof(char) * strlen(m) + 1);
  strcpy(method, m);

  fileprefix = (char *) malloc(sizeof(char) * strlen(fp) + 1);
  strcpy(fileprefix, fp);
  fps = (FILE **) malloc(sizeof(FILE *) * 102);
}
/**
 * Destructor
 */
PairWiseClusterWriter::~PairWiseClusterWriter() {
  free(fps);
  free(fileprefix);
}
/**
 * Opens and creates 102 files for storing the clusters with correlation values.
 * Each file stores a range of 1/100 Spearman correlation values.
 */
void PairWiseClusterWriter::openOutFiles() {

  char clusters_dir[50];
  char nan_dir[50];

  sprintf(clusters_dir, "./clusters-%s", method);
  sprintf(nan_dir, "%s/nan", clusters_dir);

  // Make sure the output directory exists.
  struct stat st = {0};
  if (stat(clusters_dir, &st) == -1) {
      mkdir(clusters_dir, 0700);
  }

  // Open up 102 files, one each for 100 Spearman correlation value ranges.
  // and another for those without (e.g. 'nan').
  int i =  0;
  char filename[1025];
  char dirname[1025];
  for (i = 0; i <= 100; i++) {
    sprintf(dirname, "%s/%03d", clusters_dir, i);
    if (stat(dirname, &st) == -1) {
      mkdir(dirname, 0700);
    }
    sprintf(filename, "%s/%03d/%s.clusters.%03d.%03d.txt", clusters_dir, i, fileprefix, i, id + 1);
    fps[i] = fopen(filename, "w");
  }

  if (stat(nan_dir, &st) == -1) {
    mkdir(nan_dir, 0700);
  }
  sprintf(filename, "%s/%s.clusters.nan.%03d.txt", nan_dir, fileprefix, id + 1);
  fps[i] = fopen(filename, "w");
}

/**
 * Closes the 102 files that were opened.
 */
void PairWiseClusterWriter::closeOutFiles() {
  int i =  0;
  for (i = 0; i <= 101; i++) {
    FILE * fp = fps[i];
    fprintf(fp, "#Done\n");
    fclose(fp);
  }
}


/**
 * Adds a line to the clustering file.
 *
 * The clustering file is used by KINC during pair-wise correlation analysis
 * to restrict which samples are used.  The file is tab delimited.
 * The format of the file is tab delimited with the following columns:
 *
 *   1)  gene 1 name
 *   2)  gene 2 name
 *   3)  cluster name.  A 0 indicates no clustering was performed.
 *   4)  a string of 0 and 1s indicating which samples to include when
 *       performing pair-wise comparisons.
 *
 * @param SampleCluster pws
 */
void PairWiseClusterWriter::writeClusters(PairWiseClusterSet *pwcs) {

  // The file pointer of the file to write to.
  FILE *fp;
  int cluster_id = 0;

  PairWiseCluster * curr = pwcs->head;
  while (curr != NULL) {
    // Determine which file to write the output into
    if (isnan(curr->score)) {
      fp = fps[101];
    }
    else {
      float i1 = curr->score * 100.0;
      float i2 = fabs(i1);
      int i3 = (int) i2;
      fp = fps[i3];
    }
    fprintf(fp, "%i\t%i\t%i\t%i\t", pwcs->gene1, pwcs->gene2, curr->cluster_size, cluster_id);
    for (int i = 0; i < curr->num_samples; i++) {
      fprintf(fp, "%i", curr->samples[i]);
    }
    fprintf(fp, "\t%f", curr->score);
    fprintf(fp, "\n");
    curr = curr->neighbor;
    cluster_id++;
  }
  fflush(fp);
}


/**
 * Updates the samples vector of a PairWiseClusters object.
 *
 * Because the clustering() function is recursive it is successively called
 * with smaller and smaller sample sets.  When it returns it provides  PWC
 * object with a samples array.  In the samples array, samples that are present
 * in the cluster are marked with a 1 and those not in the cluster are
 * set to 0.  The order of values corresponds to the sample order. The samples
 * array, however, only contains 1's and 0's for the samples provided to the
 * clustering() function. Therefore, the results need to be merged back into
 * the larger samples array.  This function does that.
 *
 * @param int *parent_samples
 *   The list of parent samples. It contains a list of 0's and 1's indicating
 *   which samples are to be kept.  Any samples not also found in the
 *   'new.samples' argument are set to 0.
 * @param int n
 *   The size of the parent_samples array.
 * @param new
 *   The new PWC object.
 */
/*
void update_pairwise_cluster_samples(int * parent_samples, int n, PairWiseClusters * head) {

  PairWiseClusters * curr = (PairWiseClusters *) head;
  PairWiseClusters * next = (PairWiseClusters *) head->next;

  while (curr != NULL) {
    int z;
    int w = 0;

    // Create a new kept array that will update the current samples array.
    int * new_samples = (int *) malloc(sizeof(int) * n);

    // Iterate through all of the elements of the parent samples list.
    for (z = 0; z < n; z++) {
      // If the element is 1, meaning the sample is present int he parent
      // cluster, then check the current sample to see if was preserved during
      // sub clustering.
      if (parent_samples[z] == 1) {
        if (curr->samples[w] == 0) {
          new_samples[z] = 0;
        }
        // The element is kept in the pwc so preserve the 1 it in the new array.
        else {
          new_samples[z] = 1;
        }
        w++;
      }
      // The element is not kept originally so preserve the 0 in the new array.
      else {
        new_samples[z] = 0;
      }
    }
    // Free up the old samples array and replace it with a new one.
    free(curr->samples);
    curr->samples = new_samples;
    curr->num_samples = n;

    curr = next;
    if (next != NULL) {
      next = (PairWiseClusters *) next->next;
    }
  }
}
*/
