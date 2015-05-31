#include "PairWiseClusterWriter.h"

/**
 * Constructor
 */
PairWiseClusterWriter::PairWiseClusterWriter(char * m, char * fp, int i) {
  id = i;
  method = (char *) malloc(sizeof(char) * strlen(m) + 1);
  strcpy(method, m);

  fileprefix = (char *) malloc(sizeof(char) * strlen(fp) + 1);
  strcpy(fileprefix, fp);
  fps = (ofstream **) malloc(sizeof(ofstream *) * 102);

  this->openOutFiles();
}
/**
 * Destructor
 */
PairWiseClusterWriter::~PairWiseClusterWriter() {
  this->closeOutFiles();
  free(fps);
  free(fileprefix);
  free(method);
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
    fps[i] = new ofstream;
    fps[i]->open(filename, ios::out);
  }

  if (stat(nan_dir, &st) == -1) {
    mkdir(nan_dir, 0700);
  }
  sprintf(filename, "%s/%s.clusters.nan.%03d.txt", nan_dir, fileprefix, id + 1);
  fps[i] = new ofstream;
  fps[i]->open(filename, ios::out);
}

/**
 * Closes the 102 files that were opened.
 */
void PairWiseClusterWriter::closeOutFiles() {
  int i =  0;

  for (i = 0; i <= 101; i++) {
    (*fps[i]) << "#Done" << "\n";
    fps[i]->flush();
    fps[i]->close();
  }
  free(fps);
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
void PairWiseClusterWriter::writeClusters(PairWiseClusterList *pwcl, int gene1, int gene2) {

  // The file pointer of the file to write to.
  ofstream *fp;

  PairWiseCluster * curr = pwcl->head;
  while (curr != NULL) {
    // Determine which file to write the output into
    double score = curr->pwsim->getScore();
    if (!curr->pwsim || isnan(score)) {
      fp = fps[101];
    }
    else {
      float i1 = score * 100.0;
      float i2 = fabs(i1);
      int i3 = (int) i2;
      fp = fps[i3];
    }
    (*fp) << gene1 + 1 << "\t" << gene2 + 1 << "\t" << curr->index << "\t" << pwcl->num_clusters << "\t" << curr->cluster_size << "\t" << curr->num_missing  << "\t";
    if (curr->pwsim) {
      (*fp) << score << "\t";
    }
    else {
      (*fp) << NAN << "\t";
    }
    for (int i = 0; i < curr->pwset->n_orig; i++) {
      (*fp) << curr->cluster_samples[i];
    }
    (*fp) << endl;
    curr = curr->neighbor;
    fp->flush();
  }
}
