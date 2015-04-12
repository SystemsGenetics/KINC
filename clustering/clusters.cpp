#include "clusters.h"


/**
 * Constructor.
 */
PairWiseCluster::PairWiseCluster(PairWiseSet * pwset) {

  // Initialize the class members.
  this->cluster_samples = NULL;
  this->cluster_size = 0;
  this->pwset = pwset;
  this->neighbor = NULL;
  this->pwsim = NULL;
}
/**
 * Desctructor.
 */
PairWiseCluster::~PairWiseCluster() {
  free(this->cluster_samples);
  if (this->pwsim) {
    delete this->pwsim;
  }
}

void PairWiseCluster::doSimilarity(const char * method, int min_obs) {
  // Create a new PairWiseSet for this cluster
  if (strcmp(method, "sc") == 0) {
    SpearmanSimilarity * ssim = new SpearmanSimilarity(this->pwset, this->cluster_samples, min_obs);
    ssim->run();
    this->pwsim = (PairWiseSimilarity *) ssim;
  }
  else if (strcmp(method, "pc") == 0) {
    PearsonSimilarity * psim = new PearsonSimilarity(this->pwset, this->cluster_samples, min_obs);
    psim->run();
    this->pwsim = (PairWiseSimilarity *) psim;
  }
  else if (strcmp(method, "mi") == 0) {
//    MISimilarity * msim = new MISimilarity(this->pwset, this->cluster_samples, min_obs);
//    msim->run();
//    this->pwsim = (PairWiseSimilarity *) msim;
  }
}
/**
 * Sets the PairWiseCluster samples.
 */
void PairWiseCluster::setClusterSamples(int * samples, bool from_clean) {
  this->cluster_samples = (int *) malloc(sizeof(int) * this->pwset->n_orig);

  // If the samples list is derive from the clean samples set then the size
  // of the samples is pwset->n_clean.
  int k = 0;
  if (from_clean) {
    for (int i = 0; i < this->pwset->n_orig; i++) {
      if (this->pwset->samples[i] == 1) {
        this->cluster_samples[i] = samples[k];
        k++;
      }
      else {
        //this->cluster_samples[i] = 0;
        this->cluster_samples[i] = 2;
      }
    }
  }
  // The provided samples are not derived from the clean set so just copy them
  else {
    for (int i = 0; i < this->pwset->n_orig; i++) {
      this->cluster_samples[i] = samples[i];
    }
  }
}
/**
 *
 */
void PairWiseCluster::printCluster() {
  printf("%i\t%i\t", this->pwset->gene1, this->pwset->gene2);
  if (this->cluster_samples) {
    for (int i = 0; i < this->pwset->n_orig; i++) {
      printf("%i", this->cluster_samples[i]);
    }
  }
  printf("\n");
}

PairWiseClusterList::PairWiseClusterList(PairWiseSet * pwset) {
  this->pwset = pwset;
  this->num_clusters = 0;
  this->head = NULL;
}

void PairWiseClusterList::addCluster(PairWiseCluster * pwc) {

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
PairWiseClusterWriter::PairWiseClusterWriter(char * m, char * fp, int i) {
  id = i;
  method = (char *) malloc(sizeof(char) * strlen(m) + 1);
  strcpy(method, m);

  fileprefix = (char *) malloc(sizeof(char) * strlen(fp) + 1);
  strcpy(fileprefix, fp);
  fps = (FILE **) malloc(sizeof(FILE *) * 102);

  this->openOutFiles();
}
/**
 * Destructor
 */
PairWiseClusterWriter::~PairWiseClusterWriter() {
  this->closeOutFiles();
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
void PairWiseClusterWriter::writeClusters(PairWiseClusterList *pwcl, int gene1, int gene2) {

  // The file pointer of the file to write to.
  FILE *fp;
  int cluster_id = 0;

  PairWiseCluster * curr = pwcl->head;
  while (curr != NULL) {
    // Determine which file to write the output into
    if (isnan(curr->pwsim->score)) {
      fp = fps[101];
    }
    else {
      float i1 = curr->pwsim->score * 100.0;
      float i2 = fabs(i1);
      int i3 = (int) i2;
      fp = fps[i3];
    }
    fprintf(fp, "%i\t%i\t%i\t%i\t", gene1 + 1, gene2 + 1, curr->cluster_size, cluster_id);
    for (int i = 0; i < curr->pwset->n_orig; i++) {
      fprintf(fp, "%i", curr->cluster_samples[i]);
    }
    if (curr->pwsim) {
      fprintf(fp, "\t%f", curr->pwsim->score);
    }
    else {
      fprintf(fp, "\t%f", NAN);
    }
    fprintf(fp, "\n");
    curr = curr->neighbor;
    cluster_id++;
  }
  fflush(fp);
}
