#include "PairWiseClusterWriter.h"

/**
 * Constructor
 *
 * @param char * m
 *   The method used for similarity: sc, pc, mi
 * @param char * fp
 *   The string used to prefix all files
 * @param int i
 *   The job index. This will be added to all files.
 * @param int n
 *   The number of samples.
 */
PairWiseClusterWriter::PairWiseClusterWriter(char * method, char * fileprefix, int id, int num_samples) {
  this->job_index = id;
  this->num_samples = num_samples;

  this->method = (char *) malloc(sizeof(char) * strlen(method) + 1);
  strcpy(this->method, method);

  this->fileprefix = (char *) malloc(sizeof(char) * strlen(fileprefix) + 1);
  strcpy(this->fileprefix, fileprefix);

  fps = (fstream **) malloc(sizeof(fstream *) * 102);

  last_x = (int *) malloc(sizeof(int) * 102);
  last_y = (int *) malloc(sizeof(int) * 102);
  last_seek = (int *) malloc(sizeof(int) * 102);

  // Open the files, find out what was the last coordinates used, then
  // set the position in each file to pick up where it left off.
  this->openOutFiles();
  this->findLastPositions();
}
/**
 * Destructor
 */
PairWiseClusterWriter::~PairWiseClusterWriter() {
  this->closeOutFiles();
  free(fps);
  free(fileprefix);
  free(method);
  free(last_seek);
  free(last_x);
  free(last_y);
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
    sprintf(filename, "%s/%03d/%s.clusters.%03d.%03d.txt", clusters_dir, i, fileprefix, i, job_index);
    fps[i] = new fstream;
//    fps[i]->open(filename, ios::out);
    fps[i]->open(filename, ios::in|ios::out|ios::ate);
  }

  if (stat(nan_dir, &st) == -1) {
    mkdir(nan_dir, 0700);
  }
  sprintf(filename, "%s/%s.clusters.nan.%03d.txt", nan_dir, fileprefix, job_index);
  fps[i] = new fstream;
  fps[i]->open(filename, ios::in|ios::out|ios::ate);
}

/**
 *
 */
void PairWiseClusterWriter::findLastPositions() {
  // Get the last complete line written to the file. This will contain the
  // last x,y pair-wise coordinates performed within this file.  We can
  // restart from that position.

  for (int i = 0; i < 102; i++) {
    int done = 0;
    unsigned int buffer_size = 0;
    char * buffer;
    unsigned int file_size = fps[i]->tellg();

    while(!done) {
      // If our buffer size is the same size as the file then we've
      // reached the beginning of the file and we didn't find anything
      // so quit.
      if (buffer_size == file_size) {
        done = 1;
        last_x[i] = 0;
        last_y[i] = 0;
      }
      buffer_size++;

      // Seek backwards from the end of the file
      fps[i]->seekg(file_size - buffer_size);
      buffer = (char *) malloc(sizeof(char) * buffer_size);
      fps[i]->getline(buffer, buffer_size);

      if (buffer_size == 5) {
      //      if (strcmp(buffer, "#Done") == 0) {
      //        done = 1;
      //        // Set last_x and last_y to -1 to indicate the file is completed
      //        last_x = -1;
      //        last_y = -1;
      //        break;
      //      }
      }

      // When we find a new line, add one to the file position to get past
      // the new line and check to see if the line format matches.
      if (buffer_size > 1 && strlen(buffer) == 0) {
        fps[i]->seekg(file_size - (buffer_size - 1));
        fps[i]->getline(buffer, buffer_size);

        // Get the values from this line.
        int x, y, cluster_num, num_clusters, cluster_samples, num_missing;
        float cv;
        char samples[num_samples];
        int n = sscanf(buffer, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &cv, (char *) &samples);

        // If the correct number of fields was read from the file and the
        // correct number of samples is present then this is a good line. We
        // will store the x and y coordinates for this file as well
        // as the seek position for this file.
        if (n == 8 && strlen(samples) == (unsigned int) num_samples) {
          done = 1;
          last_x[i] = x;
          last_y[i] = y;
          last_seek[i] = file_size - (buffer_size - 1);
        }
      }

      free(buffer);
    }
  }
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
  fstream *fp;

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
