#include "PairWiseClusterWriter.h"

/**
 * Constructor
 *
 * @param char ** method
 *   The list of similarity methods used: sc, pc, mi
 * @param char * fp
 *   The string used to prefix all files
 * @param int i
 *   The job index. This will be added to all files.
 * @param int n
 *   The number of samples.
 */
PairWiseClusterWriter::PairWiseClusterWriter(char ** method, int num_methods,
    char * fileprefix, int id, int num_samples) {
  this->job_index = id;
  this->num_samples = num_samples;
  // Set the recovery x and y coordinates to -1 to indicate that the
  // jobs are completed.  If the jobs are not completed then these will
  // be changed to the starting coordinate.
  this->recovery_x = -1;
  this->recovery_y = -1;

  this->method = method;
  this->num_methods = num_methods;

  this->fileprefix = (char *) malloc(sizeof(char) * strlen(fileprefix) + 1);
  strcpy(this->fileprefix, fileprefix);

  fps = (fstream **) malloc(sizeof(fstream *) * 103);

  // Open the files, find out what was the last coordinates used, then
  // set the position in each file to pick up where it left off.
  this->openOutFiles();
  this->findLastPositions();
}
/**
 * Destructor
 */
PairWiseClusterWriter::~PairWiseClusterWriter() {
  // Close and free memory for files.
  this->closeOutFiles();
  // Free strings created in constructor.
  free(fileprefix);
}
// TODO: need a function to make sure the result directory is empty or that
// the same number of job files are present.
/**
 * Opens and creates 103 files for storing the clusters with correlation values.
 * Each file stores a range of 1/100 Spearman correlation values. There is
 * a file for clusters with two few samples, and a file for skipped
 * comparisions that had too many missing values.
 */
void PairWiseClusterWriter::openOutFiles() {

  char clusters_dir[50];
  char skipped_dir[50];
  char nan_dir[50];

  // Make sure the output directory exists.
  sprintf(clusters_dir, "./clusters-%s", method[0]);
  struct stat st = {0};
  if (stat(clusters_dir, &st) == -1) {
    mkdir(clusters_dir, 0700);
  }

  // Open up 102 directories, one each for 100 Spearman/Pearson correlation value
  // ranges, and another for those without (e.g. 'nan').
  int i =  0;
  char filename[1025];
  char dirname[1025];
  for (i = 0; i <= 100; i++) {
    sprintf(dirname, "%s/%03d", clusters_dir, i);
    if (stat(dirname, &st) == -1) {
      mkdir(dirname, 0700);
    }

    // Check to make sure the file exits. If not, then create it . Otherwise,
    // open for input/output so we can check where we may have left off.
    sprintf(filename, "%s/%03d/%s.clusters.%03d.%05d.txt", clusters_dir, i, fileprefix, i, job_index);
    fps[i] = new fstream;
    if (stat(filename, &st) == -1) {
      fps[i]->open(filename, fstream::out);
      fps[i]->close();
    }
    fps[i]->open(filename, fstream::in|fstream::out|fstream::ate);
  }

  // Create the nan directory and files.
  sprintf(nan_dir, "%s/nan", clusters_dir);
  if (stat(nan_dir, &st) == -1) {
    mkdir(nan_dir, 0700);
  }
  sprintf(filename, "%s/%s.clusters.nan.%05d.txt", nan_dir, fileprefix, job_index);
  fps[101] = new fstream;
  if (stat(filename, &st) == -1) {
    fps[101]->open(filename, fstream::out);
    fps[101]->close();
  }
  fps[101]->open(filename, fstream::in|fstream::out|fstream::ate);

  // Create the skipped file that will hold all skipped pair-wise comparisions.
  sprintf(skipped_dir, "%s/skipped", clusters_dir);
  if (stat(skipped_dir, &st) == -1) {
    mkdir(skipped_dir, 0700);
  }
  sprintf(filename, "%s/%s.clusters.skipped.%05d.txt", skipped_dir, fileprefix, job_index);
  fps[102] = new fstream;
  if (stat(filename, &st) == -1) {
    fps[102]->open(filename, fstream::out);
    fps[102]->close();
  }
  fps[102]->open(filename, fstream::in|fstream::out|fstream::ate);
}

/**
 * finds the last complete line written to the file.
 *
 * This will contain the last x,y pair-wise coordinates performed within
 * this file.  We can restart from that position.
 */
void PairWiseClusterWriter::findLastPositions() {

  // Holds the last x and y values that exist in each file before
  // the job terminated without completion.
  int last_x[102];
  int last_y[102];

  // The position in the file where the last_x and last_y are found.
  int last_seek[102];

  // The maximum buffer size is the the number of samples plus the other
  // fields (estimate at 100b) * 10 (for 10 lines)
  unsigned int max_buffer = (num_samples + 100) * 10;

  for (int i = 0; i < 102; i++) {
    int done = 0;
    unsigned int buffer_size = 0;
    char * buffer;
    unsigned int file_size = fps[i]->tellg();

    // initialize the last_x and last_y
    last_x[i] = 0;
    last_y[i] = 0;

    // If we have no data in this file then there's nothing to check, just
    // set the position to 0.
    if (file_size == 0) {
      done = 1;
      last_seek[i] = 0;
    }

    while(!done) {

      // If the buffer size is zero then skip.
      if (buffer_size == 0) {
        buffer_size++;
        continue;
      }

      // If the buffer size reaches the max then we can't recover. So just quit.
      if (buffer_size > max_buffer) {
        done = 1;
        last_x[i] = 0;
        last_y[i] = 0;
        last_seek[i] = 0;
      }

      // Seek backwards from the end of the file
      fps[i]->seekg(file_size - buffer_size);
      buffer = (char *) malloc(sizeof(char) * buffer_size);
      fps[i]->getline(buffer, buffer_size);

      if (buffer_size == 6) {
        if (strcmp(buffer, "#Done") == 0) {
          done = 1;
          // Set last_x and last_y to -1 to indicate the file is completed
          last_x[i] = -1;
          last_y[i] = -1;
          last_seek[i] = file_size;
          free(buffer);
          break;
        }
      }

      // Check the line if we have reached a new line (strlen(buffer) == 0) then
      // process the previous line, or if the file only has one line or we've
      // reached the beginning process the current line.
      if (buffer_size > 1 && (strlen(buffer) == 0 || buffer_size == file_size)) {
        // If we are at the beginning of the file then process this line.
        if (buffer_size == file_size) {
          fps[i]->seekg(file_size - buffer_size);
          done = 1;
        }
        // If we are at the end of the next line then process the previous line.
        else {
          fps[i]->seekg(file_size - (buffer_size - 1));
        }
        fps[i]->getline(buffer, buffer_size);

        // Get the values from this line.
        int x, y, cluster_num, num_clusters, cluster_samples, num_missing, num_outliers, num_goutliers, num_threshold;
        float cv;
        char samples[num_samples];
        int n = sscanf(buffer, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &num_outliers, &num_goutliers, &num_threshold, &cv, (char *) &samples);

        // If the correct number of fields was read from the file and the
        // correct number of samples is present then this is a good line. We
        // will store the x and y coordinates for this file as well
        // as the seek position for this file.
        if (n == 11 && strlen(samples) == (unsigned int) num_samples) {
          done = 1;
          last_x[i] = x;
          last_y[i] = y;

          // To calculate the seek position we subtract the line size from
          // the buffer size but because the buffer size is one larger we
          // must subtract 1. Add 1 to the buffer length to account for the
          // missing \n.
          if (buffer_size == file_size) {
            last_seek[i] = file_size - (buffer_size - (strlen(buffer) + 1));
          }
          else {
            last_seek[i] = file_size - ((buffer_size - 1) - (strlen(buffer) + 1));
          }
        }
      }

      free(buffer);

      // Increment the buffer size.
      buffer_size++;

    } // end while(!done) ...
  } // end for (int i = 0; i < 102; i++) ...

  // Find the largest completed x and y coordinates. The x-coordinate takes
  // precedence, or in other words, we only consider larger values of y if
  // x is also bigger or the same.
  for (int i = 0 ; i < 102; i++) {
    if (recovery_x <= last_x[i]) {
      recovery_x = last_x[i];
      if (recovery_y < last_y[i]) {
        recovery_y = last_y[i];
      }
    }
  }

  // If all the files are done then we just need to return.
  if (recovery_x == -1 || recovery_y == -1) {
    return;
  }

  // Now iterate through the files one more time and move the file pointer
  // to the proper place. If the last_x and last_y is the same as the
  // recovery_x and the recovery_y then we want to move the file pointer to
  // the first occurance of this coordinate (there may be more than one if
  // there are multiple clusters).
  for (int i = 0; i < 102; i++) {

    // If the last x is not the same as the recovery_x then this set the
    // file pointer to the last good read line.
    if (last_x[i] < recovery_x) {
      fps[i]->seekp(last_seek[i]);
      continue;
    }

    // Here the last_x is the same as the recovery_y, but if the last y is
    // not same as the recover_y then set the file pointer to the last
    // good read line.
    if (last_y[i] < recovery_y) {
      fps[i]->seekp(last_seek[i]);
      continue;
    }

    // If we're here it's because the last_x and last_y are the same as the
    // recovery_x and recovery_y.  We want to backup the file pointer if there
    // are more than one cluster for these coordinates. This way the
    // pair-wise comparison can be re-run without having duplicates in the
    // file.
    int done = 0;
    unsigned int buffer_size = 0;
    char * buffer;
    unsigned int file_size;

    // Get the file size.
    fps[i]->seekg(0, ios_base::end);
    file_size = fps[i]->tellg();

    // If the file size is zero then this is a new file and we can skip
    // the following code
    if (file_size == 0) {
      continue;
    }

    // Set the file pointer to that of the last good read line.  This should
    // be the recovery_x and recovery_y coordinate result line.
    buffer_size = file_size - last_seek[i];

    // Read the file backwards to look for the proper line.
    while(!done) {

      // If the buffer size is zero then skip.
      if (buffer_size == 0) {
        buffer_size++;
        continue;
      }

      // Seek backwards from the end of the file
      fps[i]->seekg(file_size - buffer_size);
      buffer = (char *) malloc(sizeof(char) * buffer_size);
      fps[i]->getline(buffer, buffer_size);

      // Check the line if we have reached a new line (strlen(buffer) == 0) then
      // process the previous line, or if the file only has one line or we've
      // reached the beginning process the current line.
      if (buffer_size > 1 && (strlen(buffer) == 0 || buffer_size == file_size)) {
        // If we are at the beginning of the file then process this line.
        if (buffer_size == file_size) {
          fps[i]->seekg(file_size - buffer_size);
        }
        // If we are at the end of the next line then process the previous line.
        else {
          fps[i]->seekg(file_size - (buffer_size - 1));
        }
        fps[i]->getline(buffer, buffer_size);

        // Get the values from this line.
        int x, y, cluster_num, num_clusters, cluster_samples, num_missing;
        float cv;
        char samples[num_samples];
        int n = sscanf(buffer, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &cv, (char *) &samples);

        // If the correct number of fields was read from the file and the
        // correct number of samples is present then this is a good line. If the
        // x and y coordinates are not the same then we are done and we can
        // set the file position.
        if (n == 8 && strlen(samples) == (unsigned int) num_samples) {
          if (x != recovery_x || y != recovery_y) {
            done = 1;
            if (buffer_size == file_size) {
              last_seek[i] = file_size - (buffer_size - (strlen(buffer) + 1));
            }
            else {
              last_seek[i] = file_size - ((buffer_size - 1) - (strlen(buffer) + 1));
            }
          }
        }
      }
      free(buffer);

      // Increment the buffer size.
      buffer_size++;
    }
    fps[i]->seekp(last_seek[i]);
  }


  // Add a comment to each file to indicate where we restarted
  for (int i = 0; i < 102; i++) {

    // If the file size is zero then we are not restarting, otherwise, add
    // a comment to indicate where the job restarted.
    int file_size = fps[i]->tellg();
    if (file_size > 0) {
      (*fps[i]) << "#Restarted" << endl;
    }
  }
}
/**
 * Closes the 102 files that were opened.
 */
void PairWiseClusterWriter::closeOutFiles() {
  int i =  0;

  for (i = 0; i <= 101; i++) {
    (*fps[i]) << "#Done" << endl;
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

  // If there are no clusters then this comparison was skipped.
  if (curr == NULL) {
    fp = fps[102];
    (*fp) << gene1 + 1 << "\t" << gene2 + 1 << endl;
    fp->flush();
  }

  // Iterate through the clusters and print them.
  while (curr != NULL) {
    // Determine which file to write the output into. Use the first similarity
    // method provided for this.
    PairWiseSimilarity * pwsim = curr->getPWSimilarity(0);
    double score = pwsim->getScore();
    if (!pwsim || isnan(score)) {
      fp = fps[101];
    }
    else {
      float i1 = score * 100.0;
      float i2 = fabs(i1);
      int i3 = (int) i2;
      fp = fps[i3];
    }
    (*fp) << gene1 + 1 << "\t" << gene2 + 1 << "\t" << curr->index << "\t" << pwcl->num_clusters << "\t" << curr->cluster_size << "\t" << curr->num_missing  << "\t" << curr->num_outliers << "\t" << curr->num_goutliers << "\t" << curr->num_threshold <<  "\t";
    // Add in the scores, these should be separated by a comma.
    int i;
    char score_str[1024] = "";
    for (i = 0; i < this->num_methods; i++) {
      pwsim = curr->getPWSimilarity(i);
      score = curr->getPWSimilarity(i)->getScore();
      if (i == 0) {
        sprintf(score_str, "%f", score);
      }
      else {
        sprintf(score_str, "%s,%f", score_str, score);
      }
    }
    if (i > 0) {
      (*fp) << score_str << "\t";
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
