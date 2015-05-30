#include "mixmod.h"

MixtureModelPWClustering::MixtureModelPWClustering(int argc, char *argv[])
  : PairWiseClustering(argc, argv) {

  // Initialize some values.
  max_clusters = 5;
  strcpy(criterion, "BIC");

  // The value returned by getopt_long.
  int c;

  // Loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop.
  while(1) {
    int option_index = 0;

    // Specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below.
    const struct option long_options[] = {
      {"help",         no_argument,       0,  'h' },
      {"criterion",    required_argument, 0,  't' },
      {"max_clusters", required_argument, 0,  'l' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };
    // get the next option
    c = getopt_long(argc, argv, "e:r:c:m:o:n:f:h", long_options, &option_index);

    // if the index is -1 then we have reached the end of the options list
    // and we break out of the while loop
    if (c == -1) {
      break;
    }

    // handle the options
    switch (c) {
      case 0:
        break;
      case 't':
        strcpy(criterion, optarg);
        break;
      case 'l':
        max_clusters = atoi(optarg);
        break;
    }
  }

  // TODO: make sure means shift bandwidth arguments are numeric between 1 and 0


  // Make sure the mixture module criterion are good
  bool $mmc_is_good = false;
  if (strcmp(criterion, "BIC") == 0) {
    $mmc_is_good = true;
  }
  else if (strcmp(criterion, "ICL") == 0) {
    $mmc_is_good = true;
  }
  else if (strcmp(criterion, "NEC") == 0) {
    $mmc_is_good = true;
  }
  else if (strcmp(criterion, "CV") == 0) {
    $mmc_is_good = true;
  }
  else if (strcmp(criterion, "DCV") == 0) {
    $mmc_is_good = true;
  }
  if (!$mmc_is_good) {
    fprintf(stderr, "Error: The mixture model criterion must be one of: BIC, ICL, NEC, CV, DCV (--criterion option).\n");
    exit(-1);
  }
  if (max_clusters < 2 || max_clusters > 10 ){
    fprintf(stderr, "Error: Please select a maximum cluters between 2 and 10. (--max_clusters option).\n");
    exit(-1);
  }

  printf("  Mixture model criterion: %s\n", criterion);
  printf("  Max clusters: %d\n", max_clusters);

}

/**
 *
 */
MixtureModelPWClustering::~MixtureModelPWClustering() {

}


/**
 *
 */
void MixtureModelPWClustering::run() {
  // Register signal and signal handler
  // signal(SIGINT, signal_callback_handler);

  // variables used for timing of the code.
  time_t start_time = time(0);
  time_t now;

  // If KINC is being rung using MPI then those values take precedence.
//  if (mpi_num_procs > 1) {
//    num_jobs = mpi_num_procs;
//    job_index = mpi_id;
//  }

  // Calculate the total number of comparisons and how many will be
  // performed by this process. We subtract 1 from the first params->rows
  // because we do not calculate the diagonal.
  long long int num_rows = ematrix->getNumGenes() - 1;
  long long int total_comps  = num_rows * (num_rows + 1) / 2;
  long long int comps_per_process = total_comps / num_jobs;
  long long int comp_start = job_index * comps_per_process;
  long long int comp_stop = job_index * comps_per_process + comps_per_process;

  // If this is the last process and there are some remainder comparisons
  // then we need to add them to the stop
  if (job_index + 1 == num_jobs) {
    comp_stop = total_comps;
  }
  printf("  Job %d of %d. \n", job_index + 1, num_jobs);
  printf("  Performing comparisions %lld to %lld (%lld) of %lld\n", comp_start, comp_stop, comp_stop - comp_start, total_comps);
  fflush(stdout);

  // Creat the writer object to write out the cluters.
  PairWiseClusterWriter * pwcw = new PairWiseClusterWriter(method, ematrix->getFilePrefix(), job_index);

  // Perform the pair-wise clustering.
  int n_comps = 0;
  int my_comps = 0;

//  for (int i = 122 - 1; i < num_rows; i++) {
  for (int i = 0; i < num_rows; i++) {
    /*if (i == 50) {
      break;
    }*/
    for (int j = 0; j < num_rows; j++) {
//    for (int j = 77 - 1; j < num_rows; j++) {

      // We only need to calculate clusters in the lower triangle of the
      // full pair-wise matrix
      if (j >= i) {
        continue;
      }

      // If this computation is not meant for this process, then skip it.
      if (n_comps < comp_start || n_comps >= comp_stop) {
        n_comps++;
        continue;
      }

      // Perform pairwise clustering using mixture models
      PairWiseSet * pwset = new PairWiseSet(ematrix, i, j);
      MixtureModelPWClusters * mixmod = new MixtureModelPWClusters(pwset, min_obs, method);
      mixmod->run(criterion, max_clusters);
      pwcw->writeClusters(mixmod->pwcl, i, j);

      // Print run stats.
      if (my_comps % 100 == 0) {
        // Get the amount of memory used.
        statm_t * memory = memory_get_usage();

        // Get the percent completed.
        double percent_complete = (my_comps / (float) (comp_stop - comp_start)) * 100;

        // Calculate the time left to complete
        now = time(0);
        double seconds_passed = now - start_time;
        double seconds_per_comp = seconds_passed / my_comps;
        double total_seconds = seconds_per_comp * (comp_stop - comp_start);

        double minutes_left = (total_seconds - seconds_passed) / 60;
        double hours_left = minutes_left / 60;
        double days_left = hours_left / 24;
        double years_left = days_left / 365;

        // Write progress report.
        printf("%d. Complete: %.4f%%. Mem: %ldb. Remaining: %.2fh; %.2fd; %.2fy. Coords: %d, %d.        \n",
          job_index + 1, (float) percent_complete, memory->size, (float) hours_left, (float) days_left, (float) years_left, i, j);
        free(memory);
      }
      n_comps++;
      my_comps++;
      delete mixmod;
      delete pwset;

      // If the user cancelled the program then finish up.
      if (Close_Now == 1) {
        i = num_rows;
        j = num_rows;
      }
//      exit(1);
    }
  }
  delete pwcw;
}

/**
 * Constructor for MixMod class.
 *
 * @param EMatrix *ematrix;
 */
MixtureModelPWClusters::MixtureModelPWClusters(PairWiseSet *pwset, int min_obs, char * method) {

  // Initialize some values.
  this->pwset = pwset;
  this->min_obs = min_obs;
  this->method = method;
  labels = NULL;
  pwcl = NULL;
  data = NULL;

  // Create the PairWiseClusterList object
  pwcl = new PairWiseClusterList(pwset);

  // Make sure we have the correct number of observations before preparing
  // for the comparision.
  if (this->pwset->n_clean < this->min_obs) {
    return;
  }

  // Create the Gaussian Data object and set the dataDescription object.
  this->data = (double **) malloc(sizeof(double *) * pwset->n_clean);
  for (int i = 0; i < pwset->n_clean ; i++) {
    this->data[i] = (double *) malloc(sizeof(double) * 2);
    this->data[i][0] = pwset->x_clean[i];
    this->data[i][1] = pwset->y_clean[i];
  }

}

/**
 * Destructor for the MixMod class.
 */
MixtureModelPWClusters::~MixtureModelPWClusters() {
  // Free the data array.
  if (data) {
    for (int i = 0; i < pwset->n_clean; i++) {
      free(data[i]);
    }
    free(data);
  }
  if (labels) {
    free(labels);
  }
  delete pwcl;
}

/**
 * Executes mixture model clustering.
 */
void MixtureModelPWClusters::run(char * criterion, int max_clusters) {
  int64_t lowest_cluster_num = 9;
  int found = 0;

  // Make sure we have the correct number of observations before performing
  // the comparision.
  if (this->pwset->n_clean < this->min_obs) {
    return;
  }

  // The data we are using is qualitative, so set the data type.
  XEM::GaussianData * gdata = new XEM::GaussianData(pwset->n_clean, 2, data);
  XEM::DataDescription dataDescription(gdata);

  // Set the number of clusters to be tested to 9.
  vector<int64_t> nbCluster;
  for (int i = 1; i <= max_clusters; i++) {
    nbCluster.push_back(i);
  }

  // Create the ClusteringInput object.
  XEM::ClusteringInput cInput(nbCluster, dataDescription);

  // Set the criterion to ICL
  cInput.removeCriterion(0);
  if (strcmp(criterion, "BIC") == 0) {
    cInput.addCriterion(XEM::BIC);
  }
  else if (strcmp(criterion, "ICL") == 0) {
    cInput.addCriterion(XEM::ICL);
  }
  else if (strcmp(criterion, "NEC") == 0) {
    cInput.addCriterion(XEM::NEC);
  }
  else if (strcmp(criterion, "CV") == 0) {
    cInput.addCriterion(XEM::CV);
  }
  else if (strcmp(criterion, "DCV") == 0) {
    cInput.addCriterion(XEM::DCV);
  }

  // Finalize input: run a series of sanity checks on it
  cInput.finalize();

  // Run XEM::ClusteringMain 5 times to find the iteration with the
  // lowest number of clusters.  We do this because the MixModLib can
  // detect different numbers of clusters depending on the random starting
  // point in the algorithm.
  for (int i = 0; i < 1; i++) {

    // Find clusters.
    XEM::ClusteringMain cMain(&cInput);
    cMain.run();

    // Create a new XEM::ClusteringOutput object
    XEM::ClusteringOutput * cOutput = cMain.getOutput();

    // Order the clusters using BIC.
    if (strcmp(criterion, "BIC") == 0) {
      cOutput->sort(XEM::BIC);
    }
    else if (strcmp(criterion, "ICL") == 0) {
      cOutput->sort(XEM::ICL);
    }
    else if (strcmp(criterion, "NEC") == 0) {
      cOutput->sort(XEM::NEC);
    }
    else if (strcmp(criterion, "CV") == 0) {
      cOutput->sort(XEM::CV);
    }
    else if (strcmp(criterion, "DCV") == 0) {
      cOutput->sort(XEM::DCV);
    }

    if (cOutput->atLeastOneEstimationNoError()) {
      found = 1;

      // Get the best XEM::ClusteringModelOutput
      XEM::ClusteringModelOutput* cMOutput = cOutput->getClusteringModelOutput().front();
      XEM::LabelDescription * ldescription = cMOutput->getLabelDescription();
      XEM::Label * label = ldescription->getLabel();
      int64_t * cLabels = label->getTabLabel();

      // Find the number of clusters
      int64_t num_clusters = 0;
      for (int i = 0; i < this->pwset->n_clean; i++) {
        if (cLabels[i] > num_clusters) {
          num_clusters = cLabels[i];
        }
      }

      // If this is the fewest number of clusters we've seen then copy
      // the labels over.
      if (num_clusters < lowest_cluster_num) {
        lowest_cluster_num = num_clusters;
        if (this->labels == NULL) {
          this->labels = (int64_t *) malloc(sizeof(int64_t) * this->pwset->n_clean);
        }
        for (int i = 0; i < this->pwset->n_clean; i++) {
          this->labels[i] = cLabels[i];
        }
      }
      delete cLabels;
    }
  }
  delete gdata;
  if (found) {
//    printf("Lowest: %d\n", lowest_cluster_num);
    // Create the cluster objects using the iteration with the fewest clusters.
    int cluster_num = 1;
    int cluster_samples[this->pwset->n_clean];
    bool done = false;
    while (!done) {
      done = true;
      for (int i = 0; i < this->pwset->n_clean; i++) {
        if (this->labels[i] == cluster_num) {
          done = false;
          cluster_samples[i] = 1;
        }
        else {
          cluster_samples[i] = 0;
        }
      }
      // if we found samples with the current cluster_num then create a
      // cluster and add it to the list
      if (!done) {
        PairWiseCluster * cluster = new PairWiseCluster(this->pwset);
        cluster->setClusterSamples(cluster_samples, true);
        // Calculate the Spearman correlation for this cluster.
        cluster->doSimilarity(this->method, this->min_obs);
//        cluster->printCluster();
        this->pwcl->addCluster(cluster);
      }
      cluster_num++;
    }
  }
}
