#include "MixtureModelClustering.h"

MixtureModelClustering::MixtureModelClustering(EMatrix *ematrix, int min_obs,
    int num_jobs, int job_index, char *method, char * criterion,
    int max_clusters)
  : PairWiseClustering(ematrix, min_obs, num_jobs, job_index) {

  // Initialize some values.
  this->max_clusters = max_clusters;
  this->criterion = criterion;
  this->method = method;

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
}

/**
 *
 */
MixtureModelClustering::~MixtureModelClustering() {

}


/**
 *
 */
void MixtureModelClustering::run() {
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
  long long int comp_start = (job_index - 1) * comps_per_process;
  long long int comp_stop = (job_index - 1) * comps_per_process + comps_per_process;

  // If this is the last process and there are some remainder comparisons
  // then we need to add them to the stop
  if (job_index == num_jobs) {
    comp_stop = total_comps;
  }
  printf("  Job %d of %d. \n", job_index, num_jobs);
  printf("  Performing comparisions %lld to %lld (%lld) of %lld\n", comp_start, comp_stop, comp_stop - comp_start, total_comps);
  fflush(stdout);

  // Create the writer object to write out the cluters.
  PairWiseClusterWriter * pwcw = new PairWiseClusterWriter(method,
      ematrix->getFilePrefix(), job_index, ematrix->getNumSamples());

  // Perform the pair-wise clustering.
  int n_comps = 0;
  int my_comps = 0;

//  for (int i = 122 - 1; i < num_rows; i++) {
  for (int i = pwcw->getRecoveryX() - 1; i < num_rows; i++) {
    /*if (i == 50) {
      break;
    }*/
    for (int j = pwcw->getRecoveryY() - 1; j < num_rows; j++) {
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
      PairWiseClusterList * cluster_list = mixmod->getClusterList();
      pwcw->writeClusters(cluster_list, i, j);

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
        //if (!isnan(seconds_passed)) {
          printf("%d. Complete: %.4f%%. Mem: %ldb. Remaining: %.2fh; %.2fd; %.2fy. Coords: %d, %d.        \r",
            job_index, (float) percent_complete, memory->size, (float) hours_left, (float) days_left, (float) years_left, i, j);
        //}
        free(memory);
      }
      n_comps++;
      my_comps++;
      delete mixmod;
      delete pwset;

//      exit(1);
    }
  }
  delete pwcw;
}
