#include "PairWiseClustering.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2


// Set to 1 if Ctrl-C is pressed. Will allow for graceful termination
int Close_Now = 0;

// Define the function to be called when ctrl-c (SIGINT) signal is sent to process
void signal_callback_handler(int signum) {
   printf("Caught signal %d. Terminating program...\n",signum);
   // Cleanup and close up stuff here
   Close_Now = 1;

}
/**
 *
 */
////int do_dimreduce(int argc, char *argv[], int mpi_id, int mpi_num_procs) {
//void do_pairwise_Clustering(int argc, char *argv[], EMatrix *ematrix) {
//
//  MixtureModelPWClustering mixmodc = new MixtureModelPWClustering(argc, argv, ematrix);
//
//  mixmodc->run();
//}
/**
 * DRArgs constructor.
 */
PairWiseClustering::PairWiseClustering(EMatrix *ematrix, int min_obs, int num_jobs, int job_index) {
  this->min_obs = min_obs;
  this->num_jobs = num_jobs;
  this->job_index = job_index;
  this->ematrix = ematrix;
}

/**
 * DRArgs destructor.
 */
PairWiseClustering::~PairWiseClustering() {
  delete ematrix;
}
