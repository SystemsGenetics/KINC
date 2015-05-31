#include "kinc.h"
#include <mcheck.h>

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_usage() {
  printf("\n");
  printf("Usage: ./kinc [command]\n");
  printf("Available commands:\n");
  printf("  similarity       Performs pair-wise similarity calculations using an input expression matrix.\n");
  printf("  threshold        Identifies a threshold for cutting the similarity matrix\n");
  printf("  extract          Outputs the network edges file\n");
  printf("  help [command]   Prints these instructions. Include the command to print help\n");
  printf("                   for a specific command (e.g. kinc help similarity)\n");
  printf("\n");
}

/**
 * The main subroutine.  Parses the input parameters and executes the program
 * accordingly.
 */
int main(int argc, char *argv[]) {
  // Enable mtrace memory leak checking
  //mtrace();

  // The return value
  int retval = 0;

//  // MPI variables
//  int mpi_err, mpi_num_procs, mpi_id;
//
//  // Initialize MPI.
//  mpi_err = MPI_Init(&argc, &argv);
//
//  // Find out my process ID, and how many processes were started.
//  mpi_err |= MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
//  mpi_err |= MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);
//
//  if (mpi_err != 0) {
//    printf("MPI initialization failed\n");
//    exit(1);
//  }

  // For testing a single process... should comment out when not testing.
//  if (mpi_id + 1 != 5) {
//    mpi_err = MPI_Finalize();
//    return 1;
//  }

//  printf("Using %i out of %i processes.\n", mpi_id + 1, mpi_num_procs);

  // make sure we have at least one input argument for the command
  if (argc == 1) {
    printf("ERROR: Please provide the command to execute.\n\n");
    print_usage();
    retval = -1;
  }
//  else if (strcmp(argv[1], "dimreduce") == 0) {
////    retval =  do_dimreduce(argc, argv, mpi_id, mpi_num_procs);
//    retval =  do_dimreduce(argc, argv);
//  }
  // construct the similarity matrix
  else if (strcmp(argv[1], "similarity") == 0) {
//    retval =  do_similarity(argc, argv);
  }
  // identify the threshold for cutting the similarity matrix
  else if (strcmp(argv[1], "threshold") == 0) {
//    retval =  do_threshold(argc, argv);
  }
  // extract a given element from the matrix or a network
  else if (strcmp(argv[1], "extract") == 0) {
//    retval =  do_extract(argc, argv);
  }
  // print help documentation
  else if (strcmp(argv[1], "help") == 0) {
    if (argc == 3) {
      if (strcmp(argv[2], "similarity") == 0) {
        print_similarity_usage();
      }
      if (strcmp(argv[2], "threshold") == 0) {
        print_threshold_usage();
      }
      if (strcmp(argv[2], "extract") == 0) {
        print_extract_usage();
      }
    }
    else {
      print_usage();
    }
  }
  else {
    printf("ERROR: Unknown command.\n\n");
    print_usage();
    retval = -1;
  }

  // Wait until all other processes are completed before closing the manager 
//  MPI_Status stat;
//  char message[10];
//  if (mpi_id > 0 ) {
//    // All non master processes should report done when completed.
//    sprintf(message, "Done");
//    MPI_Send(message, strlen(message)+1, MPI_BYTE, 0,1,MPI_COMM_WORLD);
//  }
//  else {
//    // The master process should wait to get 'Done' from each process
//    // before terminating
//    int proc;
//    for (proc = 1; proc < mpi_num_procs; proc++) {
//      MPI_Recv(message, sizeof(message), MPI_BYTE, proc, 1, MPI_COMM_WORLD, &stat);
//    }
//  }
//
//  // Terminate MPI.
//  mpi_err = MPI_Finalize();

  return retval;
}

/**
 * The function for extracting results from the similarity matrix
 */
int do_extract(int argc, char *argv[]) {

  // Get the similarity matrix.
//  SimMatrixBinary * smatrix = new SimMatrixBinary(argc, argv);
  SimMatrixTabCluster * smatrix = new SimMatrixTabCluster(argc, argv);

  // If we have a threshold then we want to get the edges of the network.
  // Otherwise the user has asked to print out the similarty value for
  // two genes.
  if (smatrix->getThreshold() > 0) {
    smatrix->writeNetwork();
  }
  else {
    smatrix->getSimilarity();
  }

  return 1;
}

