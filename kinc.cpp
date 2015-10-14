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
  mtrace();

  // The return value
  int retval = 0;

  // make sure we have at least one input argument for the command
  if (argc == 1) {
    printf("ERROR: Please provide the command to execute.\n\n");
    print_usage();
    retval = -1;
  }
  // construct the similarity matrix
  else if (strcmp(argv[1], "similarity") == 0) {
    RunSimilarity * similarity = new RunSimilarity(argc, argv);
    similarity->execute();
    delete similarity;
  }
  // construct the similarity matrix
  else if (strcmp(argv[1], "index") == 0) {
    RunIndex * index = new RunIndex(argc, argv);
    index->execute();
    delete index;
  }
  // construct the similarity matrix
  else if (strcmp(argv[1], "query") == 0) {
    RunQuery * query = new RunQuery(argc, argv);
    query->execute();
    delete query;
  }
  // identify the threshold for cutting the similarity matrix
  else if (strcmp(argv[1], "threshold") == 0) {
    RunThreshold * threshold = new RunThreshold(argc, argv);
    threshold->execute();
    delete threshold;
  }
  // extract a given element from the matrix or a network
  else if (strcmp(argv[1], "extract") == 0) {
    RunExtract * extract = new RunExtract(argc, argv);
    extract->execute();
    delete extract;
  }
  // print help documentation
  else if (strcmp(argv[1], "help") == 0) {
    if (argc == 3) {
      if (strcmp(argv[2], "similarity") == 0) {
        RunSimilarity::printUsage();
      }
      if (strcmp(argv[2], "threshold") == 0) {
        RunThreshold::printUsage();
      }
      if (strcmp(argv[2], "extract") == 0) {
        RunExtract::printUsage();
      }
      if (strcmp(argv[2], "index") == 0) {
        RunIndex::printUsage();
      }
      if (strcmp(argv[2], "query") == 0) {
        RunQuery::printUsage();
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

  return retval;
}

/**
 *
 */
void start_mpi() {

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
}

/**
 *
 */
void end_mpi() {
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
}


