#include "kinc.h"

/**
 * The main subroutine.  Parses the input parameters and executes the program
 * accordingly.
 */
int main(int argc, char *argv[]) {

  // The return value
  int retval = 0;

  // MPI variables
  int mpi_err, mpi_num_procs, mpi_id;

  // Initialize MPI.
  mpi_err = MPI_Init(&argc, &argv);

  // Find out my process ID, and how many processes were started.
  mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
  mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_procs);

  printf("Using %i out of %i processes.\n", mpi_id + 1, mpi_num_procs);

  // make sure we have at least one input argument for the command
  if (argc == 1) {
    printf("ERROR: Please provide the command to execute.\n\n");
    print_usage();
    retval = -1;
  }
  else if (strcmp(argv[1], "dimreduce") == 0) {
    retval =  do_dimreduce(argc, argv, mpi_id, mpi_num_procs);
  }
  // construct the similarity matrix
  else if (strcmp(argv[1], "similarity") == 0) {
    retval =  do_similarity(argc, argv);
  }
  // identify the threshold for cutting the similarity matrix
  else if (strcmp(argv[1], "threshold") == 0) {
    retval =  do_threshold(argc, argv);
  }
  // extract a given element from the matrix or a network
  else if (strcmp(argv[1], "extract") == 0) {
    retval =  do_extract(argc, argv);
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

  // Terminate MPI.
  mpi_err = MPI_Finalize();

  return retval;
}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_usage() {
  printf("\n");
  printf("Usage: ./kinc [command]\n");
  printf("Available commands:\n");
  printf("  similarity       Constructs the similarity matrix from an input expression matrix.\n");
  printf("  threshold        Identifies a threshold using RMT for the similarity matrix\n");
  printf("  extract          Outputs the network file\n");
  printf("  help [command]   Prints these instructions. Include the command to print help\n");
  printf("                   for a specific command (e.g. kinc help similarity)\n");
  printf("\n");
}
