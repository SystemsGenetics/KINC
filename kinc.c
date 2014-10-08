#include "kinc.h"

/**
 * The main subroutine.  Parses the input parameters and executes the program
 * accordingly.
 */
int main(int argc, char *argv[]) {
  // make sure we have at least one input argument for the command
  if (argc == 1) {
    printf("ERROR: Please provide the command to execute.\n\n");
    print_usage();
    exit(-1);
  }
  if (strcmp(argv[1], "dimreduce") == 0) {
    return do_dimreduce(argc, argv);
  }
  // construct the similarity matrix
  if (strcmp(argv[1], "similarity") == 0) {
    return do_similarity(argc, argv);
  }
  // identify the threshold for cutting the similarity matrix
  else if (strcmp(argv[1], "threshold") == 0) {
    return do_threshold(argc, argv);
  }
  // extract a given element from the matrix or a network
  else if (strcmp(argv[1], "extract") == 0) {
    return do_extract(argc, argv);
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
    exit(-1);
  }
  return 0;
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
