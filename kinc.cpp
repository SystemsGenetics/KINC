#include "kinc.h"

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
