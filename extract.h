#ifndef _EXTRACT_
#define _EXTRACT_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <limits.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include "ematrix.h"

// Executes the 'extract' program of KINC.
int do_extract(int argc, char *argv[]);

// Prints the usage instructions for the 'extract' program.
void print_extract_usage();

/**
 * Similarity Matrix base class.
 *
 * Each file format (e.g binary or cluster file) must implmenet its own
 * child class for extracting results from its files.
 */
class SimilarityMatrix {

  protected:
    // The expression matrix object.
    EMatrix * ematrix;
    // Indicates if headers are present in the input EMatrix file.
    int headers;
    // The number of rows in the expression matrix.
    int rows;
    // The number of columns in the expression matrix.
    int cols;
    // Set to 1 if nothing but the sim value is shown
    int quiet;
    // The input expression matrix file name.
    char * infilename;
    // The directory where the expression matrix is found
    char * input_dir;
    // Specifies the method: cor, mi.
    char method[10];
    // The user-specified x coordinate to retrieve
    int x_coord;
    // The user-specified y coordinate to retrieve
    int y_coord;
    // The user-specified name of gene1
    char * gene1;
    // the user-specified name of gene2
    char * gene2;
    // The threshold for creating the network.
    float th;
    // Filters
    int max_missing;
    int min_cluster_size;
  public:
    // Constructor.
    SimilarityMatrix(int argc, char *argv[]);
    // Desctructor.
    ~SimilarityMatrix();

    // GETTERS
    // -------
    // Retrieves the similarity threshold specified by the user.
    float getThreshold() { return th; }

    // TO BE IMPLEMENTED BY THE CHILD CLASS
    // ------------------------------------
    // Retrieves the set of edges that match the given filtering parameters.
    // The user must have provided a threshold value.
    void writeNetwork() {};
    // Retrieves the similarity value for the given filtering paramters.
    // The user must have provided an x and y coordiante.
    void getSimilarity() {};
};

/**
 * Class for extracting a network from the original RMTGeneNet binary files.
 */
class SimMatrixBinary : public SimilarityMatrix {

  private:
    // Holds file handles for all of the binary files.
    FILE * files[50];
    // The number of lines per bin file.
    int num_lines[50];
    // Indicates the number of files in the files array.
    int num_files;

    // Opens file handles to all of the binary files.
    void openBinFiles();
    // Closes all of the file handles for the binary files.
    void closeBinFiles();

  public:
    // Constructur.
    SimMatrixBinary(int argc, char *argv[]);
    // Destructor.
    ~SimMatrixBinary();
    // Retrieves the set of edges that match the given filtering parameters.
    // The user must have provided a threshold value.
    void writeNetwork();
    // Retrieves the similarity value for the given filtering paramters.
    // The user must have provided an x and y coordiante.
    void getPosition();
};

/**
 * Class for extract a network from a tab-delimited cluster file.
 */
class SimMatrixTabCluster : public SimilarityMatrix {

  private:
    // The number of jobs that were used to construct the files.
    int num_jobs;


    // Discovers the number of jobs used to generate the clustering files.
    void getNumJobs();
  public:
    // Constructur.
    SimMatrixTabCluster(int argc, char *argv[]);
    // Destructor.
    ~SimMatrixTabCluster();
    // Retrieves the set of edges that match the given filtering parameters.
    // The user must have provided a threshold value.
    void writeNetwork();
    // Retrieves the similarity value for the given filtering paramters.
    // The user must have provided an x and y coordiante.
    void getPosition();
};

#endif
