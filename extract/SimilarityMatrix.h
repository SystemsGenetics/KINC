#ifndef _SIMMATRIX_
#define _SIMMATRIX_

#include <getopt.h>
#include <dirent.h>
#include "../ematrix/EMatrix.h"

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
    // Set to 1 if nothing but the sim value is shown
    int quiet;
    // Specifies the method: cor, mi.
    char * method;
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
  public:
    // Constructor.
    SimilarityMatrix(EMatrix *ematrix, int quiet, char * method, int x_coord,
        int y_cood, char * gene1, char * gene2, float th);
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
    void getPosition() {};
};
#endif
