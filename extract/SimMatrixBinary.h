#ifndef _EXTRACT_
#define _EXTRACT_

#include "SimilarityMatrix.h"

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
    // The directory where the expression matrix is found
    char * bin_dir;


    // Opens file handles to all of the binary files.
    void openBinFiles();
    // Closes all of the file handles for the binary files.
    void closeBinFiles();

  public:
    // Constructur.
    SimMatrixBinary(EMatrix *ematrix, int quiet, char ** method, int num_methods, char * th_method,
        int x_coord, int y_cood, char * gene1, char * gene2, float th);
    // Destructor.
    ~SimMatrixBinary();
    // Retrieves the set of edges that match the given filtering parameters.
    // The user must have provided a threshold value.
    void writeNetwork();
    // Retrieves the similarity value for the given filtering paramters.
    // The user must have provided an x and y coordiante.
    void getPosition();
};


#endif
