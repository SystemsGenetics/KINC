#ifndef _EMATRIX_
#define _EMATRIX_

#include <string.h>
#include <math.h>
#include <libgen.h>

#include "misc.h"

class EMatrix {
  private:
    // The expression matrix rows and columns.
    double ** data;
    // An array of gene names.
    char ** genes;
    // An array of sample names
    char ** samples;
    // The number of genes in the expression matrix.
    int num_genes;
    // The number of samples in the expression matrix;
    int num_samples;
    // The input filename without the prefix.
    char * file_prefix;

  public:

    // Constructor
    EMatrix(char * infilename, int rows, int cols, int headers,
        int omit_na = 0, char * na_val = NULL);
    // Destructor
    ~EMatrix();

    // GETTERS
    // Retrieves the expression matrix data array.
    double ** getMatrix() { return data; }
    // Retrieves a single row of the expression matrix.
    double * getRow(int i) { return data[i]; }
    // Retrieves the value of a single cell in the expression matrix.
    double getCell(int i, int j) { return data[i][j]; }
    // Retrieves the number of samples in the expression matrix.
    int getNumSamples() { return num_samples; }
    // Retrieves the number of genes in the expression matrix.
    int getNumGenes() { return num_genes; }
    // Retrieves the expression matrix file name without the extension
    char * getFilePrefix() { return file_prefix; }
    // Retrieves the genes array.
    char ** getGenes() { return genes; }
    // Retrieves the samples array.
    char ** getSamples() { return samples; }

    // Log transforms the values in the expression matrix.
    void logTransform();
    // Log2 transforms the values in the expression matrix.
    void log2Transform();
    // log10 transforms the values in the expression matrix.
    void log10Transform();
};

#endif
