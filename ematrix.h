#ifndef _EMATRIX_
#define _EMATRIX_

#include <string.h>
#include <math.h>

#include "misc.h"

class EMatrix {
  private:
    double ** data;
    char ** genes;
    char ** samples;
    int num_genes;
    int num_samples;
  public:
    EMatrix(char * infilename, int rows, int cols, int headers, int omit_na, char * na_val);
    ~EMatrix();

    double ** getMatrix() { return data; }
    double * getRow(int i) { return data[i]; }
    double getCell(int i, int j) { return data[i][j]; }
    int getNumSamples() { return num_samples; }
    int getNumGenes() { return num_genes; }

    void logTransform();
    void log2Transform();
    void log10Transform();
};

#endif
