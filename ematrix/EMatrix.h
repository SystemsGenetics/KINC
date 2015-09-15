#ifndef _EMATRIX_
#define _EMATRIX_

#include <string.h>
#include <math.h>
#include <libgen.h>
#include <getopt.h>
#include <unistd.h>


#include "../general/misc.h"

class EMatrix {
  private:
    // The maximum length of the sample and gene name strings.
    int max_sample_len;
    int max_gene_len;

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
    // Indicates if the expression matrix has headers.
    int headers;
    // The input file name
    char *infilename;
    // The number of rows in the input ematrix file (including the header)
    int rows;
    // The number of cols in the input ematrix file.
    int cols;
    // Indicates if missing values should be ignored in the EMatrix file.
    int omit_na;
    // Specifies the value that represents a missing value.
    char *na_val;
    // Specifies the transformation function: log2, none.
    char *func;
    // Set to 1 to perform log10 transformation.
    int do_log10;
    // Set to 1 to perform log2 transformation.
    int do_log2;
    // Set to 1 to perform log transformation.
    int do_log;

  public:


    // Constructor
    EMatrix(char * infilename, int rows, int cols, int headers, int omit_na, char *na_val, char * func);
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
    // Retrieves the input file name
    char * getInfileName() { return infilename; }
    // Indicates if missing values are omitted.
    int isMissingOmitted() { return omit_na; }

    // Return the max length of the genes and samples
    int getMaxGeneLen() { return max_gene_len; }
    int getMaxSampleLen() { return max_sample_len; }

    int getGeneCoord(char * gene);
    char * getGene(int index);

    char * getUsage();

    // Log transforms the values in the expression matrix.
    void logTransform();
    // Log2 transforms the values in the expression matrix.
    void log2Transform();
    // log10 transforms the values in the expression matrix.
    void log10Transform();
};

#endif
