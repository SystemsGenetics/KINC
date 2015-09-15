#ifndef _INDEX_
#define _INDEX_

#include <sys/stat.h>

//#include "clucene/CLuceneIndexer.h"
#include "sqlite/SQLiteIndexer.h"
#include "../ematrix/EMatrix.h"


class RunIndex {
  protected:
    // The expression matrix object.
    EMatrix * ematrix;
    // Specifies the method: sc, pc, mi.
    char * method;

    // The directory where indexes are housed.
    char indexdir[1024];

    // Variables for indexing
    // ----------------------
    // The number of samples
    int nsamples;
    // The index of this job within the total jobs.  Must be
    // between 0 and 101.
    int job_index;
    // The index of this job where indexing will begin. Must be between 0 and
    // 101.
    int job_start;

    // Variables for the expression matrix
    // -----------------------------------
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
    char func[10];

  public:
    RunIndex(int argc, char *argv[]);
    ~RunIndex();

    void execute();
    static void printUsage();
};

#endif
