#ifndef _QUERY_
#define _QUERY_


//#include "clucene/CLuceneQuery.h"
//#include "sqlite/SQLiteQuery.h"
#include "../ematrix/EMatrix.h"

class RunQuery {
  private:
    // The expression matrix object.
    EMatrix * ematrix;
    // Specifies the method: sc, pc, mi.
    char * method;

    // The directory where indexes are housed.
    char indexdir[1024];


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

    // The user-specified x coordinate to retrieve
    int x_coord;
    // The user-specified y coordinate to retrieve
    int y_coord;
    // The user-specified name of gene1
    char * gene1;
    // the user-specified name of gene2
    char * gene2;

    double score;

  public:
    RunQuery(int argc, char *argv[]);
    ~RunQuery();

    void execute();
    static void printUsage();
};

#endif
