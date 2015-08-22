#ifndef _QUERY_
#define _QUERY_


#include "clucene/CLuceneQuery.h"
#include "sqlite/SQLiteQuery.h"

class RunQuery {
  private:
    // The directory where indexes are housed.
    char * indexdir;
    // The full path to the result file.
    char * outfile;

    // The gene indices.
    int x_coord;
    int y_coord;

    double score;

  public:
    RunQuery(int argc, char *argv[]);
    ~RunQuery();

    void execute();
    static void printUsage();
};

#endif
