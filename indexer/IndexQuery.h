#ifndef _INDEXQUERY_
#define _INDEXQUERY_

#include "../ematrix/EMatrix.h"

class IndexQuery {
  protected:
    char * indexdir;
    EMatrix * ematrix;

    // The index at which to stop searching.  All indexes in directories
    // less than this value will be skipped.
    double th;

  public:
    IndexQuery(char * indexdir, EMatrix * ematrix);
    ~IndexQuery();

    void run();
};
#endif
