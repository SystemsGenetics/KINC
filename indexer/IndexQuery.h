#ifndef _INDEXQUERY_
#define _INDEXQUERY_

#include "../ematrix/EMatrix.h"

class IndexQuery {
  protected:
    char * indexdir;
    EMatrix * ematrix;

  public:
    IndexQuery(char * indexdir, EMatrix * ematrix);
    ~IndexQuery();

    void run();
};
#endif
