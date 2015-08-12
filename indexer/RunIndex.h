#ifndef _INDEX_
#define _INDEX_


#include "clucene/CLuceneIndexer.h"


class RunIndex {
  private:
    char * outdir;
    // The number of samples
    int nsamples;

  public:
    RunIndex(int argc, char *argv[]);
    ~RunIndex();

    void execute();
    static void printUsage();
};

#endif
