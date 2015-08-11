#ifndef _INDEX_
#define _INDEX_

#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>

#include <CLucene.h>
#include <CLucene/StdHeader.h>
#include <CLucene/util/CLStreams.h>
#include <CLucene/index/IndexWriter.h>

using namespace std;
using namespace lucene::index;
using namespace lucene::analysis;
using namespace lucene::util;
using namespace lucene::store;
using namespace lucene::document;

int regcomp(regex_t *, const char *, int);

class RunIndex {
  private:
    char * outdir;
    // The number of samples
    int nsamples;

    void indexFile(IndexWriter * writer, char * filepath);

  public:
    RunIndex(int argc, char *argv[]);
    ~RunIndex();
    void execute();
    static void printUsage();
};

#endif
