#ifndef _LUCENEINDEXER_
#define _LUCENEINDEXER_

#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <wchar.h>

#include "../Indexer.h"

#include <CLucene.h>
#include <CLucene/StdHeader.h>
#include <CLucene/util/CLStreams.h>
#include <CLucene/index/IndexWriter.h>

int regcomp(regex_t *, const char *, int);

using namespace std;
using namespace lucene::index;
using namespace lucene::analysis;
using namespace lucene::util;
using namespace lucene::store;
using namespace lucene::document;
using namespace lucene::search;

class CLuceneIndexer : public Indexer {
  private:
    void indexFile(IndexWriter * writer, char * filepath, int nsamples);

  public:
    CLuceneIndexer(char * indexdir);
    ~CLuceneIndexer();

    void run(int nsamples);
};
#endif
