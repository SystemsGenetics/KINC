#ifndef _LUCENEQUERY_
#define _LUCENEQUERY_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include <wchar.h>

#include "../IndexQuery.h"

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
using namespace lucene::search;
using namespace lucene::queryParser;

class CLuceneQuery : public IndexQuery {
  private:

  public:
    CLuceneQuery(char * indexdir);
    ~CLuceneQuery();

    void run(char * outfile, int x_coord, int y_coord);
};
#endif
