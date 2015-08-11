#ifndef _QUERY_
#define _QUERY_

#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <wchar.h>


#include <CLucene.h>
#include <CLucene/StdHeader.h>
#include <CLucene/util/CLStreams.h>
#include <CLucene/index/IndexWriter.h>

using namespace std;
using namespace lucene::analysis;
using namespace lucene::index;
using namespace lucene::util;
using namespace lucene::queryParser;
using namespace lucene::document;
using namespace lucene::search;

class RunQuery {
  private:
    // The directory where indexes are housed.
    char * indexes;
    // The full path to the result file.
    char * outfile;

    // The gene indicies.
    char * x_coord;
    char * y_coord;

    void indexFile(IndexWriter * writer, char * filepath);

  public:
    RunQuery(int argc, char *argv[]);
    ~RunQuery();
    void execute();
    static void printUsage();
};

#endif
