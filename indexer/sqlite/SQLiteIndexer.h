#ifndef _SQLLITEINDEXER_
#define _SQLLITEINDEXER_

#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>

#include "../Indexer.h"

#include <sqlite3.h>

class SQLiteIndexer : public Indexer {
  private:
    void createTables(sqlite3 * db);
  public:
    SQLiteIndexer(char * indexdir);
    ~SQLiteIndexer();

    void run(int nsamples);
};
#endif
