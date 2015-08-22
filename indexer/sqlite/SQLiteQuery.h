#ifndef _SQLITEQUERY_
#define _SQLITEQUERY_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include <sqlite3.h>

#include "../IndexQuery.h"


class SQLiteQuery : public IndexQuery {
  private:

  public:
    SQLiteQuery(char * indexdir);
    ~SQLiteQuery();

    void run(char * outfile, int x_coord, int y_coord, float score);
};
#endif
