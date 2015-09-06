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
#include "../../ematrix/EMatrix.h"

#include <sqlite3.h>

class SQLiteIndexer : public Indexer {
  private:
    EMatrix * ematrix;

    void createDBTables(sqlite3 * db, char * dbname);
    void insertGenes(sqlite3 *db);
    void insertSamples(sqlite3 *db);
    void insertCluster(sqlite3 *db, FILE *fp);

    void IndexFile(sqlite3 *db, char * filepath, int *mode_hist, int * num_comps);

  public:
    SQLiteIndexer(EMatrix * ematrix, char * indexdir);
    ~SQLiteIndexer();

    void run(int nsamples, int job_index);
};
#endif
