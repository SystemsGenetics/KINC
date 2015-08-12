#include "../sqlite/SQLiteIndexer.h"

/**
 * The function to call when running the 'index' command.
 */
SQLiteIndexer::SQLiteIndexer(char * indexdir) : Indexer(indexdir) {

}
/**
 * Destructor
 */
SQLiteIndexer::~SQLiteIndexer() {

}

/**
 * Performs the indexing.
 */
void SQLiteIndexer::run(int nsamples) {
  int i;

   // Iterate through the directories.
   for (i = 100; i >= 0; i--) {
     char dirname[1024];
     char dbname[1024];
     if (i > 0) {
       sprintf(dirname, "./%s/%03d", indexdir, i);
       sprintf(dbname, "%s/%03d.db", dirname, i);
     }
     else {
       sprintf(dirname, "./%s/nan", indexdir);
       sprintf(dbname, "%s/nan.db", dirname);
     }

     // Open the output directory.
     DIR * dir;
     dir = opendir(dirname);
     if (!dir) {
       fprintf(stderr, "WARNING: The output directory, %s, is missing. skipping.\n", dirname);
       continue;
     }

     // Each directory has it's own SQLite database.
     sqlite3 * db;
     if (sqlite3_open(dbname, &db)){
       fprintf(stderr, "ERROR: Can't open sqlite database: %s\n", sqlite3_errmsg(db));
       sqlite3_close(db);
       exit(-1);
     }
     createTables(db);


     // Iterate through each of the files in the directory.
     struct dirent * entry;
     while ((entry = readdir(dir)) != NULL) {
       const char * filename = entry->d_name;

       // Skip the . and .. files.
       if (strcmp(filename, ".") == 0 || strcmp(filename, "..") == 0) {
         continue;
       }

       // The file must end in a a suffix with 3 digits followed by .txt.
       // We use a regular expression to see if this is true. If not, then
       // skip this file.
       regex_t preg;
       char pattern[64] = "[0-9][0-9][0-9].txt";
       int  rc;
       size_t     nmatch = 2;
       regmatch_t pmatch[2];
       if (0 != (rc = regcomp(&preg, pattern, 0))) {
           printf("regcomp() failed, returning nonzero (%d)\n", rc);
           exit(EXIT_FAILURE);
       }
       if (0 != (rc = regexec(&preg, filename, nmatch, pmatch, 0))) {
        continue;
       }

       // Construct the full path to the file.
       char filepath[1024];
       sprintf(filepath, "%s/%s", dirname, filename);
       printf("Indexing file %s...\n", filename);
     }
     sqlite3_close(db);
   }
}

/**
 *
 */
void SQLiteIndexer::createTables(sqlite3 * db) {
  char q[2048];
  sqlite3_stmt *createStmt;

  //
  // Create the gene table.
  //
  sprintf(q, "%s", "\
    CREATE TABLE IF NOT EXISTS genes \
      (gene_id INTEGER PRIMARY KEY, \
       name CHAR(128) NOT NULL);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create %s table: %s.\n", "genes", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE UNIQUE INDEX IF NOT EXISTS genes_name ON genes(name);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "genes_name", "genes", sqlite3_errmsg(db));
    exit(-1);
  }
  //
  // Create the comps table.
  //
  sprintf(q, "%s", "\
  CREATE TABLE IF NOT EXISTS comps \
    (comp_id INTEGER PRIMARY KEY, \
     gene1 INTEGER NOT NULL REFERENCES genes(gene_id), \
     gene2 INTEGER NOT NULL REFERENCES genes(gene_id), \
     num_clusters INTEGER NOT NULL DEFAULT 0); ");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create %s table: %s.\n", "comps", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_g1 ON comps(gene1);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "comps_g1", "comps", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_g2 ON comps(gene2);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "comps_g2", "comps", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_num_clusters ON comps(num_clusters);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "comps_num_clusters", "comps", sqlite3_errmsg(db));
    exit(-1);
  }

  //
  // Create the clusters table
  //
  sprintf(q, "%s",  "\
    CREATE TABLE IF NOT EXISTS clusters \
    (cluster_id INTEGER PRIMARY KEY, \
     cluster_num INTEGER NOT NULL, \
     comp_id INTEGER NOT NULL REFERENCES comps(comp_id), \
     num_missing INTEGER NOT NULL, \
     num_samples INTEGER NOT NULL, \
     score REAL, \
     samples TEXT NOT NULL); ");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create %s table: %s.\n", "clusters", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_comp ON clusters(comp_id);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "clusters_comp", "clusters", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_score ON clusters(score);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "clusters_score", "clusters", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_missing ON clusters(num_missing);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "clusters_num_missing", "clusters", sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_samples ON clusters(num_samples);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index %s on %s table: %s.\n", "clusters_num_samples", "clusters", sqlite3_errmsg(db));
    exit(-1);
  }

}
