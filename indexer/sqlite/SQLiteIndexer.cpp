#include "../sqlite/SQLiteIndexer.h"

/**
 * The function to call when running the 'index' command.
 */
SQLiteIndexer::SQLiteIndexer(EMatrix * ematrix, char * indexdir) : Indexer(indexdir) {
  this->ematrix = ematrix;
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
     createDBTables(db, dbname);

     // Add the genes to the genes table.
     char q[2048];
     sqlite3_stmt *gene_insert_stmt;
     sqlite3_stmt *gene_select_stmt;
     int rc;

     // Prepare the SQL statements for selecting data from the SQLite database.
     sprintf(q, "%s", "SELECT * FROM genes WHERE gene_id = ?;");
     rc = sqlite3_prepare_v2(db, q, 2048, &gene_select_stmt, NULL);
     if (rc != SQLITE_OK) {
       fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }

     // Prepare the SQL statements for inserting data into the SQLite database.
     sprintf(q, "%s", "\
       INSERT INTO genes \
         (gene_id, name) VALUES (?, ?);");
     rc = sqlite3_prepare_v2(db, q, 2048, &gene_insert_stmt, NULL);
     if (rc != SQLITE_OK) {
       fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }

     // Iterate through the genes in the Expression matrix and add them to
     // the table.
     char ** genes = ematrix->getGenes();
     int num_genes = ematrix->getNumGenes();
     for (int i = 0; i < num_genes; i++) {
       // First, make sure the gene does not already exist. If it doesn't then
       // add it.
       sqlite3_bind_int(gene_select_stmt, 1, i + 1);
       if (sqlite3_step(gene_select_stmt) == SQLITE_ROW) {
         sqlite3_finalize(gene_select_stmt);
         continue;
       }
       char * gene = genes[i];
       sqlite3_bind_int(gene_insert_stmt, 1, i + 1);
       sqlite3_bind_text(gene_insert_stmt, 2, gene, strlen(gene), SQLITE_STATIC);
       printf("Inserting: %s\n", genes[i]);
       if(sqlite3_step(gene_insert_stmt) != SQLITE_ROW) {
         fprintf(stderr, "Failed to insert gene: %s\n", sqlite3_errmsg(db));
         exit(-1);
       }
       sqlite3_finalize(gene_insert_stmt);
       sqlite3_finalize(gene_select_stmt);
     }
     // Close the prepared SQL statement
     sqlite3_finalize(gene_insert_stmt);

     // Iterate through each of the files in the directory.
     struct dirent * entry;
     while ((entry = readdir(dir)) != NULL) {
       const char * filename = entry->d_name;

       // Skip the . and .. files.
       if (strcmp(filename, ".") == 0 || strcmp(filename, "..") == 0) {
         continue;
       }

       // The file must end in a suffix with 3 digits followed by .txt.
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
void SQLiteIndexer::IndexFile(sqlite3 *db, char * filepath, int nsamples) {
  // The 8 fields of the clusters tab file.
  int x, y, cluster_num, num_clusters, cluster_samples, num_missing;
  char * samples = (char *) malloc(sizeof(char) * nsamples);
  float cv;


  // Open the file for indexing.
  FILE * fp = fopen(filepath, "r");
  if (!fp) {
    fprintf(stderr, "Can not open file, %s. Cannot continue.\n", filepath);
    exit(-1);
  }

  while (!feof(fp)) {

    // Read in the fields for this line.
    int matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\n", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &cv, samples);
    // Skip lines that don't have the proper number of columns
    if (matches < 8) {
      continue;
    }



  }
  fclose(fp);

}

/**
 *
 */
void SQLiteIndexer::createDBTables(sqlite3 * db, char * dbname) {
  char q[2048];
  sqlite3_stmt *createStmt;
  int max_gene_len = ematrix->getMaxGeneLen();
  int max_sample_len = ematrix->getMaxSampleLen();

  //
  // Create the gene table.
  //
  sprintf(q, "\
    CREATE TABLE IF NOT EXISTS genes \
      (gene_id INTEGER PRIMARY KEY, \
       name CHAR(%d) NOT NULL);", max_gene_len);
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: %s.\n", "genes", dbname, sqlite3_errmsg(db));
    exit(-1);
  }

  sprintf(q, "%s", "CREATE UNIQUE INDEX IF NOT EXISTS genes_name ON genes(name);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "genes_name", "genes", dbname, sqlite3_errmsg(db));
    exit(-1);
  }

  //
  // Create the expression table
  //

  //
  // Create the sample table
  //
  sprintf(q, "\
    CREATE TABLE IF NOT EXISTS samples \
      (sample_id INTEGER PRIMARY KEY, \
       name CHAR(%d) NOT NULL);", max_sample_len);
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: '%s'.\n", "samples", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE UNIQUE INDEX IF NOT EXISTS sample_name ON samples(name);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "genes_name", "samples", dbname, sqlite3_errmsg(db));
  }

  //
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
    printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: %s.\n", "comps", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_g1 ON comps(gene1);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "comps_g1", "comps", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_g2 ON comps(gene2);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "comps_g2", "comps", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_num_clusters ON comps(num_clusters);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "comps_num_clusters", "comps", dbname, sqlite3_errmsg(db));
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
    printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: %s.\n", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_comp ON clusters(comp_id);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_comp", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_score ON clusters(score);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_score", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_missing ON clusters(num_missing);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_num_missing", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_samples ON clusters(num_samples);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_num_samples", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }

}
