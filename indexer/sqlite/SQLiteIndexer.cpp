#include "SQLiteIndexer.h"

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
  char* errorMessage;


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
       fprintf(stderr, "ERROR: Can't open SQLite database: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }

     printf("Preparing database file %s...\n", dbname);
     createDBTables(db, dbname);
     insertGenes(db);
     insertSamples(db);

     // Start a transaction to make inserts go faster.
     sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);

     // Clean up the summary and mode_hist tables
     char q[2048];
     char* errorMessage;
     sprintf(q, "%s", "DELETE FROM comps;");
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to clear comps table: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }
     sprintf(q, "%s", "DELETE FROM clusters;");
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to clear clusters table: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }
     sprintf(q, "%s", "DELETE FROM summary;");
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to clear summary table: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }
     sprintf(q, "%s", "DELETE FROM mode_hist;");
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to clear mode_hist table: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }

     // Initialize the histogram
     int max_modes = 15;
     int * mode_hist = (int *) malloc(sizeof(int) * max_modes);
     for (int i = 0; i < max_modes; i++) {
       mode_hist[i] = 0;
     }
     int total_comps;

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

       // Index the file.
       IndexFile(db, filepath, mode_hist, &total_comps);

     }

     // Now add in the stats:
     for (int i = 0; i < max_modes; i++) {
       sprintf(q, "INSERT INTO mode_hist (mode, num_comps) VALUES (%d, %d);", i, mode_hist[i]);
       if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
         fprintf(stderr, "Failed to insert mode_hist record: %s\n", sqlite3_errmsg(db));
         exit(-1);
       }
     }
     sprintf(q, "INSERT INTO summary (num_comps) VALUES (%d);", total_comps);
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to insert mode_hist record: %s\n", sqlite3_errmsg(db));
       exit(-1);
     }
     free(mode_hist);

     // Terminate the transaction
     sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);

     // Close the database.
     sqlite3_close(db);
   }
}
/**
 *
 */
void SQLiteIndexer::IndexFile(sqlite3 *db, char * filepath, int *mode_hist, int * total_comps) {
  char q[24800];
  sqlite3_stmt *comp_select_stmt;
  sqlite3_stmt *comp_insert_stmt;
  sqlite3_stmt *cluster_select_stmt;
  sqlite3_stmt *cluster_insert_stmt;
  int rc;

  // Open the file for indexing.
  FILE * fp = fopen(filepath, "r");
  if (!fp) {
    fprintf(stderr, "Can't open file, %s. Cannot continue.\n", filepath);
    exit(-1);
  }

  // Prepare the SQL statements for selecting a pair-wise comp.
  sprintf(q, "%s", "SELECT * FROM comps WHERE gene1 = ? and gene2 = ?;");
  rc = sqlite3_prepare_v2(db, q, 2048, &comp_select_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to prepare statement: %s\n", sqlite3_errmsg(db));
    exit(-1);
  }
  // Prepare the SQL statements for inserting a pair-wise comp.
  sprintf(q, "%s", "\
    INSERT INTO comps \
      (gene1, gene2, num_clusters) VALUES (?, ?, ?);");
  rc = sqlite3_prepare_v2(db, q, 2048, &comp_insert_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
    exit(-1);
  }
  // Prepare the SQL statements for selecting a cluster.
  sprintf(q, "%s", "SELECT * FROM clusters WHERE gene1 = ? and gene2 = ? and cluster_num = ?;");
  rc = sqlite3_prepare_v2(db, q, 2048, &cluster_select_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to prepare statement: %s\n", sqlite3_errmsg(db));
    exit(-1);
  }
  // Prepare the SQL statements for inserting a cluster.
  sprintf(q, "%s", "\
    INSERT INTO clusters \
      (gene1, gene2, num_clusters, cluster_num, num_missing, num_samples, score, samples) \
      VALUES (?, ?, ?, ?, ?, ?, ?, ?);");
  rc = sqlite3_prepare_v2(db, q, 2048, &cluster_insert_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
    exit(-1);
  }

  int nsamples = ematrix->getNumSamples();
  int j, k, cluster_num, num_clusters, cluster_num_samples, num_missing;
  char samples[nsamples];
  float cv;
  while (!feof(fp)) {

    // Read in the fields for this line.
    int matches = fscanf(fp, "%d\t%d\%d\t%d\%d\t%d\t%f\t%s\n", &j, &k, &cluster_num, &num_clusters, &cluster_num_samples, &num_missing, &cv, (char *)&samples);

    // Skip lines that don't have the proper number of columns
    if (matches < 8) {
      char tmp[nsamples*2];
      matches = fscanf(fp, "%s\n", (char *)&tmp);
      continue;
    }

    // Update the histogram counts
    mode_hist[num_clusters]++;

    // See if this comparision already exists.  If it doesn't then add it.
    sqlite3_bind_int(comp_select_stmt, 1, j);
    sqlite3_bind_int(comp_select_stmt, 2, k);
    if (sqlite3_step(comp_select_stmt) != SQLITE_ROW) {
      sqlite3_bind_int(comp_insert_stmt, 1, j);
      sqlite3_bind_int(comp_insert_stmt, 2, k);
      sqlite3_bind_int(comp_insert_stmt, 3, num_clusters);
      if(sqlite3_step(comp_insert_stmt) != SQLITE_DONE) {
        fprintf(stderr, "Failed to insert comparision: %s\n", sqlite3_errmsg(db));
        exit(-1);
      }
    }
    sqlite3_reset(comp_select_stmt);
    sqlite3_reset(comp_insert_stmt);

    // See if this cluster already exists. If it doesn't then add it.
    sqlite3_bind_int(cluster_select_stmt, 1, j);
    sqlite3_bind_int(cluster_select_stmt, 2, k);
    sqlite3_bind_int(cluster_select_stmt, 3, cluster_num);
    if (sqlite3_step(cluster_select_stmt) != SQLITE_ROW) {
      int sample_len = strlen(samples);
      sqlite3_bind_int(cluster_insert_stmt, 1, j);
      sqlite3_bind_int(cluster_insert_stmt, 2, k);
      sqlite3_bind_int(cluster_insert_stmt, 3, num_clusters);
      sqlite3_bind_int(cluster_insert_stmt, 4, cluster_num);
      sqlite3_bind_int(cluster_insert_stmt, 5, num_missing);
      sqlite3_bind_int(cluster_insert_stmt, 6, cluster_num_samples);
      sqlite3_bind_double(cluster_insert_stmt, 7, cv);
      sqlite3_bind_text(cluster_insert_stmt, 8, samples, sample_len, SQLITE_STATIC);
      if(sqlite3_step(cluster_insert_stmt) != SQLITE_DONE) {
        fprintf(stderr, "Failed to insert cluster: %s\n", sqlite3_errmsg(db));
        exit(-1);
      }
    }
    sqlite3_reset(cluster_select_stmt);
    sqlite3_reset(cluster_insert_stmt);


  }

}


/**
 *
 */
void SQLiteIndexer::insertGenes(sqlite3 *db) {
  // Add the genes to the genes table.
  char* errorMessage;
  char q[2048];
  sqlite3_stmt *gene_insert_stmt;
  sqlite3_stmt *gene_select_stmt;
  int rc;

  // Start a transaction to make inserts go faster.
  sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);


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
      sqlite3_reset(gene_select_stmt);
      continue;
    }
    char * gene = genes[i];
    sqlite3_bind_int(gene_insert_stmt, 1, i + 1);
    int gene_len = strlen(gene);
    sqlite3_bind_text(gene_insert_stmt, 2, gene, gene_len, SQLITE_STATIC);
    if(sqlite3_step(gene_insert_stmt) != SQLITE_DONE) {
      fprintf(stderr, "Failed to insert gene: %s\n", sqlite3_errmsg(db));
      exit(-1);
    }
    sqlite3_reset(gene_insert_stmt);
    sqlite3_reset(gene_select_stmt);
  }
  // Close the prepared SQL statement
  sqlite3_finalize(gene_select_stmt);
  sqlite3_finalize(gene_insert_stmt);

  // Terminate the transaction
  sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
}
/**
 *
 */
void SQLiteIndexer::insertSamples(sqlite3 *db) {
  // Add the genes to the genes table.
  char* errorMessage;
  char q[2048];
  sqlite3_stmt *sample_insert_stmt;
  sqlite3_stmt *sample_select_stmt;
  int rc;

  // Start a transaction to make inserts go faster.
  sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);

  // Prepare the SQL statements for selecting data from the SQLite database.
  sprintf(q, "%s", "SELECT * FROM samples WHERE sample_id = ?;");
  rc = sqlite3_prepare_v2(db, q, 2048, &sample_select_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
    exit(-1);
  }

  // Prepare the SQL statements for inserting data into the SQLite database.
  sprintf(q, "%s", "\
    INSERT INTO samples \
      (sample_id, name) VALUES (?, ?);");
  rc = sqlite3_prepare_v2(db, q, 2048, &sample_insert_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
    exit(-1);
  }

  // Iterate through the genes in the Expression matrix and add them to
  // the table.
  char ** samples = ematrix->getSamples();
  int num_samples = ematrix->getNumSamples();
  for (int i = 0; i < num_samples; i++) {
    // First, make sure the gene does not already exist. If it doesn't then
    // add it.
    sqlite3_bind_int(sample_select_stmt, 1, i + 1);
    if (sqlite3_step(sample_select_stmt) == SQLITE_ROW) {
      sqlite3_reset(sample_select_stmt);
      continue;
    }
    char * sample = samples[i];
    sqlite3_bind_int(sample_insert_stmt, 1, i + 1);
    sqlite3_bind_text(sample_insert_stmt, 2, sample, strlen(sample), SQLITE_STATIC);
    if(sqlite3_step(sample_insert_stmt) != SQLITE_DONE) {
      fprintf(stderr, "Failed to insert gene: %s\n", sqlite3_errmsg(db));
      exit(-1);
    }
    sqlite3_reset(sample_insert_stmt);
    sqlite3_reset(sample_select_stmt);
  }
  // Close the prepared SQL statement
  sqlite3_finalize(sample_select_stmt);
  sqlite3_finalize(sample_insert_stmt);

  // Terminate the transaction
  sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
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
     gene1 INTEGER NOT NULL, \
     gene2 INTEGER NOT NULL, \
     num_clusters INTEGER NOT NULL DEFAULT 0); ");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: %s.\n", "comps", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_gene1 ON comps(gene1);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "comps_gene1", "comps", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_gene2 ON comps(gene2);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "comps_gene1", "comps", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_g12 ON comps(gene1, gene2);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "comps_g12", "comps", dbname, sqlite3_errmsg(db));
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
     gene1 INTEGER NOT NULL, \
     gene2 INTEGER NOT NULL, \
     num_clusters INTEGER NOT NULL DEFAULT 0, \
     cluster_num INTEGER NOT NULL, \
     num_missing INTEGER NOT NULL, \
     num_samples INTEGER NOT NULL, \
     score REAL, \
     samples TEXT NOT NULL); ");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: %s.\n", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_gene1 ON clusters(gene1);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_gene1", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_gene2 ON clusters(gene2);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_gene1", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_g12 ON clusters(gene1, gene2);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_g12", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_g12num ON clusters(gene1, gene2, cluster_num);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_g12num", "clusters", dbname, sqlite3_errmsg(db));
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_clusters ON clusters(num_clusters);");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.  ERROR: %s.\n", "clusters_num_clusters", "clusters", dbname, sqlite3_errmsg(db));
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


  //
  // Create the mode histogram table
  //
  sprintf(q, "%s",  "\
    CREATE TABLE IF NOT EXISTS mode_hist \
    (mode INTEGER NOT NULL, \
     num_comps INTEGER NOT NULL); ");
  sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
  if (sqlite3_step(createStmt) != SQLITE_DONE) {
    printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: %s.\n", "mode_hist", dbname, sqlite3_errmsg(db));
    exit(-1);
  }

  //
  // Create the summary table
  //
  sprintf(q, "%s",  "\
      CREATE TABLE IF NOT EXISTS summary \
      (num_comps INTEGER NOT NULL); ");
    sqlite3_prepare_v2(db, q, 2048, &createStmt, NULL);
    if (sqlite3_step(createStmt) != SQLITE_DONE) {
      printf("ERROR: failed to create '%s' table for index '%s'.  ERROR: %s.\n", "summary", dbname, sqlite3_errmsg(db));
      exit(-1);
    }

}
