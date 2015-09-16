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
void SQLiteIndexer::run(int nsamples, int job_index, int job_start) {
  int i;
  char* errorMessage;

  // Convert a job_index of 101 to be -1.  This is because we want to
  // process the index files in descending order and we want the 'nan'
  // directory to be last.
  if (job_index == 101) {
    job_index = -1;
  }
  if (job_start == -2) {
    job_start = 100;
  }


   // Iterate through the directories.
   for (i = job_start; i >= -1; i--) {

     // If we have a job index and it's not -1 (which indicates no splitting
     // of the indexing into jobs).  Then only continue if i is the same
     // as the job index.
     if (job_index > -2 and i != job_index) {
       continue;
     }
     char dirname[1024];
     char dbname[1024];
     // i == -1 then we will index the 'nan' directory.  This is -1 rather
     // than 101 because we want to index it last.
     if (i == -1) {
       sprintf(dirname, "./%s/nan", indexdir);
       sprintf(dbname, "%s/nan.db", dirname);
     }
     else {
       sprintf(dirname, "./%s/%03d", indexdir, i);
       sprintf(dbname, "%s/%03d.db", dirname, i);
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

     // Turn on extended result codes:
     sqlite3_extended_result_codes(db, 1);

     printf("\nPreparing database file %s...\n", dbname);
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
       fprintf(stderr, "Failed to clear comps table\n");
       handleSQLiteError(db);
       exit(-1);
     }
     sprintf(q, "%s", "DELETE FROM clusters;");
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to clear clusters table\n");
       handleSQLiteError(db);
       exit(-1);
     }
     sprintf(q, "%s", "DELETE FROM summary;");
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to clear summary table\n");
       handleSQLiteError(db);
       exit(-1);
     }
     sprintf(q, "%s", "DELETE FROM mode_hist;");
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to clear mode_hist table\n");
       handleSQLiteError(db);
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
       regex_t re;
       char pattern[64] = "[0-9][0-9][0-9].txt";
       int  rc;
       size_t nmatch = 2;
       regmatch_t pmatch[2];
       if (0 != (rc = regcomp(&re, pattern, 0))) {
           printf("regcomp() failed, returning nonzero (%d)\n", rc);
           exit(-1);
       }
       if (regexec(&re, filename, nmatch, pmatch, 0)) {
         regfree(&re);
         continue;
       }
       regfree(&re);

       // Construct the full path to the file.
       statm_t * memory = memory_get_usage();
       char filepath[1024];
       sprintf(filepath, "%s/%s", dirname, filename);
       printf("Indexing file %s.  Mem: %ldb.\r", filename, memory->size);
       free(memory);

       // Index the file.
       IndexFile(db, filepath, mode_hist, &total_comps);

     }

     // Now add in the stats:
     for (int i = 0; i < max_modes; i++) {
       sprintf(q, "INSERT INTO mode_hist (mode, num_comps) VALUES (%d, %d);", i, mode_hist[i]);
       if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
         fprintf(stderr, "Failed to insert mode_hist record\n");
         handleSQLiteError(db);
         exit(-1);
       }
     }
     sprintf(q, "INSERT INTO summary (num_comps) VALUES (%d);", total_comps);
     if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
       fprintf(stderr, "Failed to insert mode_hist record\n");
       handleSQLiteError(db);
       exit(-1);
     }
     free(mode_hist);

     // Terminate the transaction
     sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);

     // Close the database.
     sqlite3_close(db);

     free(dir);
   }
}
/**
 *
 */
void SQLiteIndexer::IndexFile(sqlite3 *db, char * filepath, int *mode_hist, int * total_comps) {
  int cluster_len = 24800;
  char q[cluster_len];
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
  rc = sqlite3_prepare_v2(db, q, cluster_len, &comp_select_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to prepare statement\n");
    handleSQLiteError(db);
    exit(-1);
  }
  // Prepare the SQL statements for inserting a pair-wise comp.
  sprintf(q, "%s", "\
    INSERT INTO comps \
      (gene1, gene2, num_clusters) VALUES (?, ?, ?);");
  rc = sqlite3_prepare_v2(db, q, cluster_len, &comp_insert_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement\n");
    handleSQLiteError(db);
    exit(-1);
  }
  // Prepare the SQL statements for selecting a cluster.
  sprintf(q, "%s", "SELECT * FROM clusters WHERE gene1 = ? and gene2 = ? and cluster_num = ?;");
  rc = sqlite3_prepare_v2(db, q, cluster_len, &cluster_select_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to prepare statement\n");
    handleSQLiteError(db);
    exit(-1);
  }
  // Prepare the SQL statements for inserting a cluster.
  sprintf(q, "%s", "\
    INSERT INTO clusters \
      (gene1, gene2, num_clusters, cluster_num, num_missing, num_samples, score, samples) \
      VALUES (?, ?, ?, ?, ?, ?, ?, ?);");
  rc = sqlite3_prepare_v2(db, q, cluster_len, &cluster_insert_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement\n");
    handleSQLiteError(db);
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

    // See if this comparison already exists.  If it doesn't then add it.
    sqlite3_bind_int(comp_select_stmt, 1, j);
    sqlite3_bind_int(comp_select_stmt, 2, k);
    if (sqlite3_step(comp_select_stmt) != SQLITE_ROW) {
      sqlite3_bind_int(comp_insert_stmt, 1, j);
      sqlite3_bind_int(comp_insert_stmt, 2, k);
      sqlite3_bind_int(comp_insert_stmt, 3, num_clusters);
      if(sqlite3_step(comp_insert_stmt) != SQLITE_DONE) {
        fprintf(stderr, "Failed to insert comparison\n");
        handleSQLiteError(db);
        exit(-1);
      }
    }

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
        fprintf(stderr, "Failed to insert cluster.\n");
        handleSQLiteError(db);
        exit(-1);
      }
    }

    sqlite3_reset(comp_select_stmt);
    sqlite3_reset(comp_insert_stmt);
    sqlite3_reset(cluster_select_stmt);
    sqlite3_reset(cluster_insert_stmt);
  }
  sqlite3_finalize(comp_select_stmt);
  sqlite3_finalize(comp_insert_stmt);
  sqlite3_finalize(cluster_select_stmt);
  sqlite3_finalize(cluster_insert_stmt);

  fclose(fp);
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
    fprintf(stderr, "Failed to execute statement\n");
    handleSQLiteError(db);
    exit(-1);
  }

  // Prepare the SQL statements for inserting data into the SQLite database.
  sprintf(q, "%s", "\
    INSERT INTO genes \
      (gene_id, name) VALUES (?, ?);");
  rc = sqlite3_prepare_v2(db, q, 2048, &gene_insert_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement\n");
    handleSQLiteError(db);
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
      fprintf(stderr, "Failed to insert gene\n");
      handleSQLiteError(db);
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
    fprintf(stderr, "Failed to execute statement\n");
    handleSQLiteError(db);
    exit(-1);
  }

  // Prepare the SQL statements for inserting data into the SQLite database.
  sprintf(q, "%s", "\
    INSERT INTO samples \
      (sample_id, name) VALUES (?, ?);");
  rc = sqlite3_prepare_v2(db, q, 2048, &sample_insert_stmt, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Failed to execute statement.\n");
    handleSQLiteError(db);
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
    int sample_len = strlen(sample);
    sqlite3_bind_text(sample_insert_stmt, 2, sample, sample_len, SQLITE_STATIC);
    if(sqlite3_step(sample_insert_stmt) != SQLITE_DONE) {
      fprintf(stderr, "Failed to insert gene\n");
      handleSQLiteError(db);
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
  char* errorMessage;
  int max_gene_len = ematrix->getMaxGeneLen();
  int max_sample_len = ematrix->getMaxSampleLen();

  //
  // Create the gene table.
  //
  sprintf(q, "\
    CREATE TABLE IF NOT EXISTS genes \
      (gene_id INTEGER PRIMARY KEY, \
       name CHAR(%d) NOT NULL);", max_gene_len);
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create '%s' table for index '%s'.\n", "genes", dbname);
    handleSQLiteError(db);
    exit(-1);
  }

  sprintf(q, "%s", "CREATE UNIQUE INDEX IF NOT EXISTS genes_name ON genes(name);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "genes_name", "genes", dbname);
    handleSQLiteError(db);
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
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create '%s' table for index '%s'.\n", "samples", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE UNIQUE INDEX IF NOT EXISTS sample_name ON samples(name);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "genes_name", "samples", dbname);
    handleSQLiteError(db);
    exit(-1);
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
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create '%s' table for index '%s'.\n", "comps", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_gene1 ON comps(gene1);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "comps_gene1", "comps", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_gene2 ON comps(gene2);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "comps_gene1", "comps", dbname);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_g12 ON comps(gene1, gene2);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "comps_g12", "comps", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS comps_num_clusters ON comps(num_clusters);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "comps_num_clusters", "comps", dbname);
    handleSQLiteError(db);
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
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create '%s' table for index '%s'.\n", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_gene1 ON clusters(gene1);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_gene1", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_gene2 ON clusters(gene2);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_gene1", "clusters", dbname);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_g12 ON clusters(gene1, gene2);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_g12", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_g12num ON clusters(gene1, gene2, cluster_num);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_g12num", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_clusters ON clusters(num_clusters);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_num_clusters", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_score ON clusters(score);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_score", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_missing ON clusters(num_missing);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_num_missing", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }
  sprintf(q, "%s", "CREATE INDEX IF NOT EXISTS clusters_num_samples ON clusters(num_samples);");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create index '%s' on '%s' table for index '%s'.\n", "clusters_num_samples", "clusters", dbname);
    handleSQLiteError(db);
    exit(-1);
  }


  //
  // Create the mode histogram table
  //
  sprintf(q, "%s",  "\
    CREATE TABLE IF NOT EXISTS mode_hist \
    (mode INTEGER NOT NULL, \
     num_comps INTEGER NOT NULL); ");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
    printf("ERROR: failed to create '%s' table for index '%s'.\n", "mode_hist", dbname);
    handleSQLiteError(db);
    exit(-1);
  }

  //
  // Create the summary table
  //
  sprintf(q, "%s",  "\
      CREATE TABLE IF NOT EXISTS summary \
      (num_comps INTEGER NOT NULL); ");
  if(sqlite3_exec(db, q, NULL, NULL, &errorMessage) != SQLITE_OK) {
      printf("ERROR: failed to create '%s' table for index '%s'.\n", "summary", dbname);
      handleSQLiteError(db);
      exit(-1);
    }

}

/**
 *
 */
void SQLiteIndexer::handleSQLiteError(sqlite3 *db) {
  printf("ERROR: %s (%d).\n", sqlite3_errmsg(db), sqlite3_extended_errcode(db));
}
