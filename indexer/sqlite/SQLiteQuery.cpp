#include "SQLiteQuery.h"

/**
 * The function to call when running the 'index' command.
 */
SQLiteQuery::SQLiteQuery(char * indexdir, EMatrix * ematrix) :
  IndexQuery(indexdir, ematrix) {


}
/**
 * Destructor
 */
SQLiteQuery::~SQLiteQuery() {

}

/**
 * Performs the indexing.
 */
void SQLiteQuery::run(int x_coord, int y_coord, float score) {
  int i;
  char q[2048];
  sqlite3_stmt *cluster_select_stmt;
  int rc;
  int j, k, cluster_num, num_clusters, num_missing, num_samples;
  char * samples;
  float cv;

  // Iterate through the directories.
  for (i = 100; i >= 0; i--) {

    // If a threshold is set then don't examine indexes below the requested value.
    if (i / 100.0 <= score) {
      continue;
    }

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

    // Build the SQL query for filtering.
    sprintf(q, "%s", "\
      SELECT gene1, gene2, num_clusters, cluster_num, num_missing, \
        num_samples, score, samples  \
        FROM clusters C WHERE (1=1) ");
    if (x_coord) {
      sprintf(q, "%s AND (gene1 = %d OR gene2 = %d) ", q, x_coord, x_coord);
    }
    if (y_coord) {
      sprintf(q, "%s AND (gene1 = %d OR gene2 = %d) ", q, y_coord, y_coord);
    }
    if (score) {
      sprintf(q, "%s AND (score >= %f OR score <= -%f) ", q, score, score);
    }
    // sprintf(q, "%s ORDER BY score DESC;", q);
    sprintf(q, "%s;", q);
 // printf("%s\n", q);

    // Prepare the query for execution.
    rc = sqlite3_prepare_v2(db, q, 2048, &cluster_select_stmt, NULL);
    if (rc != SQLITE_OK) {
      fprintf(stderr, "Failed to prepare cluster select statement: %s\n", sqlite3_errmsg(db));
      exit(-1);
    }

    // Iterate through the results of the query
    rc = sqlite3_step(cluster_select_stmt);
    while (rc == SQLITE_ROW) {

      // Get the values for this row in the result set.
      j  = sqlite3_column_int(cluster_select_stmt, 0);
      k  = sqlite3_column_int(cluster_select_stmt, 1);
      num_clusters  = sqlite3_column_int(cluster_select_stmt, 2);
      cluster_num  = sqlite3_column_int(cluster_select_stmt, 3);
      num_missing  = sqlite3_column_int(cluster_select_stmt, 4);
      num_samples  = sqlite3_column_int(cluster_select_stmt, 5);
      cv  = sqlite3_column_double(cluster_select_stmt, 6);
      samples  = (char *) sqlite3_column_text(cluster_select_stmt, 7);

      // Write the line to output.
      printf("%s\t%s\t%0.8f\tco\t%d\t%d\t%d\t%d\t%s\n", ematrix->getGene(j), ematrix->getGene(k), cv,
          cluster_num, num_clusters, num_samples, num_missing, samples);

      // Get the next row.
      rc = sqlite3_step(cluster_select_stmt);
    }
    // Finalize the statement.
    sqlite3_finalize(cluster_select_stmt);

    // Close the database.
    sqlite3_close(db);
  }
}
