#include "CLuceneIndexer.h"

/**
 * The function to call when running the 'index' command.
 */
CLuceneIndexer::CLuceneIndexer(char * indexdir) : Indexer(indexdir) {

}
/**
 * Destructor
 */
CLuceneIndexer::~CLuceneIndexer() {

}

/**
 * Performs the indexing.
 */
void CLuceneIndexer::run(int nsamples) {
  int i;

  // Iterate through the directories.
  for (i = 100; i >= 0; i--) {
     char dirname[1024];
    if (i > 0) {
      sprintf(dirname, "./%s/%03d", indexdir, i);
    }
    else {
      sprintf(dirname, "./%s/nan", indexdir);
    }

    // Open the output directory.
    DIR * dir;
    dir = opendir(dirname);
    if (!dir) {
      fprintf(stderr, "WARNING: The output directory, %s, is missing. skipping.\n", dirname);
      continue;
    }

    // Create a new file for each directory.
    IndexWriter * writer = NULL;
    lucene::analysis::WhitespaceAnalyzer an;
    if (IndexReader::indexExists(dirname)) {
      if (IndexReader::isLocked(dirname)) {
        printf("Index was locked... unlocking it.\n");
        IndexReader::unlock(dirname);
      }
      writer = _CLNEW IndexWriter(dirname, &an, false);
    }
    else{
      writer = _CLNEW IndexWriter(dirname ,&an, true);
    }

    // LUCENE_INT32_MAX_SHOULDBE
    writer->setMaxFieldLength(0x7FFFFFFFL);
    // Turn this off to make indexing faster; we'll turn it on later before optimizing
    writer->setUseCompoundFile(false);

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
      indexFile(writer, filepath, nsamples);
    }
    writer->setUseCompoundFile(true);
    writer->optimize();

    // Close and clean up the lucene writer.
    writer->close();
    _CLLDELETE(writer);
  }
}

/**
 * Adds an individual file to the Lucene index.
 */
void CLuceneIndexer::indexFile(IndexWriter * writer, char * filepath, int nsamples) {
  // Open the file for indexing.
  FILE * fp = fopen(filepath, "r");
  if (!fp) {
    fprintf(stderr, "Can't open file, %s. Cannot continue.\n", filepath);
    exit(-1);
  }

  try {
    Document doc;

    wchar_t j[128], k[128], cluster_num[128], num_clusters[128], cluster_num_samples[128], num_missing[128];
    wchar_t samples[nsamples + 1];
    wchar_t cv[128];
    while (!feof(fp)) {

      // Read in the fields for this line.
      int matches = fwscanf(fp, L"%s\t%s\%s\t%s\%s\t%s\t%s\t%s\n", &j, &k, &cluster_num, &num_clusters, &cluster_num_samples, &num_missing, &cv, &samples);

      // Skip lines that don't have the proper number of columns
      if (matches < 8) {
        char tmp[nsamples*2];
        matches = fscanf(fp, "%s\n", (char *)&tmp);
        continue;
      }

      // Add each field as an indexable entry for lucene.
      doc.clear();
      doc.add(*_CLNEW Field(_T("gene1"), j, Field::STORE_YES | Field::INDEX_UNTOKENIZED));
      doc.add(*_CLNEW Field(_T("gene2"), k, Field::STORE_YES | Field::INDEX_UNTOKENIZED));
      doc.add(*_CLNEW Field(_T("cluster_num"), cluster_num, Field::STORE_YES | Field::INDEX_NO));
      doc.add(*_CLNEW Field(_T("num_clusters"), num_clusters, Field::STORE_YES | Field::INDEX_NO));
      doc.add(*_CLNEW Field(_T("cluster_samples"), cluster_num_samples, Field::STORE_YES | Field::INDEX_NO));
      doc.add(*_CLNEW Field(_T("num_missing"), num_missing, Field::STORE_YES | Field::INDEX_NO));
      doc.add(*_CLNEW Field(_T("similarity"), cv, Field::STORE_YES | Field::INDEX_NO));
      doc.add(*_CLNEW Field(_T("samples"), samples, Field::STORE_YES | Field::INDEX_NO));
      writer->addDocument(&doc);
    }
    fclose(fp);
  }
  catch(CLuceneError& err){
    fclose(fp);
    printf("Error: %s\n", err.what());
    exit(-1);
  }
  catch(...){
    fclose(fp);
    printf("Unknown error\n");
    exit(-1);
  }

  //clears all static memory
  //_lucene_shutdown();
}

