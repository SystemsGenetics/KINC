#include "CLuceneQuery.h"
#include <stdlib.h>

/**
 * The function to call when running the 'index' command.
 */
CLuceneQuery::CLuceneQuery(char * indexdir) : IndexQuery(indexdir) {


}
/**
 * Destructor
 */
CLuceneQuery::~CLuceneQuery() {

}

/**
 * Performs the indexing.
 */
void CLuceneQuery::run(char * outfile, int x_coord, int y_coord) {
  int i;
  standard::StandardAnalyzer analyzer;

  // Open the output file.
  FILE * fp = fopen(outfile ,"w");
  if (!fp) {
    fprintf(stderr, "The output file cannot be opened for writing. Quiting.\n");
    exit(-1);
  }

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
      fprintf(stderr, "WARNING: The indexes directory, %s, is missing. skipping.\n", dirname);
      continue;
    }

    try {
      printf("Reading index in %s...\n", dirname);
      IndexReader * reader = IndexReader::open(dirname);
      IndexSearcher searcher(reader);

      // Convert the x & y coordinates to a wide character.
      BooleanQuery bq;

      if (x_coord) {
        TCHAR temp[64];
        char xstr[64];
        sprintf(xstr, "%d", x_coord);
        mbstowcs(temp, xstr, 64);
        // true, false == MUST OCCUR
        // false, false == SHOULD OCCUR
        bq.add(QueryParser::parse(temp, _T("gene1"), &analyzer), false, false);
        bq.add(QueryParser::parse(temp, _T("gene2"), &analyzer), false, false);
      }
      if (y_coord) {
        TCHAR temp[64];
        char ystr[64];
        sprintf(ystr, "%d", y_coord);
        mbstowcs(temp, ystr, 64);
        bq.add(QueryParser::parse(temp, _T("gene1"), &analyzer), false, false);
        bq.add(QueryParser::parse(temp, _T("gene2"), &analyzer), false, false);
      }

      // Perform the search.
      Hits * h = searcher.search(&bq);
      // Iterate through the hits.
      unsigned int num_hits = h->length();
      for (unsigned int i = 0; i < num_hits; i++) {
        Document * doc = &h->doc(i);
        fwprintf(fp, L"%ls\t%ls\t%ls\t%ls\t%ls\t%ls\t%ls\t%ls\n",
            doc->get(_T("gene1")),
            doc->get(_T("gene2")),
            doc->get(_T("cluster_num")),
            doc->get(_T("num_clusters")),
            doc->get(_T("cluster_samples")),
            doc->get(_T("num_missing")),
            doc->get(_T("similarity")),
            doc->get(_T("samples")));
      }
      _CLLDELETE(reader);
    }
    catch(CLuceneError& err){
      printf("Error: %s\n", err.what());
      exit(-1);
    }
    catch(...){
      printf("Unknown error\n");
      exit(-1);
    }
  }
  fclose(fp);
}
