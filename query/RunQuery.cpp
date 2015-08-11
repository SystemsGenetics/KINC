#include "RunQuery.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunQuery::printUsage() {
  printf("\n");
  printf("Usage: ./kinc query [options]\n");
  printf("The list of required options:\n");
  printf("  --indexes|-i     The KINC directory where the indexed files are housed.\n");
  printf("  --outfile|-o     The file where results should be written.\n");
  printf("  -x               Retrieve results for gene x.\n");
  printf("  -y               Retrieve results for gene y.\n");
  printf("\n");
  printf("For Help:\n");
  printf("  --help|-h       Print these usage instructions\n");
  printf("\n");
}

/**
 * The function to call when running the 'index' command.
 */
RunQuery::RunQuery(int argc, char *argv[]) {
  //x_coord = 0;
  //y_coord = 0;

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  // The value returned by getopt_long.
  int c;
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"help",         no_argument,       0,  'h' },
      {"indexes",      required_argument, 0,  'i' },
      {"outfile",      required_argument, 0,  'o' },
      {"x",            required_argument, 0,  'x' },
      {"y",            required_argument, 0,  'y' },
      // Last element required to be all zeros.
      {0, 0, 0,  0 }
     };

     // get the next option
     c = getopt_long(argc, argv, "i:x:y:h", long_options, &option_index);

     // if the index is -1 then we have reached the end of the options list
     // and we break out of the while loop
     if (c == -1) {
       break;
     }

     // handle the options
     switch (c) {
       case 0:
         break;
       case 'i':
         indexes = optarg;
         break;
       case 'o':
         outfile = optarg;
         break;
       case 'x':
         x_coord = optarg;
         break;
       case 'y':
         y_coord = optarg;
         break;
       // Help and catch-all options.
       case 'h':
         printUsage();
         exit(-1);
         break;
       case '?':
         exit(-1);
         break;
       case ':':
         printUsage();
         exit(-1);
         break;
       default:
         printUsage();
         exit(-1);
     }
   }

  // Make sure an indexes directory is provided
  if (!indexes) {
    fprintf(stderr, "Please provide the KINC indexes directory (--indexes option).\n");
    exit(-1);
  }

  if (!outfile) {
    fprintf(stderr, "Please provide the output file name (--outfile option).\n");
    exit(-1);
  }

  if (atoi(x_coord) < 0) {
    fprintf(stderr, "Please provide a positive integer for the x genes (--x option).\n");
    exit(-1);
  }

/*  if (atoi(y_coord) < 0) {
    fprintf(stderr, "Please provide a positive integer for the x genes (--x option).\n");
    exit(-1);
  }*/

  // Make sure the output directory exists.
  struct stat st = {0};
  if (stat(indexes, &st) == -1) {
    fprintf(stderr, "The specified indexes directory ,'%s', is missing. Please check the value of the --indexes argument.\n", indexes);
    exit(-1);
  }

}
/**
 * Destructor
 */
RunQuery::~RunQuery() {

}

/**
 * Performs the indexing.
 */
void RunQuery::execute() {
  int i;
  standard::StandardAnalyzer analyzer;

  // Open the output file.
  FILE * fp = fopen(outfile ,"w");
  if (!fp) {
    fprintf(stderr, "The output file cannot be opened for writing. Quiting.\n");
    exit(-1);
  }

  // Iterate through the directories.
  for (i = 100; i >= -1; i--) {
    char dirname[1024];
    if (i > 0) {
      sprintf(dirname, "./%s/%03d", indexes, i);
    }
    else {
      sprintf(dirname, "./%s/nan", indexes);
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

      // Convert the x_coord to a wide character.
      TCHAR tmp[128];
      mbstowcs(tmp, x_coord, 128);

      Query * q = QueryParser::parse(tmp, _T("gene1"), &analyzer);
      // Perform the search.
      Hits * h = searcher.search(q);
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
