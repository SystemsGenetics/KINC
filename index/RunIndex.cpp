#include "RunIndex.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunIndex::printUsage() {
  printf("\n");
  printf("Usage: ./kinc index [options]\n");
  printf("The list of required options:\n");
  printf("  --outdir|-o     The KINC output directory from a prevous run.\n");
  printf("  --samples|-g    The number of samples in the original ematrix.\n");
  printf("                  This corresponds to the --cols argument of the\n");
  printf("                  similarity program of KINC\n");
  printf("\n");
  printf("For Help:\n");
  printf("  --help|-h       Print these usage instructions\n");
  printf("\n");
}

/**
 * The function to call when running the 'index' command.
 */
RunIndex::RunIndex(int argc, char *argv[]) {
  nsamples = 0;

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
      {"outdir",       required_argument, 0,  'o' },
      {"nsamples",     required_argument, 0,  's' },
      // Last element required to be all zeros.
      {0, 0, 0,  0 }
     };

     // get the next option
     c = getopt_long(argc, argv, "o:h", long_options, &option_index);

     // if the index is -1 then we have reached the end of the options list
     // and we break out of the while loop
     if (c == -1) {
       break;
     }

     // handle the options
     switch (c) {
       case 0:
         break;
       case 'o':
         outdir = optarg;
         break;
       case 's':
         nsamples = atoi(optarg);
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

  // Make sure an out file directory is provided
  if (!outdir) {
    fprintf(stderr, "Please provide the KINC output directory from a previous run (--outdir option).\n");
    exit(-1);
  }

  if (nsamples < 1) {
    fprintf(stderr, "Please provide a positive integer for the number of genes (--nsamples option).\n");
    exit(-1);
  }

  // Make sure the output directory exists.
  struct stat st = {0};
  if (stat(outdir, &st) == -1) {
    fprintf(stderr, "The specified output directory ,'%s', is missing. Please check the value of the --outdir argument.\n", outdir);
    exit(-1);
  }

}
/**
 * Destructor
 */
RunIndex::~RunIndex() {

}

/**
 * Performs the indexing.
 */
void RunIndex::execute() {
  int i;

   // Iterate through the directories.
   for (i = 100; i >= 0; i--) {
     char dirname[1024];
     if (i > 0) {
       sprintf(dirname, "./%s/%03d", outdir, i);
     }
     else {
       sprintf(dirname, "./%s/nan", outdir);
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
       indexFile(writer, filepath);
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
void RunIndex::indexFile(IndexWriter * writer, char * filepath) {
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
    printf("Error: %s\n", err.what());
    exit(-1);
  }
  catch(...){
    printf("Unknown error\n");
    exit(-1);
  }

  //clears all static memory
  //_lucene_shutdown();
}
