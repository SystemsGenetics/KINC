#include "SimMatrixTabCluster.h"

/**
 * Constructor
 */
SimMatrixTabCluster::SimMatrixTabCluster(EMatrix *ematrix, int quiet, char ** method, int num_methods,
    char * th_method, int x_coord,
    int y_coord, char * gene1, char * gene2, float th, int max_missing,
    int min_cluster_size, int max_modes)
  : SimilarityMatrix(ematrix, quiet, method, num_methods, th_method, x_coord, y_coord, gene1, gene2, th) {

  // Initialize the class members.
  this->max_missing = max_missing;
  this->min_cluster_size = min_cluster_size;
  this->max_modes = max_modes;
}
/**
 * Destructor
 */
SimMatrixTabCluster::~SimMatrixTabCluster() {

}

/**
 *
 */
void SimMatrixTabCluster::writeNetwork() {
   // the three network output files for edges, negative correlated edges
   // and positive correlated edges
   FILE * edges;
   FILE * edgesN;
   FILE * edgesP;

   // stores the actual name of the output files
   char edges_file[1024];
   char edgesN_file[1024];
   char edgesP_file[1024];

   // Get some information from the EMatrix object.
   char * file_prefix = ematrix->getFilePrefix();
   char ** genes = ematrix->getGenes();
   int num_samples = ematrix->getNumSamples();

   if (!quiet) {
     printf("  Creating network files...\n");
   }

   // Open the edges output file and write the headers.
   char outfile_prefix[2048];
   if (max_missing > num_samples) {
     sprintf(outfile_prefix, "%s.%s.mcs%d.md%d.mmINF.th%0.6f", file_prefix, th_method, min_cluster_size, max_modes, th);
   }
   else {
     sprintf(outfile_prefix, "%s.%s.mcs%d.md%d.mm%d.th%0.6f", file_prefix, th_method, min_cluster_size, max_modes, max_missing, th);
   }
   sprintf(edges_file, "%s.coexpnet.edges.txt", outfile_prefix);

   // Generate the header string
   char headers[2048];
   sprintf((char *) &headers, "gene1\tgene2\t");
   for (int i = 0; i < num_methods; i++) {
     sprintf((char *) (&headers) + strlen(headers), "%s\t", method[i]);
   }
   sprintf((char *) (&headers) + strlen(headers), "interaction\tcluster\tnum_clusters\tcluster_samples\tmissing_samples\tcluster_outliers\tpair_outliers\ttoo_low\tsamples\n");

   edges = fopen(edges_file, "w");
   fprintf(edges, "%s\n", headers);

   // The Spearman and Pearson correlation methods will have both negative and
   // positive values, so we want to create separate files for each one.
   if (strcmp(th_method, "pc") == 0 || strcmp(th_method, "sc") == 0) {
     sprintf(edgesN_file, "%s.neg.coexpnet.edges.txt", outfile_prefix);
     sprintf(edgesP_file, "%s.pos.coexpnet.edges.txt", outfile_prefix);
     edgesN = fopen(edgesN_file, "w");
     edgesP = fopen(edgesP_file, "w");
     fprintf(edgesN, "%s\n", headers);
     fprintf(edgesP, "%s\n", headers);
   }

   char clusterdir[100];
   char dirname[1024];
   sprintf(clusterdir, "clusters-%s", th_method);

   // Iterate in descending order through the files and
   int limit = (int) (th * 100.0);
   for (int i = 100; i >= limit; i--) {
     sprintf(dirname, "./%s/%03d", clusterdir, i);

     DIR * dir;
     dir = opendir(dirname);
     if (!dir) {
       fprintf(stderr, "The clusters sub directory, %s, is missing. Cannot continue.\n", dirname);
       continue;
     }
     struct dirent * entry;
     while ((entry = readdir(dir)) != NULL) {
       const char * filename = entry->d_name;

       // Skip the . and .. files.
       if (strcmp(filename, ".") == 0 || strcmp(filename, "..") == 0) {
         continue;
       }

       // The file must end in a suffix with 5 digits followed by .txt.
       // We use a regular expression to see if this is true. If not, then
       // skip this file.
       regex_t re;
       char pattern[64] = "[0-9][0-9][0-9][0-9][0-9].txt";
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
       char path[1024];
       sprintf(path, "%s/%s", dirname, filename);

       // Open each file and traverse the elements
       FILE * fp = fopen(path, "r");
       if (!fp) {
         fprintf(stderr, "Can't open file, %s. Cannot continue.\n", path);
         exit(-1);
       }
       int j, k, cluster_num, num_clusters, cluster_samples, num_missing, num_outliers, num_goutliers, num_threshold;
       char samples[num_samples];
       char cscores[255];
       float cv;
       while (!feof(fp)) {

         // Read in the fields for this line. We must read in 11 fields or
         // we will skip the line.
         int matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s", &j, &k, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &num_outliers, &num_goutliers, &num_threshold, (char *) &cscores, (char *) &samples);
         if (matches < 11) {
           char tmp[num_samples*2];
           matches = fscanf(fp, "%s\n", (char *)&tmp);
           continue;
         }

         // Get the score for the selected method.
         float ** scores = parseScores((char *) &cscores);
         cv = *(scores[th_method_index]);

         // Filter the records
         if (fabs(cv) >= th && cluster_samples >= min_cluster_size  &&
             num_missing <= max_missing && num_clusters <= max_modes) {

           fprintf(edges, "%s\t%s\t", genes[j-1], genes[k-1]);
           for (int p = 0; p < num_methods; p++) {
             fprintf(edges, "%0.8f\t", *scores[p]);
           }
           fprintf(edges, "co\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", cluster_num, num_clusters, cluster_samples, num_missing, num_outliers, num_goutliers, num_threshold, samples);

           // if the method id 'pc' (Pearson's correlation) then we will have
           // negative and positive values, and we'll write those to separate files
           if (strcmp(th_method, "pc") == 0 || strcmp(th_method, "sc") == 0) {
             if(cv >= 0){
               fprintf(edgesP, "%s\t%s\t", genes[j-1], genes[k-1]);
               for (int p = 0; p < num_methods; p++) {
                 fprintf(edgesP, "%0.8f\t", *scores[p]);
               }
               fprintf(edgesP, "co\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", cluster_num, num_clusters, cluster_samples, num_missing, num_outliers, num_goutliers, num_threshold, samples);
             }
             else {
               fprintf(edgesN, "%s\t%s\t", genes[j-1], genes[k-1]);
               for (int p = 0; p < num_methods; p++) {
                 fprintf(edgesN, "%0.8f\t", *scores[p]);
               }
               fprintf(edgesN, "co\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", cluster_num, num_clusters, cluster_samples, num_missing, num_outliers, num_goutliers, num_threshold, samples);

             }
           }
         }

         for (int l = 0; l < num_methods; l++) {
           free(scores[l]);
         }
         free(scores);

       }
       fclose(fp);
     }
   }

   fclose(edges);
   fclose(edgesN);
   fclose(edgesP);
}
/**
 *
 */
float ** SimMatrixTabCluster::parseScores(char * scores_str) {
  // Split the method into as many parts
  char * tmp;
  int i = 0;
  tmp = strstr(scores_str, ",");
  float ** scores = (float **) malloc(sizeof(float *) * this->num_methods);
  char tmp_score[255];

  while (tmp) {
    strncpy((char *) &tmp_score, scores_str, (tmp - scores_str));
    scores[i] = (float *) malloc(sizeof(float) * 1);
    *(scores[i]) = atof(tmp_score);
    scores_str = tmp + 1;
    tmp = strstr(scores_str, ",");
    i++;
  }
  // Get the last element of the methods_str.
  strncpy((char *) &tmp_score, scores_str, strlen(scores_str));
  scores[i] = (float *) malloc(sizeof(float) * 1);
  *(scores[i]) = atof(tmp_score);

  return scores;
}
/**
 *
 */
void SimMatrixTabCluster::getPosition() {

  // Use the largest gene index (i or j) and computing the following
  // ((n * n - 1) / 2) / jobs = comps_per_job
  // (((i* i - 1) / 2) / comps_per_job) + 1 = job_index
  // Look through all files with job_index suffix for results
}
