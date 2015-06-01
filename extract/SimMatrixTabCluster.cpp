#include "SimMatrixTabCluster.h"

/**
 * Constructor
 */
SimMatrixTabCluster::SimMatrixTabCluster(EMatrix *ematrix, int quiet, char * method, int x_coord,
    int y_coord, char * gene1, char * gene2, float th, int max_missing,
    int min_cluster_size)
  : SimilarityMatrix(ematrix, quiet, method, x_coord, y_coord, gene1, gene2, th) {

  // Initialize the class members.
  this->num_jobs = 0;
  this->max_missing = max_missing;
  this->min_cluster_size = min_cluster_size;

  getNumJobs();
}
/**
 * Destructor
 */
SimMatrixTabCluster::~SimMatrixTabCluster() {

}
/**
 *
 */
void SimMatrixTabCluster::getNumJobs() {
  // Holds the clusters directory name.
  char dirname[1024];

  // Make sure the output directory exists.
  struct stat st = {0};
  char clusterdir[100];
  sprintf(clusterdir, "clusters-%s", method);
  if (stat(clusterdir, &st) == -1) {
    fprintf(stderr, "The clusters directory is missing. Cannot continue.\n");
    exit(-1);
  }

  // Iterate through the files in the 'nan' directory to count the number of
  // jobs used to generate the clusters, we will then intitialize that many
  // file pointers for each of the 102 file arrays.
//  sprintf(dirname, "./%s/%s", clusterdir, "nan");
  sprintf(dirname, "./%s/%d", clusterdir, 100);
  DIR * dir;
  dir = opendir(dirname);
  if (!dir) {
    fprintf(stderr, "The clusters sub directory, %s, is missing. Cannot continue.\n", dirname);
    exit(-1);
  }
  struct dirent * entry;
  while ((entry = readdir(dir)) != NULL) {
    const char * filename = entry->d_name;
    // Skip the . and .. files.
    if (strcmp(filename, ".") == 0 || strcmp(filename, "..") == 0) {
      continue;
    }
    num_jobs++;
  }
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

   // The 8 fields of the clusters tab file.
   int x, y, cluster_num, num_clusters, cluster_samples, num_missing;
   char * samples = (char *) malloc(sizeof(char) * num_samples);
   float cv;

   if (!quiet) {
     printf("  Creating network files...\n");
   }

   // Open the edges output file and write the headers.
   char outfile_prefix[2048];
   if (max_missing > num_samples) {
     sprintf(outfile_prefix, "%s.%s.th%0.6f.mcs%d.mmINF", file_prefix, method, th, min_cluster_size);
   }
   else {
     sprintf(outfile_prefix, "%s.%s.th%0.6f.mcs%d.mm%d", file_prefix, method, th, min_cluster_size, max_missing);
   }
   sprintf(edges_file, "%s.coexpnet.edges.txt", outfile_prefix);
   edges = fopen(edges_file, "w");
   fprintf(edges, "gene1\tgene2\tsimilarity\tinteraction\tcluster\tnum_clusters\tcluster_samples\tmissing_samples\tsamples\n");

   // The Spearman and Pearson correlation methods will have both negative and
   // positive values, so we want to create separate files for each one.
   if (strcmp(method, "pc") == 0 || strcmp(method, "sc") == 0) {
     sprintf(edgesN_file, "%s.neg.coexpnet.edges.txt", outfile_prefix);
     sprintf(edgesP_file, "%s.pos.coexpnet.edges.txt", outfile_prefix);
     edgesN = fopen(edgesN_file, "w");
     edgesP = fopen(edgesP_file, "w");
     fprintf(edgesN, "gene1\tgene2\tsimilarity\tinteraction\tcluster\tnum_clusters\tcluster_samples\tmissing_samples\tsamples\n");
     fprintf(edgesP, "gene1\tgene2\tsimilarity\tinteraction\tcluster\tnum_clusters\tcluster_samples\tmissing_samples\tsamples\n");
   }

   char clusterdir[100];
   char dirname[1024];
   sprintf(clusterdir, "clusters-%s", method);

   // Iterate in descending order through the files and
   int limit = (int) (th * 100.0);
   for (int i = 100; i >= limit; i--) {
     sprintf(dirname, "./%s/%03d", clusterdir, i);
     for (int j = 0; j < num_jobs; j++) {

       // Open the file
       char path[1024];
       sprintf(path, "%s/%s.clusters.%03d.%03d.txt", dirname, ematrix->getFilePrefix(), i, j + 1);
       FILE * fp = fopen(path, "r");
       if (!fp) {
         fprintf(stderr, "Could not open clustering file: '%s'.\n", path);
         exit(-1);
       }

       // Get the values from the line in the file.
       int matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\n", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &cv, samples);
       while (matches == 8) {
         if (fabs(cv) >= th && cluster_samples >= min_cluster_size  && num_missing <= max_missing) {
           fprintf(edges, "%s\t%s\t%0.8f\tco\t%d\t%d\t%d\t%d\t%s\n", genes[x-1], genes[y-1], cv, cluster_num, num_clusters, cluster_samples, num_missing, samples);

           // if the method id 'pc' (Pearson's correlation) then we will have
           // negative and positive values, and we'll write those to separate files
           if (strcmp(method, "pc") == 0 || strcmp(method, "sc") == 0) {
             if(cv >= 0){
                fprintf(edgesP, "%s\t%s\t%0.8f\tco\t%d\t%d\t%d\t%d\t%s\n", genes[x-1], genes[y-1], cv, cluster_num, num_clusters, cluster_samples, num_missing, samples);
             }
             else {
               fprintf(edgesN, "%s\t%s\t%0.8f\tco\t%d\t%d\t%d\t%d\t%s\n", genes[x-1], genes[y-1], cv, cluster_num, num_clusters, cluster_samples, num_missing, samples);
             }
           }
         }
         matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\n", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &cv, samples);
       }
       fclose(fp);
     }
   }

   free(samples);
   fclose(edges);
   fclose(edgesN);
   fclose(edgesP);
}
