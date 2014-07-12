#include "extract.h"

/**
 * The function for extracting results from the similarity matrix
 */
int do_extract(int argc, char *argv[]) {

  int c;  // The value returned by getopt_long

  static NetParameters params;

  // initialize some of the program parameters
  strcpy(params.method, "pc");
  params.x_coord = 0;
  params.y_coord = 0;
  params.th = 0;
  params.num_files = 0;

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"ematrix", required_argument, 0,  'e' },
      {"method",  required_argument, 0,  'm' },
      {"th",      required_argument, 0,  't' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };

    // get the next option
    c = getopt_long(argc, argv, "e:m:x:y:t", long_options, &option_index);

    // if the index is -1 then we have reached the end of the options list
    // and we break out of the while loop
    if (c == -1) {
      break;
    }

    // handle the options
    switch (c) {
      case 0:
        break;
      case 'e':
        params.infilename = optarg;
        break;
      case 'x':
        params.x_coord = atoi(optarg);
        break;
      case 'y':
        params.y_coord = atoi(optarg);
        break;
      case 'm':
        strcpy(params.method, optarg);
        break;
      case 't':
        params.th = atof(optarg);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        print_similarity_usage();
        exit(-1);
        break;
      default:
        print_similarity_usage();
    }
  }

  // make sure the required arguments are set and appropriate
  if (!params.infilename) {
    fprintf(stderr,"Please provide an expression matrix (--ematrix option). Use the -h option for help.\n");
    exit(-1);
  }

  if (!params.method) {
    fprintf(stderr,"Please provide the method (--method option) used to construct the similarity matrix.\n");
    exit(-1);
  }

  // make sure we have a positive integer for the rows and columns of the matrix
  if ((params.x_coord >  0 && params.y_coord == 0) ||
      (params.x_coord == 0 && params.y_coord  > 0)) {
    fprintf(stderr,"Please provide a positive integer for both the x and y coordinates (-x and -y options)\n");
    exit(-1);
  }

  if (params.th > 0 && (params.x_coord > 0 || params.y_coord > 0)) {
    fprintf(stderr,"Please provide a threshold or x and y coordinates only but not both\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(params.infilename, F_OK) == -1) {
    fprintf(stderr,"The input file does not exists or is not readable.\n");
    exit(-1);
  }

  if (strcmp(params.method, "pc") != 0 && strcmp(params.method, "mi") != 0) {
    fprintf(stderr,"The method (--method option) must either be 'pc' or 'mi'.\n");
    exit(-1);
  }

  // print out some setup details
  printf("  Using method: '%s'\n", params.method);
  if (params.th > 0) {
    printf("  Using threshold of %f\n", params.th);
  }

  // remove the path and extension from the filename
  char * temp = basename(params.infilename);
  strcpy(params.fileprefix, temp);
  char * p = rindex(params.fileprefix, '.');
  if (p) {
    p[0] = 0;
  }

  if (strcmp(params.method, "mi") == 0) {
    params.inputDir = "MI";
  }
  else if (strcmp(params.method, "pc") == 0) {
    params.inputDir = "Pearson";
  }

  // open all of the bin files for easy access
  open_bin_files(params);

  // close all of the bin files
  close_bin_files(params);

  return 1;
}

/**
 * Looks in the directory where the similarity matrix is kept and opens each
 * of the bin files for reading.
 */
void open_bin_files(NetParameters params) {
  DIR * FD;                  // a directory handle
  struct dirent* curr_file;  // a file structure for the file while looping
  char * bin_pos;            // the position of the '.bin' extension in the file name
  char * method_pos;         // the position of the .mi or .pc in the file name
  char bin_num[5];           // the string representation of the number of the bin file

  // Scanning for the files of the similarity matrix
  if (NULL == (FD = opendir(params.inputDir))) {
    fprintf(stderr, "Error : Failed to open input directory - %s\n", strerror(errno));
    exit(-1);
  }
  while ((curr_file = readdir(FD))) {
    // skip the . and .. dirs
    if (strcmp(curr_file->d_name, ".") == 0 || strcmp(curr_file->d_name, "..") == 0) {
      continue;
    }

    // make sure this file has a .bin extension
    bin_pos = strstr(curr_file->d_name, ".bin");
    if(bin_pos) {
      // make sure the file has the prefix
      if(strstr(curr_file->d_name, params.fileprefix) != NULL) {
        // make sure that this file has the method name preceeded by a period
        char method[5];
        sprintf(method, ".%s", params.method);
        method_pos = strstr(curr_file->d_name, method);
        if (method_pos) {
          // get the numerical value of this bin file
          int size = (bin_pos - (method_pos + 3)) * sizeof(char);
          memcpy(bin_num, method_pos + (3 * sizeof(char)), size);
          bin_num[size] = 0; // add a terminator to the size

          // open the file and store the file handle for later use and store
          // it in the files array using the bin_num as an index.
          printf("  Opening file: %s\n", curr_file->d_name);
          params.files[atoi(bin_num)] = fopen(curr_file->d_name, "w");
          params.num_files++;
        }
      }
    }
  }
}

/**
 * Closes all of the file pointers opened by the open_bin_files
 */
void close_bin_files(NetParameters params) {
  int i;

  for (i = 0; i < params.num_files; i++) {
    fclose(params.files[i]);
  }
}
//
///**
// *
// */
//void get_edges (NetParameters params) {
//   my $t = $_[0];
//   my $files = $_[1];
//   my $psets = $_[2];
//   my ($x,$y,$j,$k,$i,$n,$m) = 0;
//   my $pos = 0;
//   my $rows_seen = 0;
//   my $r_sqr;
//   my $data;
//
//   int bin_i = 0;
//   int x,y;   // used for iterating through the n x n similarity matrix
//   int i;     // used for specifying the bin file to read
//
//   FILE * edges;
//   FILE * edgesN;
//   FILE * edgesP;
//
//   char edges_file[1024];
//   char edgesN_file[1024];
//   char edgesP_file[1024];
//
//   sprintf(edges_file, "%s.coexpnet.edges.txt", params.fileprefix);
//   sprintf(edgesN_file, "%s.neg.coexpnet.edges.txt", params.fileprefix);
//   sprintf(edgesP_file, "%s.pos.coexpnet.edges.txt", params.fileprefix);
//
//   edges = fopen(edges_file, "w");
//   edgesN = fopen(edgesN_file, "w");
//   edgesP = fopen(edgesP_file, "w");
//
//   # get the size of the matrix in one dimension (i.e. mxm)
//   $m = $files->{$bin_i}{cols}; # the matrix is square so the size is same as num of columns
//   seek($files{$bin_i}{fh},8,0); # set the file position just past the row info
//   i = 0;
//   for(x = 0; x < params.numGenes; x++){
//      if(i >= params.numLinesPerFile){
//         bin_i++;
//         $i = 0;
//         seek($files{$bin_i}{fh},8,0); # set the file position just past the row info
//      }
//      for(y = 0; y < params.numGenes; y++){
//
//         $n = read($files{$bin_i}{fh},$data,4);
//         if($n != 4){
//            print "ERROR: cannot fetch from bin file $bin_i\n";
//            exit;
//         }
//
//         last if($y >= $x);
//         $r_sqr = unpack("f",$data);
//         if(($r_sqr > 0 and  $r_sqr >= $t) or
//            ($r_sqr < 0 and -$r_sqr >= $t)){
//            print EDGES $psets->[$x]."\t".$psets->[$y]."\t". sprintf("%0.8f",$r_sqr). "\n";
//            if($r_sqr >= 0){
//               print EDGESP $psets->[$x]."\t".$psets->[$y]."\t". sprintf("%0.8f",$r_sqr). "\n";
//            } else {
//               print EDGESN $psets->[$x]."\t".$psets->[$y]."\t". sprintf("%0.8f",$r_sqr). "\n";
//            }
//         }
//      }
//      $i++;
//   }
//   fclose(edges);
//   fclose(edgesN);
//   fclose(edgesP);
//}
//#------------------------------------------------------------------------------
//sub get_position {
//   my $x = $_[0];
//   my $y = $_[1];
//   my $files = $_[2];
//   my $pos = 0;
//   my $i = 0;
//   my ($j,$k);
//   my $bin_i;
//
//   $x = $x -1;
//   $y = $y -1;
//
//   # if $y > $x then reverse the two
//   my $temp;
//   if($y > $x){
//      $temp = $x;
//      $x = $y;
//      $y = $temp;
//   }
//
//   my $max;
//   for $bin_i (sort {$a <=> $b} keys %{$files}){
//      if($x >= $i and $x < $i + $files->{$bin_i}{rows}){
//         $k = $x - $i;
//         for($j = 0; $j < $k; $j++){
//            $pos = $pos + $i + 1 + $j;
//         }
//         $pos = ($pos + $y)*4 + 8;
//         seek($files{$bin_i}{fh},$pos,0); # set the file position just past the row info
//         $n = read($files{$bin_i}{fh},$data,4);
//         if($n != 4){
//            print "ERROR: cannot fetch from bin file $bin_i\n";
//            exit;
//         }
//         $r_sqr = unpack("f",$data);
//         printf("cor(%i,%i) = %0.8f, bin = $bin_i, pos = %d\n",$x+1,$y+1,$r_sqr,$pos);
//         last;
//      }
//      $i += $files{$bin_i}{rows};
//   }
//}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_extract_usage() {
  printf("\n");
  printf("Usage: ./RMTGeneNet extract [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("                 The rows must be genes or probe sets and columns are samples\n");
  printf("  -x           When extracting a single value, the x coordinate. Must also use -y\n");
  printf("  -y           When extracting a single value, the y coordinate. Must also use -x.\n");
  printf("  --th|t       The threshold to cut the similarity matrix. Network files will be generated.\n");
  printf("  --method|-m  The correlation method to use. Supported methods include\n");
  printf("                 Pearson's correlation and Mutual Information. Provide\n");
  printf("                 either 'pc' or mi' as values respectively.\n");
  printf("\n");
}
