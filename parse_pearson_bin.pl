#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

$| = 1;

my $t;
my $x;
my $y;
my $psets;
my $bin_dir;
my $prefix;
my $help;

# get the incoming arguments
Getopt::Long::Configure ('bundling');
GetOptions (
   'x=s' => \$x,
   'y=s' => \$y,
   't=s' => \$t,
   'p|psets=s'=> \$psets,
   'b|bins=s'=>\$bin_dir,
   'o|prefix=s'=>\$prefix,
   'h|help'=>\$help,
);

if($help){
   printUsage();
   exit;
}

if(!$bin_dir){
  print "Please specify the directory where the correlation .bin files are found (-b option)\n";
  exit;
}
if(!$t and !$x and !$y){
  print "Please provide a threshold (-t option) value or x and y (-x and -y options) coordinates to extract\n";
  exit;
}
if((defined($x) and $x !~ /^\d+$/ and !$psets) or 
   (defined($y) and $y !~ /^\d+$/ and !$psets)){
  print "Please provide the file with the probe set names in the proper order (-p option)\n";
  exit;
}
if($t and ($x or $y)){
  print "please provide only a threshold or a coordiante but not both\n";
  exit;
}
if($t and !$prefix){
  print "please provide a prefix for the output file name (-o option)\n";
  exit;
}

my $bin_file;
my $n;
my $data;
my $r_sqr;
my $bin_i;
my %files;
my @psets;
my $i = 1;

# get the probeset order
if($psets){
   open(PSETS,$psets) or die "$!: $psets";
   while(<PSETS>){
      chomp($_);
      push(@psets,$_);
      # if $x and $y are not numeric then we assume they are probeset
      # names and we need to look them up
      if($x and $x eq $_){
         $x = $i;
      }
      if($y and $y eq $_){
         $y = $i;
      }
      $i++;
   }
   close(PSETS);
}

if(defined($x) and $x !~ /^\d+$/){
   print "Cannot find probeset named '$x' in the dataset (-x option)\n";
   exit;
}
if(defined($y) and $y !~ /^\d+$/){
   print "Cannot find probeset named '$y' in the dataset (-y option)\n";
   exit;
}

# first read the file and get the dimensions and file handles
opendir(BINDIR,$bin_dir) or die $!;
while($bin_file = readdir(BINDIR)){
   next if($bin_file eq '.');
   next if($bin_file eq '..');
   next if($bin_file !~ /\.bin$/);
   if($bin_file =~ /(\d+)\.bin$/){
      $bin_i = $1;
   }
   open($files{$bin_i}{fh},"$bin_dir/$bin_file") or die $!;
   binmode $files{$bin_i}{fh};
   $n = read($files{$bin_i}{fh},$data,4);
   $files{$bin_i}{cols} = unpack("i",$data);
   $n = read($files{$bin_i}{fh},$data,4);
   $files{$bin_i}{rows} = unpack("i",$data);
   if(!$t){
      printf("Bin: %d, Num Genes: %d, Lines in File: %d\n",$bin_i,$files{$bin_i}{cols},$files{$bin_i}{rows});
   }
}
close(BINDIR);

if($t and scalar(@psets) != $files{0}{cols}){
   print "ERROR: The probesets file is the incorrect size. It is expected to be ".scalar(@psets)."\n";
   exit;
}

if(defined($x) and $x > 0 and $y > 0){
   get_position($x,$y,\%files);
}
if(defined($t)){
   get_edges($t,\%files,\@psets);
}

# close the files
for $bin_file (keys %files){
   close($files{$bin_file}{fh});
}
#------------------------------------------------------------------------------
sub get_edges {
   my $t = $_[0];
   my $files = $_[1];
   my $psets = $_[2];
   my ($x,$y,$j,$k,$i,$n,$m) = 0;
   my $bin_i = 0;
   my $pos = 0;
   my $rows_seen = 0;
   my $r_sqr;
   my $data;

   # open the file
   open(EDGES,">$prefix.psetnet.edges.txt") or die "$!: $prefix.psetnet.edges.txt\n";
   open(EDGESN,">$prefix.neg.psetnet.edges.txt") or die "$!: $prefix.neg.psetnet.edges.txt\n";
   open(EDGESP,">$prefix.pos.psetnet.edges.txt") or die "$!: $prefix.pos.psetnet.edges.txt\n";

   # get the size of the matrix in one dimension (i.e. mxm)
   $m = $files->{$bin_i}{cols}; # the matrix is square so the size is same as num of columns
   seek($files{$bin_i}{fh},8,0); # set the file position just past the row info
   $i = 0;
   for($x = 0; $x < $m; $x++){
      if($i >= $files->{$bin_i}{rows}){
         $bin_i++;
         $i = 0;
         seek($files{$bin_i}{fh},8,0); # set the file position just past the row info
      }
      for($y = 0; $y < $m; $y++){

         $n = read($files{$bin_i}{fh},$data,4);
         if($n != 4){   
            print "ERROR: cannot fetch from bin file $bin_i\n";
            exit;
         }

         last if($y >= $x);
         $r_sqr = unpack("f",$data);
         if(($r_sqr > 0 and  $r_sqr >= $t) or 
            ($r_sqr < 0 and -$r_sqr >= $t)){
            print EDGES $psets->[$x]."\t".$psets->[$y]."\t". sprintf("%0.8f",$r_sqr). "\n";
            if($r_sqr >= 0){
               print EDGESP $psets->[$x]."\t".$psets->[$y]."\t". sprintf("%0.8f",$r_sqr). "\n";
            } else {
               print EDGESN $psets->[$x]."\t".$psets->[$y]."\t". sprintf("%0.8f",$r_sqr). "\n";
            }
         }
      }
      $i++;
   }
   close(EDGES);
   close(EDGESP);
   close(EDGESN);
}
#------------------------------------------------------------------------------
sub get_position {
   my $x = $_[0];
   my $y = $_[1];
   my $files = $_[2];
   my $pos = 0;
   my $i = 0;
   my ($j,$k);
   my $bin_i;

   $x = $x -1;
   $y = $y -1;

   # if $y > $x then reverse the two
   my $temp;
   if($y > $x){
      $temp = $x;
      $x = $y;
      $y = $temp;
   }
   
   my $max;
   for $bin_i (sort {$a <=> $b} keys %{$files}){
      if($x >= $i and $x < $i + $files->{$bin_i}{rows}){
         $k = $x - $i;
         for($j = 0; $j < $k; $j++){
            $pos = $pos + $i + 1 + $j;
         }    
         $pos = ($pos + $y)*4 + 8;
         seek($files{$bin_i}{fh},$pos,0); # set the file position just past the row info
         $n = read($files{$bin_i}{fh},$data,4);
         if($n != 4){   
            print "ERROR: cannot fetch from bin file $bin_i\n";
            exit;
         }
         $r_sqr = unpack("f",$data);
         printf("cor(%i,%i) = %0.8f, bin = $bin_i, pos = %d\n",$x+1,$y+1,$r_sqr,$pos);
         last;
      }
      $i += $files{$bin_i}{rows};
   }
}
# --------------------------------------------------------------------------
sub printUsage {
  print "\nUsage: parse_pearson_bin.pl <arguments>\n\n"; 
  print qq(
  This script reads the binary Pearson correlation files generated by the
  RMT ccm program.  The Pearson files can be parsed into network files
  by providing a threshold (-t option). Or the pearson files can be
  queried by providing an x and y coordinates (-x and -y respectively).
  
  If a threshold is provided three different network files are generated.
  The first contains the list of all edges after thresholding and is 
  named 'psetnet.edges.txt' with an additional user provided prefix.  The
  second file contains the list of all positively corrleated edges. This
  file is named 'pos.psetnet.edges.txt'.  The third file contains the
  list of negatively correlated edges and is named 'neg.psetnet.edges.txt'.

  When querying the correlation files using an x and y coordinate, the
  coordinates are 1 indexed meaning the first cell in the matrix is
  x=1, y=1.  Alternatively, you may provde the probeset names for 
  either the x and y coordinates rather than the numeric index. 

  The following options are available.

  -b|--bins <directory>
  (Required). Specifies the path to the location where the .bin files 
  are stroed.

  -x <integer>|<pset name>
  (Optional). Specifies the x coordinate in the correlation matrix. Must
  be used with the -y option. 

  -y <intger>|<pset name>
  (optional).  Specifies the y corrdinate in the correlation matrix. Must
  be used with the -x option.  

  -t <threshold>
  (optional).  Specifies the threshold to use for filtering the correlation
  matrix.  When this option is provided, a network file containing 3 
  columns is generated.  Each row in the file represents and edge with the
  first two columsn indicating the nodes in the edge. The thrid column
  is the Pearson correlation value, or strength of the edge.

  -p|--psets <filename>
  (Required for the -t option). A file with an ordered listing of probesets.
  The position of the probeset in the list should correspond to the index or row
  that belongs to the pset in the correlation matrix.

  -o|--prefix <string>
  (Required with -t option).  If the -t option is provided then the output
  file generated will be prefixed with the string provided to this 
  argument.

  
  );
}

