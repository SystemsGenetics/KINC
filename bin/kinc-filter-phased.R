#!/usr/bin/env Rscript

library("optparse")

option_list = list(
   make_option(c("--net"), type="character", default=NULL,
              help="The path to the KINC-derived network file. It must be in tidy format.",
              metavar="character"),
   make_option(c("--emx"), type="character", default=NULL,
              help="The path to the log2 transformed Gene Expression Matrix or Metabolite abundance matrix. This must the be same file used by KINC to create the network.", metavar="character"),
   make_option(c("--amx"), type="character", default=NULL,
              help="The path to the sample annotation matrix. This must the be same file used by KINC to perform the cond-test function.", metavar="character"),
  make_option(c("--test_col"), type="character", default=NULL,
              help="The name of the column in the sample annotation matrix that has the labels that should be tested for phased (or differential expression).", metavar="character"),
  make_option(c("--sample_col"), type="character", default=NULL,
              help="(optional). The name of the column in the sample annotation matrix that contains the sample names. The default is 'Sample'.", metavar="character"),
  make_option(c("--min_cluster_size"), type="numeric", default=10,
              help="(optional). The minimum number of samples, per label, that must exist for the label to be tested. The labels are the categories in the 'test_col' of the annotation matrix. ", metavar="numeric"),
   make_option(c("--out_prefix"), type="character", default=NULL,
              help="(optional). The file name prefix used for the ouput files. If this arugment is not provided then the original network file name is used as a prefix.", metavar="character"),
   make_option(c("--threads"), type="numeric", default=0,
              help="(optional). The number of computational threads to use for parallel processing. By default, all but 2 cores available on the local machine will be used.",
              metavar="numeric"),
   make_option(c("--th"), type="numeric", default=1e-3,
              help="(optional). The p-value threshold for performing the Hotelling test. This test checks for differential expression within the GMM cluster vs for the categorical column specifed by 'test_col'. Edges having P-values below this value are kept. The default is 1e-3.",
              metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$net)){
  print_help(opt_parser)
  stop("Please provide a network file (--net).", call.=FALSE)
}
if (is.null(opt$emx)){
  print_help(opt_parser)
  stop("Please provide the expression matrix file (--emx).", call.=FALSE)
}
if (is.null(opt$amx)){
  print_help(opt_parser)
  stop("Please provide the sample annotation matrix file (--amx).", call.=FALSE)
}
if (is.null(opt$test_col)){
    print_help(opt_parser)
    stop("Please provide the column in the annotation matrix file that should be tested for differential expression (--test_col).", call.=FALSE)
}
if (is.null(opt$sample_col)){
  opt$sample_col = 'Sample'
}
if (is.null(opt$th) || is.na(opt$th)){
  opt$th = 1e-3
}
if (is.null(opt$min_cluster_size) || is.na(opt$min_cluster_size)){
  opt$th = 10
}

net_prefix = basename(opt$net)
if (length(opt$out_prefix) > 0) {
    net_prefix = opt$out_prefix
}

# Make sure KINC.R is at the correct vresion.
if(packageVersion("KINC.R") < "1.3") {
    stop("This script requires KINC.R > 1.3")
}
suppressMessages(library("KINC.R"))

message("")
message(' __  __   ______   __  __  ____       ____        ')
message('/\\ \\/\\ \\ /\\__  _\\ /\\ \\/\\ \\/\\  _`\\    /\\  _`\\      ')
message('\\ \\ \\/\'/\'\\/_/\\ \\/ \\ \\ `\\\\ \\ \\ \\/\\_\\  \\ \\ \\L\\ \\    ')
message(' \\ \\ , <    \\ \\ \\  \\ \\ , ` \\ \\ \\/_/_  \\ \\ ,  /    ')
message('  \\ \\ \\\\`\\   \\_\\ \\__\\ \\ \\`\\ \\ \\ \\L\\ \\__\\ \\ \\\\ \\   ')
message('   \\ \\_\\ \\_\\ /\\_____\\\\ \\_\\ \\_\\ \\____/\\_\\\\ \\_\\ \\_\\ ')
message('    \\/_/\\/_/ \\/_____/ \\/_/\\/_/\\/___/\\/_/ \\/_/\\/ / ')
message("")
message("This script uses KINC.R, a companion R library for KINC")
message("https://github.com/SystemsGenetics/KINC.R")
message("-------------------------------------------------------")


# Load the expression matrix, annotation matrix and the network.
message("Loading the KINC network file...")
net = loadKINCNetwork(opt$net)

message("Loading the expression matrix file...")
emx = loadGEM(opt$emx)

message("Loading the sample annotation matrix file...")
amx = loadSampleAnnotations(opt$amx, opt$sample_col)


# Find out how many threads we should use.
threads = opt$threads
if (threads == 0) {
   threads = detectCores() - 2
   if (threads < 1) {
      threads = 1
   }
}

# Filter biased edges
message("Filtering the network for phased edges...")
message(paste0("  Num threads: max allowed - 2"))
message(paste0("  Hotellings test threshold: ", opt$th))
message(paste0("  Sample Column: ", opt$sample_col))
message(paste0("  Test Column: ", opt$test_col))
message(paste0("  Minimum Cluster Size: ", opt$min_cluster_size))
message("")

# Set this option so that the progress bar appears.
pbo = pboptions(type="txt")

message("Filtering phased edges...")
net2 = performEdgeDiff(net, emx, amx, test_column = opt$test_col,
    sample_col = opt$sample_col, min_cluster_size = opt$min_cluster_size, threads=threads)

message("Saving phased network...")
phased_net = net2[which(net2['hotelling_p_value'] < opt$th),]
net_name = paste0(net_prefix, '.phased.csGCN.txt')
saveKINCNetwork(phased_net, net_name)

message("Saving non-phased network...")
unphased_net = net2[which(net2['hotelling_p_value'] >= opt$th),]
net_name = paste0(net_prefix, '.non-phased.csGCN.txt')
saveKINCNetwork(unphased_net, net_name)
