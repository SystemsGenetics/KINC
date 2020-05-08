#!/usr/bin/env Rscript

library("optparse")

option_list = list(
   make_option(c("--net"), type="character", default=NULL,
              help="The path to the KINC-derived network file",
              metavar="character"),
   make_option(c("--emx"), type="character", default=NULL,
              help="The path to the log2 transformed Gene Expression Matrix or Metabolite abundance matrix. This must the be same file used by KINC to create the network.", metavar="character"),
   make_option(c("--threads"), type="numeric", default=0,
              help="(optional). The number of computational threads to use for parallel processing. By default, all but 2 cores available on the local machine will be used.",
              metavar="numeric"),
   make_option(c("--wa_th"), type="numeric", default=1e-3,
              help="(optional). The signficance threshold for performing the Welch Anova test. This test checks for differential expression of the GMM cluster vs non-cluster for each gene, thus, two tests are performed: one for each gene in the edge.  Edges with both tests having P-values below this value are kept. The default is 1e-3.",
              metavar="numeric"),
   make_option(c("--mtt_th"), type="numeric", default=0.1,
              help="(optional). The signficance threshold for performing the paired T-test for missingness. This test checks for signicant difference in missingness between the two genes of an edge. This is important because one gene with a high level of missginess will bias the relationship if that missigness is condition-specific. Only edges whose genes have the same pattern of missingness should be considered.  Those with p-values greater than the threshold are considered non-different. The default is 0.1.", metavar="numeric"),
   make_option(c("--out_prefix"), type="character", default=NULL,
              help="(optional). The file name prefix used for the ouput files. If this arugment is not provided then the original network file name is used as a prefix.", metavar="character"),
   make_option(c("--suffix"), type="character", default=NULL,
              help="(optional). The file ending to add to the final output files. By default this is GCN.txt"),
   make_option(c("--chunk_size"), type="numeric", default=1000000,
              help="(optional). When processing large networks, the file is divided into chunks to  prevent  overunning of local memory.  Raise or lower the size of the chunk to use or reduce memory usage.")
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
if (is.null(opt$csGCNs)) {
    opt$csGCNs = FALSE
}

# Determine the suffix for the file.
suffix = 'GCN.txt'
if (length(opt$suffix) > 0) {
   suffix = opt$suffix
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
message("Loading the expression matrix file...")
ematrix = read.table(opt$emx, header=TRUE, sep="\t")

# Load the network.
net_prefix = basename(opt$net)
if (length(opt$out_prefix) > 0) {
    net_prefix = opt$out_prefix
}

# Find out how many threads we should use.
threads = opt$threads
if (threads == 0) {
   threads = detectCores() - 2
   if (threads < 1) {
      threads = 1
   }
}

net_size = getKINCNetworkSize(opt$net)
chunk_size = opt$chunk_size
num_chunks = ceiling(net_size/chunk_size)

# Filter biased edges
message("Filtering the network for biased edges...")
message(paste0("  Num threads: max allowed - 2"))
message(paste0("  GCE Welch's Anova test threshold: ", opt$wa_th))
message(paste0("  Missigness T-test threshold: ", opt$mtt_th))
message(paste0("  Output file prefix: ", opt$out_prefix))
message(paste0("  Network Size: ", net_size))
message(paste0("  Chunk Size: ", opt$chunk_size))
message(paste0("  Number of chunks: ", num_chunks))
message("")

# Set this option so that the progress bar appears.
pbo = pboptions(type="txt")

# Iterate through the chunks and filter biased edges.
for (i in 1:num_chunks) {
    skip = chunk_size * (i - 1)
    chunk_end = chunk_size * (i - 1) + chunk_size
    nrows = chunk_size
    if (chunk_end > net_size)
      nrows = net_size - skip + 1
      chunk_end = skip + nrows

    if (skip == 0) {
        # Don't include the header in the count
        message(paste0("Working on chunk: ", i, ". Edges ", skip + 1, ' to ', chunk_end - 1))
    } else {
        message(paste0("Working on chunk: ", i, ". Edges ", skip, ' to ', chunk_end - 1))
    }

    chunk_net = loadKINCNetwork(opt$net, nrows, skip)
    chunk_net = filterBiasedEdges(chunk_net, ematrix, threads=threads, wa_th=opt$wa_th, mtt_th=opt$mtt_th)

    net_name = paste0(net_prefix, '-filtered.', suffix)
    if (i == 1) {
        saveKINCNetwork(chunk_net, net_name, append=FALSE)
    } else {
        saveKINCNetwork(chunk_net, net_name, append=TRUE)
    }
}
