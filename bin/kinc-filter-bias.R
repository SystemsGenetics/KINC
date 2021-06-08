#!/usr/bin/env Rscript

library("optparse")

option_list = list(
   make_option(c("--net"), type="character", default=NULL,
              help="The path to the KINC-derived network file",
              metavar="character"),
   make_option(c("--emx"), type="character", default=NULL,
              help="The path to the log2 transformed Gene Expression Matrix or Metabolite abundance matrix. This must the be same file used by KINC to create the network.", metavar="character"),
   make_option(c("--threads"), type="numeric", default=0,
              help="(optional). The number of computational threads to use for parallel processing. By default, all but 2 cores available on the local machine will be used.", metavar="numeric"),
   make_option(c("--wa_th"), type="numeric", default=1e-4,
              help="(optional). The p-value threshold for performing the Welch Anova test. This test checks for differential expression of the GMM cluster vs non-cluster for each gene, thus, two tests are performed: one for each gene in the edge.  Edges with both tests having P-values below this value are kept. The default is 1e-3.",
              metavar="numeric"),
   make_option(c("--amx"), type="character", default=NULL,
              help="Optional. Only needed if using the --wa_base option. The path to the sample annotation matrix. This must the be same file used by KINC to perform the cond-test function.", metavar="character"),
   make_option(c("--wa_base"), type="character", default=NULL,
              help="Optional. If the dataset consists of multiple categorical conditions (e.g. multiple treatments) then you can improve the identification of biased edges by providing the name of the categorical condition that represents the base or control condition. If this option is used then only the control samples will be used for testing the Welch Anova test. You must specify the column in the annotation matrix (specified with the --amx option), a comma, then the category representing the control. For example if the column header is \"Treatment\" and the category is \"Control\", you would specify a value of \"Treatment,Control\" for this argument. Be sure to provide proper capitalization and no spaces after the comma.", metavar="character"),
   make_option(c("--mtt_th"), type="numeric", default=0.1,
              help="(optional). The p-value threshold for performing the paired T-test for missingness. This test checks for signicant difference in missingness between the two genes of an edge. This is important because one gene with a high level of missginess will bias the relationship if that missigness is condition-specific. Only edges whose genes have the same pattern of missingness should be considered.  Those with p-values greater than the threshold are considered non-different. The default is 0.1.", metavar="numeric"),
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
if (is.null(opt$amx) & !is.null(opt$wa_base)){
  print_help(opt_parser)
  stop("If you provie the --wa_base argument you must also provide the --amx argument.", call.=FALSE)
}
if (!is.null(opt$amx) & is.null(opt$wa_base)){
  print_help(opt_parser)
  stop("If you provie the --amx argument you must also provide the --wa_base argument.", call.=FALSE)
}


# Determine the suffix for the file.
suffix = 'GCN.txt'
if (length(opt$suffix) > 0) {
   suffix = opt$suffix
}

# Make sure KINC.R is at the correct vresion.
if (packageVersion("KINC.R") < "1.3") {
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

# If a base category was provided then find the sample IDs
wa_samples = c()
if (!is.null(opt$amx) & !is.null(opt$wa_base)) {
    message("Loading the annotation matrix...")
    amx = loadSampleAnnotations(opt$amx)
    base = strsplit(opt$wa_base, ',')
    if (length(base[[1]]) != 2) {
        print_help(opt_parser)
        stop("The value provided to the --wa_base argument is malformed. It should be the column name, separated from the category with a single column.", call.=FALSE)
    }
    if (!base[[1]][1] %in% colnames(amx)) {
        print_help(opt_parser)
        stop(paste0("The value ", base[[1]][1], " provided as the first element of the --wa_base argument is not present in the annotation matrix."), call.=FALSE)
    }
    wa_samples = which(amx[base[[1]][1]] == base[[1]][2])
    if (length(wa_samples) == 0) {
        print_help(opt_parser)
        stop(paste0("The category ", base[[1]][2], " provided as the second element of the --wa_base argument cannot be found in the column ", base[[1]][1], " of the annotation amtrix."), call.=FALSE)
    }
    if (length(wa_samples) < 10) {
        print_help(opt_parser)
        stop(paste0("The number of samples with the category ", base[[1]][2], " provided as the second element of the --wa_base argument is too few."), call.=FALSE)
    }
}

# Load the expression matrix, annotation matrix and the network.
message("Loading the expression matrix file...")
ematrix = read.table(opt$emx, header=TRUE, sep="\t")

# Load the network.
message("Loading the network file...")
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
message(paste0("  Welch's Anova test threshold: ", opt$wa_th))
message(paste0("  Missigness T-test threshold: ", opt$mtt_th))
message(paste0("  Output file prefix: ", opt$out_prefix))
message(paste0("  Network Size: ", net_size))
message(paste0("  Chunk Size: ", opt$chunk_size))
message(paste0("  Number of chunks: ", num_chunks))
if (!is.null(opt$amx) & !is.null(opt$wa_base)) {
    message(paste0("  Number of base/control samples specified: ", length(wa_samples)))
}
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
    chunk_net = filterBiasedEdges(chunk_net, ematrix, threads=threads,
        wa_th=opt$wa_th, mtt_th=opt$mtt_th, wa_samples=wa_samples)

    net_name = paste0(net_prefix, '-filtered.', suffix)
    if (i == 1) {
        saveKINCNetwork(chunk_net, net_name, append=FALSE)
    } else {
        saveKINCNetwork(chunk_net, net_name, append=TRUE)
    }
}
