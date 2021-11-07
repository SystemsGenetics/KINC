#!/usr/bin/env Rscript

library("optparse")

option_list = list(
    make_option(c("--net"), type="character", default=NULL,
              help="The path to the KINC-derived network file",
              metavar="character"),
    make_option(c("--sim_col"), type="character", default='Similarity_Score',
              help="The name of the similarity column in the network file.",
              metavar="character"),
    make_option(c("--out_prefix"), type="character", default=NULL,
              help="(optional). The file name prefix used for the ouput files. If this arugment is not provided then the original network file name is used as a prefix.", metavar="character"),
    make_option(c("--no_meta"), type="logical", default=FALSE, action="store_true",
              help="(optional). Disable creation of meta modules by combinging similar modules together. Creation of meta modules is default"),
    make_option(c("--use_inverse"), type="logical", default=FALSE, action="store_true",
              help="(optional). By default inverse edges are ignored when creating clusters. Use this flag to allow inverse edges to be used in module discovery."),
    make_option(c("--th"), type="double", default=0.5,
              help="(optional). Specifies the Jaccard similarity score between two modules gene content in order for those modules to be merged. Only applies if no_meta is not used. Lower threshold results in merging being more common. Default is 0.5.",
              metavar="numeric"),
    make_option(c("--min_verticies"), type="integer", default=10,
              help="(optional). If a network is disconnected then communities will be found in each subgraph independnet of the others. This argument specifies the number of elements that must exist in a subgraph for communities to be identified. Default is 10.",
              metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$net)){
  print_help(opt_parser)
  stop("Please provide a network file (--net).", call.=FALSE)
}

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
message(paste0("KINC.R v", packageVersion("KINC.R")))
message("-------------------------------------------------------")


# Load the network.
message("Loading the network file...")
net = loadKINCNetwork(opt$net)

# For backwards compatibility, rename the "sc" columnt to "Similarity_Score"
cnames = colnames(net)
if ('sc' %in% cnames) {
  cnames[cnames == "sc"] = 'Similarity_Score'
  colnames(net) = cnames
}


net_prefix = basename(opt$net)
if (length(opt$out_prefix) > 0) {
    net_prefix = opt$out_prefix
}

meta = !opt$no_meta
ignore_inverse = !opt$use_inverse


message("Finding communities...")
if (meta == TRUE) {
    print("...will be creating meta modules.")
} else {
    print("...will NOT be creating meta modules.")
}
if (ignore_inverse == TRUE) {
    print("...inverse edges will NOT be included.")
} else {
    print("...inverse edges will be included.")
}
lcm = findLinkedCommunities(net, net_prefix, meta = meta,
        ignore_inverse = ignore_inverse, th = opt$th,
        min.vertices = opt$min_verticies)
