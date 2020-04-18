#!/usr/bin/env Rscript

library("optparse")

option_list = list(
   make_option(c("--net"), type="character", default=NULL,
              help="The path to the KINC-derived network file. It must be in tidy format.",
              metavar="character"),
      make_option(c("--out_prefix"), type="character", default=NULL,
              help="(optional). The file name prefix used for the ouput files. If this arugment is not provided then the original network file name is used as a prefix.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$net)){
  print_help(opt_parser)
  stop("Please provide a network file (--net).", call.=FALSE)
}

library("KINC.R")

message("")
message("KINC")
message("----")

# Make sure the csGCN directory exists.
if (!dir.exists(file.path('figures'))) {
  dir.create(file.path('figures'), showWarnings = FALSE)
}

# Load the network.
net_prefix = basename(opt$net)
if (length(opt$out_prefix) > 0) {
    net_prefix = opt$out_prefix
}
message("Loading the KINC network file...")
net = loadKINCNetwork(opt$net)

# Save the summary plots before filtering
message("Generating plots...")
if (!dir.exists(file.path('figures'))) {
  dir.create(file.path(location), showWarnings = FALSE)
}
setwd('figures')

# Disable warnings from the plot function
options( warn = -1 )
saveKINCplots(net, out_prefix=net_prefix)
