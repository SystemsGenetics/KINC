#!/usr/bin/env Rscript

library("optparse")

option_list = list(
   make_option(c("--net"), type="character", default=NULL,
              help="The path to the KINC-derived network file",
              metavar="character"),
   make_option(c("--out_prefix"), type="character", default=NULL,
              help="(optional). The file name prefix used for the ouput files. If this arugment is not provided then the original network file name is used as a prefix.",
              metavar="character"),
   make_option(c("--top_n"), type="numeric", default=10000,
              help="(optional). The number of edges to extract from the network file. This is the top n ranked edges per condition. If --no_conditions is provided then this is the top n ranked edges in the full network.",
              metavar="numeric"),
   make_option(c("--no_condition_rankings"), type="logical", action="store_true",
              help="(optional). By deafult the rankings are calculated independently per each condition in the network. To calculate rankings for the entire network regardless of condition set this to TRUE.",
              metavar="logical"),
   make_option(c("--save_condition_networks"), type="logical", action="store_true",
              help="(optional). If set to TRUE then a separate network file is created for each condition. Otherwise a new network containing all edges that passed the threshold will be created.",
              metavar="logical"),
  make_option(c("--unique_filter"), type="character", default=NULL,
              help="(optional). If --save_condition_networks is provided then condition-specific networks can be filtered to include an edge if it is unique in the condition as a categorical \"label\" or to the condition \"class\" itself. For example, if a condition class is named \"Treatment\" and it has several labels: \"Heat\", \"Drought\" and \"Control\" you can include edges that are only specific to a single label by providing the value \"label\" to this argument. This means the edge that are signficant for another more than one label will be excluded. Alternatively, to keep edges that are specific to the any label in a single class, (e.g. \"Treatment\"), provide the value \"class\" to this arugment.  This means that an edge can be signficant for any of the labels in the class (e.g. \"Heat\", \"Drought\" and \"Control\") but cannot be associated with another condition class. This type of filter is useful for finding edges that are specific to a label or class and nothing else (at least in the context of the experiment).  If this value is not provided then no uniqueness filter is applied.",
              metavar="character"),
   make_option(c("--score_weight"), type="numeric", default=1,
              help="(optional). The weight to apply to the Similarity Score ranking.",
              metavar="numeric"),
   make_option(c("--pval_weight"), type="numeric", default=1,
              help="(optional). The weight to apply to the p-value ranking.",
              metavar="numeric"),
   make_option(c("--rsqr_weight"), type="numeric", default=1,
              help="(optional). The weight to apply to the R-squared ranking.",
              metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$net)){
  print_help(opt_parser)
  stop("Please provide a network file (--net).", call.=FALSE)
}

if (is.null(opt$unique_filter)){
  opt$unique_filter = 'full'
}

if (is.null(opt$no_condition_rankings)){
  opt$no_condition_rankings = FALSE
}

if (is.null(opt$save_condition_networks)){
  opt$save_condition_networks = FALSE
}

# Make sure KINC.R is at the correct vresion.
if(packageVersion("KINC.R") > "1.1") {
    stop("This script requires KINC.R > 1.1")
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
message("Loading the KINC network file...")

net = loadKINCNetwork(opt$net)

# Get the edge ranks.
ranks = c()
if (opt$no_condition_rankings) {
    message("Calculating network-wide ranks...")
    ranks = getEdgeRanks(net, by_condition=FALSE, opt$score_weight, opt$pval_weight, opt$rsqr_weight)
} else {
    message("Calculating condition-specific ranks...")
    ranks = getEdgeRanks(net, by_condition=TRUE, opt$score_weight, opt$pval_weight, opt$rsqr_weight)
}
net$rank = ranks

# Filter the network to include the top n ranked edges.
netF = net[(which(ranks < opt$top_n)),]
netF = netF[order(netF$rank),]

if (opt$save_condition_networks) {
    message("Saving filtered condition-specific networks...")
    net_name = paste0(opt$out_prefix, '-th_ranked.%condition%-%filter%.csGCN.txt')
    saveConditionKINCNetwork(netF, net_name, filter=opt$unique_filter)
} else {
    message("Saving filtered network...")
    net_name = paste0(opt$out_prefix, '-th_ranked.csGCN.txt')
    saveKINCNetwork(netF, net_name)
}
