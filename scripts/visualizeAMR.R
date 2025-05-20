#!/usr/bin/env Rscript

#' @author John Juma

user <- Sys.getenv('USER')

# check version, install and load libraries
version <- as.numeric(paste0(version$major, '.', 
                       strsplit(version$minor, "\\.")[[1]][1]))
userLibrary <- paste0("/home/", user, "/R/x86_64-pc-linux-gnu-library/", version)

.libPaths(c(userLibrary, .libPaths()))

repos='http://cran.us.r-project.org'

# function to check installation of packages
ipak <- function(pkg, cran=TRUE){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (cran){
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = repos)
    sapply(pkg, require, character.only = TRUE)
  }
  else {
    if (length(new.pkg)) 
      BiocManager::install(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
}

pkgs <- c("pacman", "argparse", "colorspace", "ggplot2", "tidyverse",
          "ape", "ggnewscale", "scales", "splitstackshape")
ipak(pkg = pkgs, cran=TRUE)

pkgs <- c("ComplexHeatmap", "treeio", "ggtree", "ape", "ggtreeExtra")
ipak(pkg = pkgs, cran=FALSE)

# pacman::p_load(
#   argparse,        # parse command line options
#   colorspace,      # color palettes
#   ggplot2,         # to plot
#   tidyverse,       # general data management and visualization
#   ape,             # to import and export phylogenetic files
#   ggtree,          # to visualize phylogenetic files
#   treeio,          # to visualize phylogenetic files
#   ggtreeExtra,     # to visualize phylogenetic files
#   ggnewscale,      # to add additional layers of color schemes
#   scales,          # functions for visualization
#   splitstackshape  # stack and reshape dataframe
# )


usage <- function() {
  usage.text <- '\nUsage: visualizeAMR.R 
    --tree <path to the core SNP newick file>
    --mlst <paths to MLST reports (.tsv)> 
    --pointfinder <paths to PointFinder results (PointFinder_results.txt) files> 
    --resfinder <paths to ResFinder results (ResFinder_results_tab.txt) files> 
    --prefix <prefix to the output filename> 
    --outdir <path to the output directory name>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("--tree", 
                    default=NULL, 
                    help="paths to tree file in newick format")
parser$add_argument("--mlst", 
                    nargs="+", 
                    default=NULL, 
                    help="paths to MLST reports (.tsv)")
parser$add_argument("--pointfinder", 
                    nargs="+", 
                    default=NULL, 
                    help="paths to PointFinder results (PointFinder_results.txt) files")
parser$add_argument("--resfinder", 
                    nargs="+", 
                    default=NULL, 
                    help="paths to ResFinder results (ResFinder_results_tab.txt) files")
parser$add_argument("--prefix", 
                    default="AMR-plots", 
                    help="prefix to the output filename")
parser$add_argument("--outdir", 
                    default="./", 
                    help="output directory path")
args <- parser$parse_args()




###########################################
#######           checks           #######
###########################################

mlstFiles <- unique(unlist(strsplit(args$mlst," ")))
if (length(mlstFiles) == 0) {
  parser$print_help()
  stop("At least one MLST results file must be supplied", call.=FALSE)
}
if (!all(file.exists(mlstFiles))) {
  parser$print_help()
  stop(paste("The following mlst files don't exist:", 
             paste(mlstFiles[!file.exists(mlstFiles)], 
                   sep='', collapse=' '), sep=' '), call.=FALSE)
}


PointFinderFiles <- unique(unlist(strsplit(args$pointfinder," ")))
if (length(PointFinderFiles) == 0) {
  parser$print_help()
  stop("At least one PointFinder results file (PointFinder_results.txt) must be supplied", call.=FALSE)
}
if (!all(file.exists(PointFinderFiles))) {
  parser$print_help()
  stop(paste("The following PointFinder results files don't exist:", 
             paste(PointFinderFiles[!file.exists(PointFinderFiles)], 
                   sep='', collapse=' '), sep=' '), call.=FALSE)
}

ResFinderFiles <- unique(unlist(strsplit(args$resfinder," ")))
if (length(ResFinderFiles) == 0) {
  parser$print_help()
  stop("At least one ResFinder results file (ResFinder_results_tab.txt) must be supplied", call.=FALSE)
}
if (!all(file.exists(ResFinderFiles))) {
  parser$print_help()
  stop(paste("The following mlst files don't exist:", 
             paste(ResFinderFiles[!file.exists(ResFinderFiles)], 
                   sep='', collapse=' '), sep=' '), call.=FALSE)
}

outdir <- args$outdir
if (tail(strsplit(outdir,"")[[1]],1)!="/") {
  outdir <- paste(outdir,"/",sep='')
}
if (!file.exists(outdir)) {
  dir.create(outdir, recursive=TRUE)
}

prefix <- args$prefix

treeFile <- args$tree
if (!file.exists(treeFile)) {
  parser$print_help()
  stop("You must provide a tree file", call.=FALSE)
}
 

# baseDir <- file.path("~/trainings/ACDC_AMR2025/results/ont/klebsiella")


################################################################################
##########################      read MLST files       ##########################
################################################################################
# mlstFiles <- file.path(file.path(baseDir, "mlst"),
#                        list.files(path = file.path(baseDir, "mlst"),
#                                   recursive = TRUE,
#                                   pattern = "*.tsv"))


mlstList <- list()
for (i in 1:length(mlstFiles)){
  fn <- mlstFiles[i]
  # sample <- strsplit(basename(fn), ".", fixed=T)[[1]][1]
  df <- read.table(fn, 
                   check.names=F, 
                   sep='\t', header=F, 
                   col.names = c("path", "organism", "ST", 
                                 "gapA", "infB", "mdh", 
                                 "pgi", "phoE", "rpoB", 
                                 "tonB"))
  df$path <- sapply(strsplit(as.character(df$path), 
                             ".", fixed = TRUE), `[`, 2)
  df$Sample <- basename(df$path)
  df <- df |> dplyr::select(c(11,3,4,5,6,7,8,9,10))
  
  mlstList[[i]] <- df
}
mlst <- do.call(rbind, mlstList)

print(mlst)

################################################################################
############################ PointFinder results ###############################
################################################################################

# PointFinderFiles <- file.path(file.path(baseDir, "resfinder"),
#                             list.files(path = file.path(baseDir, "resfinder"),
#                                        recursive = TRUE,
#                                        pattern = "PointFinder_results.txt"))
PointFinderList <- list()
for (i in 1:length(PointFinderFiles)){
  fn <- PointFinderFiles[i]
  sample <- basename(dirname(fn))
  data <- data.table::fread(fn)
  data$Sample <- sample
  PointFinderList[[i]] <- data
}
PointFinderData <- do.call(rbind, PointFinderList)

print(PointFinderData)

PointFinderData$Gene <-
  sapply(strsplit(as.character(PointFinderData$Mutation),
                  " ", fixed = TRUE), `[`, 1)

PointFinderData$AA_change <-
  sapply(strsplit(as.character(PointFinderData$Mutation),
                  " ", fixed = TRUE), `[`, 2)


# summarize
PointMutations <- PointFinderData |>
  dplyr::group_by(Sample, Resistance, Gene) |>
  summarize(mutations = paste(sort(unique(Mutation)), collapse=", "))


# separate point mutations as columns
PointMutationsSplit <- cSplit_e(PointMutations, 'mutations',
                           sep= ',', type = 'character',
                           fill = 0, drop = TRUE)
names(PointMutationsSplit) <- sub('mutations_', '', names(PointMutationsSplit))

# separate point mutations as stacks
PointMutationsLong <- cSplit_l(PointMutations, 'mutations', sep= ',')

PointMutationsLongUnnest <- PointMutationsLong |> unnest(c(mutations_list))

# get unique genes
pointMutationsGenes <- sort(unique(PointFinderData$Gene))
pointMutationsResistance <- sort(unique(PointFinderData$Resistance))


pointMutationsDFs <- NULL
for (mut in pointMutationsGenes){
  for (res in pointMutationsResistance){
    df <- PointMutationsLongUnnest[
      (PointMutationsLongUnnest$Gene == mut) &
        (PointMutationsLongUnnest$Resistance == res),]

    col_names <- c("Sample", "Resistance", "Gene",
                   paste0(mut, "_", res, "_mutations"),
                   paste0(mut,"_", res, "_mutations_list")
                   )
    df <- as.data.frame(df)
    names(df) <- col_names
    if (dim(df)[1] != 0){
      df <- df |> dplyr::select(c(1,4,5))
      assign(paste0(mut,"_", res, "_mutations"), df)
      var <- paste0(mut,"_", res, "_mutations")
      pointMutationsDFs[[var]] <- df
    }
  }
}


################################################################################
############################## ResFinder results ###############################
################################################################################

# ResFinderFiles <- file.path(file.path(baseDir, "resfinder"),
#                             list.files(path = file.path(baseDir, "resfinder"),
#                                        recursive = TRUE,
#                                        pattern = "ResFinder_results_tab.txt"))
ResFinderList <- list()
for (i in 1:length(ResFinderFiles)){
  fn <- ResFinderFiles[i]
  sample <- basename(dirname(fn))
  data <- data.table::fread(fn)
  data$Sample <- sample

  # select AMR genes with 100.00 Identity and 100% coverage
  data <- data[(data$Identity == 100.00) & (data$Coverage == 100.00000),]

  # rename the columns and select columns
  data <- data |>
    dplyr::rename(
      "AMR_Gene" = "Resistance gene",
      "Length" = "Alignment Length/Gene Length"
  )
  data <- data |>
    dplyr::select(Sample, AMR_Gene, Phenotype)
  ResFinderList[[i]] <- data
}
ResFinderData <- do.call(rbind, ResFinderList)

# summarize 
ResistanceGenes <- ResFinderData |>
  dplyr::group_by(Sample) |>
  summarize(genes = paste(sort(unique(AMR_Gene)), collapse=", "))


# separate AMR genes as columns
ResistanceGenesSplit <- cSplit_e(ResistanceGenes, 'genes',
                            sep= ',', type = 'character',
                            fill = 0, drop = TRUE)
names(ResistanceGenesSplit) <- sub('genes_', '', names(ResistanceGenesSplit))

print(ResFinderData)


# separate AMR genes as stacks
ResistanceGenesLong <- cSplit_l(ResistanceGenes, 'genes', sep= ',')

ResistanceGenesLongUnnest <- ResistanceGenesLong |> unnest(c(genes_list))
print(ResistanceGenesLongUnnest)




# merge dataframes: mlst, point mutations and amr genes
sample_data <- c(list(mlst[,c(1,2)]), pointMutationsDFs)
SampleData <- sample_data %>% reduce(left_join, by='Sample')



################################################################################
################################## Tree file ###################################
################################################################################

# treeFile <- file.path(baseDir, "iqtree/core-snp.treefile")

tree <- read.tree(treeFile)

# match tip labels with sample names
all(tree$tip.label == SampleData$Sample)
SampleData <- SampleData[match(tree$tip.label, SampleData$Sample),]
all(tree$tip.label == SampleData$Sample)



# create dataframe for sequence types
sequence_type <- data.frame("ST" = SampleData[,c("ST")])
rownames(sequence_type) <- SampleData$Sample
sequence_type$ST <- as.factor(sequence_type$ST)

# create a dataframe for mutations in the acrR gene,
# which confer Fluoroquinolone resistance:
acrR_fluoroquinolone <- data.frame(
  "acrR_fluoroquinolone" = SampleData[,c("acrR_Fluoroquinolone_mutations")])
rownames(acrR_fluoroquinolone) <- SampleData$Sample


# create a dataframe for mutations in the ompK36 gene,
# which confer Carbapenem resistance:
ompK36_carbapenem <- data.frame(
  "ompK36_carbapenem" = SampleData[,c("ompK36_Carbapenem_mutations")])
rownames(ompK36_carbapenem) <- SampleData$Sample


# create a dataframe for mutations in the ompK36 gene,
# which confer Cephalosporin resistance:
ompK36_cephalosporin <- data.frame(
  "ompK36_cephalosporin" = SampleData[,c("ompK36_Cephalosporins_mutations")])
rownames(ompK36_cephalosporin) <- SampleData$Sample


# create a dataframe for mutations in the ompK37 gene,
# which confer Carbapenem resistance:
ompK37_carbapenem <- data.frame(
  "ompK37_carbapenem" = SampleData[,c("ompK37_Carbapenem_mutations")])
rownames(ompK37_carbapenem) <- SampleData$Sample


p <- ggtree(tree) %<+% SampleData +
  geom_tiplab(size = 5,
              linesize = .05,
              geom = "text",
              linetype = "dashed",
              alpha = 1,
              fontface = 2) + # labels the tips
  geom_tippoint(
    alpha = 1.0,
    size = 4,
    stroke=0.1) +
  theme_tree2()+
  xlab("genetic distance")+
  xlim(0, 0.015)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,
                                  face = "bold",
                                  hjust = 0.5,
                                  vjust = -15))
p
p <- p +
  ggnewscale::new_scale_fill()


h1 <-  gheatmap(p, sequence_type,
                offset = 0.003,
                width = 0.1,
                color="black",
                colnames = FALSE)+
  scale_fill_manual(name = "Sequence Type",
                    values = c('skyblue2', 'yellow2', 'brown4',
                               'lightblue1', 'navajowhite2',
                               'magenta3', 'purple3',
                               'green3', 'hotpink3'),
                    breaks = sort(unique(SampleData$ST)),
                    labels = sort(unique(SampleData$ST))) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical",
        legend.margin = margin())
h1

h2 <- h1 + ggnewscale::new_scale_fill()
h3 <- gheatmap(h2, acrR_fluoroquinolone,
               offset = 0.004,
               width = 0.1,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Fluoroquinolone resistance \n conferring mutation",
                    values = c("#f39544"),
                    breaks = unique(acrR_fluoroquinolone$acrR_fluoroquinolone),
                    labels = gsub("acrR ", "", unique(acrR_fluoroquinolone$acrR_fluoroquinolone)))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
h3

h4 <- h3 + ggnewscale::new_scale_fill()
h5 <- gheatmap(h4, ompK36_carbapenem,
               offset = 0.005,
               width = 0.1,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Carbapenem resistance \n conferring mutation - ompK36",
                    values = c("#fe9698", "#ea0c92"),
                    breaks = c( "ompK36 p.A217S", "ompK36 p.A217S, ompK36 p.N218H"),
                    labels = c( "p.A217S", "p.A217S, p.N218H"))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
h5


h6 <- h5 + ggnewscale::new_scale_fill()
h7 <- gheatmap(h6, ompK37_carbapenem,
               offset = 0.006,
               width = 0.1,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Carbapenem resistance \n conferring mutation - ompK37",
                    values = c("blue2", "forestgreen"),
                    breaks = unique(ompK37_carbapenem$ompK37_carbapenem),
                    labels = gsub("ompK37 ", "", unique(ompK37_carbapenem$ompK37_carbapenem)))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
h7

h8 <- h7 + ggnewscale::new_scale_fill()
h9 <- gheatmap(h8, ompK36_cephalosporin,
               offset = 0.007,
               width = 0.1,
               color = "black",
               colnames = FALSE)+
  scale_fill_manual(name = "Cephalosporin resistance \n conferring mutation",
                    values = c("brown3", "paleturquoise", "olivedrab3", "lemonchiffon4"),
                    breaks = unique(ompK36_cephalosporin$ompK36_cephalosporin),
                    labels = gsub("ompK36 ", "", unique(ompK36_cephalosporin$ompK36_cephalosporin)))+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())+
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))
h9


ggsave(file.path(outdir, paste0(prefix, ".MLST.PointMutations.pdf")),
       h9,
       bg="white",
       width = 15,
       height = 12,
       units = "in",
       limitsize = FALSE,
       dpi = 300,
       device = cairo_pdf)






# figuresDir <- file.path(baseDir, "figures")
# if (!dir.exists(figuresDir)) 
# dir.create(figuresDir, recursive = TRUE)
# 
# 
# ####################### load data ##########################
# 
# 
# 
# ############################################################
# ######### read sequence types outputs from mlst ############
# ############################################################
# mlstFile <- file.path(baseDir, "mlst/klebs-mlst.txt")
# mlstDF <- read.table(mlstFile, check.names=F, sep='\t', header=F,
#                      col.names = c("path", "organism", "ST", 
#                     "gapA", "infB", "mdh", 
#                     "pgi", "phoE", "rpoB", 
#                     "tonB"))
# mlstDF$path <- sapply(strsplit(as.character(mlstDF$path), ".", 
#                                fixed = TRUE), `[`, 2)
# mlstDF$Sample <- basename(mlstDF$path)
# mlstDF <- mlstDF |> dplyr::select(c(11,3,4,5,6,7,8,9,10))
# 
# 
# 
# # ref <- readLines("~/ResFinder_results_table.txt")
# # drugClasses <- c("Aminoglycoside", "Beta-lactam", "Colistin", 
# #                  "Fosfomycin", "Fusidic Acid",
# #                  "MLS - Macrolide, Lincosamide and Streptogramin B",
# #                  "Misc", "Nitroimidazole", "Oxazolidinone",
# #                  "Phenicol", "Fluoroquinolone", "Rifampicin",
# #                  "Sulphonamide", "Tetracycline", "Trimethoprim",
# #                  "Glycopeptide", "Pseudomonic Acid"
# #                  )
# 
# # dataList <- list()
# # for (drugClass in drugClasses){
# #   idx1 <- which(ref == drugClass)
# #   print(idx1)
# #   idx2 <- which(!nzchar(ref))
# #   print(idx2)
# #   idx3 <- idx2[idx2>idx1][1]-1
# #   print(idx3)
# #   df <- read.csv(text=gsub("\\.,", ",", sub("\\s*,$", "", ref[idx1:idx3])), 
# #                  header=FALSE, sep = "\t")
# # }
# 
# 

# # merge dataframes: mlst, point mutations and amr genes
# sample_data <- c(list(mlstDF[,c(1,2)]), pointMutationsDFs)
# SampleData <- sample_data %>% reduce(left_join, by='Sample')
# 
# 
# # ResFinder results
# ResFinderFiles <- file.path(file.path(baseDir, "resfinder"), 
#                             list.files(path = file.path(baseDir, "resfinder"),
#                                        recursive = TRUE, 
#                                        pattern = "ResFinder_results_tab.txt"))
# ResFinderList <- list()
# for (i in 1:length(ResFinderFiles)){
#   fn <- ResFinderFiles[i]
#   sample <- basename(dirname(fn))
#   data <- data.table::fread(fn)
#   data$Sample <- sample
#   
#   # select AMR genes with 100.00 Identity and 100% coverage
#   data <- data[(data$Identity == 100.00) & (data$Coverage == 100.00000),]
#   
#   # rename the columns and select columns
#   data <- data |> 
#     dplyr::rename(
#       "AMR_Gene" = "Resistance gene",
#       "Length" = "Alignment Length/Gene Length"
#   )
#   data <- data |> 
#     dplyr::select(Sample, AMR_Gene, Phenotype)
#   ResFinderList[[i]] <- data
# }
# ResFinderData <- do.call(rbind, ResFinderList)
# 
# ResistanceGenes <- ResFinderData |>
#   dplyr::group_by(Sample) |>
#   summarize(genes = paste(sort(unique(AMR_Gene)), collapse=", "))
# 
# ResistanceGenesSplit <- cSplit_e(ResistanceGenes, 'genes',
#                             sep= ',', type = 'character',
#                             fill = 0, drop = TRUE)
# 
# names(ResistanceGenesSplit) <- sub('genes_', '', names(ResistanceGenesSplit))
# 
# 
# # merge dataframes (Sequence Types, Point Mutations conferring resistance and AMR genes)
# # SampleData <- dplyr::left_join(mlstDF[,c(1,2)], PointMutations, 
# #                        by=c("name"="Sample"))
# # SampleData <- dplyr::left_join(SampleData, ResistanceGenes,
# #                                by=c("name"="Sample"))
# 
# 
# 
# # plot tree
# treeFile <- file.path(baseDir, "iqtree/core-snp.treefile")
# tree <- read.tree(treeFile)
# 
# # match tip labels with sample names
# all(tree$tip.label == SampleData$Sample)
# SampleData <- SampleData[match(tree$tip.label, SampleData$Sample),]
# all(tree$tip.label == SampleData$Sample)
# 
# 
# 
# # create dataframe for sequence types
# sequence_type <- data.frame("ST" = SampleData[,c("ST")])
# rownames(sequence_type) <- SampleData$Sample
# sequence_type$ST <- as.factor(sequence_type$ST)
# 
# # create a dataframe for mutations in the acrR gene, 
# # which confer Fluoroquinolone resistance:
# acrR_fluoroquinolone <- data.frame(
#   "acrR_fluoroquinolone" = SampleData[,c("acrR_Fluoroquinolone_mutations")])
# rownames(acrR_fluoroquinolone) <- SampleData$Sample
# 
# 
# # create a dataframe for mutations in the ompK36 gene, 
# # which confer Carbapenem resistance:
# ompK36_carbapenem <- data.frame(
#   "ompK36_carbapenem" = SampleData[,c("ompK36_Carbapenem_mutations")])
# rownames(ompK36_carbapenem) <- SampleData$Sample
# 
# 
# # create a dataframe for mutations in the ompK36 gene, 
# # which confer Cephalosporin resistance:
# ompK36_cephalosporin <- data.frame(
#   "ompK36_cephalosporin" = SampleData[,c("ompK36_Cephalosporins_mutations")])
# rownames(ompK36_cephalosporin) <- SampleData$Sample
# 
# 
# # create a dataframe for mutations in the ompK37 gene, 
# # which confer Carbapenem resistance:
# ompK37_carbapenem <- data.frame(
#   "ompK37_carbapenem" = SampleData[,c("ompK37_Carbapenem_mutations")])
# rownames(ompK37_carbapenem) <- SampleData$Sample
# 
# 
# p <- ggtree(tree) %<+% SampleData +
#   geom_tiplab(size = 3,
#               linesize = .05,
#               geom = "text",
#               linetype = "dashed",
#               alpha = 1,
#               fontface = 2) + # labels the tips
#   theme_tree2()+
#   xlab("genetic distance")+
#   xlim(0, 0.015)+
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         plot.title = element_text(size = 12, 
#                                   face = "bold",
#                                   hjust = 0.5,
#                                   vjust = -15))
# p
# p <- p + 
#   ggnewscale::new_scale_fill()
# 
# 
# h1 <-  gheatmap(p, sequence_type, 
#                 offset = 0.003,
#                 width = 0.1, 
#                 color="black", 
#                 colnames = FALSE)+
#   scale_fill_manual(name = "Sequence Type",
#                     values = c('skyblue2', 'yellow2', 'brown4', 
#                                'lightblue1', 'navajowhite2', 
#                                'magenta3', 'purple3',
#                                'green3', 'hotpink3'),
#                     breaks = sort(unique(SampleData$ST)),
#                     labels = sort(unique(SampleData$ST))) +
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.box = "vertical", 
#         legend.margin = margin())
# h1
# 
# h2 <- h1 + ggnewscale::new_scale_fill()
# h3 <- gheatmap(h2, acrR_fluoroquinolone,
#                offset = 0.004,
#                width = 0.1,
#                color = "black",
#                colnames = FALSE)+
#   scale_fill_manual(name = "Fluoroquinolone resistance \n conferring mutation",
#                     values = c("#f39544"),
#                     breaks = unique(acrR_fluoroquinolone$acrR_fluoroquinolone),
#                     labels = gsub("acrR ", "", unique(acrR_fluoroquinolone$acrR_fluoroquinolone)))+
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.box = "vertical", legend.margin = margin())+
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE))
# h3
# 
# h4 <- h3 + ggnewscale::new_scale_fill()
# h5 <- gheatmap(h4, ompK36_carbapenem,   
#                offset = 0.005, 
#                width = 0.1,
#                color = "black",
#                colnames = FALSE)+
#   scale_fill_manual(name = "Carbapenem resistance \n conferring mutation - ompK36",
#                     values = c("#fe9698", "#ea0c92"),
#                     breaks = c( "ompK36 p.A217S", "ompK36 p.A217S, ompK36 p.N218H"),
#                     labels = c( "p.A217S", "p.A217S, p.N218H"))+
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.box = "vertical", legend.margin = margin())+
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE))
# h5
# 
# 
# h6 <- h5 + ggnewscale::new_scale_fill()
# h7 <- gheatmap(h6, ompK37_carbapenem,
#                offset = 0.006,
#                width = 0.1,
#                color = "black",
#                colnames = FALSE)+
#   scale_fill_manual(name = "Carbapenem resistance \n conferring mutation - ompK37",
#                     values = c("blue2", "forestgreen"),
#                     breaks = unique(ompK37_carbapenem$ompK37_carbapenem),
#                     labels = gsub("ompK37 ", "", unique(ompK37_carbapenem$ompK37_carbapenem)))+
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.box = "vertical", legend.margin = margin())+
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE))
# h7
# 
# h8 <- h7 + ggnewscale::new_scale_fill()
# h9 <- gheatmap(h8, ompK36_cephalosporin,
#                offset = 0.007,
#                width = 0.1,
#                color = "black",
#                colnames = FALSE)+
#   scale_fill_manual(name = "Cephalosporin resistance \n conferring mutation",
#                     values = c("brown3", "paleturquoise", "olivedrab3", "lemonchiffon4"),
#                     breaks = unique(ompK36_cephalosporin$ompK36_cephalosporin),
#                     labels = gsub("ompK36 ", "", unique(ompK36_cephalosporin$ompK36_cephalosporin)))+
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.box = "vertical", legend.margin = margin())+
#   guides(fill = guide_legend(nrow = 4, byrow = TRUE))
# h9
# 
# 
# 
# 
# 
# 
# 
# p <- ggtree(tree) + 
#   xlim(0, 0.010) + # to allow more space for labels
#   geom_treescale() # adds the scale
# p
# 
# 
# p <- p %<+% mlstDF + 
#   geom_tiplab(align = TRUE,
#               linesize = .05,
#               geom = "text",
#               linetype = "dashed",
#               size = 3.0,
#               alpha = 1,
#               fontface = 2) +
#   geom_tippoint(
#     aes(color = gapA),
#     alpha = 1.0,
#     size = 4,
#     stroke=0.1) +
#   scale_color_discrete_sequential(palette = "Batlow", rev = TRUE) +
#   theme(legend.position = c(0.5,0.2)) + # no keys
#   guides(color = guide_legend(override.aes = list(size=4),
#                               title = "Sample"))
# 
# 
# sequencetype <- mlstDF |> 
#   dplyr::select(c("name", "ST")) %>% 
#   remove_rownames %>% column_to_rownames(var="name")
# sequencetype$ST <- as.factor(sequencetype$ST)
# 
# 
# 
# p <- p + 
#   ggnewscale::new_scale_fill()
# 
# # serotype annotation
# p1 <- gheatmap(p + theme(legend.position = "none"), 
#                sequencetype,
#                offset=0.001*2, 
#                width=0.25,
#                colnames_angle = 0, 
#                colnames_position = "top",
#                font.size = 8,
#                family = "Helvetica",
#                legend_title="1. Sequence Type",
#                colnames_offset_y = 0.5,
#                custom_column_labels = "1") +
#   scale_fill_discrete_divergingx(palette = "Temps", rev=FALSE) +
#   guides(fill = guide_legend(title="1. Sequence Type",
#                              order = 1,
#                              nrow = 3,
#                              override.aes = list(size = 9)))
# p1
# 
# 
# 
# # 
# 
# 
# 

DataMatrixPointMutations <- dplyr::left_join(mlst[,c(1,2)],
                                             ResistanceGenesSplit,
                                             by="Sample")



DataMatrix <- as.data.frame(DataMatrixPointMutations) |>
  dplyr::select(-c(2)) %>%
  remove_rownames %>%
  column_to_rownames(var="Sample")



mat <- as.matrix(as.data.frame(DataMatrix))
mat <- t(mat)


annotation <-  DataMatrixPointMutations |> dplyr::select(ST)

colours <- list(
  ST = c('17'='skyblue2', 
         '38'='yellow2', 
         '39'='brown4',
         '45' = 'lightblue1', 
         '48'='navajowhite2',
         '54'='magenta3', 
         '258'='purple3',
         '1427'='green3', 
         '6155'='hotpink3')
  )

colAnno <- HeatmapAnnotation(df = annotation,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),
                            simple_anno_size = unit(1.4, "cm"),
                            annotation_legend_param = list(
                              ST = list(color_bar = "discrete",
                                                  title_gp=gpar(fontsize=14,
                                                                fontface="bold"),
                                                  labels_gp=gpar(fontsize=13),
                                                  direction="vertical"))
                            )


if (dim(mat)[1] < 4) {
  height = 0.1969*nrow(mat) + (2*1.5)
  width = 0.1969*ncol(mat) + (2*1.5)
} else {
  height = 0.1969*nrow(mat) + (3*1.5)
  width = 0.35*ncol(mat) + (3.12*1.5)
}

cols <- colorRampPalette(c('#E41A1C', '#377EB8'))(2)

# heatmap
pdf(file.path(outdir, paste0(prefix, '-genes-heatmap.pdf')),
    height = height, width = width,
    bg="white", family = "Helvetica")

Heatmap(mat,
        col = cols,
        name = 'Presence/Absence',
        cluster_rows = FALSE,
        #cluster_columns = TRUE,
        #show_column_dend = TRUE,
        #show_row_dend = TRUE,
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        rect_gp = gpar(col = "white", lwd = 1.0),
        heatmap_legend_param = list(title_gp=gpar(fontsize=14,
                                                  fontface="bold"),
                                    direction="vertical",
                                    heatmap_legend_side = "topright"),
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15),
        column_order = colnames(mat),
        # row_km = 2,
        # column_km = 3,
        #row_gap = unit(c(2, 4), "mm"),
        #column_gap = unit(c(2, 4), "mm"),
        top_annotation=colAnno
)
dev.off()