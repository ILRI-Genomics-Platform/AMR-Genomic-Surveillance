
user <- Sys.getenv('USER')

# check version, install and load libraries
version <- as.numeric(paste0(version$major, '.',
                       strsplit(version$minor, "\\.")[[1]][1]))

# define library path for installing packages
userLibrary <- paste0("/home/", user, "/R/x86_64-pc-linux-gnu-library/", version)

# reorder the lipaths
.libPaths(c(userLibrary, .libPaths()))

repos='http://cran.us.r-project.org'


# function to check installation of packages
ipak <- function(pkg, cran=TRUE){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (cran){
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = repos, 
                       lib = userLibrary)
    sapply(pkg, require, character.only = TRUE)
  }
  else {
    if (length(new.pkg)) 
      BiocManager::install(new.pkg, dependencies = TRUE, lib = userLibrary)
    sapply(pkg, require, character.only = TRUE)
  }
}

pkgs <- c("pacman", "argparse", "ggplot2", "tidyverse", "ape", "openxlsx")
ipak(pkg = pkgs, cran=TRUE)

pkgs <- c("ComplexHeatmap", "ggtree", "treeio", "ggtreeExtra", "ggnewscale", "splitstackshape")
ipak(pkg = pkgs, cran=FALSE)


