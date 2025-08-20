#' @author John Juma

user <- Sys.getenv('USER')

# check version, install and load libraries
version <- as.numeric(paste0(version$major, '.',
                             strsplit(version$minor, "\\.")[[1]][1]))

# define library path for installing packages
userLibrary <- paste0("/home/", user, "/R/x86_64-pc-linux-gnu-library/", version)

# reorder the lipaths
.libPaths(c('/export/apps/R/4.4/training-libs', userLibrary, .libPaths()))

library(dplyr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)

# read tree
tree <- read.tree("~/tree.nwk")
ggtree(tree)


# read metadata
metadata <- read.csv("~/data.csv")
metadata <- metadata |>
  dplyr::select(id, everything()) |>
  dplyr::rename("taxa"="id",
         "K.type"="Best.matching.K.type..capsular.loci.",
         "K.type__color"="Best.matching.K.type..capsular.loci.__colour")
  
# match ids with tip labels and modify
metadata.match <- metadata[match(metadata$taxa, tree$tip.label),] |>
  dplyr::arrange(Country) |>
  dplyr::mutate(Country = case_when(
    Country == "USA" ~ "United States",
    TRUE ~ Country
  )) |>
  dplyr::mutate(KPC.like.variant = gsub("KPC-like KPC-", "", KPC.like.variant)) |>
  dplyr::mutate(KPC.like.variant = case_when(
    KPC.like.variant == "Not applicable" ~ "Absent",
    TRUE ~ KPC.like.variant
  )) |>
  dplyr::mutate(K.type = case_when(
    K.type == "KL156-D1" ~ "Other",
    K.type == "KL161" ~ "Other",
    K.type == "KL37" ~ "Other",
    K.type == "KL112" ~ "Other",
    K.type == "KL30" ~ "Other",
    K.type == "Unknown/failed" ~ "Other",
    K.type == "KL46" ~ "Other",
    K.type == "KL82" ~ "Other",
    K.type == "KL126" ~ "Other",
    K.type == "KL23" ~ "Other",
    K.type == "KL50" ~ "Other",
    K.type == "KL40" ~ "Other",
    K.type == "KL38" ~ "Other",
    TRUE ~ K.type
  ))
  


# extract ST information
clade258 <- metadata.match[metadata.match$ST %in% c("258", "258*"),]
clade512 <- metadata.match[metadata.match$ST %in% c("512", "512*"),]
clade868 <- metadata.match[metadata.match$ST == "868",]

metadata.match <- metadata.match |> mutate(
  clade = case_when(taxa %in% clade258$taxa ~ "ST258",
                    taxa %in% clade512$taxa ~ "ST512",
                    taxa %in% clade868$taxa ~ "ST868"
  )
)
# check if all sample ids/taxa match te tip labels in the tree
all(tree$tip.label == metadata.match$id)


# make dataframe for clade nodes
clades.df <- data.frame(clade=unique(metadata.match$clade), node=NA)

# find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  clades.df$node[i] <- MRCA(
    tree,
    metadata.match$taxa[metadata.match$clade == clades.df$clade[i]]
  )
}


# get node ids

# plot preliminary tree
gg.tree <- ggtree(tree, layout = "dendrogram") %<+% metadata.match +
  geom_highlight(data = clades.df,
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) +
  geom_tree(linewidth=0.5) +
  geom_tippoint(aes(color=Country), size=5, alpha=1.0) +
  scale_color_manual(values = metadata.match$Country__colour,
                     breaks = metadata.match$Country)


# order clades dataframe to match the tree
clades.df <- clades.df[match(gg.tree$data %>%
                               dplyr::filter(isTip == "TRUE") %>%
                               arrange(y) %>%
                               pull(clade) %>%
                               unique(),
                             clades.df$clade),]
# add a column with alternating binary value
clades.df$highlight <- rep(c(0,1), length.out=length(clades.df$clade))
clades.df$highlight <- as.factor(clades.df$highlight)
print(clades.df)

# 
KPC_colors_df <- metadata.match |> 
  dplyr::select(c(KPC.like.variant, KPC.like.variant__colour))

K_type_colors_df <- metadata.match |>
  dplyr::select(c(K.type, K.type__color))



p <- ggtree(tree, layout = "rectangular") %<+% metadata.match +
  geom_tree(linewidth=0.5) +
  geom_highlight(data = clades.df,
                 aes(node=node, fill=clade),
                 alpha = 0.1,
                 align = "right",
                 extend = 0.5,
                 show.legend = FALSE) +
  geom_cladelabel(node=652, label="ST258", color="lightblue3", 
                  barsize = 2, align=TRUE, offset=0.2, fontsize = 4.0) +
  geom_cladelabel(node=947, label="ST512", color="lightpink3", 
                  barsize = 2, align=TRUE, offset=0.2, fontsize = 4.0) +
  geom_cladelabel(node=1009, label="ST869", color="mediumpurple2", 
                  barsize = 2, align=TRUE, offset=0.2, fontsize = 4.0) +
  geom_tippoint(aes(color=Country), size=5, alpha=1.0) +
  scale_color_manual(values = metadata.match$Country__colour,
                     breaks = metadata.match$Country) +
  scale_fill_manual(values = c("lightblue3", "lightpink3", "mediumpurple2"),
                    breaks = c("ST258", "ST512", "ST868")) +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping=aes(y = taxa, fill = KPC.like.variant),
    width = 8,
    alpha=0.9,
    size=0,
    offset=0.12,
    pwidth=0.2,
    position=position_identityx()) +
  scale_fill_manual(name = "KPC like variant",
                    values = c("#d95f02", "#7570b3", "#1b9e77", "#2c7bb6", "#4d4d4d"),
                    breaks = c("2", "3", "12", "25", "Absent"),
                    labels = c("2", "3", "12", "25", "Absent"),
                    guide=guide_legend(nrow = 3,
                                       override.aes=list(size=8))) +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping=aes(y = taxa, fill = K.type),
    width = 8,
    alpha=0.9,
    size=0,
    offset=0.12,
    pwidth=0.2,
    position=position_identityx()) +
  scale_fill_manual(name = "K-type",
                    values = c("#a6cee3", "#e31a1c", "#ffff33"),
                    breaks = c("KL106", "KL107", "Other"),
                    labels = c("KL106", "KL107", "Other"),
                    guide=guide_legend(nrow = 2,
                                       override.aes=list(size=8))) +

  theme_tree() +
  theme(text = element_text(family = "Arial Narrow"),
        plot.tag = element_text(size = 30),
        legend.title = element_text(size = 24),
        legend.background = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.text =  element_text(size = 15, color="#000000"),
        plot.subtitle=element_text(size = 25, color="#000000", face = "bold")) + 
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 8)))
p
