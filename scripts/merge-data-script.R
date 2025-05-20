
install.packages("tidyverse")
install.packages("openxlsx")

library(tidyverse)
library(openxlsx)

resfinder_results <- read.csv("test-data/resfinder-results/ResFinder_results_tab.txt", sep = "\t")
amrfinder_results <- read.csv("test-data/amrfinder-results/SRR292678.tsv", sep = "\t")
mlst_results <- read.table("test-data/mlst-results/SRR292678.tsv", sep = "\t", header = FALSE)

# extract relevant columns from resfinder results
resfinder_subset <- resfinder_results %>% select(Resistance.gene, Phenotype)

# subset relevant columns from amrfinder results
amrfinder_subset <- amrfinder_results %>% select(Element.symbol, Element.name, 
                                                 Type, Subtype, 
                                                 Class, Subclass)

# subset amrfinder by type column.
amrfinder_subset_AMR <- amrfinder_subset[amrfinder_subset$Type == "AMR", ]
amrfinder_subset_virulence <- amrfinder_subset[amrfinder_subset$Type == "VIRULENCE", ]
amrfinder_subset_stress <- amrfinder_subset[amrfinder_subset$Type == "STRESS", ]


# AMR genes detected by both amrfinder and resfinder
common_genes <- resfinder_subset[resfinder_subset$Resistance.gene %in% amrfinder_subset_AMR$Element.symbol, ]

# merge the resfinder results with amrfinder results
merged_data <- merge(amrfinder_subset_AMR, resfinder_subset,
                     by.x = "Element.symbol", by.y = "Resistance.gene",
                     all.x = TRUE, all.y = TRUE)


# rename columns for AMR, Virulence and Stress data from amrfinder
merged_data <- merged_data %>% rename(AMR_gene = Element.symbol,
                                      AMR_gene_description = Element.name,
                                      Resfinder_phenotype = Phenotype)


amrfinder_subset_virulence <- amrfinder_subset_virulence  %>% 
                                      rename(Virulence_gene = Element.symbol,
                                             Virulence_gene_description = Element.name)

amrfinder_subset_stress <- amrfinder_subset_stress  %>% 
                                      rename(resistance_gene = Element.symbol,
                                             resistance_gene_description = Element.name)


# write the results to an excel file.
wb <- createWorkbook()
addWorksheet(wb, "AMRFinder-Resfinder AMR")
addWorksheet(wb, "Virulence data")
addWorksheet(wb, "Stress data")

writeData(wb, "AMRFinder-Resfinder AMR", merged_data)
writeData(wb, "Virulence data", amrfinder_subset_virulence)
writeData(wb, "Stress data", amrfinder_subset_stress)

saveWorkbook(wb, "combined-amrfinder-resfinder-results.xlsx", overwrite = TRUE)









