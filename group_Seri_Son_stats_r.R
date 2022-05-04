#Script que analiza dos tablas con los resultados de la anotación de un vfc 

library("vroom")
library("dplyr")
library("tidyr")
library("ggplot2")

#Reading table of vcf stats

stats_vcf <- vroom("group_Seri_Son_stats.csv")

#Reading table with Ensembl annotation stats

stats_vcf_ensembl <- vroom("group_Seri_Son_stats_ensembl.csv")

#Reading table with annotated vcf stats

stats_vcf_annotate <- vroom("group_Seri_Son_stats_annotated.csv")

#Merge tables by hand

VCF <- c("VCF annotate stats", "Original VCF stats")
samples <- c("4", "4")
records <- c("10771", "10817")
no_Alts <- c("0", "0")
SNPs <- c("3428", "3451")
indels <- c("6917", "6938")
MNPs <- c("0", "0")
others <- c("0", "0")
multiallelic_sites <- c("0", "0")
multiallelic_SNP_sites <- c("0", "0")
data_stats <- c( VCF, 
                 samples, 
                 records, 
                 no_Alts, 
                 SNPs, 
                 indels, 
                 MNPs,
                 others,
                 multiallelic_sites, 
                 multiallelic_SNP_sites)

data_stats2 <- VCF

stats_original_vs_annotate <- data.frame( VCF, 
                                          samples, 
                                          records, 
                                          no_Alts, 
                                          SNPs, 
                                          indels, 
                                          MNPs,
                                          others,
                                          multiallelic_sites, 
                                          multiallelic_SNP_sites
                                          )

stats_original_vs_annotate2 <- stats_original_vs_annotate %>% 
  pivot_longer(cols = data_stats )

#Plotting the difference between both stats

stats <- ggplot(data = stats_original_vs_annotate, 
                mapping = aes( x = data_stats,  
                               y = VCF, 
                               color = VCF)) +
 geom_point() +
 geom_line()

stats
