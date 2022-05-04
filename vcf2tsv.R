## This script converts a vcf file to a tab delimited file

## Charging library

library("vroom")
library("dplyr")
library("stringr")

## Up file
vcf_file <- vroom(file = "group_Seri_Son.vcf.test") #This file has the INFO column merged

## Break column. 
info_columns <- str_split_fixed(vcf_file$INFO, ";", 16) #Here we separate INFO

##changing from matrix to dataframe format 
info_columns_2 <- as.data.frame(info_columns) 


## Merging column of interest to original data 
merged_vcf <- vcf_file %>% 
  mutate(INFO = info_columns_2)

#Saving in tsv format

write.table(final_vcf, file = "Seri_Son.tsv", 
            sep = "\t", row.names = F)
