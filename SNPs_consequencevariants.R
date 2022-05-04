#script que analiza un archivo con missense variants de SNPs en Seris.
#REcuerda que en el anÃ¡lisis de la anotaciÃ³n de los Seri, te diste cuenta que el 
#65% de los SNPs estaban anotados como missense variants
#Paulina

#Libraries

library("dplyr")
library("vroom")
library("stringr")


#Leer archivo

snps_missense_variant_1.df <- vroom(file = "snp_Consequence_is_missense_variant.txt")


#Quiero hacer una lista de los genes que se muestran en la tabla

snps_missense_list <- snps_missense_variant_1.df %>% 
  select(SYMBOL) %>% 
  distinct()


#Quiero una tabla que me muestre el gen, e información de posición

mapa_variantes <- snps_missense_variant_1.df %>% 
  select(SYMBOL, Location, Allele, Existing_variation, Codons, SOURCE,
         AF) %>% 
  distinct()

#Descargando la lista

write.table(snps_missense_list, 
            file = "missense_list_snps.tsv", 
            sep = "\t")

#Descargando el mapeo de las variantes

write.table(mapa_variantes, 
            file = "map_missense_snps.tsv", 
            sep = "\t")

results <- vroom("map_missense_snps.tsv")

