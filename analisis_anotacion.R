#script para filtrar vcf con anotaciones de mutaciones
#Paulina

#Libraries

library("dplyr")
library("vroom")
library("stringr")

#En mi compu.

setwd("~/Tesis_INMEGEN/Resultados/analisis_data")

#Reading file

seri_anotado <- vroom(file = "seri_anotado.txt")

#Quiero unicamente dejar SNPs que se expresan, so vamos a filtrar

#Delimito los tipos de mutacion que quiero, en este caso, solo los que impliquen
#un alelo de nucleotido unico, o sea

SNPs <- c("A", "G", "C", "T")

#FIltro para que solo sean SNPs y que codifiquen a proteinas

seri_filtrado <- seri_anotado %>% 
  filter(Allele == SNPs) %>% 
  filter(BIOTYPE == "protein_coding")

#Guardo esta perra tabla por si acaso

write.csv(seri_filtrado, 
          file = "variaciones_seri_SNPS_exon.csv", 
          sep = ",")

#Quiero tener una lista unica de los rsIDs que encontro la herramienta de anotacion

seri_rsid <- seri_filtrado %>% 
  filter(str_detect(Existing_variation, "rs*")) %>% #para que sean solo los que tienen registro rsID
  select(Existing_variation) %>% 
 distinct()   #para que se eliminen los repetidos

#Son 1133 rsids rescatados, los voy a guardar en una tabla

write.csv(seri_rsid, 
          file = "rsid_Seri.csv", 
          sep = ",")

##Ahora quiero ver los farmacogenes

seri_citocromoCYP <- seri_filtrado %>% 
  filter(str_detect(SYMBOL, "CYP*")) ##Esto quiere decir que no hay genes de citocromos


#Guardo una tabla de citocromos

write.csv(seri_citocromoCYP2C9, 
          file = "citocromos_seri.csv", 
          sep = ",")

#Quiero ver si hay de casualidad algún gen de covi. Estas variaciones las saque
#de https://www.nature.com/articles/s41586-021-03767-x#Abs1

variaciones_covid <- c("rs1886814",
"rs2271616",
"rs10490770",
"rs11919389",
"rs72711165",
"rs912805253",
"rs10774671",
"rs1819040",
"rs77534576",
"rs2109069", 
"rs74956615",
"rs4801778",
"rs13050728")

genes_covid <- c("SLC6A20", "ABO", "OAS1", "OAS3",
                "ARHGAP27", "KANSL1", "MAPT", "SPPL2C", 
                "STH", "TYK2", "PPP1R15A", 
                "LZTFL1","ZBTB11", "RPL24", 
                "CEP97", "NXPE3",
                 "FOXP4","TMEM65",
                 "OAS2", "PLEKHM1", 
                "LINC02210", "CRHR1",
                "LRRC37A", 
                 "ARL17B", "LRRC37A2", 
                "ARL17A", "NSF",
                "WNT3",
                 "KAT7", "TAC4",
                 "DPP9","ICAM1", 
                "ICAM4", "ICAM5", 
                "ZGLP1", "FDX2",
                 "RAVER1", "ICAM3",
                "TYK2",
                 "PLEKHA4", "PPP1R15A", 
                "TULP2", "NUCB1",
                 "IFNAR2")

#Filtro de covid

seri_variantescovid <- seri_filtrado %>% 
  filter(Existing_variation %in% variaciones_covid)

seri_genescovid <- seri_filtrado %>% 
  filter(SYMBOL %in% genes_covid)

#Filtro segun algunos genes de interes medico

genes_metabolicos <- c(".")
genes_inmunes <- c(".")
genes_cancer <- c(".")
genes_renales <- c(".")
genes_infecciones <- c(".")