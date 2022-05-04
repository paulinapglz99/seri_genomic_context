#script para filtrar vcf con anotaciones de mutaciones
#Paulina

#ESPECIFICAMENTE ESTE ES EL DE LOS SNPs

#Libraries

library("dplyr")
library("vroom")
library("stringr")

#En mi compu.

setwd("~/Tesis_INMEGEN/Resultados/vcf_Seri_splits/vcf_anotado/vcf_txt")

#Reading file

vep_snps_seri <- vroom(file = "vep_snps_seri.txt")

#Los 4 genomas Seri son 

genomas <- c("SM-3MG4D", "SM-3MG4E", "SM-3MG4F", "SM-3MG4G")

#Eliminar columnas que no sirven

snps_seri_filtrado <- vep_snps_seri %>% 
  select(Location, Allele, Consequence, IMPACT, 
         SYMBOL, Gene, Feature_type, Feature, 
         BIOTYPE, Existing_variation, 
         GIVEN_REF, CLIN_SIG)

#Guardo esta perra tabla por si acaso

write.csv(snps_seri_filtrado, 
          file = "snps_seri_filtrado.csv")

#Quiero tener una lista unica de los rsIDs que encontro la herramienta de anotacion

seri_snps_rsid <- snps_seri_filtrado %>% 
  filter(str_detect(Existing_variation, "rs*")) %>% #para que sean solo los que tienen registro rsID
  select(Existing_variation) %>% 
  distinct()   #para que se eliminen los que son iguales




#Son  rsids rescatados, los voy a guardar en una tabla

write.csv(seri_snps_rsid, 
          file = "rsid_Seri.csv")

#Antes de empezar a analizar los genes y rsIDs que extrajimos de la tabla, quiero hacer una tabla
#con los genes y variantes que si tienen información en la columna de CLIN_SIG

seri_snps_clinsig <- snps_seri_filtrado %>% 
  select(Location, Allele, Consequence, IMPACT, 
         SYMBOL, Gene, Existing_variation, CLIN_SIG) %>% 
  filter(!str_detect(CLIN_SIG, "-"))

#No se encontraron filas con clinical significance patogenico

#Busqueda de genes de interes

#Lista de genes ya analizados por PCR

genes_PCR <- c("CSF1PO", "FGA", "TH01", "TPOX", "VWA", "D3S1358",
               "D5S818", "D7S820", "D8S1179", "D13S317",
               "D16S539", "D18S51", "D21S1",  "D1S1656", "D2S441", "D2S1338",
               "D10S1248", "D12S391", "D19S433", "D22S1045", "SE33", 
               "TAP1", "TAP2", "DRB1", "DQA1", "DQB1", "AKT1", "GCKR", "SOCS3", 
               "PNPLA3", "CG", "CYP2C9","D9S1120", "MTHFR", "CYP1A2", "CYP2C9",
               "CYP2C19", "CYP2D6", "CYP3A4", "CYP2D6", 
               "CSF1PO", "FGA", "TH01", "TPOX", "vWA", "D3S1358", "D5S818", 
               "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", 
               "ADRB2")

#Lista de genes ya analizados por microarreglos

genes_microarreglos <- c("SEC16B", "OLFM4", "FTO", "MC4R", 
                         "GNPDA2", "FAIM2", "FAM120AOS", "LMX1B",  
                         "HOXB5", "ADAM23", "ELP3", "RAB27B", "GPR61", 
                         "TNNI3K")

#Filtrando para variantes antes reportadas de artículos 

seri_snp_genesrepo_microarreglos <- snps_seri_filtrado %>% 
  filter(SYMBOL %in% genes_microarreglos)

seri_snp_genesrepo_pcr <- snps_seri_filtrado %>% 
  filter(SYMBOL %in% genes_PCR)

#Busqueda de rsIDs de interes

#Para covid. Quiero ver si hay de casualidad algún gen de covi. Estas variaciones las saque
#de https://www.nature.com/articles/s41586-021-03767-x#Abs1

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

#Filtrando para variantes y genes de covid en la lista rescatada del paper

seri_snp_genes_covid <- snps_seri_filtrado %>% 
  filter(Existing_variation %in% genes_covid)

seri_snp_rsID_covid <- snps_seri_filtrado %>% 
  filter(Existing_variation %in% variaciones_covid)

#Cargando algunas librerias de genes

libreria_cancer <- vroom(file = "OMIM-Gene-Map-Retrieval_Cancer.txt")

libreria_metabolicos <- read.csv(file = "OMIM-Gene-Map-Retrieval_Metabolic_Syndrome.csv")

libreria_autoinmune <- read.csv(file = "OMIM-Gene-Map-Retrieval_AutoImmune.csv")

libreria_depresion <- vroom(file = "OMIM-Gene-Map-Retrieval_Depression.csv")

libreria_renal <- vroom(file = "OMIM-Gene-Map-Retrieval_Renal.csv")
  
libreria_infecciones <- vroom(file = "OMIM-Gene-Map-Retrieval_Infection.csv")
  
#Funcion para filtrar tablas

#Solo quiero estas columnas

columnas <- c()



