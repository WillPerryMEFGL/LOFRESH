library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(microDecon)
library(ggpubr)
library(effectsize)
library(see)
library(emmeans)
library(lsmeans)
library(gratia)
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")

#________________________________________________________

#Removing 12S sequences in the taxtable that are the wrong size,then BLASTing with a broad 12S database to remove
#sequences that are assigned to other things that are not fish

#________________________________________________________

taxTab_WP2_12S <- read_table2("taxTab_WP2_12S.txt")
taxTab_WP2_12S$Kingdom <-NULL
taxTab_WP2_12S$Phylum <-NULL
taxTab_WP2_12S$Class <-NULL
taxTab_WP2_12S$Order <-NULL
taxTab_WP2_12S$Family <-NULL
taxTab_WP2_12S$Genus <-NULL
taxTab_WP2_12S$Species <-NULL

#sequence length of the ASVs
taxTab_WP2_12S_length<-aggregate(Seq~Seq_length, transform(taxTab_WP2_12S, Seq_length=Seq),
                                        FUN=function(x) nchar(unique(x)))
hist(taxTab_WP2_12S_length$Seq)

#removing ASVs that are the wrong size, over 20% larger than the average amplicon size the MiFish amplicon should be (172bp)
taxTab_WP2_12S_length<-taxTab_WP2_12S_length[!(taxTab_WP2_12S_length$Seq>206.2),]

#removing ASVs that are the wrong size, under 20% smaller than the average amplicon size the MiFish amplicon should be (172bp)
taxTab_WP2_12S_length<-taxTab_WP2_12S_length[!(taxTab_WP2_12S_length$Seq<137.6),]
hist(taxTab_WP2_12S_length$Seq)

taxTab_WP2_12S_length$Seq_number <- 1:nrow(taxTab_WP2_12S_length) 
taxTab_WP2_12S_length$Seq<-"SEQ"
taxTab_WP2_12S_length$div<-">"

taxTab_WP2_12S_length$Seq_number <-paste(taxTab_WP2_12S_length$div,taxTab_WP2_12S_length$Seq,taxTab_WP2_12S_length$Seq_number)
taxTab_WP2_12S_length$Seq_number = as.character(gsub(" ", "", taxTab_WP2_12S_length$Seq_number))

taxTab_WP2_12S_length$Seq <-NULL
taxTab_WP2_12S_length$div <-NULL

write.table(taxTab_WP2_12S_length,"length_filtered_SEQ_to_sequence_lookup_table.txt",quote = FALSE)

taxTab_WP2_12S_length$Seq <-paste(taxTab_WP2_12S_length$Seq_number,taxTab_WP2_12S_length$Seq_length)

taxTab_WP2_12S_length$Seq <- as.character(gsub(" ", "\n", taxTab_WP2_12S_length$Seq))

length_adjust_Seq_number<-taxTab_WP2_12S_length$Seq_number
length_adjust_Seq_number<-data.frame(length_adjust_Seq_number)
length_adjust_Seq_number$Seq_length<-taxTab_WP2_12S_length$Seq_length

taxTab_WP2_12S_length$Seq_number<-NULL
taxTab_WP2_12S_length$Seq_length<-NULL

#produce the sequence length adjusted BLAST input file
write.csv(taxTab_WP2_12S_length,"length_filtered_taxTab_WP2_12S_BLAST_INPUT.fasta", row.names=FALSE,quote = FALSE)

#now BLAST this on the cluster using the "12S rRNA gene" database using -perc_identity 70 -qcov_hsp_perc 80 -evalue 1

#________________________________________________________

#Read in the length adjusted BLAST output file from the first BLAST of NCBI 
#you need objects created in the code before this (^) in order for it to work
#________________________________________________________

length_filtered_BLAST_OUTPUT_12S <- read_table2("length_filtered_BLAST_OUTPUT_12S.txt", 
                                                col_names = FALSE)

TAXONOMY_filtered_BLAST_OUTPUT_12S <- read_csv("TAXONOMY_filtered_BLAST_OUTPUT_12S.csv")
length_filtered_BLAST_OUTPUT_12S<-cbind(length_filtered_BLAST_OUTPUT_12S,TAXONOMY_filtered_BLAST_OUTPUT_12S)

length_adjust_Seq_number$length_adjust_Seq_number<-as.character(gsub(">", "", length_adjust_Seq_number$length_adjust_Seq_number))
length_adjust_Seq_number<-length_adjust_Seq_number %>% rename(Seq = "length_adjust_Seq_number")
length_filtered_BLAST_OUTPUT_12S<-length_filtered_BLAST_OUTPUT_12S %>% rename(Seq = "X1")

length_adjust_Seq_number <- length_adjust_Seq_number %>% full_join(length_filtered_BLAST_OUTPUT_12S, by = c("Seq"))

#remove ASVs that could not be assigned using broad BLAST
length_adjust_Seq_number<-na.omit(length_adjust_Seq_number)

#remove ASvs that were classified as non-fish
length_adjust_Seq_number_fish<-length_adjust_Seq_number %>% filter(class == "Actinopteri")
length_adjust_Seq_number_mammals<-length_adjust_Seq_number %>% filter(class == "Mammalia")

#Format for the second BLAST

BLAST_2_mito<-length_adjust_Seq_number_fish$Seq
BLAST_2_mito<-data.frame(BLAST_2_mito)
BLAST_2_mito$Seq<-length_adjust_Seq_number_fish$Seq_length

BLAST_2_mito$Seq_BLAST <-paste(">",BLAST_2_mito$BLAST_2_mito,BLAST_2_mito$Seq)
BLAST_2_mito$Seq_BLAST<- gsub("> ", ">", BLAST_2_mito$Seq_BLAST)
BLAST_2_mito$Seq_BLAST<- as.character(gsub(" ", "\n", BLAST_2_mito$Seq_BLAST))

BLAST_2_mito$Seq<-NULL
BLAST_2_mito$BLAST_2_mito<-NULL

#produce the sequence length adjusted BLAST input file
write.csv(BLAST_2_mito,"length_filtered_taxTab_WP2_12S_BLAST_INPUT_FISH_ONLY.fasta", row.names=FALSE,quote = FALSE)
remove(BLAST_2_mito, length_adjust_Seq_number)

#now BLAST this on the cluster using the mitogenome database using -perc_identity 90 -qcov_hsp_perc 90 -evalue 0.001

#________________________________________________________

#Read in fish only BLAST output using mitogenome database

#________________________________________________________


BLAST_OUTPUT_12S <- read_table2("BLAST_OUTPUT_FISH_ONLY_12S.txt", col_names = FALSE)
BLAST_OUTPUT_12S<-BLAST_OUTPUT_12S %>% rename(rn = X2 )

#Adding the taxonomy linked with the mitogenome Ids
mitogenome_annotations <- read_csv("mitogenome_annotations.csv")
BLAST_OUTPUT_12S$rn <- gsub("ref|", "", BLAST_OUTPUT_12S$rn)
BLAST_OUTPUT_12S$rn <- gsub("\\|", "", BLAST_OUTPUT_12S$rn)

BLAST_OUTPUT_12S <- BLAST_OUTPUT_12S %>% full_join(mitogenome_annotations, by = c("rn"))
BLAST_OUTPUT_12S<-na.omit(BLAST_OUTPUT_12S)

SEQ_lookup<-length_adjust_Seq_number_fish$Seq
SEQ_lookup<-data.frame(SEQ_lookup)
SEQ_lookup$Seq<-length_adjust_Seq_number_fish$Seq_length

SEQ_lookup<-SEQ_lookup %>% rename(X1 = SEQ_lookup )

BLAST_OUTPUT_12S <- BLAST_OUTPUT_12S %>% full_join(SEQ_lookup, by = c("X1"))

BLAST_OUTPUT_12S<-na.omit(BLAST_OUTPUT_12S)

write.table(BLAST_OUTPUT_12S,"BLAST_OUTPUT_12S_IDs_tax.txt",row.names = FALSE,col.names=FALSE,quote = FALSE)
remove(mitogenome_annotations,length_adjust_Seq_number_fish,length_filtered_BLAST_OUTPUT_12S,SEQ_lookup,
       TAXONOMY_filtered_BLAST_OUTPUT_12S,taxTab_WP2_12S,taxTab_WP2_12S_length)
gc()

seqtabNoC_WP2_12S <- read_table2("seqtabNoC_WP2_12S.txt")
seqtabNoC_WP2_12S_transpose<-t(seqtabNoC_WP2_12S)
seqtabNoC_WP2_12S_transpose<-as.data.frame(seqtabNoC_WP2_12S_transpose)
seqtabNoC_WP2_12S_transpose<-setDT(seqtabNoC_WP2_12S_transpose, keep.rownames = TRUE)
names(seqtabNoC_WP2_12S_transpose)[names(seqtabNoC_WP2_12S_transpose) == "rn"] <- "Seq"
seqtabNoC_WP2_12S_transpose[] <- lapply(seqtabNoC_WP2_12S_transpose, gsub, pattern='"', replacement='')
seqtabNoC_WP2_12S_transpose <- na.omit(transform(seqtabNoC_WP2_12S_transpose, Seq = c("Seq", Seq[-nrow(seqtabNoC_WP2_12S_transpose)])))
header.true <- function(seqtabNoC_WP2_12S_transpose) {
  names(seqtabNoC_WP2_12S_transpose) <- as.character(unlist(seqtabNoC_WP2_12S_transpose[1,]))
  seqtabNoC_WP2_12S_transpose[-1,]
}
seqtabNoC_WP2_12S_transpose<-header.true(seqtabNoC_WP2_12S_transpose)
remove(seqtabNoC_WP2_12S)

#REMOVING SAMPLES WITH LOW READ DEPTH
seqtabNoC_WP2_12S_transpose<-data.frame(seqtabNoC_WP2_12S_transpose)
i <- c(2:1224) 
seqtabNoC_WP2_12S_transpose[ , i] <- apply(seqtabNoC_WP2_12S_transpose[ , i], 2,  
                                           function(x) as.numeric(as.character(x)))

rarefy_col_sum<-colSums(seqtabNoC_WP2_12S_transpose[c(2:1224)])
rarefy_col_sum<-data.frame(rarefy_col_sum)
rarefy_col_sum<-setDT(rarefy_col_sum, keep.rownames = TRUE)[]
rarefy_col_sum<-data.frame(rarefy_col_sum)

rarefy_col_sum<-rarefy_col_sum[!(rarefy_col_sum$rarefy_col_sum>=1000),]
drop<-rarefy_col_sum$rn

seqtabNoC_WP2_12S_transpose<-seqtabNoC_WP2_12S_transpose[,!(names(seqtabNoC_WP2_12S_transpose) %in% drop)]

#Replace the literal sequences from the seqtabNoc table and replace them with the corresponding SEQ Ids

seqtabNoC_WP2_12S_tax <- BLAST_OUTPUT_12S %>% full_join(seqtabNoC_WP2_12S_transpose, by = c("Seq"))
seqtabNoC_WP2_12S_tax$rn <- NULL
seqtabNoC_WP2_12S_tax$X1 <- NULL
seqtabNoC_WP2_12S_tax$X3 <- NULL
seqtabNoC_WP2_12S_tax$X4 <- NULL
seqtabNoC_WP2_12S_tax$X5 <- NULL
seqtabNoC_WP2_12S_tax$Seq<-NULL
seqtabNoC_WP2_12S_tax<-na.omit(seqtabNoC_WP2_12S_tax)

i <- c(2:992) 
seqtabNoC_WP2_12S_tax[ , i] <- apply(seqtabNoC_WP2_12S_tax[ , i], 2,  
                                     function(x) as.numeric(as.character(x)))

#CALCULATING 0.05% READ THREASHOLD FOR LATER FILTEIRNG
seqtabNoC_WP2_12S_tax_col_sum<-colSums(seqtabNoC_WP2_12S_tax[c(2:992)])
seqtabNoC_WP2_12S_tax_col_sum<-data.frame(seqtabNoC_WP2_12S_tax_col_sum)
seqtabNoC_WP2_12S_tax_col_sum$read_filter<-seqtabNoC_WP2_12S_tax_col_sum$seqtabNoC_WP2_12S_tax_col_sum *0.0005
seqtabNoC_WP2_12S_tax_col_sum<- seqtabNoC_WP2_12S_tax_col_sum %>% rownames_to_column("ID")

#REMOVE ASVS THAT CANNOT BE IDENTIFIED WITH BLAST
seqtabNoC_WP2_12S_tax<-seqtabNoC_WP2_12S_tax %>% drop_na(species)
seqtabNoC_WP2_12S_tax$X1<-NULL

#CONVERITNG ASV TABLE INTO LONG FORMAT
seqtabNoC_WP2_12S_tax_long <- gather(seqtabNoC_WP2_12S_tax, Sample, reads, 2:992, factor_key=TRUE)
seqtabNoC_WP2_12S_tax_long$reads<-   as.numeric(seqtabNoC_WP2_12S_tax_long$reads)

#READ IN METADATA
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary<-table(lofresh_metadata_WP2_3$Days,lofresh_metadata_WP2_3$SampleSite_code)
lofresh_metadata_WP2_3_summary<-data.frame(lofresh_metadata_WP2_3_summary)
write.csv(lofresh_metadata_WP2_3_summary,"lofresh_metadata_WP2_3_summary.csv")
seqtabNoC_WP2_12S_tax_long<-seqtabNoC_WP2_12S_tax_long %>% rename(ID=Sample)

#ADD METADATA TO ASV TABLE
seqtabNoC_WP2_12S_tax_long_meta <- merge(lofresh_metadata_WP2_3[, c("ID", "WP")],seqtabNoC_WP2_12S_tax_long, by="ID")
sum_reads<- sum(seqtabNoC_WP2_12S_tax_long_meta$reads)

seqtabNoC_WP2_12S_tax_long_meta_clean_filtered <- merge(seqtabNoC_WP2_12S_tax_col_sum[, c("ID", "read_filter")],seqtabNoC_WP2_12S_tax_long_meta, by="ID")

#REMOVE READS LESS THAN THE 0.005% CUT OFF
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$reads<=seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$read_filter),]
#REMOVE READS LESS THAN 20
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered$reads<20),]

write.csv(seqtabNoC_WP2_12S_tax_long_meta_clean_filtered,"WP2_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2_12S_tax_long_meta_clean_filtered <- read_csv("WP2_ASV_table_long_filtered_family.csv")

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2_12S_tax_wide_filtered<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2_12S_tax_wide_filtered$Order <- NULL 
seqtabNoC_WP2_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2_12S_tax_wide_filtered$Seq_length <- NULL 
seqtabNoC_WP2_12S_tax_wide_filtered$WP <- NULL 

seqtabNoC_WP2_12S_tax_wide_filtered<-seqtabNoC_WP2_12S_tax_wide_filtered %>% 
  group_by(ID,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(ID=ID,species=species, reads= reads,)

seqtabNoC_WP2_12S_tax_wide_filtered <- spread(seqtabNoC_WP2_12S_tax_wide_filtered,ID,reads)
seqtabNoC_WP2_12S_tax_wide_filtered[is.na(seqtabNoC_WP2_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2_12S_tax_wide_filtered))

names(seqtabNoC_WP2_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2_12S_tax_wide_filtered <- seqtabNoC_WP2_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#NORMALISING READS

#YOU MAY NEED TO CHANGE NUMBERS HERE

seqtabNoC_WP2_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2_12S_tax_wide_filtered[c(1:68)])
seqtabNoC_WP2_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2_12S_tax_col_sum2)
seqtabNoC_WP2_12S_tax_col_sum2$ID <- rownames(seqtabNoC_WP2_12S_tax_col_sum2)

seqtabNoC_WP2_12S_tax_normalised<- merge(seqtabNoC_WP2_12S_tax_col_sum2[, c("ID", "seqtabNoC_WP2_12S_tax_col_sum2")],seqtabNoC_WP2_12S_tax_long_meta_clean_filtered, by="ID")

seqtabNoC_WP2_12S_tax_normalised <- transform(seqtabNoC_WP2_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2_12S_tax_normalised , aes(x = ID, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position = "none")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2_12S_tax_normalised_wide<-seqtabNoC_WP2_12S_tax_normalised
seqtabNoC_WP2_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$SampleSite_time <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 

seqtabNoC_WP2_12S_tax_normalised_wide<-seqtabNoC_WP2_12S_tax_normalised_wide %>% 
  group_by(ID,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(ID=ID,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2_12S_tax_normalised_wide <- spread(seqtabNoC_WP2_12S_tax_normalised_wide,ID,normalised_reads)
seqtabNoC_WP2_12S_tax_normalised_wide[is.na(seqtabNoC_WP2_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2_12S_tax_normalised_wide))

names(seqtabNoC_WP2_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2_12S_tax_normalised_wide <- seqtabNoC_WP2_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2_12S_tax_normalised_wide,"NORM_WP2_ASV_table_wide_filtered_family.csv")

### WP2A  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP2A_12S_tax_long_meta_clean<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered %>% filter(WP == "WP2A")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2A_12S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP2A_12S_tax_long_meta_clean,"WP2A_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2A_12S_tax_long_meta_clean <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

### WP2B  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP2B_12S_tax_long_meta_clean<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered %>% filter(WP == "WP2B")
seqtabNoC_WP2B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2B_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2B_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2B_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2B_12S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP2B_12S_tax_long_meta_clean,"WP2B_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2B_12S_tax_long_meta_clean <- read_csv("WP2B_ASV_table_long_filtered_family.csv")

### WP3B  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP3B_12S_tax_long_meta_clean<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered %>% filter(WP == "WP3B")
seqtabNoC_WP3B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP3B_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP3B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP3B_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP3B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP3B_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP3B_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP3B_12S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP3B_12S_tax_long_meta_clean,"WP3B_ASV_table_long_filtered_family.csv")
seqtabNoC_WP3B_12S_tax_long_meta_clean <- read_csv("WP3B_ASV_table_long_filtered_family.csv")

### WP2C Splitting work packages and adding metadata ###
seqtabNoC_WP2C_12S_tax_long_meta_clean<-seqtabNoC_WP2_12S_tax_long_meta_clean_filtered %>% filter(WP == "WP2C")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_12S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Region")],seqtabNoC_WP2C_12S_tax_long_meta_clean, by="ID")

seqtabNoC_WP2C_12S_tax_long_meta_clean_Gwash<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Leicestershire")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Towy<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Carmarthenshire")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Skan<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "New_York")
seqtabNoC_WP2C_12S_tax_long_meta_clean_Glatt<-seqtabNoC_WP2C_12S_tax_long_meta_clean %>% filter(Region == "Switzerland")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Gwash,"WP2C_Gwash_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Towy,"WP2C_Towy_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Skan,"WP2C_Skan_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_12S_tax_long_meta_clean_Glatt,"WP2C_Glatt_ASV_table_long_filtered_family.csv")

#________________________________________________________

########WP2A BLAST  TAXONOMY ANALYSIS######

#________________________________________________________
library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(microDecon)
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

#sewage_sps<-c("Gadus_macrocephalus",
##"Hippoglossus_stenolepis",
#"Ammodytes_personatus",
#"Clupea_harengus",
#"Clupea_pallasii",
#"Copadichromis_virginalis",
#"Favonigobius_gymnauchen",
#"Pholis_crassispina")

#SEWAGE<-filter(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, species %in% sewage_sps)

#aggregating replicates
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
 group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
group_by(SampleSite_time,species) %>% 
summarize(reads=mean(reads)) %>%
rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Gadus_macrocephalus"),]#Pacific cod
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Hippoglossus_stenolepis"),]#Pacific halibut
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Ammodytes_personatus"),]#Pacific Sand Lance 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Clupea_harengus"),]#Atlantic herring 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Clupea_pallasii"),]#Pacific herring 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Copadichromis_virginalis"),]#Haplochromine cichlid  
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Favonigobius_gymnauchen"),]#Sand goby   
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Pholis_crassispina"),]#Rock gunnel  
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Melanogrammus_aeglefinus"),]#Haddock
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]#Starry flounder
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Scomber_scombrus"),]#Atlantic mackerel - POSITIVE CONTROL
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salmo_obtusirostris','Salmo_trutta',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Anguilla_rostrata','Anguilla_anguilla',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salvelinus_fontinalis_x_Salvelinus_malma','Salvelinus_malma',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)

#Remove negative controls
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T10"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T17"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T18"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T19"),]


#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2A_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2A_12S_tax_wide_filtered <- spread(seqtabNoC_WP2A_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2A_12S_tax_wide_filtered[is.na(seqtabNoC_WP2A_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2A_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2A_12S_tax_wide_filtered))

names(seqtabNoC_WP2A_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2A_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2A_12S_tax_wide_filtered <- seqtabNoC_WP2A_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2A_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2A_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#rarefaction_curve<-seqtabNoC_WP2A_12S_tax_wide_filtered
#rownames(rarefaction_curve) <- NULL

#i <- c(1:16) 
#rarefaction_curve[ , i] <- apply(rarefaction_curve[ , i], 2,  
                                                  # function(x) as.integer(as.character(x)))
#(raremax <- min(rowSums(rarefaction_curve)))
#rarecurve(rarefaction_curve, col = c("red", "blue", "orange", "yellow", "green"), step = 20, sample = raremax, label = FALSE)

#FILTERING BY NEGATIVE CONTROL
#seqtabNoC_WP2A_12S_tax_wide_filtered<-t(seqtabNoC_WP2A_12S_tax_wide_filtered)
#seqtabNoC_WP2A_12S_tax_wide_filtered<-data.frame(seqtabNoC_WP2A_12S_tax_wide_filtered)
#seqtabNoC_WP2A_12S_tax_wide_filtered<-setDT(seqtabNoC_WP2A_12S_tax_wide_filtered, keep.rownames = TRUE)[]
#seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_wide_filtered %>% rename(OTU_ID = rn)

#colnames_decon<-colnames(seqtabNoC_WP2A_12S_tax_wide_filtered)
#colnames_decon<-data.frame(colnames_decon)

#seqtabNoC_WP2A_12S_tax_wide_filtered_DECON<-seqtabNoC_WP2A_12S_tax_wide_filtered %>% relocate("OTU_ID","ENC_T09","ENC_T10","ENC_T17","ENC_T18","ENC_T19")
                                                                                              
#colnames_decon<-colnames(seqtabNoC_WP2A_12S_tax_wide_filtered_DECON)
#colnames_decon<-data.frame(colnames_decon)
#decontaminated <- remove.cont(data = seqtabNoC_WP2A_12S_tax_wide_filtered_DECON,numb.blanks=5,#numb.ind=c(13,#E01
                                                                                                   # 15,#E02
                                                                                                    #14,#E03
                                                                                                    #16,#E04
                                                                                                    #18,#E05
                                                                                                    #18,#E06
                                                                                                    #18,#E07
                                                                                                    #18,#E08
                                                                                                    #18,#E09
                                                                                                    #18,#E10
                                                                                                    #18,#E11
                                                                                                    #18,#E12
                                                                                                    #18,#E13
                                                                                                    #18),
#                        taxa=F)

#seqtabNoC_WP2A_12S_tax_wide_filtered<-decontaminated
#seqtabNoC_WP2A_12S_tax_wide_filtered$Mean.blank<-NULL
#seqtabNoC_WP2A_12S_tax_wide_filtered<-t(seqtabNoC_WP2A_12S_tax_wide_filtered)
#seqtabNoC_WP2A_12S_tax_wide_filtered<-data.frame(seqtabNoC_WP2A_12S_tax_wide_filtered)

#header.true <- function(seqtabNoC_WP2A_12S_tax_wide_filtered) {
 # names(seqtabNoC_WP2A_12S_tax_wide_filtered) <- as.character(unlist(seqtabNoC_WP2A_12S_tax_wide_filtered[1,]))
  #seqtabNoC_WP2A_12S_tax_wide_filtered[-1,]
#}
#seqtabNoC_WP2A_12S_tax_wide_filtered<-header.true(seqtabNoC_WP2A_12S_tax_wide_filtered)

#i <- c(1:21) 
#seqtabNoC_WP2A_12S_tax_wide_filtered[ , i] <- apply(seqtabNoC_WP2A_12S_tax_wide_filtered[ , i], 2,  
 #                                                   function(x) as.numeric(as.character(x)))

#REMOVE SAMPLES THAT NOW HAVE NO READS IN THEM
#seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_wide_filtered[rowSums(seqtabNoC_WP2A_12S_tax_wide_filtered[])>0,]

#write.csv(seqtabNoC_WP2A_12S_tax_wide_filtered,"WP2A_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2A_12S_tax_wide_filtered,"WP2A_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2A_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")


#DENSITY PLOT FOR ALL DATA IN WP2A

#seqtabNoC_WP2A_12S_tax_long_meta_clean_sum_reads<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
# group_by(species,SampleSite_time, Days,SampleSite_code) %>% 
# summarize(reads=sum(reads)) %>%
#rename(species_sum=species,SampleSite_time_sum =SampleSite_time,Days_sum=Days,SampleSite_code_sum=SampleSite_code, reads_sum= reads,)

#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_sum_reads, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  geom_density( stat = "identity",position="stack")+
#  theme_bw()+
#  theme(legend.position="none")+
#  facet_grid(~SampleSite_code_sum)

##SPLITTING BY SPECIES
#Alosa_sapidissima<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Alosa_sapidissima")
#Ammodytes_personatus<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Ammodytes_personatus")
Anguilla_anguilla<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Anguilla_anguilla")
#Anguilla_rostrata<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Anguilla_rostrata")
#Barbatula_barbatula<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Barbatula_barbatula")
#Chelon_labrosus<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Chelon_labrosus")
#Clupea_harengus<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Clupea_harengus")
#Clupea_pallasii<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Clupea_pallasii")
#Copadichromis_virginalis<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Copadichromis_virginalis")
#Cottus_perifretum<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Cottus_perifretum")
#Favonigobius_gymnauchen<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Favonigobius_gymnauchen")
#Gadus_macrocephalus<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Gadus_macrocephalus")
#Hippoglossus_stenolepis<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Hippoglossus_stenolepis")
#Leuciscus_burdigalensis<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Leuciscus_burdigalensis")
#Lota_lota<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Lota_lota")
#Myoxocephalus_scorpius<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Myoxocephalus_scorpius")
#Oncorhynchus_mykiss_x_Salmo_salar<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Oncorhynchus_mykiss_x_Salmo_salar")
#Osmerus_mordax<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Osmerus_mordax")
#Perca_fluviatilis<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Perca_fluviatilis")
#Pholis_crassispina<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Pholis_crassispina")
#Phoxinus_ujmonensis<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Phoxinus_ujmonensis")
#Platichthys_stellatus<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Platichthys_stellatus")
#Pungitius_pungitius<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Pungitius_pungitius")
#Salmo_obtusirostris<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Salmo_obtusirostris")
Salmo_trutta<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Salmo_trutta")
Salmo_salar<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Salmo_salar")
#Salvelinus_fontinalis_x_Salvelinus_malma<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Salvelinus_fontinalis_x_Salvelinus_malma")
#Salvelinus_malma<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Salvelinus_malma")
#Scomber_scombrus<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Scomber_scombrus")

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
# group_by(SampleSite_time, Days,SampleSite_code,species) %>% 
# summarize(reads=sum(reads)) %>%
# rename(SampleSite_time_sum =SampleSite_time,Days_sum=Days,SampleSite_code_sum=SampleSite_code,species_sum=species, reads_sum= reads)

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E01<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E01")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E01, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area(position = 'stack')

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E02<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E02")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E02, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area(position = 'stack')

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E03<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E03")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E03, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area(position = 'stack')

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E04<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E04")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E04, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area()

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E05<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E05")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E05, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area()

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E06<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E06")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E06, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area(position = 'stack')

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E07<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E07")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E07, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area()

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E08<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E08")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E08, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area()

#seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E09<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads %>% filter(SampleSite_code_sum == "E09")
#ggplot(data=seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered_sum_reads_E09, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  theme_bw()+geom_area()

######## Venn diagram WP2A ##########
#lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
#lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

#seqtabNoC_WP2A_12S_Venn<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$reads==0),]
#hist(seqtabNoC_WP2A_12S_Venn$reads, breaks=1000,xlim=c(0,50000))
#seqtabNoC_WP2A_12S_Venn <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "river_section")],seqtabNoC_WP2A_12S_Venn, by="SampleSite_time")

#seqtabNoC_WP2A_12S_Venn_unique<-seqtabNoC_WP2A_12S_Venn
#seqtabNoC_WP2A_12S_Venn_unique$SampleSite_time<-NULL
#seqtabNoC_WP2A_12S_Venn_unique$reads<-NULL
#seqtabNoC_WP2A_12S_Venn_unique<-unique(seqtabNoC_WP2A_12S_Venn_unique)


#upper_Venn<-seqtabNoC_WP2A_12S_Venn %>% filter(river_section == "upper")
#upper_Venn<-upper_Venn$species

#lower_Venn<-seqtabNoC_WP2A_12S_Venn %>% filter(river_section == "lower")
#lower_Venn<-lower_Venn$species

#middle_Venn<-seqtabNoC_WP2A_12S_Venn %>% filter(river_section == "middle")
#middle_Venn<-middle_Venn$species

#library(VennDiagram)
#venn.diagram(
#x = list(lower_Venn, middle_Venn,upper_Venn),
#category.names = c("lower" , "middle","upper"),
#filename = 'WP2A_venn_diagramm.png',
#output=TRUE)

#PLOT SCALED VENN DIAGRAM
#library(eulerr)
#library(limma)

#VennDiag <- euler(c("lower" = 6, "middle" = 1, "upper" = 1, "lower&middle" = 6, "lower&upper" = 0, 
 #                   "middle&upper" = 0, "lower&middle&upper" = 2))
#plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5,
 #    fill=c("darkblue", "deepskyblue2", "cornflowerblue"))


#NORMALISING READS
seqtabNoC_WP2A_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-colSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-data.frame(WP2A_12S_sps_read_summary)

seqtabNoC_WP2A_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2A_12S_tax_col_sum2)
seqtabNoC_WP2A_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2A_12S_tax_col_sum2)

seqtabNoC_WP2A_12S_tax_normalised<- merge(seqtabNoC_WP2A_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_12S_tax_col_sum2")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2A_12S_tax_col_sum3<-colSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
seqtabNoC_WP2A_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2A_12S_tax_col_sum3)

seqtabNoC_WP2A_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2A_12S_tax_col_sum3)

seqtabNoC_WP2A_12S_tax_normalised<- merge(seqtabNoC_WP2A_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_12S_tax_col_sum2")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_top<- merge(seqtabNoC_WP2A_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2A_12S_tax_col_sum3")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2A_12S_tax_normalised <- transform(seqtabNoC_WP2A_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2A_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2A_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2A_12S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Days")],seqtabNoC_WP2A_12S_tax_normalised, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_normalised, by="SampleSite_time")

########SEASONAL IMPACT ON SPECIES DETECTION########
sps_seasonal_detection<-unique(seqtabNoC_WP2A_12S_tax_normalised[c("SampleSite_time", "species")])
sps_seasonal_detection <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "season","SampleSite_code","Dist_from_lake")],sps_seasonal_detection, by="SampleSite_time")
#sps_seasonal_detection<-unique(sps_seasonal_detection[c("SampleSite_code","Dist_from_lake","season", "species")])

sps_seasonal_detection_summary<-count(sps_seasonal_detection, season, species)
sps_seasonal_detection_summary$season <- factor(sps_seasonal_detection_summary$season, levels = c("Winter", "Spring", "Summer", "Autumn"))

#sps_seasonal_detection_summary<-subset(sps_seasonal_detection_summary,duplicated(species) | duplicated(species, fromLast=TRUE))

ggplot(data=sps_seasonal_detection_summary, aes(x=season, y=n)) +
  geom_line()+
  geom_point()+facet_grid((.~species))+
  theme(axis.text.x = element_text(angle = 90))


#######SALMON CASE STUDY######

Salmo_salar_norm<-seqtabNoC_WP2A_12S_tax_normalised %>% filter(species == "Salmo_salar")
Salmo_salar_norm <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Month.Year", "season","Dist_from_lake","pH",
                                                             "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                             "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                             "temp_monthly_average")],Salmo_salar_norm, by="SampleSite_time")
Salmo_salar_norm$Days<-as.numeric(Salmo_salar_norm$Days)
Salmo_salar_norm$Month.Year<-as.factor(Salmo_salar_norm$Month.Year)

Salmo_salar_norm$Month.Year<-factor(Salmo_salar_norm$Month.Year,levels= c('4.17','5.17','6.17','7.17','8.17','9.17',
                                                    '10.17','11.17','12.17','1.18','2.18','3.18','4.18'),ordered=TRUE)                                                                                         

ggplot(Salmo_salar_norm, aes(x=Days, y=normalised_reads,color=Month.Year)) + 
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)+theme_bw()

ggplot(Salmo_salar_norm, aes(x=daily_flow_by_catchment, y=normalised_reads,color=SampleSite_code)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()

salmon_flow<-  lm(normalised_reads ~daily_flow_by_catchment*SampleSite_code, data=Salmo_salar_norm)
anova(salmon_flow)

#Salmo_salar_norm$Days<-as.numeric(Salmo_salar_norm$Days)
#ggplot(data=Salmo_salar_norm, aes(x=Days, y=normalised_reads,fill=season))+
# facet_grid(~ SampleSite_code)+theme_bw()+geom_bar(stat='identity')+ 
#  theme(axis.text.x = element_text(angle = 90))+ scale_fill_brewer(palette="Paired")+
#  geom_vline(xintercept = 156, linetype="dotted", color = "red", size=1)+
#  geom_vline(xintercept = 338, linetype="dotted", color = "red", size=1)+  
#  geom_vline(xintercept = 189,  color = "red", size=0.5)+
#  geom_vline(xintercept = 279,  color = "red", size=0.5)+
#geom_vline(xintercept = 349,  color = "blue", size=0.5)

#156=1st October
#189=2nd November
#279=31st January
#338=31st March
#349=11th April

salmon_estimates <- read_csv("salmon_estimates.csv")

Salmo_salar_norm_catch_compare<-Salmo_salar_norm %>% 
 group_by(Month.Year,SampleSite_code) %>% 
 summarize(normalised_reads=mean(normalised_reads)) %>%
rename(Month.Year=Month.Year,SampleSite_code=SampleSite_code,normalised_reads=normalised_reads)

Salmo_salar_norm_catch_compare <- merge(salmon_estimates[, c("Month.Year", "Total_LIVE_adult_biomass","DECAYING_adult_biomass",
                                                             "Total_adult_biomass",	"Total_salmon_biomass",	
                                                             "cum_live_adult_SA")],Salmo_salar_norm_catch_compare, by="Month.Year")

plive<-ggplot(Salmo_salar_norm_catch_compare, aes(x=Total_LIVE_adult_biomass, y=normalised_reads,color=SampleSite_code)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()
pdecay<-ggplot(Salmo_salar_norm_catch_compare, aes(x=DECAYING_adult_biomass, y=normalised_reads,color=SampleSite_code)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()
ptotal_adult<-ggplot(Salmo_salar_norm_catch_compare, aes(x=Total_adult_biomass, y=normalised_reads,color=SampleSite_code)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()
ptotal_all<-ggplot(Salmo_salar_norm_catch_compare, aes(x=Total_salmon_biomass, y=normalised_reads,color=SampleSite_code)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()
pSA<-ggplot(Salmo_salar_norm_catch_compare, aes(x=cum_live_adult_SA, y=normalised_reads,color=SampleSite_code)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()

ggarrange(plive, pdecay, ptotal_adult,ptotal_all,pSA, ncol = 2, nrow = 3)

Salmo_salar_norm_catch_compare_m1<-lm(normalised_reads~Total_LIVE_adult_biomass+
                                        DECAYING_adult_biomass+
                                        Total_adult_biomass+
                                        Total_salmon_biomass+
                                        cum_live_adult_SA, data =Salmo_salar_norm_catch_compare)
anova(Salmo_salar_norm_catch_compare_m1)

ggplot(Salmo_salar_norm_catch_compare, aes(x=DECAYING_adult_biomass, y=normalised_reads)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()


#EEL CASE STUDY
Anguilla_anguilla_norm<-seqtabNoC_WP2A_12S_tax_normalised %>% filter(species == "Anguilla_anguilla")
Anguilla_anguilla_norm <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "season","Dist_from_lake","pH",
                                                             "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                             "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                             "temp_monthly_average")],Anguilla_anguilla_norm, by="SampleSite_time")

ggplot(Anguilla_anguilla_norm, aes(x=daily_flow_by_catchment, y=normalised_reads,color=SampleSite_code)) + 
  geom_point()+
  geom_smooth(method=lm)+theme_bw()

Anguilla_anguilla_flow<-  lm(normalised_reads ~daily_flow_by_catchment*SampleSite_code, data=Anguilla_anguilla_norm)
anova(Anguilla_anguilla_flow)

Anguilla_anguilla_norm$Days<-as.numeric(Anguilla_anguilla_norm$Days)
ggplot(data=Anguilla_anguilla_norm, aes(x=Days, y=normalised_reads,fill=season))+
  facet_grid(~ SampleSite_code)+theme_bw()+geom_bar(stat='identity')+ 
  theme(axis.text.x = element_text(angle = 90))+ scale_fill_brewer(palette="Paired")+
  geom_vline(xintercept = 279, linetype="dotted", color = "red", size=1)+
  geom_vline(xintercept = 356, linetype="dotted", color = "red", size=1)+
  geom_vline(xintercept = 0, linetype="dotted", color = "red", size=1)+
  geom_vline(xintercept = 24, linetype="dotted", color = "red", size=1)

#0=27th April 2017
#24=31st May 2017
#279=1st Feb 2018
#356=18th April 2018

#PLOTTING TOP Taxon
seqtabNoC_WP2A_12S_tax_top <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_top, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2A_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+facet_wrap(~ SampleSite_code,scales="free",ncol = 15)

#write.csv(seqtabNoC_WP2A_12S_tax_normalised,"NORM_WP2A_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2A_12S_tax_normalised,"NORM_WP2A_ASV_table_LONG_filtered.csv")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2A_12S_tax_normalised_wide<-seqtabNoC_WP2A_12S_tax_normalised
seqtabNoC_WP2A_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2A_12S_tax_normalised_wide<-seqtabNoC_WP2A_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2A_12S_tax_normalised_wide <- spread(seqtabNoC_WP2A_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2A_12S_tax_normalised_wide[is.na(seqtabNoC_WP2A_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2A_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2A_12S_tax_normalised_wide))

names(seqtabNoC_WP2A_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2A_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2A_12S_tax_normalised_wide <- seqtabNoC_WP2A_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2A_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2A_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2A_12S_tax_normalised_wide,"NORM_WP2A_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2A_12S_tax_normalised_wide,"NORM_WP2A_ASV_table_wide_filtered.csv")

########PERMANOVA ON TAXON########
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")

NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$ID),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:16)]
adonis_data<-adonis_data[,-c(1:16)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis2(adonis_data ~ 
                            adonis_data_meta$season/adonis_data_meta$Days+
                            adonis_data_meta$Dist_from_lake+
                            adonis_data_meta$daily_flow_by_catchment+
                            adonis_data_meta$monthly_flow_by_catchment+
                            adonis_data_meta$seasonal_flow_by_catchment+
                            adonis_data_meta$week_before_flow_by_catchment+
                            adonis_data_meta$river_width_m+
                            adonis_data_meta$river_gradient_percent+
                            adonis_data_meta$rainfall_monthly_average+
                            adonis_data_meta$temp_monthly_average+
                            adonis_data_meta$gran_alk_ueq_L+
                            adonis_data_meta$conductivity_uS_cm+
                            adonis_data_meta$pH,
                            #adonis_data_meta$Arable_and_horticulture_ha	+
                            #adonis_data_meta$Broadleaf_woodland_ha	+
                            #adonis_data_meta$Coniferous_woodland_ha+	
                            #adonis_data_meta$Acid_grassland_ha	+
                            #adonis_data_meta$Neutral_grassland_ha+	
                            #adonis_data_meta$Improved_grassland_ha+	
                            #adonis_data_meta$Heather_grassland_ha	+
                            #adonis_data_meta$Heather_ha	+
                            #adonis_data_meta$Bog_ha	+
                            #adonis_data_meta$Freshwater_ha+	
                            #adonis_data_meta$Inland_rock_ha	+
                          #adonis_data_meta$Supralittoral_sediment_ha	+
                          #adonis_data_meta$Suburban_ha	+
                          #adonis_data_meta$Urban_ha,
                          method="bray",perm=999,by="margin")
print(adon.results_WP2)


#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_ID <- read.csv("NORM_WP2A_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_ID[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_ID[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake","season")],data.scores, by="SampleSite_time")

gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~season,data.scores_meta_data,mean),by="season")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=season,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))+
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")


########ALPHA DIVERSITY METRICS########
shannon_index<-diversity(seqtabNoC_WP2A_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2A_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2A_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2A<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2A<-merge(alpha_12S_WP2A,species_count,by="ID")

#write.csv(alpha_12S_WP2A,"alpha_12S_WP2A_replicates.csv")
write.csv(alpha_12S_WP2A,"alpha_12S_WP2A.csv")

######## WP2A alpha diversity modeling ##########
library(mgcv)

species_count_dataset<-species_count %>% rename(SampleSite_time = ID)
species_count_dataset<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code",
                                               "season",
                                               "Days",
                                               "Dist_from_lake",
                                               "temp_monthly_average",
                                               "conductivity_uS_cm",
                                               "gran_alk_ueq_L","daily_flow_by_catchment",
                                               "monthly_flow_by_catchment",
                                               "seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment",
                                               "river_width_m",
                                               "river_gradient_percent",
                                               "rainfall_monthly_average",
                                               "pH","river_section")],species_count_dataset, by="SampleSite_time")

species_count_dataset$Days<-as.numeric(species_count_dataset$Days)
species_count_dataset<-species_count_dataset[species_count_dataset$species_count != 0, ] 

shannon_index_dataset<-shannon_index %>% rename(SampleSite_time = ID)
shannon_index_dataset<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code",
                                                         "season",
                                                         "Days",
                                                         "Dist_from_lake",
                                                         "temp_monthly_average",
                                                         "conductivity_uS_cm",
                                                         "gran_alk_ueq_L","daily_flow_by_catchment",
                                                         "monthly_flow_by_catchment",
                                                         "seasonal_flow_by_catchment",
                                                         "week_before_flow_by_catchment",
                                                         "river_width_m",
                                                         "river_gradient_percent",
                                                         "rainfall_monthly_average",
                                                         "pH","river_section")],shannon_index_dataset, by="SampleSite_time")

shannon_index_dataset$Days<-as.numeric(shannon_index_dataset$Days)
shannon_index_dataset<-shannon_index_dataset[shannon_index_dataset$shannon_index != 0, ] 

#GAMs
write.csv(species_count_dataset,"12S_species_count_dataset.csv" )
ggplot(species_count_dataset, aes(Days, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))#+facet_wrap(~SampleSite_code)

ggplot(species_count_dataset, aes(gran_alk_ueq_L, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+facet_wrap(~SampleSite_code)

ggplot(species_count_dataset, aes(temp_monthly_average, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+facet_wrap(~SampleSite_code)

ggplot(species_count_dataset, aes(daily_flow_by_catchment, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+facet_wrap(~SampleSite_code)

ggplot(species_count_dataset, aes(Dist_from_lake, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(species_count_dataset, aes(conductivity_uS_cm, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+facet_wrap(~SampleSite_code)

ggplot(species_count_dataset, aes(pH, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+facet_wrap(~SampleSite_code)

ggplot(species_count_dataset, aes(river_width_m, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(species_count_dataset, aes(season, species_count)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+facet_wrap(~SampleSite_code)

data = species_count_dataset %>%
  mutate(river_section = factor(river_section))

mod_gam2 <- gam(species_count ~ s(Days) + season+
                  s(Dist_from_lake)+
                  s(pH)+
                  s(temp_monthly_average)+
                  s(conductivity_uS_cm)+
                  s(gran_alk_ueq_L)+
                  s(daily_flow_by_catchment)+
                  s(monthly_flow_by_catchment)+
                  s(seasonal_flow_by_catchment)+
                  s(week_before_flow_by_catchment)+
                  s(river_gradient_percent)+
                  s(rainfall_monthly_average),
                  data=data,
                  method="REML", select=TRUE)
#plot(mod_gam2)
gam.check(mod_gam2)
appraise(mod_gam2)
draw(mod_gam2)

summary(mod_gam2)

mod_gam3 <- gam(species_count ~ 
                  s(Dist_from_lake)+
                  s(temp_monthly_average)+
                  s(conductivity_uS_cm)+
                  s(gran_alk_ueq_L)+
                  s(monthly_flow_by_catchment)+
                  s(seasonal_flow_by_catchment)+
                  s(week_before_flow_by_catchment)+season,
                data=data,
                method="REML", select=TRUE)
summary(mod_gam3)
anova(mod_gam3)
plot(mod_gam3)
draw(mod_gam3)



species_count_dataset<-na.omit(species_count_dataset)

#Linear models
ggplot(species_count_dataset, aes(Dist_from_lake, species_count,colour=season)) + 
  geom_point() + 
  geom_smooth(method = lm, se = TRUE)+theme_bw()
ggplot(species_count_dataset, aes(rainfall_monthly_average, species_count)) + 
  geom_point() + 
  geom_smooth(se = TRUE)+theme_bw()
ggplot(shannon_index_dataset, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  geom_smooth(method = lm, se = TRUE)+theme_bw()

hist(species_count_dataset$species_count)

lm_12S_species_count = lm(species_count~
                            season*
                            Dist_from_lake+
                            temp_monthly_average+
                            conductivity_uS_cm+
                            gran_alk_ueq_L+
                            daily_flow_by_catchment+
                            seasonal_flow_by_catchment+
                            monthly_flow_by_catchment+
                            week_before_flow_by_catchment+
                            river_width_m+
                            river_gradient_percent+
                            rainfall_monthly_average+
                            pH, data = species_count_dataset)
anova(lm_12S_species_count)
step(lm_12S_species_count)

lm_12S_species_count_2<-lm(formula = species_count ~ season + Dist_from_lake + conductivity_uS_cm + 
                             gran_alk_ueq_L + daily_flow_by_catchment + seasonal_flow_by_catchment + 
                             river_width_m + river_gradient_percent + season:Dist_from_lake, 
                           data = species_count_dataset)
anova(lm_12S_species_count_2)

m.lst <- lstrends(lm_12S_species_count_2, "season", var="Dist_from_lake")
summary(m.lst)
plot(m.lst)
pairs(m.lst)

#effect size
cohens_f_squared_12S<-cohens_f_squared(lm_12S_species_count_2)
plot(cohens_f_squared_12S)


#pairwise comparisons 
#emmeans_season_river_section<-emmeans(lm_12S_species_count_2, pairwise ~ Dist_from_lake*season)
#plot(emmeans_season_river_section)

######## Turnover and nestedness ########
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered

SampleSite_time<-NORM_WP2A_ASV_table_wide_filtered$X1
SampleSite_time<-data.frame(SampleSite_time)
SampleSite_time <- tibble::rownames_to_column(SampleSite_time, "SiteID_1")

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

SampleSite_time<-merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","Days")],SampleSite_time, by="SampleSite_time")

adonis_data_meta<- adonis_data[c(1:16)]
adonis_data<-adonis_data[,-c(1:16)]

adonis_data<-na.omit(adonis_data)

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
#out <- nestedtemp(adonis_data)
#plot(out)
#plot(out, kind="incid")

library(betapart)
#library(qgraph)
adonis_data<-adonis_data %>% mutate_if(is.numeric, ~1 * (. != 0))
adonis_data<-adonis_data[, colSums(adonis_data != 0) > 0]
sorensen.beta<-beta.pair(adonis_data, index.family="sorensen")

nestedness <- as.matrix(sorensen.beta[[3]])
nestedness <- melt(nestedness)
names(nestedness) <- c("SiteID_1", "SiteID_2", "Beta_Distance")

###-----------------------
## get nestedness matrix (this is [[2]]) into a data frame
###-----------------------
m.beta.sorn <- as.matrix(sorensen.beta[[2]])
m.beta.sorn <- melt(m.beta.sorn)
names(m.beta.sorn) <- c("SiteID_1", "SiteID_2", "Beta_nest")
beta.measures.sor <- merge(nestedness, m.beta.sorn, by = c("SiteID_1","SiteID_2")) 

###-----------------------
## get turnover matrix (here simpson, this is [[1]]) into a data frame
###-----------------------
m.beta.sort <- as.matrix(sorensen.beta[[1]])
m.beta.sort  <- melt(m.beta.sort )
names(m.beta.sort ) <- c("SiteID_1", "SiteID_2", "Beta_turn")
beta.measures.sor <- merge(beta.measures.sor, m.beta.sort , by = c("SiteID_1","SiteID_2")) 

beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_1")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_1 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_1 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_1 = Days)

SampleSite_time<-SampleSite_time %>% rename(SiteID_2 = SiteID_1)
beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_2")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_2 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_2 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_2 = Days)

beta.measures.sor$Dist_from_lake_DIFF<-(beta.measures.sor$Dist_from_lake_2-beta.measures.sor$Dist_from_lake_1)
beta.measures.sor$Days_2<-as.numeric(beta.measures.sor$Days_2)
beta.measures.sor$Days_1<-as.numeric(beta.measures.sor$Days_1)
beta.measures.sor$Days_DIFF<-(beta.measures.sor$Days_2-beta.measures.sor$Days_1)
beta.measures.sor <- beta.measures.sor[beta.measures.sor$Days_DIFF >= 0, ]

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_turn, colour=Days_DIFF)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_nest, colour=Days_DIFF)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance, colour=Days_DIFF)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))

hist(beta.measures.sor$Days_DIFF_sqrd)

######## Sparse partial least squares analysis ##########
library(mixOmics)

NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
#adonis_data<-adonis_data[adonis_data$X1 != "E04_T07", ]  

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                               #"season",
                                               "Days",
                                               "Dist_from_lake",
                                               "temp_monthly_average",
                                               "conductivity_uS_cm",
                                               "gran_alk_ueq_L","daily_flow_by_catchment",
                                               "monthly_flow_by_catchment",
                                               "seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment",
                                               "river_width_m",
                                               "river_gradient_percent",
                                               "rainfall_monthly_average",
                                               "pH", "Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               "Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               "Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)

adonis_data_meta<- adonis_data[c(1:28)]
adonis_data<-adonis_data[,-c(1:28)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
adonis_data_meta$SampleSite_time<-NULL

ncomp = 1
#result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'regression')
result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'canonical')

design <- data.frame(samp = adonis_data_meta$sample)

tune.spls <- perf(result.spls, validation = 'Mfold', folds = 5,
                  criterion = 'all',nrepeat = 99, progressBar = TRUE)

tune.spls$R2
plot(tune.spls$Q2.total)
abline(h = 0.0975)

plot(tune.spls, criterion="R2", type = 'l')
plot(tune.spls, criterion="Q2", type = 'l')
plot(tune.spls)

cim_plot<-cim(result.spls,
              comp = 1:1,
              margins = c(15, 15))

######## Mantel test ##########
adonis_data_meta$Days<-NULL
adonis_data_meta$Dist_from_lake<-NULL
adonis_data_meta$river_width_m<-NULL

dist.abund = vegdist(adonis_data, method = "bray")
dist.meta = dist(adonis_data_meta, method = "euclidean")

abund_meta = mantel(dist.abund, dist.meta, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_meta

aa = as.vector(dist.abund)
mm = as.vector(dist.meta)
mat = data.frame(aa,mm)

ggplot(mat, aes(y = aa, x = mm)) + 
  geom_point(size = 1) + 
  theme_bw()+ geom_smooth(method = lm)+ scale_y_continuous(limits = c(0, 1))

#________________________________________________________

########WP2C_Glatt BLAST  TAXONOMY ANALYSIS########

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Glatt_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered$species=="Sardina_pilchardus"),]

#CONVERTING FILTERED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered<-seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered <- seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered,"WP2C_Glatt_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered,"WP2C_Glatt_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
Sardina_pilchardus<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Glatt_12S_Venn<-seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Glatt_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Glatt_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Glatt_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS

seqtabNoC_WP2C_Glatt_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[c(1:25)])
seqtabNoC_WP2C_Glatt_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2)
seqtabNoC_WP2C_Glatt_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2)

seqtabNoC_WP2C_Glatt_12S_tax_normalised<- merge(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Glatt_12S_tax_col_sum2")],seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Glatt_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Glatt_12S_tax_wide_filtered[c(1:25)])
seqtabNoC_WP2C_Glatt_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3)

seqtabNoC_WP2C_Glatt_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3)

seqtabNoC_WP2C_Glatt_12S_tax_normalised<- merge(seqtabNoC_WP2C_Glatt_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Glatt_12S_tax_col_sum2")],seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Glatt_12S_tax_top<- merge(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Glatt_12S_tax_col_sum3")],seqtabNoC_WP2C_Glatt_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Glatt_12S_tax_normalised <- transform(seqtabNoC_WP2C_Glatt_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Glatt_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Glatt_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")


#PLOTTING TOP Taxon
seqtabNoC_WP2C_Glatt_12S_tax_top$SampleSite_time<-factor(seqtabNoC_WP2C_Glatt_12S_tax_top$SampleSite_time,levels= c('GL01',
                         'GL13','GL02','GL03','GL04','GL05','GL06','GL07','GL08','GL09','GL10','GL11','GL12'),ordered=TRUE)

seqtabNoC_WP2C_Glatt_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Glatt_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))

#write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised,"NORM_WP2C_Glatt_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised,"NORM_WP2C_Glatt_ASV_table_LONG_filtered.csv")

#CONVERTING NORMALISED LONG FORMAT INTO SampleSite_time FORMAT
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide<-seqtabNoC_WP2C_Glatt_12S_tax_normalised
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide<-seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide <- seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide,"NORM_WP2C_Glatt_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide,"NORM_WP2C_Glatt_ASV_table_wide_filtered.csv")


#PERMANOVA ON TAXON
NORM_WP2C_Glatt_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Glatt_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(SampleSite_time = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Average_Soil_Temp_5cm_oC+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)

print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_SampleSite_time <- read.csv("NORM_WP2C_Glatt_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_SampleSite_time[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_SampleSite_time[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Glatt_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Glatt<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Glatt<-merge(alpha_12S_WP2C_Glatt,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Glatt,"alpha_12S_WP2C_Glatt_replicates.csv")
write.csv(alpha_12S_WP2C_Glatt,"alpha_12S_WP2C_Glatt.csv")

#________________________________________________________

########WP2C_Towy BLAST  TAXONOMY ANALYSIS########

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Towy_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]

#CONVERTING FILTERED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Towy_12S_tax_wide_filtered<-seqtabNoC_WP2C_Towy_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Towy_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered <- seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered,"WP2C_Towy_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered,"WP2C_Towy_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
Sardina_pilchardus<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Towy_12S_Venn<-seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Towy_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Towy_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Towy_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS
seqtabNoC_WP2C_Towy_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[c(1:10)])
seqtabNoC_WP2C_Towy_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Towy_12S_tax_col_sum2)
seqtabNoC_WP2C_Towy_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Towy_12S_tax_col_sum2)

seqtabNoC_WP2C_Towy_12S_tax_normalised<- merge(seqtabNoC_WP2C_Towy_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Towy_12S_tax_col_sum2")],seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Towy_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Towy_12S_tax_wide_filtered[c(1:10)])
seqtabNoC_WP2C_Towy_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Towy_12S_tax_col_sum3)

seqtabNoC_WP2C_Towy_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Towy_12S_tax_col_sum3)

seqtabNoC_WP2C_Towy_12S_tax_normalised<- merge(seqtabNoC_WP2C_Towy_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Towy_12S_tax_col_sum2")],seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Towy_12S_tax_top<- merge(seqtabNoC_WP2C_Towy_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Towy_12S_tax_col_sum3")],seqtabNoC_WP2C_Towy_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Towy_12S_tax_normalised <- transform(seqtabNoC_WP2C_Towy_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Towy_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Towy_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))
##scale_fill_brewer(palette="Paired")


#PLOTTING TOP Taxon
seqtabNoC_WP2C_Towy_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Towy_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))

#write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised,"NORM_WP2C_Towy_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised,"NORM_WP2C_Towy_ASV_table_LONG_filtered.csv")


#CONVERTING NORMALISED LONG FORMAT INTO SampleSite_time FORMAT
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide<-seqtabNoC_WP2C_Towy_12S_tax_normalised
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Towy_12S_tax_normalised_wide<-seqtabNoC_WP2C_Towy_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Towy_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide <- seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Towy_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide,"NORM_WP2C_Towy_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide,"NORM_WP2C_Towy_ASV_table_wide_filtered.csv")


#PERMANOVA ON TAXON
NORM_WP2C_Towy_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Towy_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Average_Soil_Temp_5cm_oC+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)
adon.results_WP2.2<-adonis(adonis_data ~ adonis_data_meta$Days+
                             adonis_data_meta$Dist_from_lake+
                             adonis_data_meta$Average_Soil_Temp_5cm_oC+
                             adonis_data_meta$Catchment_flow_m3_s+
                             #adonis_data_meta$gran_alk_ueq_L+
                             adonis_data_meta$Total_Rainfall+
                             adonis_data_meta$conductivity_uS_cm+
                             adonis_data_meta$pH,
                           method="bray",perm=999)
adon.results_WP2.3<-adonis(adonis_data ~ adonis_data_meta$Days*adonis_data_meta$Catchment_flow_m3_s+
                             adonis_data_meta$Dist_from_lake+
                             adonis_data_meta$Average_Soil_Temp_5cm_oC+
                             #adonis_data_meta$gran_alk_ueq_L+
                             adonis_data_meta$Total_Rainfall+
                             adonis_data_meta$conductivity_uS_cm+
                             adonis_data_meta$pH,
                           method="bray",perm=999)
adon.results_WP2.4<-adonis(adonis_data ~ adonis_data_meta$Days+adonis_data_meta$Catchment_flow_m3_s+
                             adonis_data_meta$Dist_from_lake+
                             adonis_data_meta$Average_Soil_Temp_5cm_oC+
                             adonis_data_meta$Total_Rainfall+
                             adonis_data_meta$conductivity_uS_cm+
                             adonis_data_meta$pH*adonis_data_meta$gran_alk_ueq_L,
                           method="bray",perm=999)
adon.results_WP2.5<-adonis(adonis_data ~ adonis_data_meta$Days+adonis_data_meta$Catchment_flow_m3_s*
                             adonis_data_meta$Dist_from_lake*
                             adonis_data_meta$Average_Soil_Temp_5cm_oC*
                             adonis_data_meta$Total_Rainfall*
                             adonis_data_meta$conductivity_uS_cm*
                             adonis_data_meta$pH*adonis_data_meta$gran_alk_ueq_L,
                           method="bray",perm=999)
print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_ID <- read.csv("NORM_WP2C_Towy_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_ID[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_ID[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Towy_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Towy<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Towy<-merge(alpha_12S_WP2C_Towy,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Towy,"alpha_12S_WP2C_Towy_replicates.csv")
write.csv(alpha_12S_WP2C_Towy,"alpha_12S_WP2C_Towy.csv")


#________________________________________________________

########WP2C_Gwash BLAST  TAXONOMY ANALYSIS#########

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Gwash_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
 summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered$species=="Sardina_pilchardus"),]

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered <- seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered,"WP2C_Gwash_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered,"WP2C_Gwash_ASV_table_wide_filtered.csv")


#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered_transpose<-t(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered)

##SPLITTING BY SPECIES
Anguilla_anguilla<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered %>% filter(species == "Anguilla_anguilla")

Sardina_pilchardus<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Gwash_12S_Venn<-seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Gwash_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Gwash_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Gwash_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS
seqtabNoC_WP2C_Gwash_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[c(1:23)])
seqtabNoC_WP2C_Gwash_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2)
seqtabNoC_WP2C_Gwash_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2)

seqtabNoC_WP2C_Gwash_12S_tax_normalised<- merge(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Gwash_12S_tax_col_sum2")],seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Gwash_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Gwash_12S_tax_wide_filtered[c(1:23)])
seqtabNoC_WP2C_Gwash_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3)

seqtabNoC_WP2C_Gwash_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3)

seqtabNoC_WP2C_Gwash_12S_tax_normalised<- merge(seqtabNoC_WP2C_Gwash_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Gwash_12S_tax_col_sum2")],seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Gwash_12S_tax_top<- merge(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Gwash_12S_tax_col_sum3")],seqtabNoC_WP2C_Gwash_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Gwash_12S_tax_normalised <- transform(seqtabNoC_WP2C_Gwash_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Gwash_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Gwash_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))
##scale_fill_brewer(palette="Paired")


#PLOTTING TOP Taxon
seqtabNoC_WP2C_Gwash_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Gwash_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))


#write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised,"NORM_WP2C_Gwash_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised,"NORM_WP2C_Gwash_ASV_table_LONG_filtered.csv")


#CONVERTING NORMALISED LONG FORMAT INTO SampleSite_time FORMAT
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide<-seqtabNoC_WP2C_Gwash_12S_tax_normalised
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide<-seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide <- seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide,"NORM_WP2C_Gwash_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide,"NORM_WP2C_Gwash_ASV_table_wide_filtered.csv")


#PERMANOVA ON TAXON
NORM_WP2C_Gwash_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Gwash_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(SampleSite_time = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Average_Soil_Temp_5cm_oC+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)

print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_SampleSite_time <- read.csv("NORM_WP2C_Gwash_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_SampleSite_time[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_SampleSite_time[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Gwash_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Gwash<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Gwash<-merge(alpha_12S_WP2C_Gwash,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Gwash,"alpha_12S_WP2C_Gwash_replicates.csv")
write.csv(alpha_12S_WP2C_Gwash,"alpha_12S_WP2C_Gwash.csv")


#________________________________________________________

#########WP2C_Skan BLAST  TAXONOMY ANALYSIS##########

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered <- read_csv("WP2C_Skan_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
 summarize(reads=mean(reads)) %>%
 rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
#seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]

#CONVERTING FILTERED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Skan_12S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Skan_12S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered))

names(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered <- seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered,"WP2C_Skan_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered,"WP2C_Skan_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
#Alosa_sapidissima<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered %>% filter(species == "Alosa_sapidissima")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Skan_12S_Venn<-seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean[!(seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Skan_12S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Skan_12S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Skan_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS

seqtabNoC_WP2C_Skan_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[c(1:22)])
seqtabNoC_WP2C_Skan_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Skan_12S_tax_col_sum2)
seqtabNoC_WP2C_Skan_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Skan_12S_tax_col_sum2)

seqtabNoC_WP2C_Skan_12S_tax_normalised<- merge(seqtabNoC_WP2C_Skan_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Skan_12S_tax_col_sum2")],seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Skan_12S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Skan_12S_tax_wide_filtered[c(1:22)])
seqtabNoC_WP2C_Skan_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Skan_12S_tax_col_sum3)

seqtabNoC_WP2C_Skan_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Skan_12S_tax_col_sum3)

seqtabNoC_WP2C_Skan_12S_tax_normalised<- merge(seqtabNoC_WP2C_Skan_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Skan_12S_tax_col_sum2")],seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Skan_12S_tax_top<- merge(seqtabNoC_WP2C_Skan_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Skan_12S_tax_col_sum3")],seqtabNoC_WP2C_Skan_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Skan_12S_tax_normalised <- transform(seqtabNoC_WP2C_Skan_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Skan_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Skan_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))

#PLOTTING TOP Taxon
seqtabNoC_WP2C_Skan_12S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Skan_12S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+ 
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(0.3, 'cm'))

#write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised,"NORM_WP2C_Skan_ASV_table_LONG_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised,"NORM_WP2C_Skan_ASV_table_LONG_filtered.csv")


#CONVERTING NORMALISED LONG FORMAT INTO WSampleSite_timeE FORMAT
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide<-seqtabNoC_WP2C_Skan_12S_tax_normalised
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide$seqtabNoC_WP2_12S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Skan_12S_tax_normalised_wide<-seqtabNoC_WP2C_Skan_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Skan_12S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide))

names(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide <- seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Skan_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

#write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide,"NORM_WP2C_Skan_ASV_table_wide_filtered_replicates.csv")
write.csv(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide,"NORM_WP2C_Skan_ASV_table_wide_filtered.csv")

#PERMANOVA ON TAXON
NORM_WP2C_Skan_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Skan_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(SampleSite_time = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Days","Catchment_flow_m3_s","Average_Soil_Temp_5cm_oC","Total_Rainfall","Dist_from_lake","pH","conductivity_uS_cm","gran_alk_ueq_L")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:10)]
adonis_data<-adonis_data[,-c(1:10)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$Days+
                           adonis_data_meta$Dist_from_lake+
                           adonis_data_meta$Average_Soil_Temp_5cm_oC+
                           adonis_data_meta$Catchment_flow_m3_s+
                           adonis_data_meta$gran_alk_ueq_L+
                           adonis_data_meta$Total_Rainfall+
                           adonis_data_meta$conductivity_uS_cm+
                           adonis_data_meta$pH,
                         method="bray",perm=999)

print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_SampleSite_time <- read.csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_SampleSite_time[,-1]
com <- abund_table_12S[,col(abund_table_12S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_12S_SampleSite_time[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")


p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Skan_12S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_12S_WP2C_Skan<-merge(shannon_index,simpson_index,by="ID")
alpha_12S_WP2C_Skan<-merge(alpha_12S_WP2C_Skan,species_count,by="ID")

#write.csv(alpha_12S_WP2C_Skan,"alpha_12S_WP2C_Skan_replicates.csv")
write.csv(alpha_12S_WP2C_Skan,"alpha_12S_WP2C_Skan.csv")


#________________________________________________________

########COMPILING WP2 CURATED NORMALISED OTU TABLES (APART FROM WP2B) FOR ALPHA AND BETA DIVERSITY #####

#________________________________________________________

alpha_12S_WP2A <- read_csv("alpha_12S_WP2A.csv")
alpha_12S_WP2C_Glatt <- read_csv("alpha_12S_WP2C_Glatt.csv")
alpha_12S_WP2C_Gwash <- read_csv("alpha_12S_WP2C_Gwash.csv")
alpha_12S_WP2C_Towy <- read_csv("alpha_12S_WP2C_Towy.csv")
#alpha_12S_WP2C_Skan <- read_csv("alpha_12S_WP2C_Skan.csv")

alpha_12S_WP2<-rbind(alpha_12S_WP2A,alpha_12S_WP2C_Glatt,alpha_12S_WP2C_Gwash,alpha_12S_WP2C_Towy)#,alpha_12S_WP2C_Skan)
alpha_12S_WP2<-alpha_12S_WP2[!(alpha_12S_WP2$shannon_index==0),]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

alpha_12S_WP2<-alpha_12S_WP2 %>% rename(SampleSite_time = ID)

alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Dist_from_lake")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "River")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Date_sampled")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "river_width_m")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "pH")],alpha_12S_WP2, by="SampleSite_time")
alpha_12S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "conductivity_uS_cm")],alpha_12S_WP2, by="SampleSite_time")

#Only keep the conwy summer samples that were coordinated with the Gwash, Glatt and Towy
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "01.09.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "02.11.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "04.01.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "11.04.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "12.06.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "12.08.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "12.10.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "14.02.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "14.12.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "18.04.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "18.05.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "22.09.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "23.11.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "25.01.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "27.04.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "28.03.18" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "29.06.17" , ] 
alpha_12S_WP2<-alpha_12S_WP2[alpha_12S_WP2$Date_sampled != "29.06.17" , ] 

ggplot(alpha_12S_WP2, aes(x=Dist_from_lake, y=species_count, color=River)) +
  geom_point() + 
  geom_smooth(method = "lm",fullrange = T)+theme_bw()+scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2","#D55E00","#999999"))

ggplot(alpha_12S_WP2, aes(x=river_width_m, y=shannon_index, color=River)) +
  geom_point() + 
  geom_smooth(method = "lm",fullrange = T)+theme_bw()+scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2","#D55E00","#999999"))


hist(alpha_12S_WP2$shannon_index)

#ggplot(alpha_12S_WP2, aes(x=Dist_from_lake, y=shannon_index, color=River)) +
 #geom_point() + 
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2))+theme_bw()+scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2","#D55E00","#999999"))

#lm_12S_WP2 = lm(shannon_index~I(Dist_from_lake^2)*River, data = alpha_12S_WP2)
#lm_12S_WP2.1 = lm(shannon_index~Dist_from_lake*River, data = alpha_12S_WP2)
#AIC(lm_12S_WP2)
#AIC(lm_12S_WP2.1)
#anova(lm_12S_WP2,lm_12S_WP2.1)

lm_12S_WP2.1 = lm(species_count~Dist_from_lake*River*river_width_m, data = alpha_12S_WP2)
anova(lm_12S_WP2.1)

step(lm_12S_WP2.1)

lm_12S_WP2.2 = lm(formula = species_count ~ Dist_from_lake * River * river_width_m, 
                  data = alpha_12S_WP2)
anova(lm_12S_WP2.2)
plot(lm_12S_WP2.2)

m.lst.2 <- lstrends(lm_12S_WP2.2, "River", var="Dist_from_lake")
summary(m.lst.2)
pairs(m.lst.2)


#BETA DIVERSITY - PERMANOVA ON TAXON
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2A_ASV_table_LONG_filtered <- read_csv("NORM_WP2A_ASV_table_LONG_filtered.csv")
NORM_WP2C_Glatt_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_LONG_filtered.csv")
NORM_WP2C_Gwash_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_LONG_filtered.csv")
NORM_WP2C_Towy_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_LONG_filtered.csv")
NORM_WP2C_Skan_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_LONG_filtered.csv")

NORM_WP2A_ASV_table_LONG_filtered$Days<-NULL
NORM_WP2A_ASV_table_LONG_filtered$SampleSite_code<-NULL

NORM_WP2A_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Glatt_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Gwash_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Towy_ASV_table_LONG_filtered[,3] <- NULL
NORM_WP2C_Skan_ASV_table_LONG_filtered[,3] <- NULL

NORM_WP2_ASV_table_LONG_filtered<-rbind(NORM_WP2A_ASV_table_LONG_filtered,NORM_WP2C_Glatt_ASV_table_LONG_filtered,
                                        NORM_WP2C_Gwash_ASV_table_LONG_filtered ,NORM_WP2C_Towy_ASV_table_LONG_filtered,NORM_WP2C_Skan_ASV_table_LONG_filtered)

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2A_12S_tax_normalised_wide<-NORM_WP2_ASV_table_LONG_filtered
seqtabNoC_WP2A_12S_tax_normalised_wide$X1<- NULL 
seqtabNoC_WP2A_12S_tax_normalised_wide$reads <- NULL 

seqtabNoC_WP2A_12S_tax_normalised_wide<-seqtabNoC_WP2A_12S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2A_12S_tax_normalised_wide <- spread(seqtabNoC_WP2A_12S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2A_12S_tax_normalised_wide[is.na(seqtabNoC_WP2A_12S_tax_normalised_wide)] <- 0
seqtabNoC_WP2A_12S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2A_12S_tax_normalised_wide))

names(seqtabNoC_WP2A_12S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2A_12S_tax_normalised_wide[1, ])
seqtabNoC_WP2A_12S_tax_normalised_wide <- seqtabNoC_WP2A_12S_tax_normalised_wide[-1, ]
seqtabNoC_WP2A_12S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2A_12S_tax_normalised_wide, function(x) type.convert(as.character(x)))

adonis_data<-seqtabNoC_WP2A_12S_tax_normalised_wide
adonis_data<-setDT(adonis_data, keep.rownames = TRUE)[]
adonis_data<-data.frame(adonis_data)
adonis_data_samples<-adonis_data$rn
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = rn)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","River","Date_sampled","Dist_from_lake","river_width_m")],adonis_data, by="SampleSite_time")

adonis_data<-adonis_data[adonis_data$Date_sampled != "01.09.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "02.11.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "04.01.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "11.04.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.06.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.08.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.10.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "14.02.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "14.12.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "18.04.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "18.05.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "22.09.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "23.11.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "25.01.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "27.04.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "28.03.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "29.06.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "29.06.17" , ] 

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:5)]
adonis_data<-adonis_data[,-c(1:5)]

adon.results_WP2<-adonis2(adonis_data ~ adonis_data_meta$River*adonis_data_meta$Dist_from_lake*adonis_data_meta$river_width_m,
                         method="bray",perm=999)
print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_12S_ID <- read.csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_12S <- abund_table_12S_ID[,-1]

com <- adonis_data[,col(adonis_data)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = adonis_data_meta[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","River")],data.scores, by="SampleSite_time")

gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~River,data.scores_meta_data,mean),by="River")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=River,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")+
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))

######## Turnover and nestedness Conwy ########
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2_ASV_table_wide_filtered

SampleSite_time<-NORM_WP2_ASV_table_wide_filtered$X1
SampleSite_time<-data.frame(SampleSite_time)
SampleSite_time <- tibble::rownames_to_column(SampleSite_time, "SiteID_1")

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average","Date_sampled"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

adonis_data<-adonis_data[adonis_data$Date_sampled != "01.09.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "02.11.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "04.01.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "11.04.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.06.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.08.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "12.10.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "14.02.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "14.12.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "18.04.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "18.05.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "22.09.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "23.11.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "25.01.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "27.04.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "28.03.18" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "29.06.17" , ] 
adonis_data<-adonis_data[adonis_data$Date_sampled != "29.06.17" , ] 

SampleSite_time<-merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","Days")],SampleSite_time, by="SampleSite_time")

adonis_data_meta<- adonis_data[c(1:17)]
adonis_data<-adonis_data[,-c(1:17)]

adonis_data<-na.omit(adonis_data)

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
#out <- nestedtemp(adonis_data)
#plot(out)
#plot(out, kind="incid")

library(betapart)
#library(qgraph)
adonis_data<-adonis_data %>% mutate_if(is.numeric, ~1 * (. != 0))
adonis_data<-adonis_data[, colSums(adonis_data != 0) > 0]
sorensen.beta<-beta.pair(adonis_data, index.family="sorensen")

nestedness <- as.matrix(sorensen.beta[[3]])
nestedness <- melt(nestedness)
names(nestedness) <- c("SiteID_1", "SiteID_2", "Beta_Distance")

###-----------------------
## get nestedness matrix (this is [[2]]) into a data frame
###-----------------------
m.beta.sorn <- as.matrix(sorensen.beta[[2]])
m.beta.sorn <- melt(m.beta.sorn)
names(m.beta.sorn) <- c("SiteID_1", "SiteID_2", "Beta_nest")
beta.measures.sor <- merge(nestedness, m.beta.sorn, by = c("SiteID_1","SiteID_2")) 

###-----------------------
## get turnover matrix (here simpson, this is [[1]]) into a data frame
###-----------------------
m.beta.sort <- as.matrix(sorensen.beta[[1]])
m.beta.sort  <- melt(m.beta.sort )
names(m.beta.sort ) <- c("SiteID_1", "SiteID_2", "Beta_turn")
beta.measures.sor <- merge(beta.measures.sor, m.beta.sort , by = c("SiteID_1","SiteID_2")) 

beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_1")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_1 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_1 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_1 = Days)

SampleSite_time<-SampleSite_time %>% rename(SiteID_2 = SiteID_1)
beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_2")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_2 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_2 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_2 = Days)

beta.measures.sor$Dist_from_lake_DIFF<-(beta.measures.sor$Dist_from_lake_2-beta.measures.sor$Dist_from_lake_1)
beta.measures.sor_Conwy<-beta.measures.sor

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_turn)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_nest)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))

######## Turnover and nestedness Glatt ########

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2_ASV_table_wide_filtered

SampleSite_time<-NORM_WP2_ASV_table_wide_filtered$X1
SampleSite_time<-data.frame(SampleSite_time)
SampleSite_time <- tibble::rownames_to_column(SampleSite_time, "SiteID_1")

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average","Date_sampled"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

SampleSite_time<-merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","Days")],SampleSite_time, by="SampleSite_time")

adonis_data_meta<- adonis_data[c(1:17)]
adonis_data<-adonis_data[,-c(1:17)]

adonis_data<-na.omit(adonis_data)

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
#out <- nestedtemp(adonis_data)
#plot(out)
#plot(out, kind="incid")

library(betapart)
#library(qgraph)
adonis_data<-adonis_data %>% mutate_if(is.numeric, ~1 * (. != 0))
adonis_data<-adonis_data[, colSums(adonis_data != 0) > 0]
sorensen.beta<-beta.pair(adonis_data, index.family="sorensen")

nestedness <- as.matrix(sorensen.beta[[3]])
nestedness <- melt(nestedness)
names(nestedness) <- c("SiteID_1", "SiteID_2", "Beta_Distance")

###-----------------------
## get nestedness matrix (this is [[2]]) into a data frame
###-----------------------
m.beta.sorn <- as.matrix(sorensen.beta[[2]])
m.beta.sorn <- melt(m.beta.sorn)
names(m.beta.sorn) <- c("SiteID_1", "SiteID_2", "Beta_nest")
beta.measures.sor <- merge(nestedness, m.beta.sorn, by = c("SiteID_1","SiteID_2")) 

###-----------------------
## get turnover matrix (here simpson, this is [[1]]) into a data frame
###-----------------------
m.beta.sort <- as.matrix(sorensen.beta[[1]])
m.beta.sort  <- melt(m.beta.sort )
names(m.beta.sort ) <- c("SiteID_1", "SiteID_2", "Beta_turn")
beta.measures.sor <- merge(beta.measures.sor, m.beta.sort , by = c("SiteID_1","SiteID_2")) 

beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_1")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_1 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_1 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_1 = Days)

SampleSite_time<-SampleSite_time %>% rename(SiteID_2 = SiteID_1)
beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_2")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_2 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_2 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_2 = Days)

beta.measures.sor$Dist_from_lake_DIFF<-(beta.measures.sor$Dist_from_lake_2-beta.measures.sor$Dist_from_lake_1)
beta.measures.sor_Glatt<-beta.measures.sor

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_turn)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_nest)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))

######## Turnover and nestedness Gwash ########

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2_ASV_table_wide_filtered

SampleSite_time<-NORM_WP2_ASV_table_wide_filtered$X1
SampleSite_time<-data.frame(SampleSite_time)
SampleSite_time <- tibble::rownames_to_column(SampleSite_time, "SiteID_1")

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average","Date_sampled"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

SampleSite_time<-merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","Days")],SampleSite_time, by="SampleSite_time")

adonis_data_meta<- adonis_data[c(1:17)]
adonis_data<-adonis_data[,-c(1:17)]

adonis_data<-na.omit(adonis_data)

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
#out <- nestedtemp(adonis_data)
#plot(out)
#plot(out, kind="incid")

library(betapart)
#library(qgraph)
adonis_data<-adonis_data %>% mutate_if(is.numeric, ~1 * (. != 0))
adonis_data<-adonis_data[, colSums(adonis_data != 0) > 0]
sorensen.beta<-beta.pair(adonis_data, index.family="sorensen")

nestedness <- as.matrix(sorensen.beta[[3]])
nestedness <- melt(nestedness)
names(nestedness) <- c("SiteID_1", "SiteID_2", "Beta_Distance")

###-----------------------
## get nestedness matrix (this is [[2]]) into a data frame
###-----------------------
m.beta.sorn <- as.matrix(sorensen.beta[[2]])
m.beta.sorn <- melt(m.beta.sorn)
names(m.beta.sorn) <- c("SiteID_1", "SiteID_2", "Beta_nest")
beta.measures.sor <- merge(nestedness, m.beta.sorn, by = c("SiteID_1","SiteID_2")) 

###-----------------------
## get turnover matrix (here simpson, this is [[1]]) into a data frame
###-----------------------
m.beta.sort <- as.matrix(sorensen.beta[[1]])
m.beta.sort  <- melt(m.beta.sort )
names(m.beta.sort ) <- c("SiteID_1", "SiteID_2", "Beta_turn")
beta.measures.sor <- merge(beta.measures.sor, m.beta.sort , by = c("SiteID_1","SiteID_2")) 

beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_1")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_1 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_1 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_1 = Days)

SampleSite_time<-SampleSite_time %>% rename(SiteID_2 = SiteID_1)
beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_2")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_2 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_2 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_2 = Days)

beta.measures.sor$Dist_from_lake_DIFF<-(beta.measures.sor$Dist_from_lake_2-beta.measures.sor$Dist_from_lake_1)
beta.measures.sor_Gwash<-beta.measures.sor

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_turn)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_nest)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))



######## Turnover and nestedness Towy ########

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2_ASV_table_wide_filtered

SampleSite_time<-NORM_WP2_ASV_table_wide_filtered$X1
SampleSite_time<-data.frame(SampleSite_time)
SampleSite_time <- tibble::rownames_to_column(SampleSite_time, "SiteID_1")

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average","Date_sampled"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

SampleSite_time<-merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","Days")],SampleSite_time, by="SampleSite_time")

adonis_data_meta<- adonis_data[c(1:17)]
adonis_data<-adonis_data[,-c(1:17)]

adonis_data<-na.omit(adonis_data)

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
#out <- nestedtemp(adonis_data)
#plot(out)
#plot(out, kind="incid")

library(betapart)
#library(qgraph)
adonis_data<-adonis_data %>% mutate_if(is.numeric, ~1 * (. != 0))
adonis_data<-adonis_data[, colSums(adonis_data != 0) > 0]
sorensen.beta<-beta.pair(adonis_data, index.family="sorensen")

nestedness <- as.matrix(sorensen.beta[[3]])
nestedness <- melt(nestedness)
names(nestedness) <- c("SiteID_1", "SiteID_2", "Beta_Distance")

###-----------------------
## get nestedness matrix (this is [[2]]) into a data frame
###-----------------------
m.beta.sorn <- as.matrix(sorensen.beta[[2]])
m.beta.sorn <- melt(m.beta.sorn)
names(m.beta.sorn) <- c("SiteID_1", "SiteID_2", "Beta_nest")
beta.measures.sor <- merge(nestedness, m.beta.sorn, by = c("SiteID_1","SiteID_2")) 

###-----------------------
## get turnover matrix (here simpson, this is [[1]]) into a data frame
###-----------------------
m.beta.sort <- as.matrix(sorensen.beta[[1]])
m.beta.sort  <- melt(m.beta.sort )
names(m.beta.sort ) <- c("SiteID_1", "SiteID_2", "Beta_turn")
beta.measures.sor <- merge(beta.measures.sor, m.beta.sort , by = c("SiteID_1","SiteID_2")) 

beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_1")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_1 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_1 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_1 = Days)

SampleSite_time<-SampleSite_time %>% rename(SiteID_2 = SiteID_1)
beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_2")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_2 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_2 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_2 = Days)

beta.measures.sor$Dist_from_lake_DIFF<-(beta.measures.sor$Dist_from_lake_2-beta.measures.sor$Dist_from_lake_1)
beta.measures.sor_Towy<-beta.measures.sor

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_turn)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_nest)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))


######## Turnover and nestedness Skan ########

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/12S_BLAST")
NORM_WP2_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2_ASV_table_wide_filtered

SampleSite_time<-NORM_WP2_ASV_table_wide_filtered$X1
SampleSite_time<-data.frame(SampleSite_time)
SampleSite_time <- tibble::rownames_to_column(SampleSite_time, "SiteID_1")

adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average","Date_sampled"
                                               #"Arable_and_horticulture_ha","Broadleaf_woodland_ha",	"Coniferous_woodland_ha",
                                               #"Acid_grassland_ha",	"Neutral_grassland_ha",	"Improved_grassland_ha",	"Heather_grassland_ha",	"Heather_ha",
                                               #"Bog_ha",	"Freshwater_ha",	"Inland_rock_ha",	"Supralittoral_sediment_ha",	"Suburban_ha",	"Urban_ha"
)],adonis_data, by="SampleSite_time")

SampleSite_time<-merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake","Days")],SampleSite_time, by="SampleSite_time")

adonis_data_meta<- adonis_data[c(1:17)]
adonis_data<-adonis_data[,-c(1:17)]

adonis_data<-na.omit(adonis_data)

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
#out <- nestedtemp(adonis_data)
#plot(out)
#plot(out, kind="incid")

library(betapart)
#library(qgraph)
adonis_data<-adonis_data %>% mutate_if(is.numeric, ~1 * (. != 0))
adonis_data<-adonis_data[, colSums(adonis_data != 0) > 0]
sorensen.beta<-beta.pair(adonis_data, index.family="sorensen")

nestedness <- as.matrix(sorensen.beta[[3]])
nestedness <- melt(nestedness)
names(nestedness) <- c("SiteID_1", "SiteID_2", "Beta_Distance")

###-----------------------
## get nestedness matrix (this is [[2]]) into a data frame
###-----------------------
m.beta.sorn <- as.matrix(sorensen.beta[[2]])
m.beta.sorn <- melt(m.beta.sorn)
names(m.beta.sorn) <- c("SiteID_1", "SiteID_2", "Beta_nest")
beta.measures.sor <- merge(nestedness, m.beta.sorn, by = c("SiteID_1","SiteID_2")) 

###-----------------------
## get turnover matrix (here simpson, this is [[1]]) into a data frame
###-----------------------
m.beta.sort <- as.matrix(sorensen.beta[[1]])
m.beta.sort  <- melt(m.beta.sort )
names(m.beta.sort ) <- c("SiteID_1", "SiteID_2", "Beta_turn")
beta.measures.sor <- merge(beta.measures.sor, m.beta.sort , by = c("SiteID_1","SiteID_2")) 

beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_1")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_1 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_1 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_1 = Days)

SampleSite_time<-SampleSite_time %>% rename(SiteID_2 = SiteID_1)
beta.measures.sor<-merge(SampleSite_time,beta.measures.sor,by="SiteID_2")
beta.measures.sor<-beta.measures.sor %>% rename(SampleSite_time_2 = SampleSite_time)
beta.measures.sor<-beta.measures.sor %>% rename(Dist_from_lake_2 = Dist_from_lake)
beta.measures.sor<-beta.measures.sor %>% rename(Days_2 = Days)

beta.measures.sor$Dist_from_lake_DIFF<-(beta.measures.sor$Dist_from_lake_2-beta.measures.sor$Dist_from_lake_1)
beta.measures.sor_Skan<-beta.measures.sor

ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_turn)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_nest)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))
ggplot(beta.measures.sor, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_color_gradientn(colours = rainbow(5))+scale_x_continuous(limits = c(0, NA))

#Compiling nestedness for different rivers
beta.measures.sor_all<-rbind(beta.measures.sor_Conwy,beta.measures.sor_Glatt,beta.measures.sor_Gwash,
                             beta.measures.sor_Towy)
                            # beta.measures.sor_Skan)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
beta.measures.sor_all<-beta.measures.sor_all %>% rename(SampleSite_time = SampleSite_time_1)

beta.measures.sor_all<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","River")],beta.measures.sor_all, by="SampleSite_time")

beta.measures.sor_all<-beta.measures.sor_all %>%
  group_by(SampleSite_time,SampleSite_time_2) %>%
  filter(all(Dist_from_lake_DIFF>=0))
#beta.measures.sor_all<-beta.measures.sor_all[beta.measures.sor_all$Dist_from_lake_DIFF != 0, ]


ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_turn,colour=River)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2"))
ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_nest,colour=River)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2"))
ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance,colour=River)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#0072B2"))

