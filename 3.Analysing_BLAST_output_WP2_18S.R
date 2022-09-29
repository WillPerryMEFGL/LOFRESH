library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(cluster)
library(pairwiseAdonis)
library(corrmorant)
library(scales)
library(gam)

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
lofresh_18s_no_cluster_ssu_otus_metazoans <- read_csv("lofresh_18s_no_cluster---ssu---otus_metazoans_filter.csv")

lofresh_18s_no_cluster_ssu_otus_metazoans$taxonomy<-paste(
      lofresh_18s_no_cluster_ssu_otus_metazoans$X10,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X11,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X12,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X13,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X14,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X15,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X16,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X17,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X18,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X19,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X20,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X21,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X22,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X23,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X24,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X25,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X26,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X27,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X28,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X29,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X30,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X31,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X32,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X33,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X34,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X35,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X36,
      lofresh_18s_no_cluster_ssu_otus_metazoans$X37,sep = ";")

#lofresh_18s_no_cluster_ssu_otus_metazoans$taxonomy <- gsub(" ", ";", lofresh_18s_no_cluster_ssu_otus_metazoans$taxonomy)
#lofresh_18s_no_cluster_ssu_otus_metazoans$taxonomy <- gsub(";NA", "", lofresh_18s_no_cluster_ssu_otus_metazoans$taxonomy)

myvars <- c("sequence", "taxonomy")
lofresh_18s_no_cluster_ssu_otus_metazoans <- lofresh_18s_no_cluster_ssu_otus_metazoans[myvars]

seqtabNoC_WP2_18S <- read_table2("seqtabNoC_WP2_18S.txt")
seqtabNoC_WP2_18S_transpose<-t(seqtabNoC_WP2_18S)
seqtabNoC_WP2_18S_transpose<-as.data.frame(seqtabNoC_WP2_18S_transpose)
seqtabNoC_WP2_18S_transpose<-setDT(seqtabNoC_WP2_18S_transpose, keep.rownames = TRUE)
names(seqtabNoC_WP2_18S_transpose)[names(seqtabNoC_WP2_18S_transpose) == "rn"] <- "Seq"
seqtabNoC_WP2_18S_transpose[] <- lapply(seqtabNoC_WP2_18S_transpose, gsub, pattern='"', replacement='')
seqtabNoC_WP2_18S_transpose <- na.omit(transform(seqtabNoC_WP2_18S_transpose, Seq = c("Seq", Seq[-nrow(seqtabNoC_WP2_18S_transpose)])))
header.true <- function(seqtabNoC_WP2_18S_transpose) {
  names(seqtabNoC_WP2_18S_transpose) <- as.character(unlist(seqtabNoC_WP2_18S_transpose[1,]))
  seqtabNoC_WP2_18S_transpose[-1,]
}
seqtabNoC_WP2_18S_transpose<-header.true(seqtabNoC_WP2_18S_transpose)

#rarefaction_curve_all<-t(seqtabNoC_WP2_18S_transpose)
#rarefaction_curve_all= rarefaction_curve_all[-1,]
#rownames(rarefaction_curve_all) <- NULL
#rarefaction_curve_all<-as.data.frame(rarefaction_curve_all)

#rarefaction_curve_all<-lapply(rarefaction_curve_all,as.integer)
#rarefaction_curve_all<-data.frame(rarefaction_curve_all)

#rarecurve_100<-rarefaction_curve_all[rowSums(rarefaction_curve_all[, -1])>100, ]
#rarecurve_100<-rarecurve_100[rowSums(rarecurve_100[, -1])<400000, ]
#rarecurve_100<-na.omit(rarecurve_100)
#(raremax <- min(rowSums(rarecurve_100)))

#rarecurve(rarecurve_100, col = c("red", "blue", "orange", "yellow", "green"), step = 20, sample = raremax, label = FALSE)

#REMOVING SAMPLES WITH LOW READ DEPTH
seqtabNoC_WP2_18S_transpose<-data.frame(seqtabNoC_WP2_18S_transpose)
i <- c(2:1219) 
seqtabNoC_WP2_18S_transpose[ , i] <- apply(seqtabNoC_WP2_18S_transpose[ , i], 2,  
                                     function(x) as.numeric(as.character(x)))

read_depth_original<-colSums(seqtabNoC_WP2_18S_transpose[c(2:1219)])
read_depth_original<-data.frame(read_depth_original)
read_depth_original <- tibble::rownames_to_column(read_depth_original, "sample")


rarefy_col_sum<-colSums(seqtabNoC_WP2_18S_transpose[c(2:1219)])
rarefy_col_sum<-data.frame(rarefy_col_sum)
rarefy_col_sum<-setDT(rarefy_col_sum, keep.rownames = TRUE)[]
rarefy_col_sum<-data.frame(rarefy_col_sum)

rarefy_col_sum<-rarefy_col_sum[!(rarefy_col_sum$rarefy_col_sum>=1000),]
drop<-rarefy_col_sum$rn

seqtabNoC_WP2_18S_transpose<-seqtabNoC_WP2_18S_transpose[,!(names(seqtabNoC_WP2_18S_transpose) %in% drop)]

read_depth_1000_filter<-colSums(seqtabNoC_WP2_18S_transpose[c(2:1143)])
read_depth_1000_filter<-data.frame(read_depth_1000_filter)
read_depth_1000_filter <- tibble::rownames_to_column(read_depth_1000_filter, "sample")
read_depth_track<-merge(read_depth_original[, c("sample", "read_depth_original")],read_depth_1000_filter, by="sample",all = TRUE)

#Add taxonomy to the seqtabNoC_WP2_18S table
lofresh_18s_no_cluster_ssu_otus_metazoans<-lofresh_18s_no_cluster_ssu_otus_metazoans %>% rename(Seq = "sequence")
seqtabNoC_WP2_18S_tax <- lofresh_18s_no_cluster_ssu_otus_metazoans %>% full_join(seqtabNoC_WP2_18S_transpose, by = c("Seq"))

#REMOVING NON-METAZOANS
seqtabNoC_WP2_18S_tax<-na.omit(seqtabNoC_WP2_18S_tax)
seqtabNoC_WP2_18S_tax<-seqtabNoC_WP2_18S_tax %>% rename(species = taxonomy )
seqtabNoC_WP2_18S_tax$Seq<-NULL

seqtabNoC_WP2_18S_tax<-setDT(seqtabNoC_WP2_18S_tax, keep.rownames = TRUE)[]
seqtabNoC_WP2_18S_tax$species<-paste(seqtabNoC_WP2_18S_tax$rn,seqtabNoC_WP2_18S_tax$species,sep=";")
seqtabNoC_WP2_18S_tax$rn<-NULL
seqtabNoC_WP2_18S_tax<-data.frame(seqtabNoC_WP2_18S_tax)
                                                          
i <- c(2:1143) 
seqtabNoC_WP2_18S_tax[ , i] <- apply(seqtabNoC_WP2_18S_tax[ , i], 2,  
                                     function(x) as.numeric(as.character(x)))

read_depth_metazoan_filter<-colSums(seqtabNoC_WP2_18S_tax[c(2:1143)])
read_depth_metazoan_filter<-data.frame(read_depth_metazoan_filter)
read_depth_metazoan_filter <- tibble::rownames_to_column(read_depth_metazoan_filter, "sample")
read_depth_track<-merge(read_depth_metazoan_filter[, c("sample", "read_depth_metazoan_filter")],read_depth_track, by="sample",all = TRUE)

seqtabNoC_WP2_18S_tax<-seqtabNoC_WP2_18S_tax[!grepl(";Primates", seqtabNoC_WP2_18S_tax$species),]

#read_depth_primate_filter<-colSums(seqtabNoC_WP2_18S_tax[c(2:1143)])
#read_depth_primate_filter<-data.frame(read_depth_primate_filter)
#read_depth_primate_filter <- tibble::rownames_to_column(read_depth_primate_filter, "sample")
#read_depth_track<-merge(read_depth_primate_filter[, c("sample", "read_depth_primate_filter")],read_depth_track, by="sample",all = TRUE)

#CALCULATING 0.05% READ THREASHOLD FOR LATER FILTEIRNG
seqtabNoC_WP2_18S_tax_col_sum<-colSums(seqtabNoC_WP2_18S_tax[c(2:1143)])
seqtabNoC_WP2_18S_tax_col_sum<-data.frame(seqtabNoC_WP2_18S_tax_col_sum)
seqtabNoC_WP2_18S_tax_col_sum$read_filter<-seqtabNoC_WP2_18S_tax_col_sum$seqtabNoC_WP2_18S_tax_col_sum *0.0005
seqtabNoC_WP2_18S_tax_col_sum<- seqtabNoC_WP2_18S_tax_col_sum %>% rownames_to_column("ID")

#CONVERITNG ASV TABLE INTO LONG FORMAT
seqtabNoC_WP2_18S_tax_long <- gather(seqtabNoC_WP2_18S_tax, Sample, reads, 2:1143, factor_key=TRUE)
seqtabNoC_WP2_18S_tax_long$reads<-   as.numeric(seqtabNoC_WP2_18S_tax_long$reads)

#READ IN METADATA
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary<-table(lofresh_metadata_WP2_3$Days,lofresh_metadata_WP2_3$SampleSite_code)
lofresh_metadata_WP2_3_summary<-data.frame(lofresh_metadata_WP2_3_summary)
write.csv(lofresh_metadata_WP2_3_summary,"lofresh_metadata_WP2_3_summary.csv")
seqtabNoC_WP2_18S_tax_long<-seqtabNoC_WP2_18S_tax_long %>% rename(ID=Sample)

#COVARIATES IN METADATA PLOT FOR WP2A
library(corrmorant)
lofresh_metadata_WP2A<-lofresh_metadata_WP2_3 %>% filter(WP == "WP2A")
lofresh_metadata_WP2A<-lofresh_metadata_WP2A %>%
  select("season","Days","Dist_from_lake","pH",
"conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
"week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
"temp_monthly_average")
corrmorant(lofresh_metadata_WP2A, style = "blue_red")

#ADD METADATA TO ASV TABLE
seqtabNoC_WP2_18S_tax_long_meta <- merge(lofresh_metadata_WP2_3[, c("ID", "WP")],seqtabNoC_WP2_18S_tax_long, by="ID")
sum_reads<- sum(seqtabNoC_WP2_18S_tax_long_meta$reads)

seqtabNoC_WP2_18S_tax_long_meta_clean_filtered <- merge(seqtabNoC_WP2_18S_tax_col_sum[, c("ID", "read_filter")],seqtabNoC_WP2_18S_tax_long_meta, by="ID")

#REMOVE READS LESS THAN THE 0.005% CUT OFF
seqtabNoC_WP2_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2_18S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2_18S_tax_long_meta_clean_filtered$reads<=seqtabNoC_WP2_18S_tax_long_meta_clean_filtered$read_filter),]

write.csv(seqtabNoC_WP2_18S_tax_long_meta_clean_filtered,"WP2_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2_18S_tax_long_meta_clean_filtered <- read_csv("WP2_ASV_table_long_filtered_family.csv")

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2_18S_tax_wide_filtered<-seqtabNoC_WP2_18S_tax_long_meta_clean_filtered
seqtabNoC_WP2_18S_tax_wide_filtered$Order <- NULL 
seqtabNoC_WP2_18S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2_18S_tax_wide_filtered$Seq_length <- NULL 
seqtabNoC_WP2_18S_tax_wide_filtered$WP <- NULL 

seqtabNoC_WP2_18S_tax_wide_filtered<-seqtabNoC_WP2_18S_tax_wide_filtered %>% 
  group_by(ID,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(ID=ID,species=species, reads= reads,)

seqtabNoC_WP2_18S_tax_wide_filtered <- spread(seqtabNoC_WP2_18S_tax_wide_filtered,ID,reads)
seqtabNoC_WP2_18S_tax_wide_filtered[is.na(seqtabNoC_WP2_18S_tax_wide_filtered)] <- 0
seqtabNoC_WP2_18S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2_18S_tax_wide_filtered))

names(seqtabNoC_WP2_18S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2_18S_tax_wide_filtered[1, ])
seqtabNoC_WP2_18S_tax_wide_filtered <- seqtabNoC_WP2_18S_tax_wide_filtered[-1, ]
seqtabNoC_WP2_18S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2_18S_tax_wide_filtered, function(x) type.convert(as.character(x)))

read_depth_0.05_filter<-colSums(seqtabNoC_WP2_18S_tax[c(2:1143)])
read_depth_0.05_filter<-data.frame(read_depth_0.05_filter)
read_depth_0.05_filter <- tibble::rownames_to_column(read_depth_0.05_filter, "sample")
read_depth_track<-merge(read_depth_0.05_filter[, c("sample", "read_depth_0.05_filter")],read_depth_track, by="sample",all = TRUE)

#tracking read depth plot
read_depth_track_long<-gather(read_depth_track, step,reads, read_depth_0.05_filter:read_depth_1000_filter, factor_key=TRUE)

read_depth_track_long$sample<-gsub('_R1', '', read_depth_track_long$sample)
read_depth_track_long$sample<-gsub('_R2', '', read_depth_track_long$sample)
read_depth_track_long$sample<-gsub('_R3', '', read_depth_track_long$sample)

read_depth_track_long<-read_depth_track_long[- grep("BLANK", read_depth_track_long$sample),]
read_depth_track_long<-read_depth_track_long[- grep("Blank", read_depth_track_long$sample),]
read_depth_track_long<-read_depth_track_long[- grep("C0", read_depth_track_long$sample),]
read_depth_track_long<-read_depth_track_long[- grep("C1", read_depth_track_long$sample),]
read_depth_track_long<-read_depth_track_long[- grep("PC", read_depth_track_long$sample),]


read_depth_track_long<-read_depth_track_long %>% 
  group_by(sample,step) %>% 
  summarize(reads=sum(reads)) %>%
  rename(sample=sample,step=step, reads= reads)

read_depth_track_long$step <- factor(read_depth_track_long$step, levels = c("read_depth_original", "read_depth_1000_filter", "read_depth_metazoan_filter", "read_depth_0.05_filter"))


ggplot(data=read_depth_track_long, aes(x=step, y=reads,fill=step)) +
  geom_bar(stat="identity")+scale_y_continuous(labels = comma)+theme_bw()+facet_wrap(.~sample,nrow=4)+theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks=element_blank(),
    axis.text.x=element_blank(),
    panel.grid=element_blank(),
    panel.spacing.x = unit(0, "null")) +
  labs(x='')


#NORMALISING READS

#YOU MAY NEED TO CHANGE NUMBERS HERE
seqtabNoC_WP2_18S_tax_col_sum2<-rowSums(seqtabNoC_WP2_18S_tax_wide_filtered[c(1:4267)])
seqtabNoC_WP2_18S_tax_col_sum2<-data.frame(seqtabNoC_WP2_18S_tax_col_sum2)
seqtabNoC_WP2_18S_tax_col_sum2$ID <- rownames(seqtabNoC_WP2_18S_tax_col_sum2)

seqtabNoC_WP2_18S_tax_normalised<- merge(seqtabNoC_WP2_18S_tax_col_sum2[, c("ID", "seqtabNoC_WP2_18S_tax_col_sum2")],seqtabNoC_WP2_18S_tax_long_meta_clean_filtered, by="ID")

seqtabNoC_WP2_18S_tax_normalised <- transform(seqtabNoC_WP2_18S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2_18S_tax_col_sum2)

ggplot(seqtabNoC_WP2_18S_tax_normalised , aes(x = ID, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position = "none")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2_18S_tax_normalised_wide<-seqtabNoC_WP2_18S_tax_normalised
seqtabNoC_WP2_18S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$SampleSite_time <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2_18S_tax_normalised_wide$seqtabNoC_WP2_18S_tax_col_sum2 <- NULL 

seqtabNoC_WP2_18S_tax_normalised_wide<-seqtabNoC_WP2_18S_tax_normalised_wide %>% 
  group_by(ID,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(ID=ID,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2_18S_tax_normalised_wide <- spread(seqtabNoC_WP2_18S_tax_normalised_wide,ID,normalised_reads)
seqtabNoC_WP2_18S_tax_normalised_wide[is.na(seqtabNoC_WP2_18S_tax_normalised_wide)] <- 0
seqtabNoC_WP2_18S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2_18S_tax_normalised_wide))

names(seqtabNoC_WP2_18S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2_18S_tax_normalised_wide[1, ])
seqtabNoC_WP2_18S_tax_normalised_wide <- seqtabNoC_WP2_18S_tax_normalised_wide[-1, ]
seqtabNoC_WP2_18S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2_18S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2_18S_tax_normalised_wide,"NORM_WP2_ASV_table_wide_filtered_family.csv")

### WP2A  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP2A_18S_tax_long_meta_clean<-seqtabNoC_WP2_18S_tax_long_meta_clean_filtered %>% filter(WP == "WP2A")
seqtabNoC_WP2A_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2A_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2A_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2A_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2A_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2A_18S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP2A_18S_tax_long_meta_clean,"WP2A_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2A_18S_tax_long_meta_clean <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

### WP2B  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP2B_18S_tax_long_meta_clean<-seqtabNoC_WP2_18S_tax_long_meta_clean_filtered %>% filter(WP == "WP2B")
seqtabNoC_WP2B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2B_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2B_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2B_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2B_18S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP2B_18S_tax_long_meta_clean,"WP2B_ASV_table_long_filtered_family.csv")
seqtabNoC_WP2B_18S_tax_long_meta_clean <- read_csv("WP2B_ASV_table_long_filtered_family.csv")

### WP3B  Splitting work packages and adding metadata ###
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
seqtabNoC_WP3B_18S_tax_long_meta_clean<-seqtabNoC_WP2_18S_tax_long_meta_clean_filtered %>% filter(WP == "WP3B")
seqtabNoC_WP3B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP3B_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP3B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP3B_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP3B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP3B_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP3B_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP3B_18S_tax_long_meta_clean, by="ID")

write.csv(seqtabNoC_WP3B_18S_tax_long_meta_clean,"WP3B_ASV_table_long_filtered_family.csv")
seqtabNoC_WP3B_18S_tax_long_meta_clean <- read_csv("WP3B_ASV_table_long_filtered_family.csv")

### WP2C Splitting work packages and adding metadata ###
seqtabNoC_WP2C_18S_tax_long_meta_clean<-seqtabNoC_WP2_18S_tax_long_meta_clean_filtered %>% filter(WP == "WP2C")
seqtabNoC_WP2C_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_time")],seqtabNoC_WP2C_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "SampleSite_code")],seqtabNoC_WP2C_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Date_sampled")],seqtabNoC_WP2C_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Days")],seqtabNoC_WP2C_18S_tax_long_meta_clean, by="ID")
seqtabNoC_WP2C_18S_tax_long_meta_clean <- merge(lofresh_metadata_WP2_3[, c("ID", "Region")],seqtabNoC_WP2C_18S_tax_long_meta_clean, by="ID")

seqtabNoC_WP2C_18S_tax_long_meta_clean_Gwash<-seqtabNoC_WP2C_18S_tax_long_meta_clean %>% filter(Region == "Leicestershire")
seqtabNoC_WP2C_18S_tax_long_meta_clean_Towy<-seqtabNoC_WP2C_18S_tax_long_meta_clean %>% filter(Region == "Carmarthenshire")
seqtabNoC_WP2C_18S_tax_long_meta_clean_Skan<-seqtabNoC_WP2C_18S_tax_long_meta_clean %>% filter(Region == "New_York")
seqtabNoC_WP2C_18S_tax_long_meta_clean_Glatt<-seqtabNoC_WP2C_18S_tax_long_meta_clean %>% filter(Region == "Switzerland")

write.csv(seqtabNoC_WP2C_18S_tax_long_meta_clean_Gwash,"WP2C_Gwash_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_18S_tax_long_meta_clean_Towy,"WP2C_Towy_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_18S_tax_long_meta_clean_Skan,"WP2C_Skan_ASV_table_long_filtered_family.csv")

write.csv(seqtabNoC_WP2C_18S_tax_long_meta_clean_Glatt,"WP2C_Glatt_ASV_table_long_filtered_family.csv")

#________________________________________________________

#######WP2A BLAST  TAXONOMY ANALYSIS#######

#________________________________________________________

library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(microDecon)

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

#aggregating time points for lentic investigation
#seqtabNoC_WP2A_18S_tax_lentic<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered
#lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
#lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
#seqtabNoC_WP2A_18S_tax_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_18S_tax_lentic, by="SampleSite_time")
#seqtabNoC_WP2A_18S_tax_lentic$SampleSite_time<-NULL
#seqtabNoC_WP2A_18S_tax_lentic<-seqtabNoC_WP2A_18S_tax_lentic %>% 
#  group_by(SampleSite_code,species) %>% 
#  summarize(reads=mean(reads)) %>%
#  rename(SampleSite_code=SampleSite_code,species=species, reads= reads)
#seqtabNoC_WP2A_18S_tax_lentic <- spread(seqtabNoC_WP2A_18S_tax_lentic,SampleSite_code,reads)
#seqtabNoC_WP2A_18S_tax_lentic[is.na(seqtabNoC_WP2A_18S_tax_lentic)] <- 0
#seqtabNoC_WP2A_18S_tax_lentic <- data.frame(t(seqtabNoC_WP2A_18S_tax_lentic))
#names(seqtabNoC_WP2A_18S_tax_lentic) <- as.matrix(seqtabNoC_WP2A_18S_tax_lentic[1, ])
#seqtabNoC_WP2A_18S_tax_lentic <- seqtabNoC_WP2A_18S_tax_lentic[-1, ]
#seqtabNoC_WP2A_18S_tax_lentic[] <- lapply(seqtabNoC_WP2A_18S_tax_lentic, function(x) type.convert(as.character(x)))
#seqtabNoC_WP2A_18S_tax_lentic<-t(seqtabNoC_WP2A_18S_tax_lentic)
#write.csv(seqtabNoC_WP2A_18S_tax_lentic,"seqtabNoC_WP2A_18S_tax_lentic.csv")


#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2A_18S_tax_wide_filtered<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered
seqtabNoC_WP2A_18S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2A_18S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2A_18S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2A_18S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2A_18S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2A_18S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2A_18S_tax_wide_filtered<-seqtabNoC_WP2A_18S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2A_18S_tax_wide_filtered <- spread(seqtabNoC_WP2A_18S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2A_18S_tax_wide_filtered[is.na(seqtabNoC_WP2A_18S_tax_wide_filtered)] <- 0
seqtabNoC_WP2A_18S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2A_18S_tax_wide_filtered))

names(seqtabNoC_WP2A_18S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2A_18S_tax_wide_filtered[1, ])
seqtabNoC_WP2A_18S_tax_wide_filtered <- seqtabNoC_WP2A_18S_tax_wide_filtered[-1, ]
seqtabNoC_WP2A_18S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2A_18S_tax_wide_filtered, function(x) type.convert(as.character(x)))

rarefaction_curve<-seqtabNoC_WP2A_18S_tax_wide_filtered
rownames(rarefaction_curve) <- NULL

#i <- c(1:2587) 
#rarefaction_curve[ , i] <- apply(rarefaction_curve[ , i], 2,  
#                                 function(x) as.integer(as.character(x)))
#(raremax <- min(rowSums(rarefaction_curve)))

#rarecurve(rarefaction_curve, col = c("red", "blue", "orange", "yellow", "green"), step = 20, sample = raremax, label = FALSE)


#FILTERING BY NEGATIVE CONTROL
seqtabNoC_WP2A_18S_tax_wide_filtered<-t(seqtabNoC_WP2A_18S_tax_wide_filtered)
seqtabNoC_WP2A_18S_tax_wide_filtered<-data.frame(seqtabNoC_WP2A_18S_tax_wide_filtered)
seqtabNoC_WP2A_18S_tax_wide_filtered<-setDT(seqtabNoC_WP2A_18S_tax_wide_filtered, keep.rownames = TRUE)[]
seqtabNoC_WP2A_18S_tax_wide_filtered<-seqtabNoC_WP2A_18S_tax_wide_filtered %>% rename(OTU_ID = rn)

colnames_decon<-colnames(seqtabNoC_WP2A_18S_tax_wide_filtered)
colnames_decon<-data.frame(colnames_decon)

seqtabNoC_WP2A_18S_tax_wide_filtered_DECON<-seqtabNoC_WP2A_18S_tax_wide_filtered %>% relocate("OTU_ID","EBlank_T02",
                "EBlank_T14","ENC_T03","ENC_T04","ENC_T05","ENC_T07","ENC_T08","ENC_T09","ENC_T10","ENC_T11","ENC_T12",
                "ENC_T13","ENC_T17", "ENC_T18","ENC_T19")
                
colnames_decon<-colnames(seqtabNoC_WP2A_18S_tax_wide_filtered_DECON)
                
decontaminated <- decon(data = seqtabNoC_WP2A_18S_tax_wide_filtered_DECON,numb.blanks=15,numb.ind=c(17,#E01
                                                                                                    15,#E02
                                                                                                    16,#E03
                                                                                                    18,#E04
                                                                                                    18,#E05
                                                                                                    18,#E06
                                                                                                    18,#E07
                                                                                                    18,#E08
                                                                                                    18,#E09
                                                                                                    18,#E10
                                                                                                    18,#E11
                                                                                                    18,#E12
                                                                                                    18,#E13
                                                                                                    18#E14
                                                                                                    ),taxa=F)
reads_removed<-decontaminated$reads.removed



seqtabNoC_WP2A_18S_tax_wide_filtered<-decontaminated$decon.table
seqtabNoC_WP2A_18S_tax_wide_filtered$Mean.blank<-NULL
seqtabNoC_WP2A_18S_tax_wide_filtered<-t(seqtabNoC_WP2A_18S_tax_wide_filtered)
seqtabNoC_WP2A_18S_tax_wide_filtered<-data.frame(seqtabNoC_WP2A_18S_tax_wide_filtered)

header.true <- function(seqtabNoC_WP2A_18S_tax_wide_filtered) {
  names(seqtabNoC_WP2A_18S_tax_wide_filtered) <- as.character(unlist(seqtabNoC_WP2A_18S_tax_wide_filtered[1,]))
  seqtabNoC_WP2A_18S_tax_wide_filtered[-1,]
}
seqtabNoC_WP2A_18S_tax_wide_filtered<-header.true(seqtabNoC_WP2A_18S_tax_wide_filtered)

i <- c(1:2484) 
seqtabNoC_WP2A_18S_tax_wide_filtered[ , i] <- apply(seqtabNoC_WP2A_18S_tax_wide_filtered[ , i], 2,  
                                           function(x) as.numeric(as.character(x)))

#REMOVE SAMPLES THAT NOW HAVE NO READS IN THEM
seqtabNoC_WP2A_18S_tax_wide_filtered<-seqtabNoC_WP2A_18S_tax_wide_filtered[rowSums(seqtabNoC_WP2A_18S_tax_wide_filtered[])>0,]

#CONVERTING INTO LONG FORMAT
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_18S_tax_wide_filtered
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered <- tibble::rownames_to_column(seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered, "rn")
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered <- gather(seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered, species, reads, 2:2485, factor_key=TRUE)
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered$reads<-   as.numeric(seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered$reads)
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered %>% rename(SampleSite_time = rn)

write.csv(seqtabNoC_WP2A_18S_tax_wide_filtered,"WP2A_ASV_table_wide_filtered.csv")

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2A_18S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

#DENSITY PLOT FOR ALL DATA IN WP2A

#seqtabNoC_WP2A_18S_tax_long_meta_clean_sum_reads<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered %>% 
# group_by(species,SampleSite_time, Days,SampleSite_code) %>% 
# summarize(reads=sum(reads)) %>%
#rename(species_sum=species,SampleSite_time_sum =SampleSite_time,Days_sum=Days,SampleSite_code_sum=SampleSite_code, reads_sum= reads,)

#ggplot(data=seqtabNoC_WP2A_18S_tax_long_meta_clean_sum_reads, aes(x=Days_sum, y=reads_sum,fill=species_sum))+
#  geom_density( stat = "identity",position="stack")+
#  theme_bw()+
#  theme(legend.position="none")+
#  facet_grid(~SampleSite_code_sum)

##SPLITTING BY TAXONOMY
Arthropoda_count <- seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered[grep("Arthropoda", seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered$species), ]
Arthropoda_count<-sum(Arthropoda_count$reads)


######## Venn diagram WP2A ##########
#lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
#lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

#seqtabNoC_WP2A_18S_Venn<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered$reads==0),]
#hist(seqtabNoC_WP2A_18S_Venn$reads, breaks=1000,xlim=c(0,50000))
#seqtabNoC_WP2A_18S_Venn <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "river_section")],seqtabNoC_WP2A_18S_Venn, by="SampleSite_time")

#seqtabNoC_WP2A_18S_Venn_unique<-seqtabNoC_WP2A_18S_Venn
#seqtabNoC_WP2A_18S_Venn_unique$SampleSite_time<-NULL
#seqtabNoC_WP2A_18S_Venn_unique$reads<-NULL
#seqtabNoC_WP2A_18S_Venn_unique<-unique(seqtabNoC_WP2A_18S_Venn_unique)


#upper_Venn<-seqtabNoC_WP2A_18S_Venn %>% filter(river_section == "upper")
#upper_Venn<-upper_Venn$species

#lower_Venn<-seqtabNoC_WP2A_18S_Venn %>% filter(river_section == "lower")
#lower_Venn<-lower_Venn$species

#middle_Venn<-seqtabNoC_WP2A_18S_Venn %>% filter(river_section == "middle")
#middle_Venn<-middle_Venn$species

#library(VennDiagram)
#venn.diagram(
#  x = list(lower_Venn, middle_Venn,upper_Venn),
#  category.names = c("lower" , "middle","upper"),
#  filename = 'WP2A_venn_diagramm.png',
#  output=TRUE)

#PLOT SCALED VENN DIAGRAM
#library(eulerr)
#library(limma)

#VennDiag <- euler(c("lower" = 653, "middle" = 622, "upper" = 331, "lower&middle" = 307, "lower&upper" = 55, 
 #                   "middle&upper" = 155, "lower&middle&upper" = 236))
#plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5,
#     fill=c("darkblue", "deepskyblue2", "cornflowerblue"))

#venn diagram season
#seqtabNoC_WP2A_18S_Venn<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered$reads==0),]
#seqtabNoC_WP2A_18S_Venn <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "season")],seqtabNoC_WP2A_18S_Venn, by="SampleSite_time")

#Autumn_Venn<-seqtabNoC_WP2A_18S_Venn %>% filter(season == "Autumn")
#Autumn_Venn<-Autumn_Venn$species

#Winter_Venn<-seqtabNoC_WP2A_18S_Venn %>% filter(season == "Winter")
#Winter_Venn<-Winter_Venn$species

#Spring_Venn<-seqtabNoC_WP2A_18S_Venn %>% filter(season == "Spring")
#Spring_Venn<-Spring_Venn$species

#Summer_Venn<-seqtabNoC_WP2A_18S_Venn %>% filter(season == "Summer")
#Summer_Venn<-Summer_Venn$species

#library(VennDiagram)
#venn.diagram(
 # x = list(Autumn_Venn, Winter_Venn,Spring_Venn,Summer_Venn),
#  category.names = c("Autumn" , "Winter","Spring","Summer"),
#  filename = 'WP2A_venn_diagramm.png',
#  output=TRUE)

#PLOT SCALED VENN DIAGRAM
#library(eulerr)
#library(limma)

#VennDiag <- euler(c("lower" = 653, "middle" = 622, "upper" = 331, "lower&middle" = 307, "lower&upper" = 55, 
 #                   "middle&upper" = 155, "lower&middle&upper" = 236))
#plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5,
 #    fill=c("darkblue", "deepskyblue2", "cornflowerblue"))

#VennDiag <- euler(c("Autumn"=335, "Winter"=552,"Spring"=301,"Summer"=328,
#                    "Autumn&Summer"=91,"Autumn&Winter"=93,
#                    "Spring&Autumn"=46,"Winter&Summer"=45,"Spring&Summer"=51,
#                    "Spring&Winter"=122,"Autumn&Spring&Summer"=43,"Summer&Winter&Spring"=43,
#                    "Autumn&Summer&Winter"=44,"Winter&Spring&Autumn"=74,
#                    "Autumn&Summer&Winter&Spring"=191))
#plot(VennDiag, counts = TRUE, font=1, cex=1, alpha=0.5,
#     fill=c("darkblue", "deepskyblue2", "cornflowerblue","steelblue4"))



#NORMALISING READS
seqtabNoC_WP2A_18S_tax_col_sum2<-rowSums(seqtabNoC_WP2A_18S_tax_wide_filtered[c(1:2484)])
seqtabNoC_WP2A_18S_tax_col_sum2<-data.frame(seqtabNoC_WP2A_18S_tax_col_sum2)
seqtabNoC_WP2A_18S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2A_18S_tax_col_sum2)

seqtabNoC_WP2A_18S_tax_normalised<- merge(seqtabNoC_WP2A_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_18S_tax_col_sum2")],seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2A_18S_tax_col_sum3<-colSums(seqtabNoC_WP2A_18S_tax_wide_filtered[c(1:2484)])
seqtabNoC_WP2A_18S_tax_col_sum3<-data.frame(seqtabNoC_WP2A_18S_tax_col_sum3)

seqtabNoC_WP2A_18S_tax_col_sum3$species <- rownames(seqtabNoC_WP2A_18S_tax_col_sum3)

seqtabNoC_WP2A_18S_tax_normalised<- merge(seqtabNoC_WP2A_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_18S_tax_col_sum2")],seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered, by="SampleSite_time")

seqtabNoC_WP2A_18S_tax_top<- merge(seqtabNoC_WP2A_18S_tax_col_sum3[, c("species", "seqtabNoC_WP2A_18S_tax_col_sum3")],seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2A_18S_tax_normalised <- transform(seqtabNoC_WP2A_18S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2A_18S_tax_col_sum2)

ggplot(seqtabNoC_WP2A_18S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2A_18S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Days")],seqtabNoC_WP2A_18S_tax_normalised, by="SampleSite_time")
seqtabNoC_WP2A_18S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_18S_tax_normalised, by="SampleSite_time")

#Looking at individual species
#Salmo_salar_norm<-seqtabNoC_WP2A_18S_tax_normalised %>% filter(species == "Salmo_salar")

#Salmo_salar_norm_sum_reads<-Salmo_salar_norm %>% 
 # group_by(SampleSite_time, Days,SampleSite_code) %>% 
  #summarize(normalised_reads=sum(normalised_reads)) %>%
  #rename(SampleSite_time_sum =SampleSite_time,Days_sum=Days,SampleSite_code_sum=SampleSite_code,normalised_reads_sum= normalised_reads)

#Salmo_salar_norm_sum_reads$Days_sum<-as.numeric(Salmo_salar_norm_sum_reads$Days_sum)

#ggplot(data=Salmo_salar_norm_sum_reads, aes(x=Days_sum, y=normalised_reads_sum))+
 # facet_grid(~ SampleSite_code_sum)+theme_bw()+geom_bar(stat='identity')+ 
  #theme(axis.text.x = element_text(angle = 90))

#PLOTTING TOP Taxon
#seqtabNoC_WP2A_18S_tax_top %>%
#  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2A_18S_tax_col_sum3))) %>%
#  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
#  geom_bar(position="fill",stat = "identity", colour = "black")+
#  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+scale_fill_brewer(palette="Paired")

write.csv(seqtabNoC_WP2A_18S_tax_normalised,"NORM_WP2A_ASV_table_LONG_filtered.csv")

#PLOTTING PHYLA
seqtabNoC_WP2A_18S_tax_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2A_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2A_18S_tax_phyla<-cbind(seqtabNoC_WP2A_18S_tax_top,seqtabNoC_WP2A_18S_tax_phyla)

seqtabNoC_WP2A_18S_tax_phyla<-seqtabNoC_WP2A_18S_tax_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2A_18S_tax_col_sum3=sum(seqtabNoC_WP2A_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2A_18S_tax_col_sum3=seqtabNoC_WP2A_18S_tax_col_sum3)

seqtabNoC_WP2A_18S_tax_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_18S_tax_phyla, by="SampleSite_time")
unique(seqtabNoC_WP2A_18S_tax_phyla[c("X2")])
seqtabNoC_WP2A_18S_tax_phyla$X2<-factor(seqtabNoC_WP2A_18S_tax_phyla$X2,levels= c('Arthropoda','Rotifera','Mollusca','Chordata',
                                                                                  'Annelida','Nematoda','Gastrotricha', 'Tardigrada','Cnidaria','Porifera','Echinodermata','Nemertea',
                                                                                  'Xenacoelomorpha',
                                                                                  'Nematomorpha', 'Acanthocephala <thorny-headed worm>','Chaetognatha','Hemichordata',
                                                                                  'Ctenophora <comb jellyfish phylum>','Bryozoa','Platyhelminthes', 'Micrognathozoa'),ordered=TRUE)                                                                                         

                                                                                                        
ggplot(data=seqtabNoC_WP2A_18S_tax_phyla,aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90))+facet_wrap(~ SampleSite_code,scales="free",ncol = 15)+scale_fill_manual(values = c("#A6CEE3",#Arthropods
                                                                             "#1F78B4", #Rotifera
                                                                             "#B2DF8A",#Molluscs
                                                                             "#33A02C",#Chordates
                                                                             "#FB9A99",#annelids
                                                                             "#E31A1C",#nematodes
                                                                             "#FDBF6F",#gastrotricha
                                                                             "#FF7F00",#tardigrada
                                                                             "#CAB2D6",#Cnideria
                                                                             "#6A3D9A",#porifera
                                                                             "#FFFF99",#echinodermata
                                                                             "#B15928",#nemertea,
                                                                             "white",
                                                                             "white","white","white","white","white","white","white","white"))
#PLOTTING TOP PHYLA
seqtabNoC_WP2A_18S_tax_top_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2A_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2A_18S_tax_top_phyla<-cbind(seqtabNoC_WP2A_18S_tax_top,seqtabNoC_WP2A_18S_tax_top_phyla)

seqtabNoC_WP2A_18S_tax_top_phyla<-seqtabNoC_WP2A_18S_tax_top_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2A_18S_tax_col_sum3=sum(seqtabNoC_WP2A_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2A_18S_tax_col_sum3=seqtabNoC_WP2A_18S_tax_col_sum3)

seqtabNoC_WP2A_18S_tax_top_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2A_18S_tax_top_phyla, by="SampleSite_time")

#seqtabNoC_WP2A_18S_tax_top_phyla %>%
  #mutate(X2 = fct_reorder(X2, desc(seqtabNoC_WP2A_18S_tax_col_sum3))) %>%
#  ggplot(aes(x = SampleSite_time, fill = X2, y = reads)) + 
#  geom_bar(position="fill",stat = "identity", colour = "black")+
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank())+
 # theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")+facet_wrap(~ SampleSite_code,scales="free",ncol = 15)

seqtabNoC_WP2A_18S_tax_top_phyla<-seqtabNoC_WP2A_18S_tax_top_phyla %>% 
  group_by(X2) %>% 
  summarize(reads=sum(reads)) %>%
  rename(X2=X2, reads= reads)



#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2A_18S_tax_normalised_wide<-seqtabNoC_WP2A_18S_tax_normalised
seqtabNoC_WP2A_18S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$seqtabNoC_WP2_18S_tax_col_sum2 <- NULL 

seqtabNoC_WP2A_18S_tax_normalised_wide<-seqtabNoC_WP2A_18S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2A_18S_tax_normalised_wide <- spread(seqtabNoC_WP2A_18S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2A_18S_tax_normalised_wide[is.na(seqtabNoC_WP2A_18S_tax_normalised_wide)] <- 0
seqtabNoC_WP2A_18S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2A_18S_tax_normalised_wide))

names(seqtabNoC_WP2A_18S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2A_18S_tax_normalised_wide[1, ])
seqtabNoC_WP2A_18S_tax_normalised_wide <- seqtabNoC_WP2A_18S_tax_normalised_wide[-1, ]
seqtabNoC_WP2A_18S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2A_18S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2A_18S_tax_normalised_wide,"NORM_WP2A_ASV_table_wide_filtered.csv")

#Investigating lentic species by time point
#seqtabNoC_WP2A_18S_tax_normalised_wide_transposed<-t(seqtabNoC_WP2A_18S_tax_normalised_wide)
#seqtabNoC_WP2A_18S_tax_normalised_wide_transposed<-data.frame(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed)
#lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),] 

#library(ggpubr)
#T01
#T01_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T01==0),]
#T01_lentic <- T01_lentic %>% select(contains("T01"))
#T01_lentic<-setDT(T01_lentic, keep.rownames = TRUE)[]
#T01_lentic <- gather(T01_lentic, SampleSite_time, normalised_reads, 2:14, factor_key=TRUE)
#T01_lentic<-data.frame(T01_lentic)
#T01_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T01_lentic, by="SampleSite_time")

#p01<-ggplot(data=T01_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")
#ggplotly(p01)
  
#T02  
#T02_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T02==0),]
#T02_lentic <- T02_lentic %>% select(contains("T02"))
#T02_lentic<-setDT(T02_lentic, keep.rownames = TRUE)[]
#T02_lentic <- gather(T02_lentic, SampleSite_time, normalised_reads, 2:14, factor_key=TRUE)
#T02_lentic<-data.frame(T02_lentic)
#T02_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T02_lentic, by="SampleSite_time")
#T02_lentic<-na.omit(T02_lentic)

#p02<-ggplot(data=T02_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T03  
#T03_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T03==0),]
#T03_lentic <- T03_lentic %>% select(contains("T03"))
#T03_lentic<-setDT(T03_lentic, keep.rownames = TRUE)[]
#T03_lentic <- gather(T03_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T03_lentic<-data.frame(T03_lentic)
#T03_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T03_lentic, by="SampleSite_time")
#T03_lentic<-na.omit(T03_lentic)

#p03<-ggplot(data=T03_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T04  
#T04_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T04==0),]
#T04_lentic <- T04_lentic %>% select(contains("T04"))
#T04_lentic<-setDT(T04_lentic, keep.rownames = TRUE)[]
#T04_lentic <- gather(T04_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T04_lentic<-data.frame(T04_lentic)
#T04_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T04_lentic, by="SampleSite_time")
#T04_lentic<-na.omit(T04_lentic)

#p04<-ggplot(data=T04_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T05  
#T05_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T05==0),]
#T05_lentic <- T05_lentic %>% select(contains("T05"))
#T05_lentic<-setDT(T05_lentic, keep.rownames = TRUE)[]
#T05_lentic <- gather(T05_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T05_lentic<-data.frame(T05_lentic)
#T05_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T05_lentic, by="SampleSite_time")
#T05_lentic<-na.omit(T05_lentic)

#p05<-ggplot(data=T05_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T06  
#T06_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T06==0),]
#T06_lentic <- T06_lentic %>% select(contains("T06"))
#T06_lentic<-setDT(T06_lentic, keep.rownames = TRUE)[]
#T06_lentic <- gather(T06_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T06_lentic<-data.frame(T06_lentic)
#T06_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T06_lentic, by="SampleSite_time")

#p06<-ggplot(data=T06_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#NO DATA FOR T07

#T08  
#T08_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T08==0),]
#T08_lentic <- T08_lentic %>% select(contains("T08"))
#T08_lentic<-setDT(T08_lentic, keep.rownames = TRUE)[]
#T08_lentic <- gather(T08_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T08_lentic<-data.frame(T08_lentic)
#T08_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T08_lentic, by="SampleSite_time")
#T08_lentic<-na.omit(T08_lentic)

#p08<-ggplot(data=T08_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T09  
#T09_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T09==0),]
#T09_lentic <- T09_lentic %>% select(contains("T09"))
#T09_lentic<-setDT(T09_lentic, keep.rownames = TRUE)[]
#T09_lentic <- gather(T09_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T09_lentic<-data.frame(T09_lentic)
#T09_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T09_lentic, by="SampleSite_time")
#T09_lentic<-na.omit(T09_lentic)

#p09<-ggplot(data=T09_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T10  
#T10_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T10==0),]
#T10_lentic <- T10_lentic %>% select(contains("T10"))
#T10_lentic<-setDT(T10_lentic, keep.rownames = TRUE)[]
#T10_lentic <- gather(T10_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T10_lentic<-data.frame(T10_lentic)
#T10_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T10_lentic, by="SampleSite_time")
#T10_lentic<-na.omit(T10_lentic)

#p10<-ggplot(data=T10_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
  #geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T11  
#T11_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T11==0),]
#T11_lentic <- T11_lentic %>% select(contains("T11"))
#T11_lentic<-setDT(T11_lentic, keep.rownames = TRUE)[]
#T11_lentic <- gather(T11_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T11_lentic<-data.frame(T11_lentic)
#T11_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T11_lentic, by="SampleSite_time")
#T11_lentic<-na.omit(T11_lentic)

#p11<-ggplot(data=T11_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
 # geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T12  
#T12_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T12==0),]
#T12_lentic <- T12_lentic %>% select(contains("T12"))
#T12_lentic<-setDT(T12_lentic, keep.rownames = TRUE)[]
#T12_lentic <- gather(T12_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T12_lentic<-data.frame(T12_lentic)
#T12_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T12_lentic, by="SampleSite_time")
#T12_lentic<-na.omit(T12_lentic)

#p12<-ggplot(data=T12_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T13  
#T13_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T13==0),]
#T13_lentic <- T13_lentic %>% select(contains("T13"))
#T13_lentic<-setDT(T13_lentic, keep.rownames = TRUE)[]
#T13_lentic <- gather(T13_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T13_lentic<-data.frame(T13_lentic)
#T13_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T13_lentic, by="SampleSite_time")
#T13_lentic<-na.omit(T13_lentic)

#p13<-ggplot(data=T13_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T14  
#T14_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T14==0),]
#T14_lentic <- T14_lentic %>% select(contains("T14"))
#T14_lentic<-setDT(T14_lentic, keep.rownames = TRUE)[]
#T14_lentic <- gather(T14_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T14_lentic<-data.frame(T14_lentic)
#T14_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T14_lentic, by="SampleSite_time")
#T14_lentic<-na.omit(T14_lentic)

#p14<-ggplot(data=T14_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
 # geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T15  
#T15_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T15==0),]
#T15_lentic <- T15_lentic %>% select(contains("T15"))
#T15_lentic<-setDT(T15_lentic, keep.rownames = TRUE)[]
#T15_lentic <- gather(T15_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T15_lentic<-data.frame(T15_lentic)
#T15_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T15_lentic, by="SampleSite_time")

#p15<-ggplot(data=T15_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T17  
#T17_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T17==0),]
#T17_lentic <- T17_lentic %>% select(contains("T17"))
#T17_lentic<-setDT(T17_lentic, keep.rownames = TRUE)[]
#T17_lentic <- gather(T17_lentic, SampleSite_time, normalised_reads, 2:15, factor_key=TRUE)
#T17_lentic<-data.frame(T17_lentic)
#T17_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T17_lentic, by="SampleSite_time")
#T17_lentic<-na.omit(T17_lentic)

#p17<-ggplot(data=T17_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
 # geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T18  
#T18_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T18==0),]
#T18_lentic <- T18_lentic %>% select(contains("T18"))
#T18_lentic<-setDT(T18_lentic, keep.rownames = TRUE)[]
#T18_lentic <- gather(T18_lentic, SampleSite_time, normalised_reads, 2:13, factor_key=TRUE)
#T18_lentic<-data.frame(T18_lentic)
#T18_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T18_lentic, by="SampleSite_time")
#T18_lentic<-na.omit(T18_lentic)

#p18<-ggplot(data=T18_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")

#T19  
#T19_lentic<-seqtabNoC_WP2A_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2A_18S_tax_normalised_wide_transposed$E01_T19==0),]
#T19_lentic <- T19_lentic %>% select(contains("T19"))
#T19_lentic<-setDT(T19_lentic, keep.rownames = TRUE)[]
#T19_lentic <- gather(T19_lentic, SampleSite_time, normalised_reads, 2:14, factor_key=TRUE)
#T19_lentic<-data.frame(T19_lentic)
#T19_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake","timepoint")],T19_lentic, by="SampleSite_time")
#T19_lentic<-na.omit(T19_lentic)

#p19<-ggplot(data=T19_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
#  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")


#ggarrange(p01,p02,p03,p04,p05,p06,p08,p09,p10,p11,p12,p13,p14,p15,p17,p18,p19, #17plots
#          ncol = 3, nrow = 6)

#lentic_all<-rbind(T01_lentic,T02_lentic,T03_lentic,T04_lentic,T05_lentic,T06_lentic,T08_lentic,T09_lentic,T10_lentic,T11_lentic,T12_lentic,T13_lentic,T14_lentic,T15_lentic,T17_lentic,T18_lentic,T19_lentic)
#write.csv(lentic_all,"lentic_conwy.csv")

###### PERMANOVA ON ALL TAXON ######
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data<-adonis_data[adonis_data$X1 != "E04_T07", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T19", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E02_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E14_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T09", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T18", ]  

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

adonis_data<-na.omit(adonis_data)

adonis_data_meta<- adonis_data[c(1:16)]
adonis_data<-adonis_data[,-c(1:16)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)

adon.results_WP2<-adonis2(adonis_data ~ adonis_data_meta$Dist_from_lake+ 
                            adonis_data_meta$season+
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
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Blank","River","SampleSite_code","Days","Dist_from_lake","season")],data.scores, by="SampleSite_time")


p<-
  ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)

gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~season,data.scores_meta_data,mean),by="season")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=season,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
  geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))+
  geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")

#Canonical Correspondence Analysis (CCA)
ccamodel_adonis_data_meta<-adonis_data_meta
ccamodel_adonis_data_meta$SampleSite_time<-NULL

ccamodel<-cca(adonis_data ~.,ccamodel_adonis_data_meta)

ccamodel_adonis_data_meta$season<-as.factor(ccamodel_adonis_data_meta$season)
ccamodel_adonis_data_meta$SampleSite_code<-as.factor(ccamodel_adonis_data_meta$SampleSite_code)
ccamodel_adonis_data_meta$river_section<-as.factor(ccamodel_adonis_data_meta$river_section)

plot(ccamodel,scaling="sites")
points(ccamodel, pch=21, cex=1.2, scaling="sites")
with(ccamodel_adonis_data_meta, points(ccamodel, scaling="sites", pch=21, col=1, bg=river_section))
with(ccamodel_adonis_data_meta, legend("bottomleft", levels(river_section), pch=21, pt.bg=1:12, bty="n"))


###### PERMANOVA ON ARTHROPODS ######
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data_Arthropoda<-NORM_WP2A_ASV_table_wide_filtered

#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E12_T03", ]  
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E05_T10", ]  
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E04_T05", ]  
adonis_data_Arthropoda<-adonis_data_Arthropoda[adonis_data_Arthropoda$X1 != "E06_T10", ]  

adonis_data_Arthropoda<-adonis_data_Arthropoda %>% rename(SampleSite_time_Arthropoda = X1)
adonis_data_Arthropoda <- adonis_data_Arthropoda[,grepl("Arthropoda", colnames(adonis_data_Arthropoda))]
adonis_data_Arthropoda<-adonis_data_Arthropoda[, colSums(adonis_data_Arthropoda != 0) > 0]
adonis_data_Arthropoda<-adonis_data_Arthropoda[apply(adonis_data_Arthropoda[,-1], 1, function(x) !all(x==0)),]

adonis_data_Arthropoda_samples<-adonis_data_Arthropoda$SampleSite_time_Arthropoda
adonis_data_Arthropoda_samples<-data.frame(adonis_data_Arthropoda_samples)
adonis_data_Arthropoda_samples<-adonis_data_Arthropoda_samples %>% rename(ID = adonis_data_Arthropoda_samples)
adonis_data_Arthropoda<-adonis_data_Arthropoda %>% rename(SampleSite_time = SampleSite_time_Arthropoda)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data_Arthropoda<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                               "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                               "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                               "temp_monthly_average")],adonis_data_Arthropoda, by="SampleSite_time")
adonis_data_Arthropoda<-na.omit(adonis_data_Arthropoda)

adonis_data_Arthropoda_meta<- adonis_data_Arthropoda[c(1:16)]
adonis_data_Arthropoda<-adonis_data_Arthropoda[,-c(1:16)]

adonis_data_Arthropoda_meta$Days<-as.numeric(adonis_data_Arthropoda_meta$Days)

adon.results_WP2_Arthropoda<-adonis2(adonis_data_Arthropoda ~ 
                            adonis_data_Arthropoda_meta$season/adonis_data_Arthropoda_meta$Days+
                            adonis_data_Arthropoda_meta$Dist_from_lake+
                            adonis_data_Arthropoda_meta$daily_flow_by_catchment+
                            adonis_data_Arthropoda_meta$monthly_flow_by_catchment+
                            adonis_data_Arthropoda_meta$seasonal_flow_by_catchment+
                            adonis_data_Arthropoda_meta$week_before_flow_by_catchment+
                            adonis_data_Arthropoda_meta$river_width_m+
                            adonis_data_Arthropoda_meta$river_gradient_percent+
                            adonis_data_Arthropoda_meta$rainfall_monthly_average+
                            adonis_data_Arthropoda_meta$temp_monthly_average+
                            adonis_data_Arthropoda_meta$gran_alk_ueq_L+
                            adonis_data_Arthropoda_meta$conductivity_uS_cm+
                            adonis_data_Arthropoda_meta$pH,
                          method="bray",perm=999,by="margin")
print(adon.results_WP2_Arthropoda)

#PLOTTING BETA DIVERSITY NMDS
com_Arthropoda <- adonis_data_Arthropoda[,col(adonis_data_Arthropoda)]
m_com_Arthropoda <- as.matrix(com_Arthropoda)
set.seed(123)
nmds_Arthropoda = metaMDS(m_com_Arthropoda, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores_Arthropoda = as.data.frame(scores(nmds_Arthropoda))
#add columns to data frame 
data.scores_Arthropoda$SampleSite_time = adonis_data_Arthropoda_meta[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_Arthropoda_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Blank","SampleSite_code","Days","Dist_from_lake")],data.scores_Arthropoda, by="SampleSite_time")

p<-
  ggplot(data.scores_Arthropoda_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#Canonical Correspondence Analysis (CCA) ARTHROPODA
ccamodel_adonis_data_Arthropoda_meta<-adonis_data_Arthropoda_meta
ccamodel_adonis_data_Arthropoda_meta$SampleSite_time<-NULL

ccamodel_Arthropoda<-cca(adonis_data_Arthropoda ~.,ccamodel_adonis_data_Arthropoda_meta)

ccamodel_adonis_data_Arthropoda_meta$season<-as.factor(ccamodel_adonis_data_Arthropoda_meta$season)
ccamodel_adonis_data_Arthropoda_meta$SampleSite_code<-as.factor(ccamodel_adonis_data_Arthropoda_meta$SampleSite_code)
ccamodel_adonis_data_Arthropoda_meta$river_section<-as.factor(ccamodel_adonis_data_Arthropoda_meta$river_section)

plot(ccamodel_Arthropoda,scaling="sites")
points(ccamodel_Arthropoda, pch=21, cex=1.2, scaling="sites")
with(ccamodel_adonis_data_Arthropoda_meta, points(ccamodel_Arthropoda, scaling="sites", pch=21, col=1, bg=river_section))
with(ccamodel_adonis_data_Arthropoda_meta, legend("bottomleft", levels(river_section), pch=21, pt.bg=1:12, bty="n"))

###### PERMANOVA ON ROTIFERA #######
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")

NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data_Rotifera<-NORM_WP2A_ASV_table_wide_filtered

#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data_Rotifera<-adonis_data_Rotifera[adonis_data_Rotifera$X1 != "E02_T01", ]  
adonis_data_Rotifera<-adonis_data_Rotifera[adonis_data_Rotifera$X1 != "E01_T19", ]  

adonis_data_Rotifera<-adonis_data_Rotifera %>% rename(SampleSite_time_Rotifera = X1)
adonis_data_Rotifera <- adonis_data_Rotifera[,grepl("Rotifera", colnames(adonis_data_Rotifera))]
adonis_data_Rotifera<-adonis_data_Rotifera[, colSums(adonis_data_Rotifera != 0) > 0]
adonis_data_Rotifera<-adonis_data_Rotifera[apply(adonis_data_Rotifera[,-1], 1, function(x) !all(x==0)),]

adonis_data_Rotifera_samples<-adonis_data_Rotifera$SampleSite_time_Rotifera
adonis_data_Rotifera_samples<-data.frame(adonis_data_Rotifera_samples)
adonis_data_Rotifera_samples<-adonis_data_Rotifera_samples %>% rename(ID = adonis_data_Rotifera_samples)
adonis_data_Rotifera<-adonis_data_Rotifera %>% rename(SampleSite_time = SampleSite_time_Rotifera)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data_Rotifera<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","season","Days","Dist_from_lake","pH",
                                                          "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                          "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                          "temp_monthly_average","river_section")],adonis_data_Rotifera, by="SampleSite_time")
adonis_data_Rotifera<-na.omit(adonis_data_Rotifera)

adonis_data_Rotifera_meta<- adonis_data_Rotifera[c(1:16)]
adonis_data_Rotifera<-adonis_data_Rotifera[,-c(1:16)]

adonis_data_Rotifera_meta$Days<-as.numeric(adonis_data_Rotifera_meta$Days)

adon.results_WP2_Rotifera<-adonis2(adonis_data_Rotifera ~ adonis_data_Rotifera_meta$Dist_from_lake+ adonis_data_Rotifera_meta$Days+
                                       adonis_data_Rotifera_meta$season/adonis_data_Rotifera_meta$Days+
                                       adonis_data_Rotifera_meta$Dist_from_lake+
                                       adonis_data_Rotifera_meta$daily_flow_by_catchment+
                                       adonis_data_Rotifera_meta$monthly_flow_by_catchment+
                                       adonis_data_Rotifera_meta$seasonal_flow_by_catchment+
                                       adonis_data_Rotifera_meta$week_before_flow_by_catchment+
                                       adonis_data_Rotifera_meta$river_width_m+
                                       adonis_data_Rotifera_meta$river_gradient_percent+
                                       adonis_data_Rotifera_meta$rainfall_monthly_average+
                                       adonis_data_Rotifera_meta$temp_monthly_average+
                                       adonis_data_Rotifera_meta$gran_alk_ueq_L+
                                       adonis_data_Rotifera_meta$conductivity_uS_cm+
                                       adonis_data_Rotifera_meta$pH,
                                     method="bray",perm=999,by="margin")
print(adon.results_WP2_Rotifera)

#PLOTTING BETA DIVERSITY NMDS
com_Rotifera <- adonis_data_Rotifera[,col(adonis_data_Rotifera)]
m_com_Rotifera <- as.matrix(com_Rotifera)
set.seed(123)
nmds_Rotifera = metaMDS(m_com_Rotifera, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores_Rotifera = as.data.frame(scores(nmds_Rotifera))
#add columns to data frame 
data.scores_Rotifera$SampleSite_time = adonis_data_Rotifera_meta[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_Rotifera_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Blank","SampleSite_code","Days","Dist_from_lake")],data.scores_Rotifera, by="SampleSite_time")

p<-
  ggplot(data.scores_Rotifera_meta_data, aes(x = NMDS1, y = NMDS2,colour = SampleSite_time)) + 
  geom_point(size = 2, aes(colour = SampleSite_time))+
  theme_bw()
ggplotly(p)

p<-
  ggplot(data.scores_Rotifera_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)

#Canonical Correspondence Analysis (CCA) ROTIFERA
ccamodel_adonis_data_Rotifera_meta<-adonis_data_Rotifera_meta
ccamodel_adonis_data_Rotifera_meta$SampleSite_time<-NULL

ccamodel_Rotifera<-cca(adonis_data_Rotifera ~.,ccamodel_adonis_data_Rotifera_meta)

ccamodel_adonis_data_Rotifera_meta$season<-as.factor(ccamodel_adonis_data_Rotifera_meta$season)
ccamodel_adonis_data_Rotifera_meta$SampleSite_code<-as.factor(ccamodel_adonis_data_Rotifera_meta$SampleSite_code)
ccamodel_adonis_data_Rotifera_meta$river_section<-as.factor(ccamodel_adonis_data_Rotifera_meta$river_section)

plot(ccamodel_Rotifera,scaling="sites")
points(ccamodel_Rotifera, pch=21, cex=1.2, scaling="sites")
with(ccamodel_adonis_data_Rotifera_meta, points(ccamodel_Rotifera, scaling="sites", pch=21, col=1, bg=river_section))
with(ccamodel_adonis_data_Rotifera_meta, legend("bottomleft", levels(river_section), pch=21, pt.bg=1:12, bty="n"))

###### PERMANOVA ON Nematoda #######

NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data_Nematoda<-NORM_WP2A_ASV_table_wide_filtered

#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E08_T03", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E10_T04", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E02_T03", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E12_T04", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E12_T08", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E11_T10", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E09_T04", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E08_T02", ]  
adonis_data_Nematoda<-adonis_data_Nematoda[adonis_data_Nematoda$X1 != "E06_T10", ]  

adonis_data_Nematoda<-adonis_data_Nematoda %>% rename(SampleSite_time_Nematoda = X1)
adonis_data_Nematoda <- adonis_data_Nematoda[,grepl("Nematoda", colnames(adonis_data_Nematoda))]
adonis_data_Nematoda<-adonis_data_Nematoda[, colSums(adonis_data_Nematoda != 0) > 0]
adonis_data_Nematoda<-adonis_data_Nematoda[apply(adonis_data_Nematoda[,-1], 1, function(x) !all(x==0)),]

adonis_data_Nematoda_samples<-adonis_data_Nematoda$SampleSite_time_Nematoda
adonis_data_Nematoda_samples<-data.frame(adonis_data_Nematoda_samples)
adonis_data_Nematoda_samples<-adonis_data_Nematoda_samples %>% rename(ID = adonis_data_Nematoda_samples)
adonis_data_Nematoda<-adonis_data_Nematoda %>% rename(SampleSite_time = SampleSite_time_Nematoda)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data_Nematoda<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                                        "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                        "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                        "temp_monthly_average")],adonis_data_Nematoda, by="SampleSite_time")
adonis_data_Nematoda<-na.omit(adonis_data_Nematoda)

adonis_data_Nematoda_meta<- adonis_data_Nematoda[c(1:16)]
adonis_data_Nematoda<-adonis_data_Nematoda[,-c(1:16)]

adonis_data_Nematoda_meta$Days<-as.numeric(adonis_data_Nematoda_meta$Days)

adon.results_WP2_Nematoda<-adonis2(adonis_data_Nematoda ~ adonis_data_Nematoda_meta$Dist_from_lake+ adonis_data_Nematoda_meta$Days+
                                     adonis_data_Nematoda_meta$season/adonis_data_Nematoda_meta$Days+
                                     adonis_data_Nematoda_meta$Dist_from_lake+
                                     adonis_data_Nematoda_meta$daily_flow_by_catchment+
                                     adonis_data_Nematoda_meta$monthly_flow_by_catchment+
                                     adonis_data_Nematoda_meta$seasonal_flow_by_catchment+
                                     adonis_data_Nematoda_meta$week_before_flow_by_catchment+
                                     adonis_data_Nematoda_meta$river_width_m+
                                     adonis_data_Nematoda_meta$river_gradient_percent+
                                     adonis_data_Nematoda_meta$rainfall_monthly_average+
                                     adonis_data_Nematoda_meta$temp_monthly_average+
                                     adonis_data_Nematoda_meta$gran_alk_ueq_L+
                                     adonis_data_Nematoda_meta$conductivity_uS_cm+
                                     adonis_data_Nematoda_meta$pH,
                                   method="bray",perm=999,by="margin")
print(adon.results_WP2_Nematoda)

#Canonical Correspondence Analysis (CCA) NEMATODA
ccamodel_adonis_data_Nematoda_meta<-adonis_data_Nematoda_meta
ccamodel_adonis_data_Nematoda_meta$SampleSite_time<-NULL

ccamodel_Nematoda<-cca(adonis_data_Nematoda ~.,ccamodel_adonis_data_Nematoda_meta)

ccamodel_adonis_data_Nematoda_meta$season<-as.factor(ccamodel_adonis_data_Nematoda_meta$season)
ccamodel_adonis_data_Nematoda_meta$SampleSite_code<-as.factor(ccamodel_adonis_data_Nematoda_meta$SampleSite_code)
ccamodel_adonis_data_Nematoda_meta$river_section<-as.factor(ccamodel_adonis_data_Nematoda_meta$river_section)

plot(ccamodel_Nematoda,scaling="sites")
points(ccamodel_Nematoda, pch=21, cex=1.2, scaling="sites")
with(ccamodel_adonis_data_Nematoda_meta, points(ccamodel_Nematoda, scaling="sites", pch=21, col=1, bg=river_section))
with(ccamodel_adonis_data_Nematoda_meta, legend("bottomleft", levels(river_section), pch=21, pt.bg=1:12, bty="n"))

###### PERMANOVA ON Annelida #######

NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data_Annelida<-NORM_WP2A_ASV_table_wide_filtered

#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data_Annelida<-adonis_data_Annelida %>% rename(SampleSite_time_Annelida = X1)
adonis_data_Annelida <- adonis_data_Annelida[,grepl("Annelida", colnames(adonis_data_Annelida))]
adonis_data_Annelida<-adonis_data_Annelida[, colSums(adonis_data_Annelida != 0) > 0]
adonis_data_Annelida<-adonis_data_Annelida[apply(adonis_data_Annelida[,-1], 1, function(x) !all(x==0)),]

adonis_data_Annelida_samples<-adonis_data_Annelida$SampleSite_time_Annelida
adonis_data_Annelida_samples<-data.frame(adonis_data_Annelida_samples)
adonis_data_Annelida_samples<-adonis_data_Annelida_samples %>% rename(ID = adonis_data_Annelida_samples)
adonis_data_Annelida<-adonis_data_Annelida %>% rename(SampleSite_time = SampleSite_time_Annelida)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data_Annelida<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","season","Days","Dist_from_lake","pH",
                                                        "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                        "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                        "temp_monthly_average","river_section")],adonis_data_Annelida, by="SampleSite_time")
adonis_data_Annelida<-na.omit(adonis_data_Annelida)

adonis_data_Annelida_meta<- adonis_data_Annelida[c(1:16)]
adonis_data_Annelida<-adonis_data_Annelida[,-c(1:16)]

adonis_data_Annelida_meta$Days<-as.numeric(adonis_data_Annelida_meta$Days)

adon.results_WP2_Annelida<-adonis2(adonis_data_Annelida ~ adonis_data_Annelida_meta$Dist_from_lake+ adonis_data_Annelida_meta$Days+
                                     adonis_data_Annelida_meta$season/adonis_data_Annelida_meta$Days+
                                     adonis_data_Annelida_meta$Dist_from_lake+
                                     adonis_data_Annelida_meta$daily_flow_by_catchment+
                                     adonis_data_Annelida_meta$monthly_flow_by_catchment+
                                     adonis_data_Annelida_meta$seasonal_flow_by_catchment+
                                     adonis_data_Annelida_meta$week_before_flow_by_catchment+
                                     adonis_data_Annelida_meta$river_width_m+
                                     adonis_data_Annelida_meta$river_gradient_percent+
                                     adonis_data_Annelida_meta$rainfall_monthly_average+
                                     adonis_data_Annelida_meta$temp_monthly_average+
                                     adonis_data_Annelida_meta$gran_alk_ueq_L+
                                     adonis_data_Annelida_meta$conductivity_uS_cm+
                                     adonis_data_Annelida_meta$pH,
                                   method="bray",perm=999,by="margin")
print(adon.results_WP2_Annelida)

#PLOTTING BETA DIVERSITY NMDS
com_Annelida <- adonis_data_Annelida[,col(adonis_data_Annelida)]
m_com_Annelida <- as.matrix(com_Annelida)
set.seed(123)
nmds_Annelida = metaMDS(m_com_Annelida, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores_Annelida = as.data.frame(scores(nmds_Annelida))
#add columns to data frame 
data.scores_Annelida$SampleSite_time = adonis_data_Annelida_meta[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_Annelida_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Blank","SampleSite_code","Days","Dist_from_lake")],data.scores_Annelida, by="SampleSite_time")

p<-
  ggplot(data.scores_Annelida_meta_data, aes(x = NMDS1, y = NMDS2,colour = SampleSite_time)) + 
  geom_point(size = 2, aes(colour = SampleSite_time))+
  theme_bw()
ggplotly(p)

p<-
  ggplot(data.scores_Annelida_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)

#Canonical Correspondence Analysis (CCA) NEMATODA
ccamodel_adonis_data_Annelida_meta<-adonis_data_Annelida_meta
ccamodel_adonis_data_Annelida_meta$SampleSite_time<-NULL

ccamodel_Annelida<-cca(adonis_data_Annelida ~.,ccamodel_adonis_data_Annelida_meta)

ccamodel_adonis_data_Annelida_meta$season<-as.factor(ccamodel_adonis_data_Annelida_meta$season)
ccamodel_adonis_data_Annelida_meta$SampleSite_code<-as.factor(ccamodel_adonis_data_Annelida_meta$SampleSite_code)
ccamodel_adonis_data_Annelida_meta$river_section<-as.factor(ccamodel_adonis_data_Annelida_meta$river_section)

plot(ccamodel_Annelida,scaling="sites")
points(ccamodel_Annelida, pch=21, cex=1.2, scaling="sites")
with(ccamodel_adonis_data_Annelida_meta, points(ccamodel_Annelida, scaling="sites", pch=21, col=1, bg=river_section))
with(ccamodel_adonis_data_Annelida_meta, legend("bottomleft", levels(river_section), pch=21, pt.bg=1:12, bty="n"))


###### PERMANOVA ON Mollusca ######
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data_Mollusca<-NORM_WP2A_ASV_table_wide_filtered

#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data_Mollusca<-adonis_data_Mollusca[adonis_data_Mollusca$X1 != "E12_T03", ]  
adonis_data_Mollusca<-adonis_data_Mollusca[adonis_data_Mollusca$X1 != "E05_T10", ]  
adonis_data_Mollusca<-adonis_data_Mollusca[adonis_data_Mollusca$X1 != "E04_T05", ]  
adonis_data_Mollusca<-adonis_data_Mollusca[adonis_data_Mollusca$X1 != "E06_T10", ]  

adonis_data_Mollusca<-adonis_data_Mollusca %>% rename(SampleSite_time_Mollusca = X1)
adonis_data_Mollusca <- adonis_data_Mollusca[,grepl("Mollusca", colnames(adonis_data_Mollusca))]
adonis_data_Mollusca<-adonis_data_Mollusca[, colSums(adonis_data_Mollusca != 0) > 0]
adonis_data_Mollusca<-adonis_data_Mollusca[apply(adonis_data_Mollusca[,-1], 1, function(x) !all(x==0)),]

adonis_data_Mollusca_samples<-adonis_data_Mollusca$SampleSite_time_Mollusca
adonis_data_Mollusca_samples<-data.frame(adonis_data_Mollusca_samples)
adonis_data_Mollusca_samples<-adonis_data_Mollusca_samples %>% rename(ID = adonis_data_Mollusca_samples)
adonis_data_Mollusca<-adonis_data_Mollusca %>% rename(SampleSite_time = SampleSite_time_Mollusca)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data_Mollusca<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","river_section","season","Days","Dist_from_lake","pH",
                                                          "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                          "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                          "temp_monthly_average")],adonis_data_Mollusca, by="SampleSite_time")
adonis_data_Mollusca<-na.omit(adonis_data_Mollusca)

adonis_data_Mollusca_meta<- adonis_data_Mollusca[c(1:16)]
adonis_data_Mollusca<-adonis_data_Mollusca[,-c(1:16)]

adonis_data_Mollusca_meta$Days<-as.numeric(adonis_data_Mollusca_meta$Days)

adon.results_WP2_Mollusca<-adonis2(adonis_data_Mollusca ~ adonis_data_Mollusca_meta$Dist_from_lake+ adonis_data_Mollusca_meta$Days+
                                       adonis_data_Mollusca_meta$season/adonis_data_Mollusca_meta$Days+
                                       adonis_data_Mollusca_meta$Dist_from_lake+
                                       adonis_data_Mollusca_meta$daily_flow_by_catchment+
                                       adonis_data_Mollusca_meta$monthly_flow_by_catchment+
                                       adonis_data_Mollusca_meta$seasonal_flow_by_catchment+
                                       adonis_data_Mollusca_meta$week_before_flow_by_catchment+
                                       adonis_data_Mollusca_meta$river_width_m+
                                       adonis_data_Mollusca_meta$river_gradient_percent+
                                       adonis_data_Mollusca_meta$rainfall_monthly_average+
                                       adonis_data_Mollusca_meta$temp_monthly_average+
                                       adonis_data_Mollusca_meta$gran_alk_ueq_L+
                                       adonis_data_Mollusca_meta$conductivity_uS_cm+
                                       adonis_data_Mollusca_meta$pH,
                                     method="bray",perm=999,by="margin")
print(adon.results_WP2_Mollusca)

##########   ALPHA DIVERSITY METRICS ################
shannon_index<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2A_18S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_18S_WP2A<-merge(shannon_index,simpson_index,by="ID")
alpha_18S_WP2A<-merge(alpha_18S_WP2A,species_count,by="ID")

write.csv(alpha_18S_WP2A,"alpha_18S_WP2A.csv")

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

alpha_18S_WP2A_meta_data <- alpha_18S_WP2A %>% rename(SampleSite_time= ID)
alpha_18S_WP2A_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days","season")],alpha_18S_WP2A_meta_data, by="SampleSite_time")

  ggplot(alpha_18S_WP2A_meta_data, aes(x=Dist_from_lake, y=shannon_index, color=season)) +
  geom_point() + 
    geom_smooth(method = "gam", formula = y ~ s(x))+
    theme_bw()


#ALPHA DIVERSITY Arthropoda
seqtabNoC_WP2A_18S_tax_normalised_wide_Arthropoda <- seqtabNoC_WP2A_18S_tax_normalised_wide[ , grepl( "Arthropoda" , names( seqtabNoC_WP2A_18S_tax_normalised_wide ) ) ]
shannon_index_Arthropoda<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide_Arthropoda, index = "shannon")
shannon_index_Arthropoda<-data.frame(shannon_index_Arthropoda)
ID_Arthropoda <- rownames(shannon_index_Arthropoda)
shannon_index_Arthropoda <- cbind(shannon_index_Arthropoda,ID_Arthropoda)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Arthropoda_meta_data <- shannon_index_Arthropoda %>% rename(SampleSite_time= ID_Arthropoda)
shannon_index_Arthropoda_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Arthropoda_meta_data, by="SampleSite_time")
shannon_index_Arthropoda_meta_data$phyla <- paste0("Arthropoda", shannon_index_Arthropoda_meta_data$phyla)
shannon_index_Arthropoda_meta_data<-shannon_index_Arthropoda_meta_data %>% rename(shannon_index = shannon_index_Arthropoda)

p_Arthropoda<-ggplot(shannon_index_Arthropoda_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

shannon_index_Arthropoda_meta_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                                                  "season",
                                                                  "temp_monthly_average",
                                                                  "conductivity_uS_cm",
                                                                  "gran_alk_ueq_L","daily_flow_by_catchment",
                                                                  "monthly_flow_by_catchment",
                                                                  "seasonal_flow_by_catchment",
                                                                  "week_before_flow_by_catchment",
                                                                  "river_width_m",
                                                                  "river_gradient_percent",
                                                                  "rainfall_monthly_average",
                                                                  "pH")],shannon_index_Arthropoda_meta_data, by="SampleSite_time")

shannon_index_Arthropoda_meta_data$Days<-as.numeric(shannon_index_Arthropoda_meta_data$Days)
shannon_index_Arthropoda_meta_data<-shannon_index_Arthropoda_meta_data[shannon_index_Arthropoda_meta_data$shannon_index != 0, ] 

shannon_index_Arthropoda_meta_data$season<-as.factor(shannon_index_Arthropoda_meta_data$season)

gam_Arthropoda_plot<-ggplot(shannon_index_Arthropoda_meta_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

ggplot(shannon_index_Arthropoda_meta_data, aes(daily_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

ggplot(shannon_index_Arthropoda_meta_data, aes(week_before_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

ggplot(shannon_index_Arthropoda_meta_data, aes(rainfall_monthly_average, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")



mod_gam <- gam(shannon_index ~ 
                 season+
                 s(Dist_from_lake)+
                 s(pH)+
                 s(daily_flow_by_catchment)+
                 s(monthly_flow_by_catchment)+
                 s(seasonal_flow_by_catchment)+
                 s(week_before_flow_by_catchment)+
                 s(river_gradient_percent)+
                 s(rainfall_monthly_average)+
                 s(temp_monthly_average)+
                 s(gran_alk_ueq_L)+
                 s(conductivity_uS_cm), 
               data=shannon_index_Arthropoda_meta_data)
summary(mod_gam)


#ALPHA DIVERSITY Rotifera
seqtabNoC_WP2A_18S_tax_normalised_wide_Rotifera <- seqtabNoC_WP2A_18S_tax_normalised_wide[ , grepl( "Rotifera" , names( seqtabNoC_WP2A_18S_tax_normalised_wide ) ) ]
shannon_index_Rotifera<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide_Rotifera, index = "shannon")
shannon_index_Rotifera<-data.frame(shannon_index_Rotifera)
ID_Rotifera <- rownames(shannon_index_Rotifera)
shannon_index_Rotifera <- cbind(shannon_index_Rotifera,ID_Rotifera)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Rotifera_meta_data <- shannon_index_Rotifera %>% rename(SampleSite_time= ID_Rotifera)
shannon_index_Rotifera_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Rotifera_meta_data, by="SampleSite_time")
shannon_index_Rotifera_meta_data$phyla <- paste0("Rotifera", shannon_index_Rotifera_meta_data$phyla)
shannon_index_Rotifera_meta_data<-shannon_index_Rotifera_meta_data %>% rename(shannon_index = shannon_index_Rotifera)

p_Rotifera<-ggplot(shannon_index_Rotifera_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

shannon_index_Rotifera_meta_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                                                      "season",
                                                                      "temp_monthly_average",
                                                                      "conductivity_uS_cm",
                                                                      "gran_alk_ueq_L","daily_flow_by_catchment",
                                                                      "monthly_flow_by_catchment",
                                                                      "seasonal_flow_by_catchment",
                                                                      "week_before_flow_by_catchment",
                                                                      "river_width_m",
                                                                      "river_gradient_percent",
                                                                      "rainfall_monthly_average",
                                                                      "pH")],shannon_index_Rotifera_meta_data, by="SampleSite_time")

shannon_index_Rotifera_meta_data$Days<-as.numeric(shannon_index_Rotifera_meta_data$Days)
shannon_index_Rotifera_meta_data<-shannon_index_Rotifera_meta_data[shannon_index_Rotifera_meta_data$shannon_index != 0, ] 

shannon_index_Rotifera_meta_data$season<-as.factor(shannon_index_Rotifera_meta_data$season)

gam_Rotifera_plot<-ggplot(shannon_index_Rotifera_meta_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

mod_gam <- gam(shannon_index ~ 
                 season+
                 s(Dist_from_lake)+
                 s(pH)+
                 s(daily_flow_by_catchment)+
                 s(monthly_flow_by_catchment)+
                 s(seasonal_flow_by_catchment)+
                 s(week_before_flow_by_catchment)+
                 s(river_gradient_percent)+
                 s(rainfall_monthly_average)+
                 s(temp_monthly_average)+
                 s(gran_alk_ueq_L)+
                 s(conductivity_uS_cm), 
               data=shannon_index_Rotifera_meta_data)
summary(mod_gam)

#ALPHA DIVERSITY Nematoda
seqtabNoC_WP2A_18S_tax_normalised_wide_Nematoda <- seqtabNoC_WP2A_18S_tax_normalised_wide[ , grepl( "Nematoda" , names( seqtabNoC_WP2A_18S_tax_normalised_wide ) ) ]
shannon_index_Nematoda<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide_Nematoda, index = "shannon")
shannon_index_Nematoda<-data.frame(shannon_index_Nematoda)
ID_Nematoda <- rownames(shannon_index_Nematoda)
shannon_index_Nematoda <- cbind(shannon_index_Nematoda,ID_Nematoda)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Nematoda_meta_data <- shannon_index_Nematoda %>% rename(SampleSite_time= ID_Nematoda)
shannon_index_Nematoda_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Nematoda_meta_data, by="SampleSite_time")
shannon_index_Nematoda_meta_data$Days<-as.numeric(shannon_index_Nematoda_meta_data$Days)
shannon_index_Nematoda_meta_data$phyla <- paste0("Nematoda", shannon_index_Nematoda_meta_data$phyla)
shannon_index_Nematoda_meta_data<-shannon_index_Nematoda_meta_data %>% rename(shannon_index = shannon_index_Nematoda)

p_Nematoda<-ggplot(shannon_index_Nematoda_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

shannon_index_Nematoda_meta_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                                                    "season",
                                                                    "temp_monthly_average",
                                                                    "conductivity_uS_cm",
                                                                    "gran_alk_ueq_L","daily_flow_by_catchment",
                                                                    "monthly_flow_by_catchment",
                                                                    "seasonal_flow_by_catchment",
                                                                    "week_before_flow_by_catchment",
                                                                    "river_width_m",
                                                                    "river_gradient_percent",
                                                                    "rainfall_monthly_average",
                                                                    "pH")],shannon_index_Nematoda_meta_data, by="SampleSite_time")

shannon_index_Nematoda_meta_data$Days<-as.numeric(shannon_index_Nematoda_meta_data$Days)
shannon_index_Nematoda_meta_data<-shannon_index_Nematoda_meta_data[shannon_index_Nematoda_meta_data$shannon_index != 0, ] 

shannon_index_Nematoda_meta_data$season<-as.factor(shannon_index_Nematoda_meta_data$season)

gam_Nematoda_plot<-ggplot(shannon_index_Nematoda_meta_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

mod_gam <- gam(shannon_index ~ 
                 season+
                 s(Dist_from_lake)+
                 s(pH)+
                 s(daily_flow_by_catchment)+
                 s(monthly_flow_by_catchment)+
                 s(seasonal_flow_by_catchment)+
                 s(week_before_flow_by_catchment)+
                 s(river_gradient_percent)+
                 s(rainfall_monthly_average)+
                 s(temp_monthly_average)+
                 s(gran_alk_ueq_L)+
                 s(conductivity_uS_cm), 
               data=shannon_index_Nematoda_meta_data)
summary(mod_gam)


#ALPHA DIVERSITY Mollusca
seqtabNoC_WP2A_18S_tax_normalised_wide_Mollusca <- seqtabNoC_WP2A_18S_tax_normalised_wide[ , grepl( "Mollusca" , names( seqtabNoC_WP2A_18S_tax_normalised_wide ) ) ]
shannon_index_Mollusca<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide_Mollusca, index = "shannon")
shannon_index_Mollusca<-data.frame(shannon_index_Mollusca)
ID_Mollusca <- rownames(shannon_index_Mollusca)
shannon_index_Mollusca <- cbind(shannon_index_Mollusca,ID_Mollusca)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Mollusca_meta_data <- shannon_index_Mollusca %>% rename(SampleSite_time= ID_Mollusca)
shannon_index_Mollusca_meta_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Dist_from_lake",
                                                                    "season",
                                                                    "temp_monthly_average",
                                                                    "conductivity_uS_cm",
                                                                    "gran_alk_ueq_L","daily_flow_by_catchment",
                                                                    "monthly_flow_by_catchment",
                                                                    "seasonal_flow_by_catchment",
                                                                    "week_before_flow_by_catchment",
                                                                    "river_width_m",
                                                                    "river_gradient_percent",
                                                                    "rainfall_monthly_average",
                                                                    "pH")],shannon_index_Mollusca_meta_data, by="SampleSite_time")
shannon_index_Mollusca_meta_data$phyla <- paste0("Mollusca", shannon_index_Mollusca_meta_data$phyla)
shannon_index_Mollusca_meta_data<-shannon_index_Mollusca_meta_data %>% rename(shannon_index = shannon_index_Mollusca)



p_Mollusca<-ggplot(shannon_index_Mollusca_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

mod_gam <- gam(shannon_index ~ season+
                 s(Dist_from_lake)+
                 s(pH)+
                 s(daily_flow_by_catchment)+
                 s(monthly_flow_by_catchment)+
                 s(seasonal_flow_by_catchment)+
                 s(week_before_flow_by_catchment)+
                 s(river_gradient_percent)+
                 s(rainfall_monthly_average)+
                 s(temp_monthly_average)+
                 s(gran_alk_ueq_L)+
                 s(conductivity_uS_cm), 
               data=shannon_index_Mollusca_meta_data)
summary(mod_gam)

gam_Mollusca_plot<-ggplot(shannon_index_Mollusca_meta_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")


#ALPHA DIVERSITY Porifera
seqtabNoC_WP2A_18S_tax_normalised_wide_Porifera <- seqtabNoC_WP2A_18S_tax_normalised_wide[ , grepl( "Porifera" , names( seqtabNoC_WP2A_18S_tax_normalised_wide ) ) ]
shannon_index_Porifera<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide_Porifera, index = "shannon")
shannon_index_Porifera<-data.frame(shannon_index_Porifera)
ID_Porifera <- rownames(shannon_index_Porifera)
shannon_index_Porifera <- cbind(shannon_index_Porifera,ID_Porifera)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Porifera_meta_data <- shannon_index_Porifera %>% rename(SampleSite_time= ID_Porifera)
shannon_index_Porifera_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Porifera_meta_data, by="SampleSite_time")
shannon_index_Porifera_meta_data$phyla <- paste0("Porifera", shannon_index_Porifera_meta_data$phyla)
shannon_index_Porifera_meta_data<-shannon_index_Porifera_meta_data %>% rename(shannon_index = shannon_index_Porifera)

p_Porifera<-ggplot(shannon_index_Porifera_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

#ALPHA DIVERSITY Chordata
seqtabNoC_WP2A_18S_tax_normalised_wide_Chordata <- seqtabNoC_WP2A_18S_tax_normalised_wide[ , grepl( "Chordata" , names( seqtabNoC_WP2A_18S_tax_normalised_wide ) ) ]
shannon_index_Chordata<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide_Chordata, index = "shannon")
shannon_index_Chordata<-data.frame(shannon_index_Chordata)
ID_Chordata <- rownames(shannon_index_Chordata)
shannon_index_Chordata <- cbind(shannon_index_Chordata,ID_Chordata)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Chordata_meta_data <- shannon_index_Chordata %>% rename(SampleSite_time= ID_Chordata)
shannon_index_Chordata_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Chordata_meta_data, by="SampleSite_time")
shannon_index_Chordata_meta_data$phyla <- paste0("Chordata", shannon_index_Chordata_meta_data$phyla)
shannon_index_Chordata_meta_data<-shannon_index_Chordata_meta_data %>% rename(shannon_index = shannon_index_Chordata)

p_Chordata<-ggplot(shannon_index_Chordata_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

#ALPHA DIVERSITY Annelida
seqtabNoC_WP2A_18S_tax_normalised_wide_Annelida <- seqtabNoC_WP2A_18S_tax_normalised_wide[ , grepl( "Annelida" , names( seqtabNoC_WP2A_18S_tax_normalised_wide ) ) ]
shannon_index_Annelida<-diversity(seqtabNoC_WP2A_18S_tax_normalised_wide_Annelida, index = "shannon")
shannon_index_Annelida<-data.frame(shannon_index_Annelida)
ID_Annelida <- rownames(shannon_index_Annelida)
shannon_index_Annelida <- cbind(shannon_index_Annelida,ID_Annelida)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

shannon_index_Annelida_meta_data <- shannon_index_Annelida %>% rename(SampleSite_time= ID_Annelida)
shannon_index_Annelida_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Date_sampled","Dist_from_lake","Days")],shannon_index_Annelida_meta_data, by="SampleSite_time")
shannon_index_Annelida_meta_data$phyla <- paste0("Annelida", shannon_index_Annelida_meta_data$phyla)
shannon_index_Annelida_meta_data<-shannon_index_Annelida_meta_data %>% rename(shannon_index = shannon_index_Annelida)

p_Annelida<-ggplot(shannon_index_Annelida_meta_data, aes(x=SampleSite_code, y=shannon_index,colour=Dist_from_lake)) +
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter()+ 
  theme_bw()

shannon_index_Annelida_meta_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                                                    "season",
                                                                    "temp_monthly_average",
                                                                    "conductivity_uS_cm",
                                                                    "gran_alk_ueq_L","daily_flow_by_catchment",
                                                                    "monthly_flow_by_catchment",
                                                                    "seasonal_flow_by_catchment",
                                                                    "week_before_flow_by_catchment",
                                                                    "river_width_m",
                                                                    "river_gradient_percent",
                                                                    "rainfall_monthly_average",
                                                                    "pH")],shannon_index_Annelida_meta_data, by="SampleSite_time")

shannon_index_Annelida_meta_data$Days<-as.numeric(shannon_index_Annelida_meta_data$Days)
shannon_index_Annelida_meta_data<-shannon_index_Annelida_meta_data[shannon_index_Annelida_meta_data$shannon_index != 0, ] 

shannon_index_Annelida_meta_data$season<-as.factor(shannon_index_Annelida_meta_data$season)

gam_Annelida_plot<-ggplot(shannon_index_Annelida_meta_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+theme_bw()+
  theme(legend.position="none")

mod_gam <- gam(shannon_index ~ season+
                 s(Dist_from_lake)+
                 s(pH)+
                 s(daily_flow_by_catchment)+
                 s(monthly_flow_by_catchment)+
                 s(seasonal_flow_by_catchment)+
                 s(week_before_flow_by_catchment)+
                 s(river_gradient_percent)+
                 s(rainfall_monthly_average)+
                 s(temp_monthly_average)+
                 s(gran_alk_ueq_L)+
                 s(conductivity_uS_cm), 
               data=shannon_index_Annelida_meta_data)
summary(mod_gam)

library(ggpubr)
ggarrange(p_Annelida,p_Arthropoda,p_Mollusca,p_Nematoda,p_Rotifera,ncol = 2,nrow=3)
ggarrange(gam_Arthropoda_plot,gam_Rotifera_plot,gam_Nematoda_plot,gam_Annelida_plot,gam_Mollusca_plot,ncol = 2,nrow=3)


#combine alpha diversity from different phyla
shannon_index_phlya_split_data<-rbind(shannon_index_Nematoda_meta_data,
                                     # shannon_index_Mollusca_meta_data,
                                      shannon_index_Annelida_meta_data,
                                      shannon_index_Rotifera_meta_data,
                                      shannon_index_Arthropoda_meta_data)

######## WP2A alpha diversity modeling ##########
library(mgcv)

shannon_index_phlya_split_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time",
                                                         "season",
                                                         "temp_monthly_average",
                                                         "conductivity_uS_cm",
                                                         "gran_alk_ueq_L","daily_flow_by_catchment",
                                                         "monthly_flow_by_catchment",
                                                         "seasonal_flow_by_catchment",
                                                         "week_before_flow_by_catchment",
                                                         "river_width_m",
                                                         "river_gradient_percent",
                                                         "rainfall_monthly_average",
                                                         "pH")],shannon_index_phlya_split_data, by="SampleSite_time")

shannon_index_phlya_split_data$Days<-as.numeric(shannon_index_phlya_split_data$Days)
shannon_index_phlya_split_data<-shannon_index_phlya_split_data[shannon_index_phlya_split_data$shannon_index != 0, ] 

#GAMs
ggplot(shannon_index_phlya_split_data, aes(Days, shannon_index,colour=phyla)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(temp_monthly_average, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(daily_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(conductivity_uS_cm, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(pH, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(gran_alk_ueq_L, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(daily_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(rainfall_monthly_average, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(monthly_flow_by_catchment, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(river_width_m, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(shannon_index_phlya_split_data, aes(river_gradient_percent, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))

ggplot(shannon_index_phlya_split_data, aes(season, shannon_index)) + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~SampleSite_code)

ggplot(shannon_index_phlya_split_data, aes(Dist_from_lake, shannon_index,colour=season)) + 
  geom_point() + 
  #geom_smooth()+
  geom_smooth(method = "gam", formula = y ~ s(x))+
  facet_wrap(~phyla)+theme_bw()

shannon_index_phlya_split_data$phyla<-as.factor(shannon_index_phlya_split_data$phyla)
shannon_index_phlya_split_data$season<-as.factor(shannon_index_phlya_split_data$season)

mod_gam <- gam(shannon_index ~ phyla+
                 s(Days)+
                 season+
                 Dist_from_lake+
                 s(Dist_from_lake, by=season)+
                 s(Dist_from_lake, by=phyla)+
                 s(Days, by=phyla)+
                 s(Dist_from_lake, by=pH)+
                 daily_flow_by_catchment+
                 monthly_flow_by_catchment+
                 seasonal_flow_by_catchment+
                 week_before_flow_by_catchment+
                 river_width_m+
                 river_gradient_percent+
                 rainfall_monthly_average+
                 temp_monthly_average+
                 s(gran_alk_ueq_L)+
                 s(conductivity_uS_cm), 
  data=shannon_index_phlya_split_data)
summary(mod_gam)

plot(mod_gam, shade = TRUE, pages = 1, scale = 0)

plot_smooths(
  model = mod_gam2,
  series = Dist_from_lake,
  comparison = Days,
  facet_terms = phyla,
  split = list(Daysphyla = c("Days", "phyla"))
)

######## Turnover and nestedness ########
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
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
library(ggdendro)
library(plotly)

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
NORM_WP2A_ASV_table_wide_filtered <- read_csv("NORM_WP2A_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2A_ASV_table_wide_filtered
#REMOVE OUTLIERS IDENTIFIED BY THE NMDS PLOT
adonis_data<-adonis_data[adonis_data$X1 != "E04_T07", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T19", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E02_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E14_T01", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T09", ]  
adonis_data<-adonis_data[adonis_data$X1 != "E01_T18", ]  

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

adonis_data_meta<- adonis_data[c(1:29)]
adonis_data<-adonis_data[,-c(1:29)]

adonis_data_meta$Days<-as.numeric(adonis_data_meta$Days)
adonis_data_meta$SampleSite_time<-NULL

ncomp = 1
result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'canonical')

design <- data.frame(samp = adonis_data_meta$sample)

#result.spls <- spls(adonis_data, adonis_data_meta, ncomp = ncomp,  mode = 'regression')
#tune.spls <- perf(result.spls, validation = 'Mfold', folds = 5,
 #                 criterion = 'all',nrepeat = 99, progressBar = TRUE)

#tune.spls$R2
##plot(tune.spls$Q2.total)
#abline(h = 0.0975)

#plot(tune.spls, criterion="R2", type = 'l')
#plot(tune.spls, criterion="Q2", type = 'l')
#plot(tune.spls)

cim_plot<-cim(result.spls,
              comp = 1:1,
             margins = c(15, 5))

cl<-colnames(adonis_data)
cl<-data.frame(cl)
cl<- data.frame(do.call('rbind', strsplit(as.character(cl$cl),';',fixed=TRUE)))
cl<-subset(cl,select=c(X2))

unique(cl[c("X2")])


#gene.col<-c(                               'Annelida'='coral3',
#                                         'Arthropoda'='chartreuse4',
#                                           'Chordata'='cadeblue3',
#                                       'Gastrotricha'='burlywood3',
#                                         'Mollusca'='brown1',
#                                          'Nematoda'='blueviolet',
#                                   'Platyhelminthes'='blue4',
#                                          'Porifera'='blue',
#                                          'Rotifera'='black',
#                                        'Tardigrada'='azure4',
#                                         'Bryozoa'='bisque1',
#                                    'Chaetognatha'='aquamarine2',
#                                        'Cnidaria'='deeppink',
#              'Ctenophora <comb jellyfish phylum>'='darkorange',
#                                       'Hemichordata'='darkgreen',
#              'Acanthocephala <thorny-headed worm>'='cyan',
#                                     'Nematomorpha'='gold',
#                                        'Nemertea'='yellow',
#                                'Xenacoelomorpha'='red')


stim.col<-c("darkblue","purple", "green4","red3","orange", "cyan","darkviolet","yellow",
            "darkorange","red","green","black","darkcyan","darkred","darkorchid","darksalmon","gray","violet","springgreen4")#19
cl$X2<-as.factor(cl$X2)
stim.col <- stim.col[as.numeric(cl$X2)]

cim_plot<-cim(result.spls,row.sideColors = stim.col,
              comp = 1:1,
    margins = c(15, 5))

#extracting values from cim plot

#ggdendro_plot<-ggdendrogram(dendrogram)+ coord_flip() +  scale_y_reverse(expand = c(0.2, 0))+
 # theme_dendro() 


spls_correlation<-cim_plot[["mat"]]
spls_correlation<-data.frame(spls_correlation)
spls_correlation <- tibble::rownames_to_column(spls_correlation, "taxa")

#spls_correlation_ordered<-spls_correlation[(cim_plot$rowInd), ]
#spls_correlation<-tibble::rowid_to_column(spls_correlation, "row_number")
#row_order_heatmap<-cim_plot[["rowInd"]]
#row_order_heatmap<-data.frame(row_order_heatmap)

#Calculating 95% confidence interval on correlations
conf_95_spls_correlation<-spls_correlation
conf_95_spls_correlation$row_number<-NULL
conf_95_spls_correlation$taxa<-NULL
conf_95_spls_correlation<-data.frame(x=unlist(conf_95_spls_correlation[]))
hist(conf_95_spls_correlation$x)
boxplot(conf_95_spls_correlation$x)

quantile(conf_95_spls_correlation$x, probs = seq(0, 1, 0.01))
mean(conf_95_spls_correlation$x)

#investigating the 95% percentile of correlations 
spls_correlation_upper <- gather(spls_correlation, env_variable, correlation, Inland_rock_ha:rainfall_monthly_average, factor_key=TRUE)
spls_correlation_upper<-spls_correlation_upper %>%
  filter(correlation > 0.2695705083)
spls_correlation_upper$env_variable<-NULL
spls_correlation_upper$correlation<-NULL
spls_correlation_upper<-unique(spls_correlation_upper$taxa)
spls_correlation_upper<-data.frame(spls_correlation_upper)
spls_correlation_upper$keep<-"keep"
spls_correlation_upper<-spls_correlation_upper %>% rename(species = spls_correlation_upper)


setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

spls_correlation_upper_abundance_plot<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered
spls_correlation_upper_abundance_plot<-merge(seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered,spls_correlation_upper,by="species")

spls_correlation_upper_abundance_plot$species<-gsub("NA;", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(";NA", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(" <Metazoa>", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(" <vertebrate>", "", spls_correlation_upper_abundance_plot$species)
spls_correlation_upper_abundance_plot$species<-gsub(" <chordata>", "", spls_correlation_upper_abundance_plot$species)

ggplot(spls_correlation_upper_abundance_plot, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=1))

#EPT INVESTIGATION

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

EPT<-seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered[grep("Trichoptera|Ephemeroptera|Plecoptera",
                                seqtabNoC_WP2A_18S_tax_long_meta_clean_filtered$species), ]

EPT$species<-gsub(".*Trichoptera","Trichoptera",EPT$species)
#EPT$species<-gsub("Trichoptera.*","Trichoptera",EPT$species)
EPT$species<-gsub(".*Ephemeroptera","Ephemeroptera",EPT$species)
#EPT$species<-gsub("Ephemeroptera.*","Ephemeroptera",EPT$species)
EPT$species<-gsub(".*Plecoptera","Plecoptera",EPT$species)
#EPT$species<-gsub("Plecoptera.*","Plecoptera",EPT$species)

ggplot(EPT, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=1))


EPT_correlation<-spls_correlation[grep("Trichoptera|Ephemeroptera|Plecoptera",
                                       spls_correlation$taxa), ]

EPT_correlation<- gather(EPT_correlation, env_variable, correlation, Inland_rock_ha:rainfall_monthly_average, factor_key=TRUE)
ggplot(EPT_correlation, aes(x = env_variable, y = correlation)) + geom_boxplot()
 

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

#####WP2C_Glatt BLAST  TAXONOMY ANALYSIS ######

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered <- read_csv("WP2C_Glatt_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered<-seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered<-seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered))

names(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered <- seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered,"WP2C_Glatt_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
p<-ggplot(seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

ggplotly(p)

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
Sardina_pilchardus<-seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Glatt_18S_Venn<-seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean[!(seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Glatt_18S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Glatt_18S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Glatt_venn_diagramm.png',
#output=TRUE)

#NORMALISING READS
seqtabNoC_WP2C_Glatt_18S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered[c(1:257)])
seqtabNoC_WP2C_Glatt_18S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Glatt_18S_tax_col_sum2)
seqtabNoC_WP2C_Glatt_18S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Glatt_18S_tax_col_sum2)

seqtabNoC_WP2C_Glatt_18S_tax_normalised<- merge(seqtabNoC_WP2C_Glatt_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Glatt_18S_tax_col_sum2")],seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Glatt_18S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Glatt_18S_tax_wide_filtered[c(1:257)])
seqtabNoC_WP2C_Glatt_18S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Glatt_18S_tax_col_sum3)

seqtabNoC_WP2C_Glatt_18S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Glatt_18S_tax_col_sum3)

seqtabNoC_WP2C_Glatt_18S_tax_normalised<- merge(seqtabNoC_WP2C_Glatt_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Glatt_18S_tax_col_sum2")],seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Glatt_18S_tax_top<- merge(seqtabNoC_WP2C_Glatt_18S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Glatt_18S_tax_col_sum3")],seqtabNoC_WP2C_Glatt_18S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Glatt_18S_tax_normalised <- transform(seqtabNoC_WP2C_Glatt_18S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Glatt_18S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Glatt_18S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

#PLOTTING TOP Taxon
#seqtabNoC_WP2C_Glatt_18S_tax_top %>%
#  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Glatt_18S_tax_col_sum3))) %>%
#  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
#  geom_bar(position="fill",stat = "identity", colour = "black")+
#  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+scale_fill_brewer(palette="Paired")

#write.csv(seqtabNoC_WP2C_Glatt_18S_tax_normalised,"NORM_WP2C_Glatt_ASV_table_LONG_filtered.csv")

#PLOTTING PHYLA
seqtabNoC_WP2C_Glatt_18S_tax_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Glatt_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Glatt_18S_tax_phyla<-cbind(seqtabNoC_WP2C_Glatt_18S_tax_top,seqtabNoC_WP2C_Glatt_18S_tax_phyla)

seqtabNoC_WP2C_Glatt_18S_tax_phyla<-seqtabNoC_WP2C_Glatt_18S_tax_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Glatt_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Glatt_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Glatt_18S_tax_col_sum3=seqtabNoC_WP2C_Glatt_18S_tax_col_sum3)

seqtabNoC_WP2C_Glatt_18S_tax_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Glatt_18S_tax_phyla, by="SampleSite_time")
unique(seqtabNoC_WP2C_Glatt_18S_tax_phyla[c("X2")])
seqtabNoC_WP2C_Glatt_18S_tax_phyla$X2<-factor(seqtabNoC_WP2C_Glatt_18S_tax_phyla$X2,levels= c('Arthropoda','Rotifera','Mollusca','Chordata',
                                                                                  'Annelida','Nematoda','Gastrotricha', 'Tardigrada','Cnidaria','Porifera','Echinodermata','Xenacoelomorpha',
                                                                                  'Nemertea','Nematomorpha', 'Acanthocephala <thorny-headed worm>','Chaetognatha','Hemichordata',
                                                                                  'Ctenophora <comb jellyfish phylum>','Bryozoa','Platyhelminthes', 'Micrognathozoa'),ordered=TRUE)                                                                                         

ggplot(data=seqtabNoC_WP2C_Glatt_18S_tax_phyla,aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_manual(values = c("#A6CEE3",#Arthropods
                                                                             "#1F78B4", #Rotifera
                                                                             "#B2DF8A",#Molluscs
                                                                             "#33A02C",#Chordates
                                                                             "#FB9A99",#annelids
                                                                             "#E31A1C",#nematodes
                                                                             "#FDBF6F",#gastrotricha
                                                                             "#FF7F00",#tardigrada
                                                                             "#CAB2D6",#Cnideria
                                                                             "#6A3D9A",#porifera
                                                                             "#FFFF99",#echinodermata
                                                                             "white",
                                                                             "white","white","white","white"))

#PLOTTING TOP PHYLA
seqtabNoC_WP2C_Glatt_18S_tax_top_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Glatt_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Glatt_18S_tax_top_phyla<-cbind(seqtabNoC_WP2C_Glatt_18S_tax_top,seqtabNoC_WP2C_Glatt_18S_tax_top_phyla)

seqtabNoC_WP2C_Glatt_18S_tax_top_phyla<-seqtabNoC_WP2C_Glatt_18S_tax_top_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Glatt_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Glatt_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Glatt_18S_tax_col_sum3=seqtabNoC_WP2C_Glatt_18S_tax_col_sum3)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2C_Glatt_18S_tax_top_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Glatt_18S_tax_top_phyla, by="SampleSite_time")

seqtabNoC_WP2C_Glatt_18S_tax_top_phyla$SampleSite_time<-factor(seqtabNoC_WP2C_Glatt_18S_tax_top_phyla$SampleSite_time,levels= c('GL01',
                        'GL13','GL02','GL03','GL04','GL05','GL06','GL07','GL08','GL09','GL10','GL11','GL12'),ordered=TRUE)

seqtabNoC_WP2C_Glatt_18S_tax_top_phyla %>%
  mutate(X2 = fct_reorder(X2, desc(seqtabNoC_WP2C_Glatt_18S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide<-seqtabNoC_WP2C_Glatt_18S_tax_normalised
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide$seqtabNoC_WP2_18S_tax_col_sum2 <- NULL 

seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide<-seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide))

names(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide <- seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide,"NORM_WP2C_Glatt_ASV_table_wide_filtered.csv")

#Investigating lentic species
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide_transposed<-t(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide)
seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide_transposed<-data.frame(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide_transposed)
lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

#GL01
GL01_lentic<-seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide_transposed$GL01==0),]
GL01_lentic<-setDT(GL01_lentic, keep.rownames = TRUE)[]
GL01_lentic <- gather(GL01_lentic, SampleSite_time, normalised_reads, 2:14, factor_key=TRUE)
GL01_lentic<-data.frame(GL01_lentic)
GL01_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake")],GL01_lentic, by="SampleSite_time")
GL01_lentic<-na.omit(GL01_lentic)

pGL_01<-ggplot(data=GL01_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")
#ggplotly(pGL_01)

write.csv(GL01_lentic,"lentic_Glatt_18S.csv")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
GL01_lentic_wide<-GL01_lentic
GL01_lentic_wide$Dist_from_lake <- NULL 

GL01_lentic_wide<-GL01_lentic_wide %>% rename(species = rn)

GL01_lentic_wide<-GL01_lentic_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

GL01_lentic_wide <- spread(GL01_lentic_wide,SampleSite_time,normalised_reads)
GL01_lentic_wide[is.na(GL01_lentic_wide)] <- 0
GL01_lentic_wide <- data.frame(t(GL01_lentic_wide))

names(GL01_lentic_wide) <- as.matrix(GL01_lentic_wide[1, ])
GL01_lentic_wide <- GL01_lentic_wide[-1, ]
GL01_lentic_wide[] <- lapply(GL01_lentic_wide, function(x) type.convert(as.character(x)))

lentic_shannon_index<-diversity(GL01_lentic_wide, index = "shannon")
lentic_shannon_index<-data.frame(lentic_shannon_index)
ID <- rownames(lentic_shannon_index)
lentic_shannon_index_Glatt <- cbind(lentic_shannon_index,ID)
write.csv(lentic_shannon_index_Glatt,"lentic_shannon_index_Glatt.csv")

ggplot(data=lentic_shannon_index_Glatt, aes(x=ID, y=lentic_shannon_index, group=1)) +
  geom_line()+
  geom_point()

#PERMANOVA ON TAXON
NORM_WP2C_Glatt_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Glatt_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","river_width_m","Dist_from_lake")],adonis_data, by="SampleSite_time")

#adonis_data<-na.omit(adonis_data)#remove this if you want to keep in negative controls
adonis_data_meta<- adonis_data[c(1:4)]
adonis_data<-adonis_data[,-c(1:4)]

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$river_width_m*
                           adonis_data_meta$Dist_from_lake,
                         method="bray",perm=999)

print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_18S_ID <- read.csv("NORM_WP2C_Glatt_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_18S <- abund_table_18S_ID[,-1]
com <- abund_table_18S[,col(abund_table_18S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_18S_ID[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")

p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)

#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Glatt_18S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_18S_WP2C_Glatt<-merge(shannon_index,simpson_index,by="ID")
alpha_18S_WP2C_Glatt<-merge(alpha_18S_WP2C_Glatt,species_count,by="ID")

write.csv(alpha_18S_WP2C_Glatt,"alpha_18S_WP2C_Glatt.csv")

#________________________________________________________

#WP2C_Towy BLAST  TAXONOMY ANALYSIS

#________________________________________________________
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered <- read_csv("WP2C_Towy_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered<-seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Towy_18S_tax_wide_filtered<-seqtabNoC_WP2C_Towy_18S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Towy_18S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered))

names(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered <- seqtabNoC_WP2C_Towy_18S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Towy_18S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered,"WP2C_Towy_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
#Sardina_pilchardus<-seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

#Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
#  group_by(SampleSite_time) %>% 
#  summarize(reads=sum(reads)) %>%
#  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

#ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
#  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Towy_18S_Venn<-seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean[!(seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Towy_18S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Towy_18S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Towy_18S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Towy_18S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Towy_18S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Towy_18S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Towy_venn_diagramm.png',
#output=TRUE)


#NORMALISING READS
seqtabNoC_WP2C_Towy_18S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered[c(1:241)])
seqtabNoC_WP2C_Towy_18S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Towy_18S_tax_col_sum2)
seqtabNoC_WP2C_Towy_18S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Towy_18S_tax_col_sum2)

seqtabNoC_WP2C_Towy_18S_tax_normalised<- merge(seqtabNoC_WP2C_Towy_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Towy_18S_tax_col_sum2")],seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Towy_18S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Towy_18S_tax_wide_filtered[c(1:241)])
seqtabNoC_WP2C_Towy_18S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Towy_18S_tax_col_sum3)

seqtabNoC_WP2C_Towy_18S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Towy_18S_tax_col_sum3)

seqtabNoC_WP2C_Towy_18S_tax_normalised<- merge(seqtabNoC_WP2C_Towy_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Towy_18S_tax_col_sum2")],seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Towy_18S_tax_top<- merge(seqtabNoC_WP2C_Towy_18S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Towy_18S_tax_col_sum3")],seqtabNoC_WP2C_Towy_18S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Towy_18S_tax_normalised <- transform(seqtabNoC_WP2C_Towy_18S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Towy_18S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Towy_18S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

#PLOTTING TOP Taxon
seqtabNoC_WP2C_Towy_18S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Towy_18S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+scale_fill_brewer(palette="Paired")

write.csv(seqtabNoC_WP2C_Towy_18S_tax_normalised,"NORM_WP2C_Towy_ASV_table_LONG_filtered.csv")

#PLOTTING PHYLA
seqtabNoC_WP2C_Towy_18S_tax_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Towy_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Towy_18S_tax_phyla<-cbind(seqtabNoC_WP2C_Towy_18S_tax_top,seqtabNoC_WP2C_Towy_18S_tax_phyla)

seqtabNoC_WP2C_Towy_18S_tax_phyla<-seqtabNoC_WP2C_Towy_18S_tax_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Towy_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Towy_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Towy_18S_tax_col_sum3=seqtabNoC_WP2C_Towy_18S_tax_col_sum3)

seqtabNoC_WP2C_Towy_18S_tax_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Towy_18S_tax_phyla, by="SampleSite_time")
unique(seqtabNoC_WP2C_Towy_18S_tax_phyla[c("X2")])
seqtabNoC_WP2C_Towy_18S_tax_phyla$X2<-factor(seqtabNoC_WP2C_Towy_18S_tax_phyla$X2,levels= c('Arthropoda','Rotifera','Mollusca','Chordata',
                                                                                            'Annelida','Nematoda','Gastrotricha', 'Tardigrada','Cnidaria','Porifera','Echinodermata','Xenacoelomorpha',
                                                                                            'Nemertea','Nematomorpha', 'Acanthocephala <thorny-headed worm>','Chaetognatha','Hemichordata',
                                                                                            'Ctenophora <comb jellyfish phylum>','Bryozoa','Platyhelminthes', 'Micrognathozoa'),ordered=TRUE)                                                                                        

ggplot(data=seqtabNoC_WP2C_Towy_18S_tax_phyla,aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_manual(values = c("#A6CEE3",#Arthropods
                                                                             "#1F78B4", #Rotifera
                                                                             "#B2DF8A",#Molluscs
                                                                             "#33A02C",#Chordates
                                                                             "#FB9A99",#annelids
                                                                             "#E31A1C",#nematodes
                                                                             "#FDBF6F",#gastrotricha
                                                                             "#FF7F00",#tardigrada
                                                                             "#CAB2D6",#Cnideria
                                                                             "#6A3D9A",#porifera
                                                                             "white",#echinodermata
                                                                             "white",
                                                                             "white","white","white","white"))
#PLOTTING TOP PHYLA
seqtabNoC_WP2C_Towy_18S_tax_top_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Towy_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Towy_18S_tax_top_phyla<-cbind(seqtabNoC_WP2C_Towy_18S_tax_top,seqtabNoC_WP2C_Towy_18S_tax_top_phyla)

seqtabNoC_WP2C_Towy_18S_tax_top_phyla<-seqtabNoC_WP2C_Towy_18S_tax_top_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Towy_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Towy_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Towy_18S_tax_col_sum3=seqtabNoC_WP2C_Towy_18S_tax_col_sum3)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2C_Towy_18S_tax_top_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Towy_18S_tax_top_phyla, by="SampleSite_time")

seqtabNoC_WP2C_Towy_18S_tax_top_phyla %>%
  mutate(X2 = fct_reorder(X2, desc(seqtabNoC_WP2C_Towy_18S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide<-seqtabNoC_WP2C_Towy_18S_tax_normalised
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide$seqtabNoC_WP2_18S_tax_col_sum2 <- NULL 

seqtabNoC_WP2C_Towy_18S_tax_normalised_wide<-seqtabNoC_WP2C_Towy_18S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Towy_18S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide))

names(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide <- seqtabNoC_WP2C_Towy_18S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide,"NORM_WP2C_Towy_ASV_table_wide_filtered.csv")

#Investigating lentic species
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide_transposed<-t(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide)
seqtabNoC_WP2C_Towy_18S_tax_normalised_wide_transposed<-data.frame(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide_transposed)
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

#TW01
TW02_lentic<-seqtabNoC_WP2C_Towy_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide_transposed$TW02==0),]
TW02_lentic<-setDT(TW02_lentic, keep.rownames = TRUE)[]
TW02_lentic <- gather(TW02_lentic, SampleSite_time, normalised_reads, 2:13, factor_key=TRUE)
TW02_lentic<-data.frame(TW02_lentic)
TW02_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake")],TW02_lentic, by="SampleSite_time")
TW02_lentic<-na.omit(TW02_lentic)

pGL_01<-ggplot(data=TW02_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")
#ggplotly(pGL_01)

write.csv(TW02_lentic,"lentic_Towy_18S.csv")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
TW02_lentic_wide<-TW02_lentic
TW02_lentic_wide$Dist_from_lake <- NULL 

TW02_lentic_wide<-TW02_lentic_wide %>% rename(species = rn)

TW02_lentic_wide<-TW02_lentic_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

TW02_lentic_wide <- spread(TW02_lentic_wide,SampleSite_time,normalised_reads)
TW02_lentic_wide[is.na(TW02_lentic_wide)] <- 0
TW02_lentic_wide <- data.frame(t(TW02_lentic_wide))

names(TW02_lentic_wide) <- as.matrix(TW02_lentic_wide[1, ])
TW02_lentic_wide <- TW02_lentic_wide[-1, ]
TW02_lentic_wide[] <- lapply(TW02_lentic_wide, function(x) type.convert(as.character(x)))

lentic_shannon_index<-diversity(TW02_lentic_wide, index = "shannon")
lentic_shannon_index<-data.frame(lentic_shannon_index)
ID <- rownames(lentic_shannon_index)
lentic_shannon_index_Towy <- cbind(lentic_shannon_index,ID)
write.csv(lentic_shannon_index_Towy,"lentic_shannon_index_Towy.csv")

ggplot(data=lentic_shannon_index_Towy, aes(x=ID, y=lentic_shannon_index, group=1)) +
  geom_line()+
  geom_point()

#PERMANOVA ON TAXON
NORM_WP2C_Towy_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Towy_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Dist_from_lake","river_width_m")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)#remove this if you want to keep in negative controls
adonis_data_meta<- adonis_data[c(1:4)]
adonis_data<-adonis_data[,-c(1:4)]

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$river_width_m*
                           adonis_data_meta$Dist_from_lake,
                         method="bray",perm=999)
print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_18S_ID <- read.csv("NORM_WP2C_Towy_ASV_table_wide_filtered.csv", header = TRUE)
abund_table_18S_ID<-abund_table_18S_ID[!(abund_table_18S_ID$X=="TWNC"),]

abund_table_18S <- abund_table_18S_ID[,-1]
com <- abund_table_18S[,col(abund_table_18S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_18S_ID[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Blank","Days","Dist_from_lake")],data.scores, by="SampleSite_time")

p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Towy_18S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_18S_WP2C_Towy<-merge(shannon_index,simpson_index,by="ID")
alpha_18S_WP2C_Towy<-merge(alpha_18S_WP2C_Towy,species_count,by="ID")

write.csv(alpha_18S_WP2C_Towy,"alpha_18S_WP2C_Towy.csv")

#________________________________________________________

#WP2C_Gwash BLAST  TAXONOMY ANALYSIS

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered <- read_csv("WP2C_Gwash_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered$species=="Sardina_pilchardus"),]

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered))

names(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered <- seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#FILTERING BY NEGATIVE CONTROL
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-t(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-data.frame(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-setDT(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered, keep.rownames = TRUE)[]
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered %>% rename(OTU_ID = rn)

colnames_decon<-colnames(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)
colnames_decon<-data.frame(colnames_decon)

seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered_DECON<-seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered %>% relocate("OTU_ID","GWNC","GW01",
                                                          "GW02","GW03","GW04","GW05","GW06","GW08","GW09","GW10","GW11")
                                                                                                        
colnames_decon<-colnames(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered_DECON)

decontaminated <- decon(data = seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered_DECON,numb.blanks=1,numb.ind=c(10),taxa=F)
reads_removed<-decontaminated$reads.removed

seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-decontaminated$decon.table
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered$GWNC<-NULL
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-t(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-data.frame(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)

header.true <- function(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered) {
  names(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered) <- as.character(unlist(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[1,]))
  seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[-1,]
}
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-header.true(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)

i <- c(1:447) 
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[ , i] <- apply(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[ , i], 2,  
                                                         function(x) as.numeric(as.character(x)))

#REMOVE SAMPLES THAT NOW HAVE NO READS IN THEM
seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[rowSums(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[])>0,]

#CONVERTING INTO LONG FORMAT
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered <- tibble::rownames_to_column(seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered, "rn")
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered <- gather(seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered, species, reads, 2:448, factor_key=TRUE)
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered$reads<-   as.numeric(seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered$reads)
seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered %>% rename(SampleSite_time = rn)

write.csv(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered,"WP2C_Gwash_ASV_table_wide_filtered.csv")

seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered_transpose<-t(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered)

##SPLITTING BY SPECIES
Anguilla_anguilla<-seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered %>% filter(species == "Anguilla_anguilla")
Sardina_pilchardus<-seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered %>% filter(species == "Sardina_pilchardus")

#Venn diagram
#seqtabNoC_WP2C_Gwash_18S_Venn<-seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean[!(seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Gwash_18S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Gwash_18S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Gwash_venn_diagramm.png',
#output=TRUE)

#NORMALISING READS
seqtabNoC_WP2C_Gwash_18S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[c(1:447)])
seqtabNoC_WP2C_Gwash_18S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Gwash_18S_tax_col_sum2)
seqtabNoC_WP2C_Gwash_18S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Gwash_18S_tax_col_sum2)

seqtabNoC_WP2C_Gwash_18S_tax_normalised<- merge(seqtabNoC_WP2C_Gwash_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Gwash_18S_tax_col_sum2")],seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Gwash_18S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Gwash_18S_tax_wide_filtered[c(1:447)])
seqtabNoC_WP2C_Gwash_18S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Gwash_18S_tax_col_sum3)

seqtabNoC_WP2C_Gwash_18S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Gwash_18S_tax_col_sum3)

seqtabNoC_WP2C_Gwash_18S_tax_normalised<- merge(seqtabNoC_WP2C_Gwash_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Gwash_18S_tax_col_sum2")],seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Gwash_18S_tax_top<- merge(seqtabNoC_WP2C_Gwash_18S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Gwash_18S_tax_col_sum3")],seqtabNoC_WP2C_Gwash_18S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Gwash_18S_tax_normalised <- transform(seqtabNoC_WP2C_Gwash_18S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Gwash_18S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Gwash_18S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

#PLOTTING TOP Taxon
seqtabNoC_WP2C_Gwash_18S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Gwash_18S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+scale_fill_brewer(palette="Paired")

write.csv(seqtabNoC_WP2C_Gwash_18S_tax_normalised,"NORM_WP2C_Gwash_ASV_table_LONG_filtered.csv")

#PLOTTING PHYLA
seqtabNoC_WP2C_Gwash_18S_tax_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Gwash_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Gwash_18S_tax_phyla<-cbind(seqtabNoC_WP2C_Gwash_18S_tax_top,seqtabNoC_WP2C_Gwash_18S_tax_phyla)

seqtabNoC_WP2C_Gwash_18S_tax_phyla<-seqtabNoC_WP2C_Gwash_18S_tax_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Gwash_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Gwash_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Gwash_18S_tax_col_sum3=seqtabNoC_WP2C_Gwash_18S_tax_col_sum3)

seqtabNoC_WP2C_Gwash_18S_tax_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Gwash_18S_tax_phyla, by="SampleSite_time")
unique(seqtabNoC_WP2C_Gwash_18S_tax_phyla[c("X2")])
seqtabNoC_WP2C_Gwash_18S_tax_phyla$X2<-factor(seqtabNoC_WP2C_Gwash_18S_tax_phyla$X2,levels= c('Arthropoda','Rotifera','Mollusca','Chordata',
                                                                                              'Annelida','Nematoda','Gastrotricha', 'Tardigrada','Cnidaria','Porifera','Echinodermata','Xenacoelomorpha',
                                                                                              'Nemertea','Nematomorpha', 'Acanthocephala <thorny-headed worm>','Chaetognatha','Hemichordata',
                                                                                              'Ctenophora <comb jellyfish phylum>','Bryozoa','Platyhelminthes', 'Micrognathozoa'),ordered=TRUE)                                                                                          

ggplot(data=seqtabNoC_WP2C_Gwash_18S_tax_phyla,aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_manual(values = c("#A6CEE3",#Arthropods
                                                                             "#1F78B4", #Rotifera
                                                                             "#B2DF8A",#Molluscs
                                                                             "#33A02C",#Chordates
                                                                             "#FB9A99",#annelids
                                                                             "#E31A1C",#nematodes
                                                                             "#FDBF6F",#gastrotricha
                                                                             "#FF7F00",#tardigrada
                                                                             "#CAB2D6",#Cnideria
                                                                             "#6A3D9A",#porifera
                                                                             "white",#echinodermata
                                                                             "white",
                                                                             "white","white","white","white"))

#PLOTTING TOP PHYLA
seqtabNoC_WP2C_Gwash_18S_tax_top_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Gwash_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Gwash_18S_tax_top_phyla<-cbind(seqtabNoC_WP2C_Gwash_18S_tax_top,seqtabNoC_WP2C_Gwash_18S_tax_top_phyla)

seqtabNoC_WP2C_Gwash_18S_tax_top_phyla<-seqtabNoC_WP2C_Gwash_18S_tax_top_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Gwash_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Gwash_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Gwash_18S_tax_col_sum3=seqtabNoC_WP2C_Gwash_18S_tax_col_sum3)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2C_Gwash_18S_tax_top_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Gwash_18S_tax_top_phyla, by="SampleSite_time")

seqtabNoC_WP2C_Gwash_18S_tax_top_phyla %>%
  mutate(X2 = fct_reorder(X2, desc(seqtabNoC_WP2C_Gwash_18S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")

#CONVERTING NORMALISED LONG FORMAT INTO SampleSite_time FORMAT
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide<-seqtabNoC_WP2C_Gwash_18S_tax_normalised
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide$seqtabNoC_WP2_18S_tax_col_sum2 <- NULL 

seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide<-seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide))

names(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide <- seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide,"NORM_WP2C_Gwash_ASV_table_wide_filtered.csv")

#Investigating lentic species
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide_transposed<-t(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide)
seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide_transposed<-data.frame(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide_transposed)
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

#GW01
GW01_lentic<-seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide_transposed$GW01==0),]
GW01_lentic<-setDT(GW01_lentic, keep.rownames = TRUE)[]
GW01_lentic <- gather(GW01_lentic, SampleSite_time, normalised_reads, 2:12, factor_key=TRUE)
GW01_lentic<-data.frame(GW01_lentic)
GW01_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake")],GW01_lentic, by="SampleSite_time")
GW01_lentic<-na.omit(GW01_lentic)

pGW_01<-ggplot(data=GW01_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")
#ggplotly(pGL_01)

write.csv(GW01_lentic,"lentic_Gwash_18S.csv")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
GW01_lentic_wide<-GW01_lentic
GW01_lentic_wide$Dist_from_lake <- NULL 

GW01_lentic_wide<-GW01_lentic_wide %>% rename(species = rn)

GW01_lentic_wide<-GW01_lentic_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

GW01_lentic_wide <- spread(GW01_lentic_wide,SampleSite_time,normalised_reads)
GW01_lentic_wide[is.na(GW01_lentic_wide)] <- 0
GW01_lentic_wide <- data.frame(t(GW01_lentic_wide))

names(GW01_lentic_wide) <- as.matrix(GW01_lentic_wide[1, ])
GW01_lentic_wide <- GW01_lentic_wide[-1, ]
GW01_lentic_wide[] <- lapply(GW01_lentic_wide, function(x) type.convert(as.character(x)))

lentic_shannon_index<-diversity(GW01_lentic_wide, index = "shannon")
lentic_shannon_index<-data.frame(lentic_shannon_index)
ID <- rownames(lentic_shannon_index)
lentic_shannon_index_Gwash <- cbind(lentic_shannon_index,ID)
write.csv(lentic_shannon_index_Gwash,"lentic_shannon_index_Gwash.csv")

ggplot(data=lentic_shannon_index_Gwash, aes(x=ID, y=lentic_shannon_index, group=1)) +
  geom_line()+
  geom_point()

#PERMANOVA ON TAXON
NORM_WP2C_Gwash_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Gwash_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Dist_from_lake","river_width_m")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:4)]
adonis_data<-adonis_data[,-c(1:4)]

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$river_width_m*
                           adonis_data_meta$Dist_from_lake,
                           method="bray",perm=999)
print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_18S_ID <- read.csv("NORM_WP2C_Gwash_ASV_table_wide_filtered.csv", header = TRUE)
#abund_table_18S_ID<-abund_table_18S_ID[!(abund_table_18S_ID$X=="GWNC"),]

abund_table_18S <- abund_table_18S_ID[,-1]
com <- abund_table_18S[,col(abund_table_18S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_18S_ID[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")

p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Gwash_18S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_18S_WP2C_Gwash<-merge(shannon_index,simpson_index,by="ID")
alpha_18S_WP2C_Gwash<-merge(alpha_18S_WP2C_Gwash,species_count,by="ID")

write.csv(alpha_18S_WP2C_Gwash,"alpha_18S_WP2C_Gwash.csv")

#________________________________________________________

#WP2C_Skan BLAST  TAXONOMY ANALYSIS

#________________________________________________________
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered <- read_csv("WP2C_Skan_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_18S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2C_Skan_18S_tax_wide_filtered <- spread(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[is.na(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered)] <- 0
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered))

names(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[1, ])
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered <- seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[-1, ]
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered, function(x) type.convert(as.character(x)))

#FILTERING BY NEGATIVE CONTROL
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-t(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered)
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-data.frame(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered)
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-setDT(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered, keep.rownames = TRUE)[]
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_18S_tax_wide_filtered %>% rename(OTU_ID = rn)

colnames_decon<-colnames(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered)
colnames_decon<-data.frame(colnames_decon)

seqtabNoC_WP2C_Skan_18S_tax_wide_filtered_DECON<-seqtabNoC_WP2C_Skan_18S_tax_wide_filtered %>% relocate("OTU_ID","SKNC","SK01",
                                                                                                          "SK02","SK03","SK04","SK05","SK06","SK07","SK08","SK09","SK10","SK11")

colnames_decon<-colnames(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered_DECON)

decontaminated <- decon(data = seqtabNoC_WP2C_Skan_18S_tax_wide_filtered_DECON,numb.blanks=1,numb.ind=c(11),taxa=F)
reads_removed<-decontaminated$reads.removed

seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-decontaminated$decon.table
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered$SKNC<-NULL
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-t(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered)
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-data.frame(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered)

header.true <- function(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered) {
  names(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered) <- as.character(unlist(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[1,]))
  seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[-1,]
}
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-header.true(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered)

i <- c(1:173) 
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[ , i] <- apply(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[ , i], 2,  
                                                          function(x) as.numeric(as.character(x)))

#REMOVE SAMPLES THAT NOW HAVE NO READS IN THEM
seqtabNoC_WP2C_Skan_18S_tax_wide_filtered<-seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[rowSums(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[])>0,]

#CONVERTING INTO LONG FORMAT
seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Skan_18S_tax_wide_filtered
seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered <- tibble::rownames_to_column(seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered, "rn")
seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered <- gather(seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered, species, reads, 2:174, factor_key=TRUE)
seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered$reads<-   as.numeric(seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered$reads)
seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered<-seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered %>% rename(SampleSite_time = rn)

write.csv(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered,"WP2C_Skan_ASV_table_wide_filtered.csv")

#RELATIVE AMPLICON FREQUENCY ABUNDANCE
ggplot(seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered, aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

#ABSOLUTE AMPLICON ABUNDANCE
#ggplot(seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean , aes(x = SampleSite_time, fill = species, y = reads)) + 
# geom_bar(stat = "identity", colour = "black")+
#theme(axis.text.x = element_text(angle = 90),legend.position = "none")
##scale_fill_brewer(palette="Paired")

##SPLITTING BY SPECIES
#Alosa_sapidissima<-seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered %>% filter(species == "Alosa_sapidissima")

Anguilla_anguilla_sum_reads<-Anguilla_anguilla %>% 
  group_by(SampleSite_time) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time_sum =SampleSite_time,reads_sum= reads)

ggplot(data=Anguilla_anguilla_sum_reads, aes(x=SampleSite_time_sum, y=reads_sum))+
  geom_col()

#Venn diagram
#seqtabNoC_WP2C_Skan_18S_Venn<-seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean[!(seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean$reads==0),]
#hist(seqtabNoC_WP2C_Skan_18S_Venn$reads, breaks=1000,xlim=c(0,50000))

#E01_Venn<-seqtabNoC_WP2C_Skan_18S_Venn %>% filter(SampleSite_code == "E01")
#E01_Venn<-E01_Venn$species

#E02_Venn<-seqtabNoC_WP2C_Skan_18S_Venn %>% filter(SampleSite_code == "E02")
#E02_Venn<-E02_Venn$species

#E03_Venn<-seqtabNoC_WP2C_Skan_18S_Venn %>% filter(SampleSite_code == "E03")
#E03_Venn<-E03_Venn$species

#E04_Venn<-seqtabNoC_WP2C_Skan_18S_Venn %>% filter(SampleSite_code == "E04")
#E04_Venn<-E04_Venn$species

#E05_Venn<-seqtabNoC_WP2C_Skan_18S_Venn %>% filter(SampleSite_code == "E05")
#E05_Venn<-E05_Venn$species


#venn.diagram(
# x = list(E01_Venn, E02_Venn,E03_Venn,E04_Venn,E05_Venn),
#category.names = c("E01" , "E02","E03","E04","E05"),
#filename = 'WP2C_Skan_venn_diagramm.png',
#output=TRUE)

#NORMALISING READS
seqtabNoC_WP2C_Skan_18S_tax_col_sum2<-rowSums(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[c(1:173)])
seqtabNoC_WP2C_Skan_18S_tax_col_sum2<-data.frame(seqtabNoC_WP2C_Skan_18S_tax_col_sum2)
seqtabNoC_WP2C_Skan_18S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2C_Skan_18S_tax_col_sum2)

seqtabNoC_WP2C_Skan_18S_tax_normalised<- merge(seqtabNoC_WP2C_Skan_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Skan_18S_tax_col_sum2")],seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2C_Skan_18S_tax_col_sum3<-colSums(seqtabNoC_WP2C_Skan_18S_tax_wide_filtered[c(1:173)])
seqtabNoC_WP2C_Skan_18S_tax_col_sum3<-data.frame(seqtabNoC_WP2C_Skan_18S_tax_col_sum3)

seqtabNoC_WP2C_Skan_18S_tax_col_sum3$species <- rownames(seqtabNoC_WP2C_Skan_18S_tax_col_sum3)

seqtabNoC_WP2C_Skan_18S_tax_normalised<- merge(seqtabNoC_WP2C_Skan_18S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2C_Skan_18S_tax_col_sum2")],seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2C_Skan_18S_tax_top<- merge(seqtabNoC_WP2C_Skan_18S_tax_col_sum3[, c("species", "seqtabNoC_WP2C_Skan_18S_tax_col_sum3")],seqtabNoC_WP2C_Skan_18S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2C_Skan_18S_tax_normalised <- transform(seqtabNoC_WP2C_Skan_18S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2C_Skan_18S_tax_col_sum2)

ggplot(seqtabNoC_WP2C_Skan_18S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

#PLOTTING TOP Taxon
seqtabNoC_WP2C_Skan_18S_tax_top %>%
  mutate(species = fct_reorder(species, desc(seqtabNoC_WP2C_Skan_18S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = species, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")+scale_fill_brewer(palette="Paired")

write.csv(seqtabNoC_WP2C_Skan_18S_tax_normalised,"NORM_WP2C_Skan_ASV_table_LONG_filtered.csv")

#PLOTTING PHYLA
seqtabNoC_WP2C_Skan_18S_tax_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Skan_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Skan_18S_tax_phyla<-cbind(seqtabNoC_WP2C_Skan_18S_tax_top,seqtabNoC_WP2C_Skan_18S_tax_phyla)

seqtabNoC_WP2C_Skan_18S_tax_phyla<-seqtabNoC_WP2C_Skan_18S_tax_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Skan_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Skan_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Skan_18S_tax_col_sum3=seqtabNoC_WP2C_Skan_18S_tax_col_sum3)

seqtabNoC_WP2C_Skan_18S_tax_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Skan_18S_tax_phyla, by="SampleSite_time")
unique(seqtabNoC_WP2C_Skan_18S_tax_phyla[c("X2")])
seqtabNoC_WP2C_Skan_18S_tax_phyla$X2<-factor(seqtabNoC_WP2C_Skan_18S_tax_phyla$X2,levels= c('Arthropoda','Rotifera','Mollusca','Chordata',
                                                                                            'Annelida','Nematoda','Gastrotricha', 'Tardigrada','Cnidaria','Porifera','Echinodermata','Xenacoelomorpha',
                                                                                            'Nemertea','Nematomorpha', 'Acanthocephala <thorny-headed worm>','Chaetognatha','Hemichordata',
                                                                                            'Ctenophora <comb jellyfish phylum>','Bryozoa','Platyhelminthes', 'Micrognathozoa'),ordered=TRUE)                                                                                          

ggplot(data=seqtabNoC_WP2C_Skan_18S_tax_phyla,aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_manual(values = c("#A6CEE3",#Arthropods
                                                                             "#1F78B4", #Rotifera
                                                                             "#B2DF8A",#Molluscs
                                                                             "#33A02C",#Chordates
                                                                             "#FB9A99",#annelids
                                                                             "#E31A1C",#nematodes
                                                                             "#FDBF6F",#gastrotricha
                                                                             "#FF7F00",#tardigrada
                                                                             "#CAB2D6",#Cnideria
                                                                             "#6A3D9A",#porifera
                                                                             "#FFFF99",#echinodermata
                                                                             "white",
                                                                             "white","white","white","white"))

#PLOTTING TOP PHYLA
seqtabNoC_WP2C_Skan_18S_tax_top_phyla <- data.frame(do.call('rbind', strsplit(as.character(seqtabNoC_WP2C_Skan_18S_tax_top$species),';',fixed=TRUE)))
seqtabNoC_WP2C_Skan_18S_tax_top_phyla<-cbind(seqtabNoC_WP2C_Skan_18S_tax_top,seqtabNoC_WP2C_Skan_18S_tax_top_phyla)

seqtabNoC_WP2C_Skan_18S_tax_top_phyla<-seqtabNoC_WP2C_Skan_18S_tax_top_phyla %>% 
  group_by(SampleSite_time,X2) %>% 
  summarize(reads=sum(reads),seqtabNoC_WP2C_Skan_18S_tax_col_sum3=sum(seqtabNoC_WP2C_Skan_18S_tax_col_sum3)) %>%
  rename(SampleSite_time=SampleSite_time,X2=X2, reads= reads,seqtabNoC_WP2C_Skan_18S_tax_col_sum3=seqtabNoC_WP2C_Skan_18S_tax_col_sum3)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

seqtabNoC_WP2C_Skan_18S_tax_top_phyla <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "SampleSite_code")],seqtabNoC_WP2C_Skan_18S_tax_top_phyla, by="SampleSite_time")

seqtabNoC_WP2C_Skan_18S_tax_top_phyla %>%
  mutate(X2 = fct_reorder(X2, desc(seqtabNoC_WP2C_Skan_18S_tax_col_sum3))) %>%
  ggplot(aes(x = SampleSite_time, fill = X2, y = reads)) + 
  geom_bar(position="fill",stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_brewer(palette="Paired")



#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide<-seqtabNoC_WP2C_Skan_18S_tax_normalised
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$read_filter <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$Days <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$Date_sampled <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$SampleSite_code <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$WP <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$Seq_length <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$reads <- NULL 
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide$seqtabNoC_WP2_18S_tax_col_sum2 <- NULL 


seqtabNoC_WP2C_Skan_18S_tax_normalised_wide<-seqtabNoC_WP2C_Skan_18S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2C_Skan_18S_tax_normalised_wide <- spread(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide[is.na(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide)] <- 0
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide))

names(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide[1, ])
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide <- seqtabNoC_WP2C_Skan_18S_tax_normalised_wide[-1, ]
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide,"NORM_WP2C_Skan_ASV_table_wide_filtered.csv")

#Investigating lentic species
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide_transposed<-t(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide)
seqtabNoC_WP2C_Skan_18S_tax_normalised_wide_transposed<-data.frame(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide_transposed)
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

#SK01
SK01_lentic<-seqtabNoC_WP2C_Skan_18S_tax_normalised_wide_transposed[!(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide_transposed$SK01==0),]
SK01_lentic<-setDT(SK01_lentic, keep.rownames = TRUE)[]
SK01_lentic <- gather(SK01_lentic, SampleSite_time, normalised_reads, 2:13, factor_key=TRUE)
SK01_lentic<-data.frame(SK01_lentic)
SK01_lentic <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Dist_from_lake")],SK01_lentic, by="SampleSite_time")
SK01_lentic<-na.omit(SK01_lentic)

pSK_01<-ggplot(data=SK01_lentic, aes(x=Dist_from_lake, y=normalised_reads, group=rn)) +
  geom_line(aes(color=rn))+geom_point(aes(color=rn))+theme_bw()+theme(legend.position = "none")
#ggplotly(pGL_01)

write.csv(SK01_lentic,"lentic_Skan_18S.csv")

#CONVERTING NORMALISED LONG FORMAT INTO WIDE FORMAT
SK01_lentic_wide<-SK01_lentic
SK01_lentic_wide$Dist_from_lake <- NULL 

SK01_lentic_wide<-SK01_lentic_wide %>% rename(species = rn)

SK01_lentic_wide<-SK01_lentic_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

SK01_lentic_wide <- spread(SK01_lentic_wide,SampleSite_time,normalised_reads)
SK01_lentic_wide[is.na(SK01_lentic_wide)] <- 0
SK01_lentic_wide <- data.frame(t(SK01_lentic_wide))

names(SK01_lentic_wide) <- as.matrix(SK01_lentic_wide[1, ])
SK01_lentic_wide <- SK01_lentic_wide[-1, ]
SK01_lentic_wide[] <- lapply(SK01_lentic_wide, function(x) type.convert(as.character(x)))

lentic_shannon_index<-diversity(SK01_lentic_wide, index = "shannon")
lentic_shannon_index<-data.frame(lentic_shannon_index)
ID <- rownames(lentic_shannon_index)
lentic_shannon_index_Skan <- cbind(lentic_shannon_index,ID)
write.csv(lentic_shannon_index_Skan,"lentic_shannon_index_Skan.csv")

ggplot(data=lentic_shannon_index_Skan, aes(x=ID, y=lentic_shannon_index, group=1)) +
  geom_line()+
  geom_point()

#PERMANOVA ON TAXON
NORM_WP2C_Skan_ASV_table_wide_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv")
adonis_data<-NORM_WP2C_Skan_ASV_table_wide_filtered
adonis_data_samples<-adonis_data$X1
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = X1)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","WP","Dist_from_lake","river_width_m")],adonis_data, by="SampleSite_time")

adonis_data<-na.omit(adonis_data)
adonis_data_meta<- adonis_data[c(1:4)]
adonis_data<-adonis_data[,-c(1:4)]

adon.results_WP2<-adonis(adonis_data ~ adonis_data_meta$river_width_m*
                           adonis_data_meta$Dist_from_lake,
                         method="bray",perm=999)
print(adon.results_WP2)

#PLOTTING BETA DIVERSITY NMDS
abund_table_18S_ID <- read.csv("NORM_WP2C_Skan_ASV_table_wide_filtered.csv", header = TRUE)
#abund_table_18S_ID<-abund_table_18S_ID[!(abund_table_18S_ID$X=="SKNC"),]

abund_table_18S <- abund_table_18S_ID[,-1]
com <- abund_table_18S[,col(abund_table_18S)]
m_com <- as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", try = 1000, trymax = 1000, k = 3)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds))
#add columns to data frame 
data.scores$SampleSite_time = abund_table_18S_ID[,1]

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","SampleSite_code","Days","Dist_from_lake")],data.scores, by="SampleSite_time")

p<-ggplot(data.scores_meta_data, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + 
  geom_point(size = 2, aes(colour = Dist_from_lake))+
  theme_bw()
ggplotly(p)


#ALPHA DIVERSITY METRICS
shannon_index<-diversity(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide, index = "shannon")
shannon_index<-data.frame(shannon_index)
ID <- rownames(shannon_index)
shannon_index <- cbind(shannon_index,ID)

simpson_index<-diversity(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide, index = "simpson")
simpson_index<-data.frame(simpson_index)
ID <- rownames(simpson_index)
simpson_index <- cbind(simpson_index,ID)

species_count <- rowSums(seqtabNoC_WP2C_Skan_18S_tax_normalised_wide >0)
species_count<-data.frame(species_count)
species_count$ID <- rownames(species_count)

alpha_18S_WP2C_Skan<-merge(shannon_index,simpson_index,by="ID")
alpha_18S_WP2C_Skan<-merge(alpha_18S_WP2C_Skan,species_count,by="ID")

write.csv(alpha_18S_WP2C_Skan,"alpha_18S_WP2C_Skan.csv")

#________________________________________________________

#COMPILING WP2 CURATED NORMALISED OTU TABLES (APART FROM WP2B) FOR ALPHA AND BETA DIVERSITY 

#________________________________________________________

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")

alpha_18S_WP2A <- read_csv("alpha_18S_WP2A.csv")
alpha_18S_WP2C_Glatt <- read_csv("alpha_18S_WP2C_Glatt.csv")
alpha_18S_WP2C_Gwash <- read_csv("alpha_18S_WP2C_Gwash.csv")
alpha_18S_WP2C_Towy <- read_csv("alpha_18S_WP2C_Towy.csv")
alpha_18S_WP2C_Skan <- read_csv("alpha_18S_WP2C_Skan.csv")


alpha_18S_WP2<-rbind(alpha_18S_WP2A,alpha_18S_WP2C_Glatt,alpha_18S_WP2C_Gwash,alpha_18S_WP2C_Towy,alpha_18S_WP2C_Skan)
alpha_18S_WP2<-alpha_18S_WP2 %>% rename(SampleSite_time = "ID")

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]  

alpha_18S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Dist_from_lake")],alpha_18S_WP2, by="SampleSite_time")
alpha_18S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "River")],alpha_18S_WP2, by="SampleSite_time")
alpha_18S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Date_sampled")],alpha_18S_WP2, by="SampleSite_time")
alpha_18S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "river_width_m")],alpha_18S_WP2, by="SampleSite_time")
alpha_18S_WP2 <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "pH")],alpha_18S_WP2, by="SampleSite_time")

#Only keep the conwy summer samples that were coordinated with the Gwash, Glatt and Towy
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "01.09.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "02.11.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "04.01.18" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "11.04.18" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "12.06.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "12.08.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "12.10.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "14.02.18" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "14.12.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "18.04.18" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "18.05.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "22.09.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "23.11.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "25.01.18" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "27.04.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "28.03.18" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "29.06.17" , ] 
alpha_18S_WP2<-alpha_18S_WP2[alpha_18S_WP2$Date_sampled != "29.06.17" , ] 

ggplot(alpha_18S_WP2, aes(x=Dist_from_lake, y=shannon_index, color=River)) +
  geom_point() + 
  geom_smooth(method = "lm",fullrange = T)+theme_bw()+scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#D55E00","#0072B2"))

#PLOTTING A POLYNOMIAL FIT, RATHER THAN LINEAR
#ggplot(alpha_18S_WP2, aes(x=Dist_from_lake, y=shannon_index, color=River)) +
# geom_point() + 
#geom_smooth(method = "lm", formula = y ~ poly(x, 2))+theme_bw()

#TESTING A POLYNOMIAL FIT, RATHER THAN LINEAR
#lm_18S_WP2 = lm(shannon_index~I(Dist_from_lake^2)*River, data = alpha_18S_WP2)
#lm_18S_WP2.1 = lm(shannon_index~Dist_from_lake*River, data = alpha_18S_WP2)

#AIC(lm_18S_WP2)
#AIC(lm_18S_WP2.1)
#anova(lm_18S_WP2.1)

lm_18S_WP2.1 = lm(shannon_index~Dist_from_lake*River*river_width_m, data = alpha_18S_WP2)
anova(lm_18S_WP2.1)

step(lm_18S_WP2.1)

lm_18S_WP2.2<-lm(formula = shannon_index ~ Dist_from_lake * River * river_width_m, 
   data = alpha_18S_WP2)
anova(lm_18S_WP2.2)

#BETA DIVERSITY - PERMANOVA ON TAXON
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
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
seqtabNoC_WP2A_18S_tax_normalised_wide<-NORM_WP2_ASV_table_LONG_filtered
seqtabNoC_WP2A_18S_tax_normalised_wide$X1<- NULL 
seqtabNoC_WP2A_18S_tax_normalised_wide$reads <- NULL 

seqtabNoC_WP2A_18S_tax_normalised_wide<-seqtabNoC_WP2A_18S_tax_normalised_wide %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(normalised_reads=sum(normalised_reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, normalised_reads= normalised_reads,)

seqtabNoC_WP2A_18S_tax_normalised_wide <- spread(seqtabNoC_WP2A_18S_tax_normalised_wide,SampleSite_time,normalised_reads)
seqtabNoC_WP2A_18S_tax_normalised_wide[is.na(seqtabNoC_WP2A_18S_tax_normalised_wide)] <- 0
seqtabNoC_WP2A_18S_tax_normalised_wide <- data.frame(t(seqtabNoC_WP2A_18S_tax_normalised_wide))

names(seqtabNoC_WP2A_18S_tax_normalised_wide) <- as.matrix(seqtabNoC_WP2A_18S_tax_normalised_wide[1, ])
seqtabNoC_WP2A_18S_tax_normalised_wide <- seqtabNoC_WP2A_18S_tax_normalised_wide[-1, ]
seqtabNoC_WP2A_18S_tax_normalised_wide[] <- lapply(seqtabNoC_WP2A_18S_tax_normalised_wide, function(x) type.convert(as.character(x)))

write.csv(seqtabNoC_WP2A_18S_tax_normalised_wide,"NORM_WP2_ASV_table_wide_filtered.csv")

adonis_data<-seqtabNoC_WP2A_18S_tax_normalised_wide
adonis_data<-setDT(adonis_data, keep.rownames = TRUE)[]
adonis_data<-data.frame(adonis_data)
adonis_data_samples<-adonis_data$rn
adonis_data_samples<-data.frame(adonis_data_samples)
adonis_data_samples<-adonis_data_samples %>% rename(ID = adonis_data_samples)
adonis_data<-adonis_data %>% rename(SampleSite_time = rn)

lofresh_metadata_WP2_3 <- read_csv("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/metadata/lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3 = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]
adonis_data<- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","River","Blank","Date_sampled","Dist_from_lake","river_width_m")],adonis_data, by="SampleSite_time")
adonis_data<-na.omit(adonis_data)

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

adonis_data<-adonis_data[!is.na(adonis_data$SampleSite_time), ]
adonis_data_meta<- adonis_data[c(1:6)]
adonis_data<-adonis_data[,-c(1:6)]
adonis_data<-data.frame(adonis_data)
adonis_data_meta<-data.frame(adonis_data_meta)

adon.results_WP2<-adonis2(adonis_data ~ adonis_data_meta$River*adonis_data_meta$Dist_from_lake*adonis_data_meta$river_width_m,
                         method="bray",perm=999)
print(adon.results_WP2)

pairwise.adonis_River<-pairwise.adonis(adonis_data,adonis_data_meta$River)


#PLOTTING BETA DIVERSITY NMDS
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
data.scores_meta_data <- merge(lofresh_metadata_WP2_3[, c("SampleSite_time","Blank","Dist_from_lake","River")],data.scores, by="SampleSite_time")

gg <- merge(data.scores_meta_data,aggregate(cbind(NMDS1.x=NMDS1,NMDS2.y=NMDS2)~River,data.scores_meta_data,mean),by="River")
ggplot(gg, aes(x = NMDS1, y = NMDS2,colour = Dist_from_lake)) + scale_shape_manual(values=c(15,16,17,8,18))+
  geom_point(size = 4, aes(shape=River,colour = Dist_from_lake))+
  theme_bw()+  theme(axis.title.x=element_blank())+  theme(axis.title.y=element_blank())+ 
 geom_point(aes(x=NMDS1.x,y=NMDS2.y),size=3,colour="red")+
 geom_segment(aes(x=NMDS1.x, y=NMDS2.y, xend=NMDS1, yend=NMDS2))

NORM_WP2A_ASV_table_LONG_filtered <- read_csv("NORM_WP2A_ASV_table_LONG_filtered.csv")
NORM_WP2C_Glatt_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Glatt_ASV_table_LONG_filtered.csv")
NORM_WP2C_Gwash_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Gwash_ASV_table_LONG_filtered.csv")
NORM_WP2C_Towy_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Towy_ASV_table_LONG_filtered.csv")
NORM_WP2C_Skan_ASV_table_LONG_filtered <- read_csv("NORM_WP2C_Skan_ASV_table_LONG_filtered.csv")

######## Turnover and nestedness Conwy ########
setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
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

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
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

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
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

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
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

setwd("C:/Users/bsp81d/OneDrive - Bangor University/LOFRESH/WP2_analysis/18S_SILVAngs")
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
                             beta.measures.sor_Towy,beta.measures.sor_Skan)

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
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#D55E00","#0072B2"))
ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_nest,colour=River)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#D55E00","#0072B2"))
ggplot(beta.measures.sor_all, aes(x=(Dist_from_lake_DIFF), y=Beta_Distance,colour=River)) +
  geom_point(size=2, shape=23)+theme_bw()+
  geom_smooth()+scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=c("#56B4E9","#F0E442","#E69F00","#D55E00","#0072B2"))


#COMPILING LENTIC SHANNON INDEXES
lentic_shannon_index_Glatt<- read_csv("lentic_shannon_index_Glatt.csv")
lentic_shannon_index_Gwash <- read_csv("lentic_shannon_index_Gwash.csv")
lentic_shannon_index_Skan <- read_csv("lentic_shannon_index_Skan.csv")
lentic_shannon_index_Towy <- read_csv("lentic_shannon_index_Towy.csv")

lentic_shannon_index_ALL<-rbind(lentic_shannon_index_Glatt,lentic_shannon_index_Gwash,lentic_shannon_index_Skan,lentic_shannon_index_Towy)
lentic_shannon_index_ALL<-lentic_shannon_index_ALL %>% rename(SampleSite_code = ID)
lentic_shannon_index_ALL <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_code","Dist_from_lake","River")],lentic_shannon_index_ALL, by="SampleSite_code")

ggplot(data=lentic_shannon_index_ALL, aes(x=Dist_from_lake, y=lentic_shannon_index, group=River)) +
  geom_line(aes(color=River))+
  geom_point(aes(color=River))+
  theme_bw()