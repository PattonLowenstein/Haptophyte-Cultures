# Kyle Mayers's Cultures from 2017(?) 
# Screened a LOBset of the first 30 or so samples to get a normal phase rt database, and have now run that through. 

setwd("C:/Users/Daniel Lowenstein/Desktop/Haptophyte Paper/")

library(tidyverse)

sample_details <- read.csv("Mayers_Cultures_2017_Sample_List_Combined.csv")

load("Mayers_Combined_Long.RData")

Mayers_Normal_QC_Corrected <- read.csv("Mayers_Normal_6IPLs_QC_Corrected.csv")
Mayers_Normal_Not_6IPL <- read.csv("Mayers_Normal_Cultures_Sorted.csv")
# load("Mayers_Hummel_Long.RData")
LOB_viewdata(Mayers_Normal_Not_6IPL)

#################
## prelim input and variable addition
#################
Mayers_Normal_long_Not_6IPL <- Mayers_Normal_Not_6IPL %>% 
  filter(!grepl(pattern = "FALSE", final_cut)) %>% 
  select(-elem_formula, -peakgroup_mz, 
       -LOBdbase_ppm_match,
       -peakgroup_mzmin, -peakgroup_mzmax, -peakgroup_rtmin, -peakgroup_rt,
       -peakgroup_rtmax, -xcms_peakgroup, -CAMERA_pseudospectrum, -LOBdbase_frag_ID, 
       -LOBdbase_exact_parent_neutral_mass, 
       -major_adduct, -C1, -C1x, -C2a, -C2b, -C3c, 
       -C3f, -C3r, -C4, -C5, -C6a, -C6b, -casecodes, -iso_C3r.match_ID, 
       -iso_C3f.match_ID, -iso_C3c.match_ID, 
       -DNPPE_Factor, 
       -Flag, -DBase_DNPPE_RF, -Flag, -DBase_DNPPE_RF, -code, -ms2_file, -X1, -yesno) %>% 
  gather(key = Sample_ID,
         value = peak_area,
         -compound_name,
         -lipid_class, -FA_total_no_C, -FA_total_no_DB, -LOBdbase_mz, 
         -degree_oxidation, -species, -match_ID, -final_cut) %>% 
  filter(species != "MGDG", 
         species != "DGDG", 
         species != "SQDG", 
         species != "PG", 
         species != "PC", 
         species != "PE",
         species != "DNPPE", 
         lipid_class!= "pigment")

Mayers_Normal_long_QC_Corrected <- Mayers_Normal_QC_Corrected %>% 
  filter(!grepl(pattern = "FALSE", final_cut)) %>% 
  select(-elem_formula, -peakgroup_mz, 
         -LOBdbase_ppm_match,
         -peakgroup_mzmin, -peakgroup_mzmax, -peakgroup_rtmin, -peakgroup_rt,
         -peakgroup_rtmax, -xcms_peakgroup, -CAMERA_pseudospectrum, -LOBdbase_frag_ID, 
         -LOBdbase_exact_parent_neutral_mass, 
         -major_adduct, -C1, -C1x, -C2a, -C2b, -C3c, 
         -C3f, -C3r, -C4, -C5, -C6a, -C6b, -casecodes, -iso_C3r.match_ID, 
         -iso_C3f.match_ID, -iso_C3c.match_ID, 
         -DNPPE_Factor, 
         -Flag, -DBase_DNPPE_RF, -Flag, -DBase_DNPPE_RF, -code, -ms2_file, -X1) %>% 
  gather(key = Sample_ID,
         value = peak_area,
         -compound_name,
         -lipid_class, -FA_total_no_C, -FA_total_no_DB, -LOBdbase_mz, 
         -degree_oxidation, -species, -match_ID, -final_cut)

Mayers_Normal_long <- rbind(Mayers_Normal_long_Not_6IPL, Mayers_Normal_long_QC_Corrected)

save(Mayers_Normal_long, file = "Mayers_Normal_long.RData")
###########################
## Quick DNPPE Response Check
###########################

QC_DNPPE <- Mayers_Normal_long %>% 
  filter(Experiment == "QC", species == "DNPPE") %>% 
  mutate(Middle = median(peak_area),
         Average = mean(peak_area),
         SD = sd(peak_area),
         SE = SD/Average,
         difference_multiple = peak_area/Average)


all_DNPPE <- Mayers_Normal_long %>% filter(compound_name == "DNPPE")%>% 
  mutate(Middle = median(peak_area),
         Average = mean(peak_area),
         SD = sd(peak_area),
         SE = SD/Average,
         difference_multiple = peak_area/Middle)
  
QC_PC <- Mayers_Normal_long %>% filter(Experiment == "QC", compound_name == "PC 32:0" | species == "DNPPE")

ggplot(all_DNPPE, aes(x = Sample_ID, y = difference_multiple))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))

###############
# Adding in the Hummel Set (and saving it and loading it into the envt)
###############
# quickly save Mayers_Normal_Long to reload and recombine

sample_details <- read.csv("Mayers_Cultures_2017_Sample_List_Combined.csv")
load("Mayers_Normal_long.RData")
load("Mayers_Hummel_Long.RData")

 

Mayers_Normal_long <- Mayers_Normal_long %>% 
  mutate(QE_Number = as.numeric(str_extract(string = Sample_ID, pattern = "(?<=QE00)\\d+")),
         Treatment = NA,
         Blank_Corrected = peak_area,
         Max_Blank_Peak = NA,
         Average_Blank_Peak = NA,
         Method = paste("Normal")) %>% 
  filter(species != "dLCB_GSL_No_FA_OH",
         compound_name != "Chl_c2", 
         compound_name != "DivinylChl_a", 
         compound_name != "Fucoxanthin", 
         compound_name != "Pheophytin_a", 
         compound_name != "Pheophytin_b2",
         species != "hGSL") %>% #and filter out the compounds found in both sets--taking the values from Hummel
  select(-final_cut, -match_ID)

Mayers_Normal_long$Eluent_Sequence <- sapply(Mayers_Normal_long$Sample_ID, function(x) paste(sample_details$Normal_Eluents[which(sample_details$File_Name == x)]))

Mayers_Hummel_Long <- Mayers_Hummel_Long %>% 
  dplyr::rename(peak_area = Peak_Area) %>% 
  filter(species != "DGCC") # now that I've done the comparison, I'm filtering it out

# check again
intersection <- intersect(Mayers_Hummel_Long$compound_name, Mayers_Normal_long$compound_name)

Mayers_Combined_Long <- rbind(Mayers_Hummel_Long, Mayers_Normal_long)

Mayers_Combined_Long$G_species <- sapply(Mayers_Combined_Long$Sample_ID, function(x) paste(sample_details$G_species[which(sample_details$File_Name == x)]))

Mayers_Combined_Long$Experiment <- sapply(Mayers_Combined_Long$Sample_ID, function(x) paste(sample_details$Experiment[which(sample_details$File_Name == x)]))

Mayers_Combined_Long$Treatment <- sapply(Mayers_Combined_Long$Sample_ID, function(x) paste(sample_details$Treatment[which(sample_details$File_Name == x)]))

Mayers_Combined_Long$Time <- sapply(Mayers_Combined_Long$Sample_ID, function(x) paste(sample_details$Time[which(sample_details$File_Name == x)]))

Mayers_Combined_Long$Replicate <- sapply(Mayers_Combined_Long$Sample_ID, function(x) paste(sample_details$Replicate[which(sample_details$File_Name == x)]))

Mayers_Combined_Long$Cells_per_sample <- sapply(Mayers_Combined_Long$Sample_ID, function(x) as.numeric(paste(sample_details$Cells_per_sample[which(sample_details$File_Name == x)])))

Mayers_Combined_Long$E_hux_Strain <- sapply(Mayers_Combined_Long$Sample_ID, function(x) paste(sample_details$E_hux_Strain[which(sample_details$File_Name == x)]))

Mayers_Combined_Long$Light_Dark <- sapply(Mayers_Combined_Long$Sample_ID, function(x) paste(sample_details$Light_Dark[which(sample_details$File_Name == x)]))

Mayers_Combined_Long$Blank_Corrected[Mayers_Combined_Long$Blank_Corrected < 0] <- 0
Mayers_Combined_Long$Blank_Corrected[is.na(Mayers_Combined_Long$Blank_Corrected)] <- 0



Mayers_Combined_Long <- Mayers_Combined_Long %>%
  filter(Sample_ID != "QE003414") %>% # removing an air blank I forgot
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Experiment, "_", 
                                     Treatment, "_", 
                                     Time, "_", 
                                     Replicate),
         QE_Number = as.numeric(str_extract(string = Sample_ID, pattern = "(?<=QE00)\\d+")),
         Peak_Per_Cell = Blank_Corrected/Cells_per_sample)

save(Mayers_Combined_Long, file = "Mayers_Combined_Long.RData")
#load("Mayers_Combined_Long.RData")


########################
## Looking at Blanks
########################

Blanks <- Mayers_Combined_Long %>% 
  filter(Treatment == "blank")

Blank_FFAs <- Blanks %>% filter(species == "FFA")%>% 
  group_by(Experiment) %>% 
  mutate(Total_FFAs = sum(peak_area)) %>% 
  ungroup()

ggplot(Blank_FFAs, aes(x = compound_name, y = peak_area, color = G_species))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))

FFAs <- Mayers_Combined_Long %>% 
  filter(compound_name == "FFA 16:0" |
           compound_name == "FFA 16:1" |
           compound_name == "FFA 18:0" |
           compound_name == "FFA 18:1") %>% 
  group_by(Sample_ID, Treatment, G_species, Experiment) %>% 
  mutate(Total_FFAs = sum(peak_area)) %>% 
  ungroup()

ggplot(FFAs)+
  geom_point(position = position_dodge(width = 0.75), aes(x = Sample_ID, y = peak_area, color = Experiment, shape = Treatment))+
  geom_point(position = position_dodge(width = 0.75), aes(x = Sample_ID, y = Blank_Corrected, color = Experiment, shape = Treatment))+
  facet_grid(rows = vars(compound_name))+
  theme(axis.text.x = element_text(angle = 90))+
  #scale_y_continuous(breaks = seq(0, 1e10, 1e9), limits = c(0, 1e10))+
  ggtitle(label = "All Experiments - FFAs - blank check - zoomed out")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(0, 16, 2, 3, 4, 15, 5, 7, 8))



TAGs <- Mayers_Combined_Long %>% 
  filter(compound_name == "TAG 46:1"|
           compound_name == "TAG 46:2"|
           compound_name == "TAG 48:1"|
           compound_name == "TAG 48:2"|
           compound_name == "TAG 50:1" |
           compound_name == "TAG 50:2") %>% 
  group_by(Sample_ID, Treatment, G_species, Experiment) %>% 
  mutate(Total_TAGs = sum(peak_area)) %>% 
  ungroup()

ggplot(TAGs)+
  geom_point(position = position_dodge(width = 0.75), aes(x = Sample_ID, y = peak_area, color = Eluent_Sequence, shape = Treatment))+
  geom_point(position = position_dodge(width = 0.75), aes(x = Sample_ID, y = Blank_Corrected, color = "blue", shape = Treatment))+  facet_grid(rows = vars(compound_name))+
  theme(axis.text.x = element_text(angle = 90))+
  scale_shape_manual(values = c(0, 1, 16, 2, 3, 4, 15, 5, 7, 8))+
  scale_y_continuous(breaks = seq(0, 5e7, 1e7), limits = c(0, 5e7))+
  scale_color_manual(values = c("1"="#66CD00", "2"="#FF3030","3"="#0000FF", "4"="#FF1493","5"="#228B22","blue"="#000000"))

Light_Dark_TAGS <- TAGs %>% filter(Experiment == "Light_Dark", FA_total_no_C > 44, FA_total_no_C < 55)

ggplot(Light_Dark_TAGS, aes(x = compound_name, y = peak_area, color = Treatment, shape = Light_Dark))+
  geom_point(position = position_dodge(width = 0.5))+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"))+
  #scale_y_continuous(breaks = seq(0, 1e10, 1e9), limits = c(0, 1e10))+
  ggtitle(label = "L/D Experiment - TAGs")


test <- Mayers_Combined_Long %>% 
  group_by(Sample_ID, Method) %>% 
  summarize(Total_peak_area = sum(peak_area)) %>% 
  ungroup()

ggplot(test, aes(x = Sample_ID, y = Total_peak_area, color = Method))+
  geom_point(position = position_dodge(width = 0.5), size = 2)+
  #facet_grid(rows = vars(compound_name), scales = "free")+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"))+
  ggtitle(label = "peak area across samples")

PG_QC <- Mayers_Combined_Long %>% 
  filter(species == "PG", Experiment == "QC")

ggplot(PG_QC, aes(Sample_ID, peak_area, color = Method))+
  geom_point()+
  ylim(0, 1.5e9)


###############
## DNPPE Correction (no standards run with Normal Phase samples)
###############
dnppe_averages <- Mayers_Combined_Long %>% 
  filter(compound_name == "DNPPE", 
         Treatment != "QC", 
         Treatment != "blank", 
         Treatment != "Chloroplasts",
         Treatment != "Whole_Cells",
         Treatment != "Plastid_Supernatant") %>% 
  mutate(DNPPE_Average = mean(peak_area), 
         DNPPE_Median = median(peak_area),
         DNPPE_Factor = peak_area/DNPPE_Average, 
         Standard_Deviation = sd(peak_area), 
         Sample_Deviation = abs(DNPPE_Median - peak_area)/Standard_Deviation)

dnppe_averages_Normal <- Mayers_Combined_Long %>% 
  filter(compound_name == "DGCC", 
         Treatment != "QC", 
         Treatment != "blank", 
         Treatment != "Chloroplasts",
         Treatment != "Whole_Cells",
         Treatment != "Plastid_Supernatant",
         Method == "Normal") %>% 
  mutate(DNPPE_Average = mean(peak_area), 
         DNPPE_Median = median(peak_area),
         DNPPE_Factor = peak_area/DNPPE_Average, 
         Standard_Deviation = sd(peak_area), 
         Sample_Deviation = abs(DNPPE_Median - peak_area)/Standard_Deviation)

dnppe_averages_Hummel <- Mayers_Combined_Long %>% 
  filter(compound_name == "DNPPE", 
         Treatment != "QC", 
         Treatment != "blank", 
         Treatment != "Chloroplasts",
         Treatment != "Whole_Cells",
         Treatment != "Plastid_Supernatant",
         Method == "Hummel") %>% 
  mutate(DNPPE_Average = mean(peak_area), 
         DNPPE_Median = median(peak_area),
         DNPPE_Factor = peak_area/DNPPE_Average, 
         Standard_Deviation = sd(peak_area), 
         Sample_Deviation = abs(DNPPE_Median - peak_area)/Standard_Deviation)


ggplot(dnppe_averages, aes(Sample_ID, peak_area, color = Method))+
  geom_point()+
  geom_text(aes(label = QE_Number))

ggplot(dnppe_averages_Hummel, aes(x = peak_area))+
  geom_histogram(bins = 50)

DNPPE_Corrected <- Mayers_Combined_Long

DNPPE_Corrected$DNPPE_Factor <- sapply(DNPPE_Corrected$Sample_ID, function(x) as.numeric(paste(dnppe_averages$DNPPE_Factor[which(dnppe_averages$Sample_ID == x)])))

DNPPE_Corrected <- DNPPE_Corrected %>% 
  mutate(peak_cutoff = if_else(Blank_Corrected > 1e6, Blank_Corrected, 0), # cutoff peaks below 1 million
          Corrected_Peak_Area = as.numeric(peak_cutoff) * as.numeric(DNPPE_Factor),
         Corrected_Per_Cell = Corrected_Peak_Area / as.numeric(Cells_per_sample))

#save(DNPPE_Corrected, file = "Mayers_DNPPE_Corrected.RData")
load("Mayers_DNPPE_Corrected.RData")
###############
# NMDS Experiment
###############
library(vegan)

NMDS_data <- read.csv("Mayers_NMDS_experiment_30March2020.csv")
NMDA_limited <- NMDS_data %>% filter(G_species == "E_huxleyi" |
                                       G_species == "I_galbana" |
                                       G_species == "P_parvum" |
                                       G_species == "P_gyrans")

NMDS_1 <- NMDA_limited[, 8:565]
NMDS_2 <- NMDA_limited[, 1:7]
NMDS <- metaMDS(NMDS_1, distance = "bray", k = 2)

co <- c("red", "blue", "green")

shape <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

plot(NMDS$points)

## trying the ggplot way....
data.scores <- as.data.frame(scores(NMDS))
data.scores$G_species <- NMDS_2$G_species
data.scores$Treatment <- NMDS_2$Treatment
data.scores$Time <- NMDS_2$Time


species.scores <- as.data.frame(scores(NMDS, "species"))
species.scores$species <- rownames(species.scores)

ggplot()+
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=G_species,colour=Treatment),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=Time),size=6,hjust=1, vjust = 1.5) +  # add the site labels
  coord_equal() +
  theme_bw()+
  ylim(-0.4, 0.4)+
  xlim(-0.6, 0.5)

######################
# PCA Analysis
######################
library(ggfortify)
library(ggbiplot)

Species_Sums <- aggregate(formula = Corrected_Peak_Area_Per_Cell ~ species + Sample_ID + G_species + Experiment + Treatment + Time, data = DNPPE_Corrected, FUN = sum)
Species_Sums_spread <- Species_Sums %>% spread("species", Corrected_Peak_Area_Per_Cell)
PCA_data <- Species_Sums_spread %>% filter(Treatment == "Replete", G_species != "P_antarctica")# |
                                        #Treatment == "N_Limited" |
                                       # Treatment == "P_Limited",
                                         # G_species == "E_huxleyi"|
                                         #   G_species == "I_galbana" |
                                         #   G_species == "P_parvum" |
                                         #   G_species == "P_gyrans" | 
                                         #  G_species == "C_ericina")

Experiment_PCA <- prcomp(PCA_data[,c(6:10, 12:26)], center = TRUE, scale. = TRUE)
autoplot(Experiment_PCA, data = PCA_data, colour = "G_species", size = 5, label = TRUE, loadings = TRUE, loadings.colour = "blue",
         loadings.label = TRUE)

PCA_data <- DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A") %>% 
  select(Corrected_Per_Cell, Variable_Treatment, species, E_hux_Strain, G_species) %>% 
  group_by(Variable_Treatment, species, E_hux_Strain, G_species) %>% 
  summarize(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  spread("species", Total_Corrected_Per_Cell)

write.csv(PCA_data, file = "Mayers_PCA_by_compound.csv")
PCA_data <- read.csv("Mayers_PCA_by_compound.csv")

PCA_data <- PCA_data[, colSums(PCA_data !=0) > 0]

Experiment_PCA <- prcomp(PCA_data[, 4:51], scale. = TRUE)

autoplot(Experiment_PCA, data = PCA_data, shape = "G_species", colour = "E_hux_Strain", loadings = TRUE, loadings.label = TRUE)

######################
## cluster analysis
######################

library(pheatmap)
library(dendextend)

Cultures_Cluster <-  DNPPE_Corrected %>% 
  filter(Treatment == "Replete",
         Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         G_species != "C_leptoporus") %>% # get rid of single samples (i.e. not run in dup or triplicate)
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% # peak proportion per sample
  ungroup() %>%
  group_by(compound_name, Variable_Treatment) %>% # group by compound name so that double peaks are aggregated within each sample
  mutate(Total_by_species = sum(Corrected_Per_Cell), 
         Proportion_by_species = Total_by_species/Total_Corrected_Per_Cell) %>% # combine double peaks in each sample
  ungroup() %>% 
  select(-LOBdbase_mz,
         -lipid_class,
         -species,
         -FA_total_no_C,
         -FA_total_no_DB,
         -degree_oxidation,
         -peak_area,
         -G_species,
         -Treatment,
         -Experiment,
         -Time,
         -Replicate,
         #-Sample_ID,
         -Cells_per_sample,
         -Total_Corrected_Per_Cell,
         -QE_Number,
         -Eluent_Sequence,
         -Max_Blank_Peak,
         -Average_Blank_Peak,
         -Blank_Corrected,
         -Peak_Proportion,
         -Light_Dark, 
         -E_hux_Strain,
         -Total_by_species,
         -Peak_Per_Cell) %>%# pause here and make sure there are matching sample ID's and
# variable_treatments--there are (2 Samp_ID's per variable treatment- Hummel and Normal
# now summarize and spread
  group_by(compound_name, Variable_Treatment) %>% 
  summarise(Proportion = mean(Proportion_by_species)) %>% 
  spread("Variable_Treatment", "Proportion")

#write.csv(Cultures_Cluster, file = "Mayers_Culture_Cluster.csv")
Cultures_Cluster <- read.csv("Mayers_Culture_Cluster.csv") #import after changing NA's to 0

Cluster <- as.data.frame(Cultures_Cluster)

# transfer column with compound name to row names
rownames(Cluster) <- Cluster[,1] 
Cluster <- Cluster[ ,-1]

# scale it, though we're not actually using this one because log cluster makes it way more obvious
set.seed(3)
scaled_df <- scale(Cluster)

pheatmap(scaled_df, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 18, 
         main = "Mayers Scaled Cluster - 'ward D2' clustering", 
         treeheight_row = 333, treeheight_col = 100,
         legend = TRUE)

###################
## Log Cluster
###################
# attemped to log transform 
log_df <- log10(Cluster)
#test <- replace(log_df, is.nan, values = -20) # this was a test to check out DAG pair clustering

#write.csv(log_df, "Mayers_Log_Cluster.csv") #write and remove -Inf values....
log_df <- read.csv("Mayers_Log_Cluster.csv")
log_df <- as.data.frame(log_df)

rownames(log_df) <- log_df[,1]
log_df <- log_df[ ,-1]

pheatmap(log_df, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 18, 
         main = "Haptophyte Lipid Cluster \n Log Scale - 'ward D2' clustering", 
         treeheight_row = 333, treeheight_col = 100,
         legend = TRUE,
         show_rownames = FALSE)

# testing separating clusters by kmeans - determining correct no of clusters with pamk in fpc packages
library(fpc)
test_kmeans <- pamk(log_df, usepam = FALSE, criterion = "asw")
j <- c()
# testing removing low significance molecules then 
set.seed(3)

for(i in 3:3){
  
  j <- c(j, length(test$molecule))
}
##################
## Variance cluster test
##################
test <- Cluster
# test$compound_name <- rownames(test) # test with only plastid lipids
# test <- test[which(grepl(pattern = "MGDG|DGDG|SQDG|PG", x = test$compound_name)),-53]
test$variance <- apply(test, 1, var)
test$row_median <- apply(test[1:52], 1, median)
test$row_mean <- apply(test[1:52], 1, mean)
test$row_min <- apply(test[1:52], 1, min)
test <- rownames_to_column(test, var = "molecule")


test <- test %>%
  filter(row_min > 0) 


row_median > summary(row_median)[5],
variance < summary(variance)[3],
row_mean > summary(row_mean)[4]

test <- test %>%
  select(-row_median, -variance, -row_mean, -row_min) 
rownames(test) <- test[,1]
test <- test[ ,-1]



summary(test$row_mean)
ggplot(test, aes(x = log(variance)))+
  geom_histogram(bins = 100)




log_test <- log(test)
write.csv(log_test, file = "Plastid_Lipids_Log_Cluster.csv")
log_test <- read.csv("Plastid_Lipids_Log_Cluster.csv")
rownames(log_test) <- log_test[,1]
log_test <- log_test[ ,-1]

# color coding clusters w/annotations, along with 
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
log_cluster <- log_df
hclust_mayers <- hclust(dist(log_cluster), method = "ward.D2")
as.dendrogram(hclust_mayers) %>% 
   plot(horiz = TRUE)

lipids <- cutree(tree = as.dendrogram(hclust_mayers), k = 8)
lipids_column <- data.frame(cluster = lipids)
lipids_column$Cluster <- sapply(lipids_column$cluster, function(x) 
  paste0("Cluster ", x))
lipids_column <- lipids_column[,-1]
lipids_column <- data.frame(Cluster = lipids_column)
row.names(lipids_column) <- rownames(log_cluster)

ggplot(lipids_column, aes(x = Cluster))+
  geom_histogram(stat = "count", binwidth = 1, bins = unique(Cluster))

species_column <- data.frame(Species = rep(NA, length(colnames(log_cluster))))
row.names(species_column) <- colnames(log_cluster)
species_column$row_name <- rownames(species_column)
species_column$Species <- sapply(species_column$row_name, 
                                 function(x) str_extract(string = x, 
                                                         pattern = "^[:upper:]_[:lower:]+(?=_)"))
species_column <- species_column[,-2]
species_column <- data.frame(Species = species_column)
row.names(species_column) <- colnames(log_cluster)

pheatmap(log_cluster, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         #main = "Haptophyte Lipid Cluster \n Log Scale - 'ward D2' clustering", 
         treeheight_row = 50, treeheight_col = 75,
         legend = TRUE,
         annotation_row = lipids_column,
         annotation_col = species_column , 
         show_rownames = TRUE,
         fontsize = 13,
         angle_col = 90
         )



# Testing bootstrap values
# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(log_df, method.hclust="ward.D2",
               method.dist="euclidean")
fit10k <- pvclust(log_df, method.hclust="ward.D2",
               method.dist="euclidean", nboot = 10000, parallel = TRUE)
#compound_cluster <- fit # ran it first with a flipped cluster, so it calc'ed bootstrap values for the compound cluster--stored in separate variable, compound_cluster--flipping it back to do species/strain
seplot(fit10k, identify = TRUE)
test <- plot(fit10k)   # dendogram with p values
# add rectangles around groups highly supported by the data
test%>% pvrect(fit10k, alpha=0.95)




# Dissimilarity matrix
d <- dist(log_df, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "ward.D2" )

# Plot the obtained dendrogram
plot(hc1)

######################
## Quick plots
######################

Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         # E_hux_Strain == "E_huxleyi_374" |
           # E_hux_Strain == "I_galbana" |
           # E_hux_Strain == "P_parvum" | E_hux_Strain == "C_ericina"
           # E_hux_Strain == "P_gyrans"
         ) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment, E_hux_Strain, Replicate, species) %>% 
  summarize(Proportion_by_lipid_class = sum(Peak_Proportion)) %>% 
  #filter(
    # compound_name == "dLCB_GSL_No_FA_OH 37:4 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 38:2 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 38:3 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 38:4 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 38:5 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 39:4 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 39:4 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 39:5 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 39:6 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 41:5 +1O"|
    #        compound_name == "dLCB_GSL_No_FA_OH 41:4 +1O",species == "PDPT"
    #        )  %>% 
  ungroup() %>% 
  filter(species == "PE", Treatment == "Replete")

#%>% 
 # mutate(DAG = paste(str_extract(string = compound_name, pattern = "(?<= )\\d+:\\d+")))


# %>% 
#   group_by(Variable_Treatment, E_hux_Strain, species, FA_total_no_DB, lipid_class, Treatment) %>% 
#   summarize(Proportion_by_species = sum(Peak_Proportion)) %>% 
#   ungroup() %>% 
#   group_by(Variable_Treatment, E_hux_Strain, lipid_class, Treatment) %>% 
#   mutate(Mean_DB_VT = weighted.mean(FA_total_no_DB, Proportion_by_species))
# 
# %>% 
#   ungroup() %>% 
#   group_by(E_hux_Strain, species, FA_total_no_DB) %>% 
#   summarize(Mean_DB = mean(Mean_DB_VT), 
#             SD_DB = sd(Mean_DB_VT))



Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N_Limited", "P_Limited"))

Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("P_carterae",
                                      "C_leptoporus",
                                      "E_huxleyi_1N", 
                                      "E_huxleyi_F", 
                                      "E_huxleyi_370", 
                                      "E_huxleyi_374", 
                                      "E_huxleyi_379", 
                                      "I_galbana", 
                                      "P_parvum", 
                                      "C_ericina", 
                                      "P_globosa", 
                                      "P_globosa_ST", 
                                      "P_antarctica", 
                                      "P_antarctica_ST", 
                                      "P_gyrans"))
# 
# test <- Cultures_GSL[order(Cultures_GSL$compound_name),]
# 
# 
# agg_GSL <- aggregate(formula = Peak_Proportion ~ E_hux_Strain + species + Treatment + Time + Replicate + Variable_Treatment + Method + FA_total_no_DB, 
#                       data = Cultures_GSL,
#                       FUN = sum)


ggplot(Cultures_GSL, aes(x = E_hux_Strain, y = Proportion_by_lipid_class, fill = Replicate))+
  # facet_wrap(~E_hux_Strain)+
  ggtitle("PE")+
  geom_bar(position = position_dodge(width = 0.5), stat = "identity")+
  geom_point(position = position_dodge(width = 0.5))+
  theme(axis.text.x = element_text(angle = 90, size = 10), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 10))+
  xlab("Species")+
  ylab("Peak Proportion Per Cell")+
  #ylim(0,0.4)+
  labs(color = "Species")

+
  geom_text(aes(label = QE_Number))






#########################
## Checking Normal QCs to see if standards are ok to use
## Ended up manually response-correcting eluent sections of the Normal run so that their QCs lined up
########################

Mayers_Normal_6IPL_QC_Peaks_Long <- Mayers_Normal_6IPL_QC_Peaks %>% 
  gather(Species, Peak_Area, -Sample_ID, -Eluent_Chunk) %>%
  filter(Sample_ID != "QE004041" ,
         Sample_ID != "QE004043",
         Sample_ID != "QE004058",
         Sample_ID != "QE004073",
         Sample_ID != "QE004088",
         Sample_ID != "QE004100",
         Sample_ID != "QE004113",
         Sample_ID != "QE003938",
         Sample_ID != "QE003939",
         Sample_ID != "QE003950",
         Sample_ID != "QE004120",
         Sample_ID != "QE004121",
         Sample_ID != "QE004124",
         Sample_ID != "QE004125",
         Sample_ID != "QE004116",
         Sample_ID != "QE004117",
         Sample_ID != "QE004118",
         Sample_ID != "QE004119")

Mayers_Normal_6IPL_QC_Peaks_Long$Eluent_Chunk <- sapply(Mayers_Normal_6IPL_QC_Peaks_Long$Eluent_Chunk, as.factor)

a <- ggplot(Mayers_Normal_6IPL_QC_Peaks_Long)+
  geom_point(aes(x = Sample_ID, y = Peak_Area, color = Eluent_Chunk))+
  facet_grid(rows = vars(Species), scales = "free")+
  theme(axis.text.x = element_text(angle = 90), 
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.75, linetype = "solid"))+
  scale_color_manual(values = c("1"="#66CD00", "2"="#FF3030","3"="#0000FF", "4"="#FF1493","5"="#228B22","6"="#000000"))



Mayers_Normal_6IPL_QC_Props_Long <- Mayers_Normal_6IPL_QC_Proportions %>% 
  gather(Species, Peak_Area, -Sample_ID, -Eluent_Chunk)%>%
  filter(Sample_ID != "QE004041" ,
         Sample_ID != "QE004043",
         Sample_ID != "QE004058",
         Sample_ID != "QE004073",
         Sample_ID != "QE004088",
         Sample_ID != "QE004100",
         Sample_ID != "QE004113",
         Sample_ID != "QE003938",
         Sample_ID != "QE003939",
         Sample_ID != "QE003950",
         Sample_ID != "QE004120",
         Sample_ID != "QE004121",
         Sample_ID != "QE004124",
         Sample_ID != "QE004125")

Mayers_Normal_6IPL_QC_Props_Long$Eluent_Chunk <- sapply(Mayers_Normal_6IPL_QC_Props_Long$Eluent_Chunk, as.factor)

b <- ggplot(Mayers_Normal_6IPL_QC_Props_Long)+
  geom_point(aes(x = Sample_ID, y = Peak_Area, color = Eluent_Chunk))+
  facet_grid(rows = vars(Species), scales = "free")+
  theme(axis.text.x = element_text(angle = 90), 
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray"), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.75, linetype = "solid"))+
  scale_color_manual(values = c("1"="#66CD00", "2"="#FF3030","3"="#0000FF", "4"="#FF1493","5"="#228B22","6"="#000000"))

grid_arrange_shared_legend(a, b, ncol = 2, position = "right")

# samples between Kyle's samples and standards
%>%
  filter(Sample_ID != "QE004041" ,
         Sample_ID != "QE004043",
         Sample_ID != "QE004058",
         Sample_ID != "QE004073",
         Sample_ID != "QE004088",
         Sample_ID != "QE004100",
         Sample_ID != "QE004113")


########################
## DAG pairs cluster 
########################

DAG_Pairs_Cluster <- DNPPE_Corrected %>% 
  filter(lipid_class == "IP_DAG", Experiment == "Nutrient_Limitation", Treatment == "Replete", degree_oxidation == 0) %>% 
  mutate(DAG = paste(str_extract(string = compound_name, pattern = "(?<= )\\d+:\\d+"))) %>% 
  group_by(Variable_Treatment, DAG) %>% 
  mutate(DAG_sums = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(DAG_Proportion = DAG_sums/sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  select(-LOBdbase_mz,
         -lipid_class,
         -species,
         -FA_total_no_C,
         -FA_total_no_DB,
         -degree_oxidation,
         -peak_area,
         -G_species,
         -Treatment,
         -Experiment,
         -Time,
         -Replicate,
         -Sample_ID,
         -Cells_per_sample,
         -QE_Number,
         -Eluent_Sequence,
         -Max_Blank_Peak,
         -Average_Blank_Peak,
         -Blank_Corrected,
         -Light_Dark, 
         -E_hux_Strain,
         -Peak_Per_Cell,
         -Method,
         -DNPPE_Factor,
         -Corrected_Peak_Area,
         -Corrected_Per_Cell,
         -DAG_sums,
         -compound_name) %>% 
  group_by(DAG, Variable_Treatment) %>% 
  summarise(DAG_Proportion = sum(DAG_Proportion)) %>% 
  mutate(Log10_DAGs = log10(DAG_Proportion))

DAG_Pairs_Cluster$Log10_DAGs <- sapply(DAG_Pairs_Cluster$Log10_DAGs, function(x) replace(x, is.infinite(x), -7.5))

DAG_Pairs_Cluster <- DAG_Pairs_Cluster %>% 
  select(-DAG_Proportion)%>% 
  spread("Variable_Treatment", "Log10_DAGs")

DAG_Pairs_Cluster <- as.data.frame(DAG_Pairs_Cluster)

DAG_Pairs_Cluster$DAG <- sapply(DAG_Pairs_Cluster$DAG, as.character)
rownames(DAG_Pairs_Cluster) <- DAG_Pairs_Cluster[,1]
DAG_Pairs_Cluster <- DAG_Pairs_Cluster[ ,-1]

scaled_dags <- scale(DAG_Pairs_Cluster)

pheatmap(scaled_dags, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 18, 
         main = "Mayers Log10 Scaled DAG Pair Cluster - 'ward D2' clustering", 
         treeheight_row = 333, treeheight_col = 100,
         legend = TRUE)

######################
## DAG Pairs
######################

DAG_Pairs <- DNPPE_Corrected %>% 
  filter(lipid_class == "IP_DAG", 
         Experiment == "Nutrient_Limitation", 
         Treatment == "Replete", 
         degree_oxidation == 0,
         E_hux_Strain != "C_leptoporus") %>% 
  mutate(DAG = paste(str_extract(string = compound_name, pattern = "(?<= )\\d+:\\d+"))) %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_per_sample = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment, FA_total_no_DB, E_hux_Strain, species, Total_per_sample) %>% 
  summarize(DB_Total = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(FA_total_no_DB, Variable_Treatment, species, E_hux_Strain) %>% 
  summarise(DB_Proportion = DB_Total/Total_per_sample) %>% 
  ungroup() %>% 
  group_by(FA_total_no_DB, species, E_hux_Strain) %>% 
  summarise(DB_Proportion = mean(DB_Proportion))

DAG_Pairs$DAG <- factor(DAG_Pairs$DAG, levels = c("26:0",  "27:0",  "27:1",  "28:0",  "28:1",  "29:0",  "29:1",  "30:0",  "30:1",  "30:2",  "30:3", "30:4",  "30:5",  "31:0",  "31:1",  "31:2",  "31:3",  "32:0",  "32:1",  "32:2",  "32:3",  "32:4", "32:5",  "32:6",  "32:7",  "32:8",  "33:1",  "33:2",  "33:3",  "33:4",  "33:5",  "33:7",  "34:0", "34:1",  "34:2",  "34:3",  "34:4",  "34:5",  "34:6",  "34:7",  "34:8",  "34:9",  "35:3",  "35:4", "35:5",  "35:6",  "35:7",  "35:8",  "36:1",  "36:2",  "36:3",  "36:4",  "36:5",  "36:6", "36:7",  "36:8",  "36:9", "36:10",  "37:2",  "37:3",  "37:5",  "37:6",  "37:7",  "37:8",  "38:1", "38:3",  "38:4",  "38:5",  "38:6",  "38:7",  "38:8",  "38:9",  "38:10",  "39:5",  "39:6",  "39:7",  "39:8", "40:5",  "40:6",  "40:7",  "40:8",  "40:9", "40:10", "40:11",  "41:7",  "41:8", "42:6",  "42:7",  "42:8",  "42:9",  "42:10", "42:11",  "44:11", "44:12"))

DAG_Pairs$FA_total

ggplot(DAG_Pairs, aes(x = FA_total_no_DB, y = DB_Proportion, fill = species))+
  geom_bar(stat = "identity", position = "stack")+
  facet_grid(rows = vars(E_hux_Strain))+
  theme(axis.text.x = element_text(angle = 90))


########################
## Rearranging Variable Treatment to look at Hummel and Normal DGCC side by side
## General result is that it looks....ok....but really reiterates our need for thorough and frequent standards
########################

hummel_eluent_grab <- DNPPE_Corrected %>% 
  filter(species == "DGCC", Experiment == "Nutrient_Limitation", Method == "Hummel") %>% 
  select(Variable_Treatment, Method, compound_name, Corrected_Peak_Area, FA_total_no_C, FA_total_no_DB, E_hux_Strain, Eluent_Sequence) %>% 
  group_by(Variable_Treatment, Method, compound_name, FA_total_no_C, FA_total_no_DB, E_hux_Strain, Eluent_Sequence) %>% 
  summarize(Total_Peak_Area = sum(Corrected_Peak_Area))%>% 
  spread("Method", "Total_Peak_Area") 

normal_eluent_grab <- DNPPE_Corrected %>% 
  filter(species == "DGCC", Experiment == "Nutrient_Limitation", Method == "Normal") %>% 
  select(Variable_Treatment, Method, compound_name, Corrected_Peak_Area, FA_total_no_C, FA_total_no_DB, E_hux_Strain, Eluent_Sequence) %>% 
  group_by(Variable_Treatment, Method, compound_name, FA_total_no_C, FA_total_no_DB, E_hux_Strain, Eluent_Sequence) %>% 
  summarize(Total_Peak_Area = sum(Corrected_Peak_Area))%>% 
  spread("Method", "Total_Peak_Area") 

hummel_eluent_grab$Hummel_Eluents <- hummel_eluent_grab$Eluent_Sequence
hummel_eluent_grab <- hummel_eluent_grab %>% select(-Eluent_Sequence)

normal_eluent_grab$Normal_Eluents <- normal_eluent_grab$Eluent_Sequence
normal_eluent_grab <- normal_eluent_grab %>% select(-Eluent_Sequence)



Comparo_DGCC <- DNPPE_Corrected %>% 
  filter(species == "DGCC", Experiment == "Nutrient_Limitation") %>% 
  select(Variable_Treatment, Method, compound_name, Corrected_Peak_Area, FA_total_no_C, FA_total_no_DB, E_hux_Strain) %>% 
  group_by(Variable_Treatment, Method, compound_name, FA_total_no_C, FA_total_no_DB, E_hux_Strain) %>% 
  summarize(Total_Peak_Area = sum(Corrected_Peak_Area)) %>% 
  spread("Method", "Total_Peak_Area") 

Comparo_DGCC <- merge(hummel_eluent_grab, Comparo_DGCC, by = c("Variable_Treatment", "compound_name", "FA_total_no_C", "FA_total_no_DB", "E_hux_Strain", "Hummel"))
Comparo_DGCC <- merge(normal_eluent_grab, Comparo_DGCC, by = c("Variable_Treatment", "compound_name", "FA_total_no_C", "FA_total_no_DB", "E_hux_Strain", "Normal"))

Comparo_DGCC[is.na(Comparo_DGCC)] <- 0
Comparo_DGCC$FA_total_no_DB <- sapply(Comparo_DGCC$FA_total_no_DB, as.factor)


Comparo_DGCC <- Comparo_DGCC %>% 
  mutate(Ratio = Hummel/Normal)

ggplot(Comparo_DGCC, aes(x = (Hummel), y = (Normal), color = Hummel_Eluents))+
  geom_point()+
  geom_text(aes(label = paste0(FA_total_no_C, ":", FA_total_no_DB)))+
  geom_abline(aes(intercept = 0, slope = 1))+
  facet_wrap(~E_hux_Strain, scales = "free")+
  ggtitle("Free axis limits, Hummel Eluent Chunks")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"))

#+
  xlim(0, 2.5e8)+
  ylim(0, 2.5e8)
# +
#   geom_abline(aes(intercept = 0, slope = 1))+
#   geom_abline(aes(intercept = -1e8, slope = 1))


#

#+
# facet_grid(rows = vars(FA_total_no_DB), cols = vars(FA_total_no_C))

########################
# Diff stats by lipid species----don't use this one--go down to where I said "Redoing diff stats by species"
########################


No_double_peaks <- DNPPE_Corrected %>% 
    filter(Experiment == "Nutrient_Limitation", 
           species != "DNPPE", 
           Replicate != "Box4_A",
           Replicate != "Box4_B",
           Replicate != "A",
           G_species == "I_galbana" |
           G_species == "P_parvum" |
           G_species == "P_gyrans" |
           E_hux_Strain == "E_huxleyi_374",
           Time == 0) %>% 
    group_by(
      #compound_name,
             Variable_Treatment, 
             species,
             # FA_total_no_C, 
             # FA_total_no_DB,
             # degree_oxidation,
             Sample_ID,
             Time,
             Experiment,
             Treatment, 
             E_hux_Strain, lipid_class
             ) %>% 
    summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
    ungroup() %>% 
    group_by(Variable_Treatment) %>% 
    mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
           Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
    ungroup()
    

No_double_peaks_majors <- No_double_peaks %>% 
    filter(lipid_class != "pigment" ,
           lipid_class != "ubiquinone" ,
           !grepl("plast", lipid_class),
           lipid_class != "IP_MAG",
           lipid_class != "FFA", 
           !grepl("GSL", lipid_class))

# separate repletes from treatment samples, 
Repletes <- No_double_peaks %>% # changed to nodoublepeaksMAJORS for some analyses
  filter(Treatment == "Replete") %>% 
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
    group_by(
      #compound_name,
             Variable_Treatment, 
             E_hux_Strain, 
             species,
             #,
             # FA_total_no_C, 
             # FA_total_no_DB,
             # degree_oxidation
             lipid_class
             ) %>% 
    summarise(Replete_Proportions = list(Peak_Proportion),
              Replete_Average = mean(Peak_Proportion),
              Replete_St_Dev = sd(Peak_Proportion)) %>% 
    ungroup() %>% 
    select(-Variable_Treatment)


  
diff_tests <- No_double_peaks %>% # changed to nodoublepeaksMAJORS for some analyses
  filter(Time == 0, Treatment != "Replete")%>%
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    #compound_name,
           Variable_Treatment, 
           E_hux_Strain, 
           species,
           # FA_total_no_C, 
           # FA_total_no_DB,
           # degree_oxidation,
           lipid_class,
           Treatment) %>%
  summarise(Treatment_Proportions = list(Peak_Proportion),
            Treatment_Average = mean(Peak_Proportion),
            Treatment_St_Dev = sd(Peak_Proportion))
  

diff_stats <- merge(Repletes, diff_tests, by = c("E_hux_Strain", "species", "lipid_class")) 
# "FA_total_no_C", "FA_total_no_DB", "degree_oxidation", "lipid_class" "species", , "compound_name"


diff_stats <- diff_stats %>% 
  group_by(Variable_Treatment, species, lipid_class) %>% #species, compound_name
  mutate(Fold_Change = Treatment_Average / Replete_Average,
         Log2_FC = log2(Fold_Change),
         Relative_Difference = Treatment_Average - Replete_Average,
         p_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$p.value,
         t_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$statistic,
         significant = if_else(p_value <= 0.05, true = paste("Significant"), false = paste("Not Significant"))) %>% 
  ungroup() 

diff_stats_all <- rbind(diff_stats, minor_diff_stats)
diff_stats_no_pigs <- diff_stats %>% filter(lipid_class == "TAG" |
                                              lipid_class == "IP_DAG")

diff_stats_no_pigs$species <- sapply(diff_stats_no_pigs$species, function (x) str_replace(x, "dLCB_GSL_No_FA_OH", "GSL"))
diff_stats_no_pigs$lipid_class <- sapply(diff_stats_no_pigs$lipid_class, function (x) str_replace(x, "IP_DAG", "Polar Lipids"))
species_labels <- c("E. huxleyi 374", "I. galbana", "P. gyrans", "P. parvum")
names(species_labels) <- c("E_huxleyi_374", "I_galbana", "P_gyrans", "P_parvum")
treatments <- c("N-Limited", "P-Limited")
names(treatments) <- c("N_Limited", "P_Limited")

only_sigs <- diff_stats %>% 
  filter(significant == "Significant")

ggplot(diff_stats_no_pigs, aes(x = lipid_class, y = Relative_Difference, fill = species))+
  geom_bar(stat = "identity", color = "black", position = "stack")+
  facet_grid(rows = vars(E_hux_Strain),
             cols = vars(Treatment), 
             labeller = labeller(E_hux_Strain = species_labels, Treatment = treatments))+
  # ggtitle("Treatment Proportion Minus Replete Proportion \n Under Nutrient Stress")+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid", 
                                        color = "gray"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))+
  ylab("Mean Difference in Peak Proportion")+
  xlab("Lipid Class")+
  geom_abline(slope = 0, intercept = 0, size = 1)




ggplot(diff_stats, aes(x = Replete_Average, y = Log2_FC, color = species))+
  geom_point()+
  facet_grid(rows = vars(E_hux_Strain),
             cols = vars(Treatment))+
  ggtitle("Fold Changes Under Nutrient Stress")+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"))+
  ylab("Log 2 Fold Change")+
  xlab("Mean Replete Peak Proportion")+
  geom_abline(slope = 0, intercept = 0)

#### test to make sure stats are working
test <- diff_stats %>% 
  filter(E_hux_Strain == "E_huxleyi_374", species == "dLCB_GSL_No_FA_OH")

a <- c(0.052210181845263, 0.0338025263816172, 0.0400582)
b <- c(0.188987774794788, 0.203721099878221, 0.187577755)
test <- wilcox.test(a, b)
###############################
# Diff stats by compound name
###############################


No_double_peaks <- DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         G_species == "I_galbana" |
           G_species == "P_parvum" |
           G_species == "P_gyrans" |
           E_hux_Strain == "E_huxleyi_374",
         Time == 0) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class
  ) %>% 
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup()


# separate repletes from treatment samples, 
Repletes <- No_double_peaks %>% 
  filter(Treatment == "Replete") %>% 
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    E_hux_Strain, 
    species,
    #,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation
    lipid_class
  ) %>% 
  summarise(Replete_Proportions = list(Peak_Proportion),
            Replete_Average = mean(Peak_Proportion),
            Replete_St_Dev = sd(Peak_Proportion)) %>% 
  ungroup() %>% 
  select(-Variable_Treatment)



Replete_Counts <- No_double_peaks %>% 
  filter(Corrected_Per_Cell != 0, Treatment == "Replete") %>% 
  count(Variable_Treatment, E_hux_Strain) %>% 
  group_by(E_hux_Strain) %>% 
  summarise(All_Counts = list(n),
            Average_Count = mean(n),
            SD = sd(n))

Replete_Counts_TAGS <- No_double_peaks %>% 
  filter(Corrected_Per_Cell != 0, Treatment == "Replete", species == "TAG") %>% 
  count(Variable_Treatment, E_hux_Strain) %>% 
  group_by(E_hux_Strain) %>% 
  summarise(All_Counts = list(n),
            Average_Count = mean(n),
            SD = sd(n))

# Replete_Counts <- Repletes %>% filter(Replete_Average != 0) %>% # easy first try to get total number of lipids per species in replete cultures, going to do it again to get stdev as well
#   count(E_hux_Strain)
Unique_Replete_Count <- No_double_peaks %>% # gets the same value as the above "easy first try" to get replete counts, which really counted the number of values where the average peak value was greater than 0 (i.e. any of the replicates had a positive value)
  filter(Corrected_Per_Cell != 0, Treatment == "Replete") %>% 
  group_by(E_hux_Strain) %>% 
  summarize(total = length(unique(compound_name)))

Unique_Replete_TAGs <- No_double_peaks %>% # gets the same value as the above "easy first try" to get replete counts, which really counted the number of values where the average peak value was greater than 0 (i.e. any of the replicates had a positive value)
  filter(Corrected_Per_Cell != 0, Treatment == "Replete", species == "TAG") %>% 
  group_by(E_hux_Strain) %>% 
  summarize(total = length(unique(compound_name)))

diff_tests <- No_double_peaks %>% 
  filter(Time == 0, Treatment != "Replete")%>%
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    E_hux_Strain, 
    species,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation,
    lipid_class,
    Treatment) %>%
  summarise(Treatment_Proportions = list(Peak_Proportion),
            Treatment_Average = mean(Peak_Proportion),
            Treatment_St_Dev = sd(Peak_Proportion))


diff_stats <- merge(Repletes, diff_tests, by = c("E_hux_Strain", "species", "compound_name", "lipid_class")) 
# "FA_total_no_C", "FA_total_no_DB", "degree_oxidation", "lipid_class" "species", , "compound_name"


diff_stats <- diff_stats %>% 
  group_by(Variable_Treatment, species, lipid_class, compound_name) %>% #species, compound_name
  mutate(Fold_Change = Treatment_Average / Replete_Average,
         Log2_FC = log2(Fold_Change),
         Relative_Difference = Treatment_Average - Replete_Average,
         p_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$p.value,
         t_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$statistic,
         significant = if_else(p_value <= 0.05, true = paste("Significant"), false = paste("Not Significant"))) %>% 
  ungroup() 

diff_stats_to_write <- diff_stats %>% select(-Replete_Proportions, -Treatment_Proportions)
#write.csv(diff_stats_to_write, "diff_stats_by_compound_name.csv")

diff_stats_no_pigs <- diff_stats %>% filter(lipid_class == "TAG" |
                                              lipid_class == "IP_DAG")

diff_stats_no_pigs$species <- sapply(diff_stats_no_pigs$species, function (x) str_replace(x, "dLCB_GSL_No_FA_OH", "GSL"))
diff_stats_no_pigs$lipid_class <- sapply(diff_stats_no_pigs$lipid_class, function (x) str_replace(x, "IP_DAG", "Polar Lipids"))
species_labels <- c("E. huxleyi 374", "I. galbana", "P. gyrans", "P. parvum")
names(species_labels) <- c("E_huxleyi_374", "I_galbana", "P_gyrans", "P_parvum")
treatments <- c("N-Limited", "P-Limited")
names(treatments) <- c("N_Limited", "P_Limited")

# only_sigs <- diff_stats %>%
#   filter(significant == "Significant") %>%
#   mutate(PosOrNeg = ifelse(Log2_FC>1, "Positive_Log2FC", "Negative_Log2FC"))
#   ## I did Benjamini-Hochberg Screening by hand--in Mayers Cultures Tables xlsx
#  loading it back in now:

# need to redo only sigs if I"m going to use it again--got rid of the p parvum replicates in this version
only_sigs <- read.csv("only_sigs_benj_hoch.csv")
only_sigs$Log2_FC[which(only_sigs$Log2_FC == "#NAME?")] <- as.numeric(-Inf)
only_sigs$Log2_FC <- as.numeric(only_sigs$Log2_FC)
only_sigs <- only_sigs %>% 
    mutate(PosOrNeg = ifelse(Log2_FC>1, "Positive_Log2FC", "Negative_Log2FC"))

only_sig_counts <- only_sigs %>% 
  count(Variable_Treatment, Treatment, PosOrNeg, sort = TRUE) %>% 
  pivot_wider(names_from = PosOrNeg, values_from = n) %>% 
  mutate(Total_Count = Positive_Log2FC + Negative_Log2FC)

#useless plot
# ggplot(only_sig_counts)+
#   facet_wrap(~Variable_Treatment)+
#   geom_bar(aes(x = Treatment, y = Negative_Log2FC), stat = "identity", width = 1)+
#   coord_polar("y", start = 0)+
#   theme_void() +
#   theme(legend.position="none") +
#   geom_text(aes(x = Treatment, y = Negative_Log2FC, label = Variable_Treatment)) +
#   scale_fill_brewer(palette="Set1")

only_sigs_TAGs <- only_sigs %>% 
  filter(species == "TAG")
only_sig_counts_TAGs <- only_sigs_TAGs %>% 
  count(Variable_Treatment, Treatment, PosOrNeg, sort = TRUE) %>% 
  pivot_wider(names_from = PosOrNeg, values_from = n) %>% 
  mutate(Total_Count = Positive_Log2FC + Negative_Log2FC)

levels(diff_stats$species) <- c(levels(diff_stats$species), "GSL")
diff_stats$species[diff_stats$species == "dLCB_GSL_No_FA_OH"] <- "GSL"
levels(only_sigs$species) <- c(levels(only_sigs$species), "GSL")
only_sigs$species[only_sigs$species == "dLCB_GSL_No_FA_OH"] <- "GSL"

diff_stats_no_pigs <- diff_stats %>% 
  filter(lipid_class != "pigment" ,
         lipid_class != "ubiquinone" ,
         !grepl("plast", lipid_class))

only_sigs_no_pigs <- only_sigs %>% 
  filter(lipid_class != "pigment" ,
         lipid_class != "ubiquinone" ,
         !grepl("plast", lipid_class))

ggplot(diff_stats_no_pigs)+
  geom_point(data = diff_stats, aes(x = log2(Replete_Average), y = log2(Treatment_Average)), shape = 1)+
  geom_point(data = only_sigs_no_pigs, aes(x = log2(Replete_Average), y = log2(Treatment_Average), color = species), shape = 19, size = 3)+
  facet_grid(rows = vars(E_hux_Strain),
             cols = vars(Treatment),
             labeller = labeller(Treatment = treatments, E_hux_Strain = species_labels))+
  geom_abline(aes(slope = 1, intercept = 0))+
  ggtitle("")+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid", 
                                        color = "gray"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))+
  ylab(expression("Log "[2]*" (Peak Proportion)"[italic("Treatment")]))+
  xlab(expression("Log "[2]*" (Peak Proportion)"[italic("Replete")]))+
  # geom_text(data = diff_stats_of_note, aes(label = species, size = 0.1))+
  labs(color = "Species")

### test to make sure stats are working
test <- diff_stats %>% 
  filter(E_hux_Strain == "E_huxleyi_374", species == "dLCB_GSL_No_FA_OH")

a <- c(0.052210181845263, 0.0338025263816172, 0.0400582)
b <- c(0.188987774794788, 0.203721099878221, 0.187577755)
test <- wilcox.test(a, b)

length(setdiff(unique(Repletes$compound_name), unique(diff_tests$compound_name)))
length(setdiff(unique(diff_tests$compound_name), unique(Repletes$compound_name)))

#######
# Aggregate by species to summarize
#######
agg_diff_stats <- diff_stats

agg_diff_stats$species <- ifelse(grepl("PQ", agg_diff_stats$species), 
                                            as.character("Plastoquinones"),
                                            as.character(agg_diff_stats$species))

agg_diff_stats$species <- ifelse(grepl("UQ", agg_diff_stats$species), 
                                            as.character("Ubiquinones"),
                                            as.character(agg_diff_stats$species))


agg_diff_stats$species <- ifelse(grepl("Chl_a|Pheophytin_a|Chlide_a|DivinylChl_a|Pheophytin_a2", agg_diff_stats$species), 
                                            as.character("Chl_a_etc"),
                                            as.character(agg_diff_stats$species))

agg_diff_stats$species <- ifelse(grepl("Chl_b|Pheophytin_b2", agg_diff_stats$species), 
                                            as.character("Chl_b_etc"),
                                            as.character(agg_diff_stats$species))

agg_diff_stats$species <- ifelse(grepl("Chl_c1|Chl_c2|Chl_c3|Chl_c2 MGDG 14:0_14:0|Chl_c2 MGDG 18:4_14:0", agg_diff_stats$species), 
                                            as.character("Chl_c_etc"),
                                            as.character(agg_diff_stats$species))

agg_diff_stats$species <- ifelse(grepl("19prime_but_fuco|19prime_hex_fuco|Alloxanthin|Crocoxanthin|Ddc_dd|Diatoxanthin|Echinenone|Lutein|Neo_Nosxanthin|Prasinoxanthin|Violaxanthin|Zeaxanthin", agg_diff_stats$species), 
                                            as.character("Minor Xanthophylls"),
                                            as.character(agg_diff_stats$species))

agg_diff_stats$species <- ifelse(grepl("Fucoxanthin", agg_diff_stats$species), 
                                            as.character("Fucoxanthin"),
                                            as.character(agg_diff_stats$species))

agg_diff_stats$species <- ifelse(grepl("Carotene", agg_diff_stats$species), 
                                            as.character("Carotenes"),
                                            as.character(agg_diff_stats$species))
agg_diff_stats <- agg_diff_stats %>% 
  group_by(Variable_Treatment, E_hux_Strain, Treatment, species) %>% 
  summarize(Replete_proportion_by_headgroup = sum(Replete_Average),
            Treatment_proportion_by_headgroup = sum(Treatment_Average),
            Differential_Response = Treatment_proportion_by_headgroup - Replete_proportion_by_headgroup,
            Agg_Fold_Change = Treatment_proportion_by_headgroup / Replete_proportion_by_headgroup,
            Log2_Agg_FC = log2(Agg_Fold_Change)) 

ggplot(agg_diff_stats, aes(x = species, 
                           y = Differential_Response, 
                           color = species, 
                           fill = species))+
  geom_bar(stat = "identity", position = "stack")+
  facet_grid(cols = vars(Treatment),
             rows = vars(E_hux_Strain))+
  ggtitle("Treatment Minus Replete Response \n Under Nutrient Stress by DB number")+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"))


################
#### Formattable - N- and P-limitation tables by species - not useful anymore
################
library(data.table)
library(formattable)
library(tidyverse)
library(gridExtra)
library(formattable)

library("htmltools")
library("webshot")    

export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

sign_formatter <- formatter("span", 
                            style = x ~ style(color = ifelse(x > 0, "green", 
                                                             ifelse(x < 0, "red", "black"))))
# Up and down arrow with greater than comparison from the vignette

significance_formatter <- formatter("span", 
                                    style = x ~ style(color = ifelse(as.numeric(x) < 0.05, "green", "black")), 
                                    x ~ icontext(ifelse(as.numeric(x) < 0.05, "star", ""), x))


# E hux
E_hux_n_limited <- diff_stats %>% 
  filter(Treatment == "N_Limited", E_hux_Strain == "E_huxleyi_374") %>% 
  select(E_hux_Strain, species, Replete_Average, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment)%>% 
  mutate(Replete_Average = formatC(Replete_Average, format = "e", digits = 2),
         Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("Replete Average" = Replete_Average,
         "N-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

E_hux_p_limited <- diff_stats %>% 
  filter(Treatment == "P_Limited", E_hux_Strain == "E_huxleyi_374") %>% 
  select(E_hux_Strain, species, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment) %>% 
  mutate(Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("P-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

eh_tab_1 <- formattable(E_hux_n_limited,
            align = c("c", "c", "c", "c", "c", "c"),
            list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter,
                 "Replete Average" = formatter("span", 
                                               style = x ~ style(color = "blue",
                                                                 font.weight = ifelse(as.numeric(x) > 0.001, "bold", "")))))

eh_tab_2 <- formattable(E_hux_p_limited,
                        align = c("c", "c", "c", "c", "c", "c"),
                        list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter))

export_formattable(eh_tab_1, "Ehux_n_limited.png")
export_formattable(eh_tab_2, "Ehux_p_limited.png")


# I galb
I_galbana_n_limited <- diff_stats %>% 
  filter(Treatment == "N_Limited", E_hux_Strain == "I_galbana") %>% 
  select(E_hux_Strain, species, Replete_Average, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment)%>% 
  mutate(Replete_Average = formatC(Replete_Average, format = "e", digits = 2),
         Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("Replete Average" = Replete_Average,
         "N-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

I_galbana_p_limited <- diff_stats %>% 
  filter(Treatment == "P_Limited", E_hux_Strain == "I_galbana") %>% 
  select(E_hux_Strain, species, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment) %>% 
  mutate(Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3))) %>% 
  rename("P-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

I_gal_tab_1 <- formattable(I_galbana_n_limited,
            align = c("c", "c", "c", "c", "c", "c"),
            list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter,
            "Replete Average" = formatter("span", 
                                          style = x ~ style(color = "blue",
                                                            font.weight = ifelse(as.numeric(x) > 0.001, "bold", "")))))
I_gal_tab_2 <- formattable(I_galbana_p_limited,
            align = c("c", "c", "c", "c", "c", "c"),
            list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter))

export_formattable(I_gal_tab_1, "I_gal_n_lim.png")
export_formattable(I_gal_tab_2, "I_gal_p_lim.png")

## P parvum
P_parvum_n_limited <- diff_stats %>% 
  filter(Treatment == "N_Limited", E_hux_Strain == "P_parvum") %>% 
  select(E_hux_Strain, species, Replete_Average, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment) %>% 
  mutate(Replete_Average = formatC(Replete_Average, format = "e", digits = 2),
         Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("Replete Average" = Replete_Average,
         "N-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

P_parvum_p_limited <- diff_stats %>% 
  filter(Treatment == "P_Limited", E_hux_Strain == "P_parvum") %>% 
  select(E_hux_Strain, species, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment)%>% 
  mutate(Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("P-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

P_parv_tab1 <- formattable(P_parvum_n_limited,
            align = c("c", "c", "c", "c", "c", "c"),
            list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter,
            "Replete Average" = formatter("span", 
                                          style = x ~ style(color = "blue",
                                                            font.weight = ifelse(as.numeric(x) > 0.001, "bold", "")))))
P_parv_tab2 <- formattable(P_parvum_p_limited,
            align = c("c", "c", "c", "c", "c", "c"),
            list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter))

export_formattable(P_parv_tab1, "P_parv_n_lim.png")
export_formattable(P_parv_tab2, "P_parv_p_lim.png")

# P gyr
P_gyrans_n_limited <- diff_stats %>% 
  filter(Treatment == "N_Limited", E_hux_Strain == "P_gyrans") %>% 
  select(E_hux_Strain, species, Replete_Average, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment)%>% 
  mutate(Replete_Average = formatC(Replete_Average, format = "e", digits = 2),
         Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("Replete Average" = Replete_Average,
         "N-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

P_gyrans_p_limited <- diff_stats %>% 
  filter(Treatment == "P_Limited", E_hux_Strain == "P_gyrans") %>% 
  select(E_hux_Strain, species, Treatment_Average, Treatment, Log2_FC, p_value, -E_hux_Strain, -Treatment) %>% 
  mutate(Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("P-Limited Average" = Treatment_Average,
         "Log 2 Fold Change" = Log2_FC)

P_gyr_tab1 <- formattable(P_gyrans_n_limited,
            align = c("l", "c", "c", "c", "c", "c"),
            list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter,
            "Replete Average" = formatter("span", 
                                          style = x ~ style(color = "blue",
                                                            font.weight = ifelse(as.numeric(x) > 0.001, "bold", "")))))
P_gyr_tab2 <- formattable(P_gyrans_p_limited,
            align = c("l", "c", "c", "c", "c", "c"),
            list("Log 2 Fold Change" = sign_formatter,
                 "p_value" = significance_formatter))

export_formattable(P_gyr_tab1, "P_gyr_n_lim.png")
export_formattable(P_gyr_tab2, "P_gyr_p_lim.png")

################
# Export nut lim tables split by class
################

export_n_lim <- diff_stats %>% 
  filter(Treatment == "N_Limited") %>% 
  select(E_hux_Strain, 
         species, 
         Replete_Average, 
         Replete_St_Dev, 
         Treatment_Average, 
         Treatment_St_Dev, 
         Treatment, 
         Log2_FC, 
         p_value, 
         -Treatment)%>% 
  mutate(Replete_Average = formatC(Replete_Average, format = "e", digits = 2),
         Replete_St_Dev = paste0("±", formatC(Replete_St_Dev, format = "e", digits = 2)),
         
         Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         Treatment_St_Dev = paste0("±", formatC(Treatment_St_Dev, format = "e", digits = 2)),
         Log2_FC = round(Log2_FC, digits = 3),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3))) %>% 
  rename("Replete Average" = Replete_Average,
         "R_SD" = Replete_St_Dev,
         "N-Limited Average" = Treatment_Average,
         "N-lim SD" = Treatment_St_Dev,
         "N-Lim Log 2 Fold Change" = Log2_FC,
         "Lipid" = species,
         "Species" = E_hux_Strain, 
         "N-lim p value" = p_value)

export_p_lim <- diff_stats %>% 
  filter(Treatment == "P_Limited") %>% 
  select(E_hux_Strain, species, Treatment_Average, Treatment_St_Dev, Treatment, Log2_FC, p_value, -Treatment)%>% 
  mutate(Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         Treatment_St_Dev = paste0("±", formatC(Treatment_St_Dev, format = "e", digits = 2)),
         Log2_FC = round(Log2_FC, digits = 3),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("P-Limited Average" = Treatment_Average,
         "P-lim_SD" = Treatment_St_Dev,
         "P-lim Log 2 Fold Change" = Log2_FC,
         "Lipid" = species,
         "Species" = E_hux_Strain, 
         "P-lim p value" = p_value)


###################
## Pigments Nut Lim Table Export
###################
pigments_n_lim <- export_n_lim %>% filter( 
Lipid == "Carotene" | 
  Lipid == "Chl_a" | 
  Lipid == "Chl_b" | 
  Lipid == "Chl_c1" | 
  Lipid == "Chl_c2" | 
  Lipid == "Chl_c2 MGDG 14:0_14:0" | 
  Lipid == "Chl_c2 MGDG 18:4_14:0" | 
  Lipid == "Chl_c3" | 
  Lipid == "Chlide_a" | 
  Lipid == "DivinylChl_a" |
  Lipid == "Pheophytin_a" | 
  Lipid == "Pheophytin_a2" | 
  Lipid == "Pheophytin_b2")

pigments_p_lim <- export_p_lim %>% filter( 
  Lipid == "Carotene" | 
    Lipid == "Chl_a" | 
    Lipid == "Chl_b" | 
    Lipid == "Chl_c1" | 
    Lipid == "Chl_c2" | 
    Lipid == "Chl_c2 MGDG 14:0_14:0" | 
    Lipid == "Chl_c2 MGDG 18:4_14:0" | 
    Lipid == "Chl_c3" | 
    Lipid == "Chlide_a" | 
    Lipid == "DivinylChl_a" |
    Lipid == "Pheophytin_a" | 
    Lipid == "Pheophytin_a2" | 
    Lipid == "Pheophytin_b2")

pigments_table_export <- merge(pigments_n_lim, pigments_p_lim, by = c("Species", "Lipid"))
write.csv(pigments_table_export, file = "pigments_table_export.csv")
#############
## pigment separation/export details, probably not useful anymore
#############
Pigments_tab1 <- formattable(E_hux_pigs_n_lim,
                          align = c("l", "c", "c", "c", "c", "c"),
                          list("Log 2 Fold Change" = sign_formatter,
                               "p_value" = significance_formatter,
                               "Replete Average" = formatter("span", 
                                                             style = x ~ style(color = "blue",
                                                                               font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                    "bold", 
                                                                                                    "")))))

Pigments_tab2 <- formattable(E_hux_pigs_p_lim,
                          align = c("c", "c", "c", "c", "c", "c"),
                          list("Log 2 Fold Change" = sign_formatter,
                               "p_value" = significance_formatter))

Pigments_tab3 <- formattable(I_gal_pigs_n_lim,
                             align = c("l", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter,
                                  "Replete Average" = formatter("span", 
                                                                style = x ~ style(color = "blue",
                                                                                  font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                       "bold", 
                                                                                                       "")))))

Pigments_tab4 <- formattable(I_gal_pigs_p_lim,
                             align = c("c", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter))

Pigments_tab5 <- formattable(P_gyr_pigs_n_lim,
                             align = c("l", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter,
                                  "Replete Average" = formatter("span", 
                                                                style = x ~ style(color = "blue",
                                                                                  font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                       "bold", 
                                                                                                       "")))))

Pigments_tab6 <- formattable(P_gyr_pigs_p_lim,
                             align = c("c", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter))

Pigments_tab7 <- formattable(P_parv_pigs_n_lim,
                             align = c("l", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter,
                                  "Replete Average" = formatter("span", 
                                                                style = x ~ style(color = "blue",
                                                                                  font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                       "bold", 
                                                                                                       "")))))

Pigments_tab8 <- formattable(P_parv_pigs_p_lim,
                             align = c("c", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter))

export_formattable(Pigments_tab1, "E_hux_pigs_n_lim.png")
export_formattable(Pigments_tab2, "E_hux_pigs_p_lim.png")
export_formattable(Pigments_tab3, "I_gal_pigs_n_lim.png")
export_formattable(Pigments_tab4, "I_gal_pigs_p_lim.png")
export_formattable(Pigments_tab5, "P_gyr_pigs_n_lim.png")
export_formattable(Pigments_tab6, "P_gyr_pigs_p_lim.png")
export_formattable(Pigments_tab7, "P_parv_pigs_n_lim.png")
export_formattable(Pigments_tab8, "P_parv_pigs_p_lim.png")

##################
## Xanthophylls Nut Lim Table Export
##################

xanths_n_lim <- export_n_lim %>% 
  filter(Lipid == "19prime_but_fuco" | 
           Lipid == "19prime_hex_fuco" | 
           Lipid == "Alloxanthin" | 
           Lipid == "Crocoxanthin" | 
           Lipid == "Ddc_dd" | 
           Lipid == "Diatoxanthin" | 
           Lipid == "Echinenone" | 
           Lipid == "Fucoxanthin" | 
           Lipid == "Lutein" | 
           Lipid == "Neo_nosxanthin" | 
           Lipid == "Prasinoxanthin" | 
           Lipid == "Violaxanthin" | 
           Lipid == "Zeaxanthin") 

xanths_p_lim <- export_p_lim %>% 
  filter(Lipid == "19prime_but_fuco" | 
           Lipid == "19prime_hex_fuco" | 
           Lipid == "Alloxanthin" | 
           Lipid == "Crocoxanthin" | 
           Lipid == "Ddc_dd" | 
           Lipid == "Diatoxanthin" | 
           Lipid == "Echinenone" | 
           Lipid == "Fucoxanthin" | 
           Lipid == "Lutein" | 
           Lipid == "Neo_nosxanthin" | 
           Lipid == "Prasinoxanthin" | 
           Lipid == "Violaxanthin" | 
           Lipid == "Zeaxanthin") 

xanth_table_export <- merge(xanths_n_lim, xanths_p_lim, by = c("Species", "Lipid"))
write.csv(xanth_table_export, file = "xanth_table_export.csv")


###########
## Xanth separation details, probably not useful anymore
###########
E_hux_xans_n_lim <- xanths_n_lim %>% filter(Species == "E_huxleyi_374")
E_hux_xans_p_lim <- xanths_p_lim %>% filter(Species == "E_huxleyi_374") %>% select(-Species, -Lipid)
I_gal_xans_n_lim <- xanths_n_lim %>% filter(Species == "I_galbana")
I_gal_xans_p_lim <- xanths_p_lim %>% filter(Species == "I_galbana")%>% select(-Species, -Lipid)
P_gyr_xans_n_lim <- xanths_n_lim %>% filter(Species == "P_gyrans")
P_gyr_xans_p_lim <- xanths_p_lim %>% filter(Species == "P_gyrans")%>% select(-Species, -Lipid)
P_parv_xans_n_lim <- xanths_n_lim %>% filter(Species == "P_parvum")
P_parv_xans_p_lim <- xanths_p_lim %>% filter(Species == "P_parvum")%>% select(-Species, -Lipid)




xanths_tab1 <- formattable(E_hux_xans_n_lim,
                             align = c("l", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter,
                                  "Replete Average" = formatter("span", 
                                                                style = x ~ style(color = "blue",
                                                                                  font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                       "bold", 
                                                                                                       "")))))

xanths_tab2 <- formattable(E_hux_xans_p_lim,
                             align = c("c", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter))

xanths_tab3 <- formattable(I_gal_xans_n_lim, 
                             align = c("l", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter,
                                  "Replete Average" = formatter("span", 
                                                                style = x ~ style(color = "blue",
                                                                                  font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                       "bold", 
                                                                                                       "")))))

xanths_tab4 <- formattable(I_gal_xans_p_lim,
                             align = c("c", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter))

xanths_tab5 <- formattable(P_gyr_xans_n_lim,
                             align = c("l", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter,
                                  "Replete Average" = formatter("span", 
                                                                style = x ~ style(color = "blue",
                                                                                  font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                       "bold", 
                                                                                                       "")))))

xanths_tab6 <- formattable(P_gyr_xans_p_lim,
                             align = c("c", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter))

xanths_tab7 <- formattable(P_parv_xans_n_lim,
                             align = c("l", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter,
                                  "Replete Average" = formatter("span", 
                                                                style = x ~ style(color = "blue",
                                                                                  font.weight = ifelse(as.numeric(x) > 0.001, 
                                                                                                       "bold", 
                                                                                                       "")))))

xanths_tab8 <- formattable(P_parv_xans_p_lim,
                             align = c("c", "c", "c", "c", "c", "c"),
                             list("Log 2 Fold Change" = sign_formatter,
                                  "p_value" = significance_formatter))

export_formattable(xanths_tab1, "E_hux_xans_n_lim.png")
export_formattable(xanths_tab2, "E_hux_xans_p_lim.png")
export_formattable(xanths_tab3, "I_gal_xans_n_lim.png")
export_formattable(xanths_tab4, "I_gal_xans_p_lim.png")
export_formattable(xanths_tab5, "P_gyr_xans_n_lim.png")
export_formattable(xanths_tab6, "P_gyr_xans_p_lim.png")
export_formattable(xanths_tab7, "P_parv_xans_n_lim.png")
export_formattable(xanths_tab8, "P_parv_xans_p_lim.png")




##############
# Glycolipids Nut Lim Table Export
##############

glycos_n_lim <- export_n_lim %>% filter(Lipid == "MGDG" | 
                                          Lipid == "DGDG" | 
                                          Lipid == "SQDG" |
                                          Lipid == "LMGDG" |
                                          Lipid == "GADG")
glycos_p_lim <- export_p_lim %>% filter(Lipid == "MGDG" | 
                                          Lipid == "DGDG" | 
                                          Lipid == "SQDG"|
                                          Lipid == "LMGDG" |
                                          Lipid == "GADG")
glycos_table_export <- merge(glycos_n_lim, glycos_p_lim, by = c("Species", "Lipid"))
write.csv(glycos_table_export, file = "glycos_table_export.csv")

###############
## Phosophos Nut Lim Table Export
###############
phosphos_n_lim <- export_n_lim %>% filter(Lipid == "PC" | 
                                          Lipid == "PE" | 
                                          Lipid == "PG")
phosphos_p_lim <- export_p_lim %>% filter(Lipid == "PC" | 
                                          Lipid == "PE" | 
                                          Lipid == "PG")
phosphos_table_export <- merge(phosphos_n_lim, phosphos_p_lim, by = c("Species", "Lipid"))
write.csv(phosphos_table_export, file = "phosphos_table_export.csv")

##################
## Betaine Lipids Nut Lim Table Export
##################
betaines_n_lim <- export_n_lim %>% filter(Lipid == "BLL" | 
                                            Lipid == "DGCC" | 
                                            Lipid == "LDGCC" |
                                            Lipid == "DGTS_DGTA")
betaines_p_lim <- export_p_lim %>% filter(Lipid == "BLL" | 
                                            Lipid == "DGCC" | 
                                            Lipid == "LDGCC" |
                                            Lipid == "DGTS_DGTA")
betaines_table_export <- merge(betaines_n_lim, betaines_p_lim, by = c("Species", "Lipid"))
write.csv(betaines_table_export, file = "betaines_table_export.csv")

#################
# Special Lipids Nut Lim Table Export
#################
odds_n_lim <- export_n_lim %>% filter(Lipid == "dLCB_GSL_No_FA_OH" |
                                            Lipid == "hGSL" |
                                            Lipid == "sGSL" | 
                                            Lipid == "PDPT")
odds_p_lim <- export_p_lim %>% filter(Lipid == "dLCB_GSL_No_FA_OH" |
                                            Lipid == "hGSL" |
                                            Lipid == "sGSL" | 
                                            Lipid == "PDPT")
odds_table_export <- merge(odds_n_lim, odds_p_lim, by = c("Species", "Lipid"))
write.csv(odds_table_export, file = "odds_table_export.csv")

#################
# Quinones Nut Lim Table Export
#################
quinones_n_lim <- export_n_lim %>% filter(Lipid == "UQ10:10" | 
                                            Lipid == "UQ9:9" |
                                            Lipid == "UQ8:8" |
                                            Lipid == "PQ9" |
                                            Lipid == "PQ9OH" | 
                                            Lipid == "PQ9OH2")
quinones_p_lim <- export_p_lim %>% filter(Lipid == "UQ10:10" | 
                                            Lipid == "UQ9:9" |
                                            Lipid == "UQ8:8" |
                                            Lipid == "PQ9" |
                                            Lipid == "PQ9OH" | 
                                            Lipid == "PQ9OH2")
quinones_table_export <- merge(quinones_n_lim, quinones_p_lim, by = c("Species", "Lipid"))
write.csv(quinones_table_export, file = "quinones_table_export.csv")

#################
# Glyceros Nut Lim Table Export
#################
glyceros_n_lim <- export_n_lim %>% filter(Lipid == "FFA" |
                                            Lipid == "DAG" |
                                            Lipid == "TAG")
glyceros_p_lim <- export_p_lim %>% filter(Lipid == "FFA" |
                                            Lipid == "DAG" |
                                            Lipid == "TAG")
glyceros_table_export <- merge(glyceros_n_lim, glyceros_p_lim, by = c("Species", "Lipid"))
write.csv(glyceros_table_export, file = "glyceros_table_export.csv")






##################
# Replete table creation
# Plan: Sum by lipid/headgroup, separate by Strain/Species to get averages and st devs,
# relabel columns then merge by lipid/headgroup to separate for export
##################

Replete_Table <-  DNPPE_Corrected %>% 
  filter(Treatment == "Replete",
         Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         E_hux_Strain != "C_leptoporus") %>% # get rid of single samples (i.e. not run in dup or triplicate)
  group_by(Variable_Treatment, species, E_hux_Strain) %>% 
  summarize(Total_per_species_per_VT = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_per_VT = sum(Total_per_species_per_VT))%>% # peak proportion per sample
  ungroup() %>%
  group_by(species, Variable_Treatment) %>% # group by compound name so that double peaks are aggregated within each sample
  mutate(Proportion = Total_per_species_per_VT/Total_per_VT) %>% 
  ungroup() 

Proportions_by_G_species <- Replete_Table %>% 
  group_by(E_hux_Strain, 
            species) %>%
  summarise(Treatment_Proportions = list(Proportion),
            Treatment_Average = formatC(mean(Proportion), 
                                        format = "e", 
                                        digits = 2),
            Treatment_St_Dev = paste0("±", formatC(sd(Proportion), 
                                                   format = "e", 
                                                   digits = 2))) %>% 
  rename("Lipid" = species,
         "G_species" = E_hux_Strain) %>% 
  ungroup()


Totals_by_G_species <- Replete_Table %>% 
  group_by(E_hux_Strain) %>% 
  summarise(Totals = list(Total_per_VT),
         Total_Average = formatC(mean(Total_per_VT), 
                                 format = "e", 
                                 digits = 2),
         Total_St_Dev = paste0("±", formatC(sd(Total_per_VT), 
                                            format = "e", 
                                            digits = 2)))



Proportions_by_G_species$Treatment_n <- sapply(Proportions_by_G_species$Treatment_Proportions, 
                                    function(x) {length(x)})
#####################
# Separation by type
#####################
C_ericina <- Proportions_by_G_species %>% filter(G_species == "C_ericina") %>% 
  select(-Treatment_Proportions) %>% 
  rename("C_ericina_Average" = Treatment_Average,
         "C_ericina_SD" = Treatment_St_Dev,
         "C_ericina_n" = Treatment_n) %>% 
  select(-G_species)
# 
# C_leptoporus <- Proportions_by_G_species %>% filter(G_species == "C_leptoporus")%>% 
#   select(-Treatment_Proportions)%>% 
#   rename("C_leptoporus_Average" = Treatment_Average,
#          "C_leptoporus_SD" = Treatment_St_Dev,
#          "C_leptoporus_n" = Treatment_n)%>% 
#   select(-G_species)

E_huxleyi_1N <- Proportions_by_G_species %>% filter(G_species == "E_huxleyi_1N")%>% 
  select(-Treatment_Proportions)%>% 
  rename("E_huxleyi_1N_Average" = Treatment_Average,
         "E_huxleyi_1N_SD" = Treatment_St_Dev,
         "E_huxleyi_1N_n" = Treatment_n)%>% 
  select(-G_species)

E_huxleyi_370 <- Proportions_by_G_species %>% filter(G_species == "E_huxleyi_370")%>% 
  select(-Treatment_Proportions)%>% 
  rename("E_huxleyi_370_Average" = Treatment_Average,
         "E_huxleyi_370_SD" = Treatment_St_Dev,
         "E_huxleyi_370_n" = Treatment_n)%>% 
  select(-G_species)

E_huxleyi_374    <- Proportions_by_G_species %>% filter(G_species == "E_huxleyi_374")%>% 
  select(-Treatment_Proportions)%>% 
  rename("E_huxleyi_374_Average" = Treatment_Average,
         "E_huxleyi_374_SD" = Treatment_St_Dev,
         "E_huxleyi_374_n" = Treatment_n)%>% 
  select(-G_species)

E_huxleyi_379    <- Proportions_by_G_species %>% filter(G_species == "E_huxleyi_379")%>% 
  select(-Treatment_Proportions)%>% 
  rename("E_huxleyi_379_Average" = Treatment_Average,
         "E_huxleyi_379_SD" = Treatment_St_Dev,
         "E_huxleyi_379_n" = Treatment_n)%>% 
  select(-G_species)

E_huxleyi_F      <- Proportions_by_G_species %>% filter(G_species == "E_huxleyi_F")%>% 
  select(-Treatment_Proportions)%>% 
  rename("E_huxleyi_3266_Average" = Treatment_Average,
         "E_huxleyi_3266_SD" = Treatment_St_Dev,
         "E_huxleyi_3266_n" = Treatment_n)%>% 
  select(-G_species)

I_galbana       <- Proportions_by_G_species %>% filter(G_species == "I_galbana")%>% 
  select(-Treatment_Proportions)%>% 
  rename("I_galbana_Average" = Treatment_Average,
         "I_galbana_SD" = Treatment_St_Dev,
         "I_galbana_n" = Treatment_n)%>% 
  select(-G_species)

P_antarctica     <- Proportions_by_G_species %>% filter(G_species == "P_antarctica")%>% 
  select(-Treatment_Proportions)%>% 
  rename("P_antarctica_Average" = Treatment_Average,
         "P_antarctica_SD" = Treatment_St_Dev,
         "P_antarctica_n" = Treatment_n)%>% 
  select(-G_species)

P_antarctica_ST  <- Proportions_by_G_species %>% filter(G_species == "P_antarctica_ST")%>% 
  select(-Treatment_Proportions)%>% 
  rename("P_antarctica_ST_Average" = Treatment_Average,
         "P_antarctica_ST_SD" = Treatment_St_Dev,
         "P_antarctica_ST_n" = Treatment_n)%>% 
  select(-G_species)

P_carterae       <- Proportions_by_G_species %>% filter(G_species == "P_carterae")%>% 
  select(-Treatment_Proportions)%>% 
  rename("P_carterae_Average" = Treatment_Average,
         "P_carterae_SD" = Treatment_St_Dev,
         "P_carterae_n" = Treatment_n)%>% 
  select(-G_species)

P_globosa       <- Proportions_by_G_species %>% filter(G_species == "P_globosa")%>% 
  select(-Treatment_Proportions)%>% 
  rename("P_globosa_Average" = Treatment_Average,
         "P_globosa_SD" = Treatment_St_Dev,
         "P_globosa_n" = Treatment_n)%>% 
  select(-G_species)

P_globosa_ST     <- Proportions_by_G_species %>% filter(G_species == "P_globosa_ST")%>% 
  select(-Treatment_Proportions)%>% 
  rename("P_globosa_ST_Average" = Treatment_Average,
         "P_globosa_ST_SD" = Treatment_St_Dev,
         "P_globosa_ST_n" = Treatment_n)%>% 
  select(-G_species)

P_gyrans         <- Proportions_by_G_species %>% filter(G_species == "P_gyrans")%>% 
  select(-Treatment_Proportions)%>% 
  rename("P_gyrans_Average" = Treatment_Average,
         "P_gyrans_SD" = Treatment_St_Dev,
         "P_gyrans_n" = Treatment_n)%>% 
  select(-G_species)

P_parvum        <- Proportions_by_G_species %>% filter(G_species == "P_parvum")%>% 
  select(-Treatment_Proportions)%>% 
  rename("P_parvum_Average" = Treatment_Average,
         "P_parvum_SD" = Treatment_St_Dev,
         "P_parvum_n" = Treatment_n)%>% 
  select(-G_species)

#Clade_Table <- merge(C_ericina, C_leptoporus, all = TRUE, by = "Lipid")
Clade_Table <- merge(C_ericina, E_huxleyi_1N, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, E_huxleyi_370, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, E_huxleyi_374, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, E_huxleyi_379, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, E_huxleyi_F, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, I_galbana, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, P_antarctica, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, P_antarctica_ST, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, P_carterae, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, P_globosa, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, P_globosa_ST, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, P_gyrans, all = TRUE, by = "Lipid")
Clade_Table <- merge(Clade_Table, P_parvum, all = TRUE, by = "Lipid")






Clade_Pigments <- Clade_Table %>% 
  filter( 
    Lipid == "Carotene" | 
      Lipid == "Chl_a" | 
      Lipid == "Chl_c1" | 
      Lipid == "Chl_c2" | 
      Lipid == "Chl_c2 MGDG 14:0_14:0" | 
      Lipid == "Chl_c2 MGDG 18:4_14:0" | 
      Lipid == "Chl_c3" | 
      Lipid == "Chlide_a" | 
      Lipid == "Pheophytin_a" | 
      Lipid == "Pheophytin_a2")

Clade_Xanths <- Clade_Table %>% 
  filter(Lipid == "19prime_but_fuco" | 
           Lipid == "19prime_hex_fuco" | 
           Lipid == "Alloxanthin" | 
           Lipid == "Crocoxanthin" | 
           Lipid == "Ddc_dd" | 
           Lipid == "Diatoxanthin" | 
           Lipid == "Echinenone" | 
           Lipid == "Fucoxanthin" | 
           Lipid == "Lutein" | 
           Lipid == "Neo_nosxanthin" | 
           Lipid == "Prasinoxanthin" | 
           Lipid == "Violaxanthin" | 
           Lipid == "Zeaxanthin")

Clade_Glyceros <- Clade_Table %>% 
  filter(Lipid == "FFA" |
           Lipid == "DAG" |
           Lipid == "TAG")

Clade_Glycos<- Clade_Table %>% 
  filter(Lipid == "MGDG" | 
           Lipid == "DGDG" | 
           Lipid == "SQDG" |
           Lipid == "LMGDG" |
           Lipid == "GADG")

Clade_Phosphos <- Clade_Table %>% 
  filter(Lipid == "PC" | 
           Lipid == "PE" | 
           Lipid == "PG")

Clade_Betaines<- Clade_Table %>% 
  filter(Lipid == "BLL" | 
           Lipid == "DGCC" | 
           Lipid == "LDGCC" |
           Lipid == "DGTS_DGTA")

Clade_Odds<- Clade_Table %>% 
  filter(Lipid == "dLCB_GSL_No_FA_OH" |
           Lipid == "hGSL" |
           Lipid == "sGSL" | 
           Lipid == "PDPT")

Clade_Quinones<- Clade_Table %>% 
  filter(Lipid == "UQ10:10" | 
           Lipid == "UQ9:9" |
           Lipid == "UQ8:8" |
           Lipid == "PQ9" |
           Lipid == "PQ9OH" | 
           Lipid == "PQ9OH2")

Organized_Clade_Table <- rbind(Clade_Pigments, Clade_Glyceros, Clade_Glycos, Clade_Phosphos, Clade_Betaines, Clade_Odds, Clade_Quinones, Clade_Xanths)

write.csv(Organized_Clade_Table, file = "Clade_Table.csv")


####################
# Tukey HSD of Repletes
####################
library(multcompView)
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}







# 
# Pairwise_tests <-  DNPPE_Corrected %>% 
#   filter(species != "DNPPE",
#          Experiment == "Nutrient_Limitation", 
#          Replicate != "Box4_A",
#          Replicate != "Box4_B",
#          Replicate != "A") %>% 
#   group_by(Variable_Treatment) %>% 
#   mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
#          Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) 
# 
# Pairwise_tests$Treatment <- factor(Pairwise_tests$Treatment, levels = c("Replete", "N_Limited", "P_Limited"))
# 
# 
# 
# agg_species <- aggregate(formula = Peak_Proportion ~ species + E_hux_Strain + Treatment + Time + Replicate + Variable_Treatment + Method, 
#                      data = Pairwise_tests,
#                      FUN = sum)

Replete_Table <- Replete_Table %>% 
  mutate(Lipid = str_replace_all(species, ":", "_"),
         Species_Strain = E_hux_Strain)
lipid_list <- unique(Replete_Table$Lipid)


Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "C_ericina")] <- paste("H. ericina")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_carterae")] <- paste("P. carterae")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 1N")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "I_galbana")] <- paste("I. galbana")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_parvum")] <- paste("P. parvum")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_globosa")] <- paste("P. globosa")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_globosa_ST")] <- paste("P. globosa ST")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_antarctica")] <- paste("P. antarctica")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_antarctica_ST")] <- paste("P. antarctica ST")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_globosa")] <- paste("P. globosa")
Replete_Table$Species_Strain[which(Replete_Table$Species_Strain == "P_gyrans")] <- paste("P. gyrans")
# test$E_hux_Strain <- sapply(test$E_hux_Strain, as.factor)

lipid_list <- c("DGCC", "BLL", "PC") #these are just the lipids I'm including in the supp. materials, and I just went through all the above busy work to redo the species names

for(i in 1:length(lipid_list)){
  
  lipid <- lipid_list[i]
  test <- Replete_Table %>% filter(Lipid == paste(lipid))
  test.lm <- lm(Proportion ~ Species_Strain, data = test)
  test.av <- aov(test.lm)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=test.av, conf.level=0.95)
  
  labels<-generate_label_df(TUKEY, "Species_Strain")#generate labels using function
  
  names(labels)<-c('Letters','Species_Strain')#rename columns for merging
  
  yvalue<-aggregate(Proportion~Species_Strain, data=test, mean)# obtain letter position for y axis using means
  
  final<-merge(labels,yvalue) #merge dataframes
  
  print(ggplot(test, aes(x = Species_Strain, y = Proportion)) +
    geom_blank() +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = 'Species/Strain', y = 'Peak Proportion') +
    ggtitle(paste(test$Lipid[1])) +
    theme(plot.title = element_text(hjust = 0.5, face='bold')) +
    #annotate(geom = "rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf, alpha = 0.6, fill = "grey90") +
    geom_boxplot(stat = "boxplot") +
    geom_point(color = "red")+
    geom_text(data = final, aes(x = Species_Strain, y = Proportion, label = Letters),vjust=-3.5,hjust=-.5, size = 5) +
    #geom_vline(aes(xintercept=4.5), linetype="dashed") +
    theme(axis.text.x = element_text(angle = 90, size = 20),
      axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 23),
    axis.title.y = element_text(size = 23),
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                    linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15))+
  ggsave(filename = paste0(lipid_list[i], "_HaptoTukeyHSD.jpg"),
         plot = last_plot(),
         device = "jpg",
         width = 11, height = 8.5))
  
}

#########################
# De-epoxidation Ratio (Dt/[Dt+dd)])
#########################

Diato_cycle <- No_double_peaks %>% 
  filter(species == "Ddc_dd" | species == "Diatoxanthin") %>%
  select(-Corrected_Per_Cell, -Total_Corrected_Per_Cell) %>% 
  spread("species", "Peak_Proportion") %>% 
  mutate(De_epoxidation_State = Diatoxanthin/(Ddc_dd+Diatoxanthin))
  
ggplot(Diato_cycle)+
  geom_point(aes(x = E_hux_Strain, y = De_epoxidation_State))+
  facet_grid(cols = vars(Treatment))
  
#######################
# Aggregating minor species to reduce number of statistical tests
#######################

minor_pre_aggregation <- No_double_peaks %>% filter(lipid_class == "pigment" |
                                     lipid_class == "ubiquinone" |
                                     grepl("plast", lipid_class), 
                                     Time == 0)

minor_pre_aggregation$lipid_class <- ifelse(grepl("PQ", minor_pre_aggregation$species), 
                           as.character("Plastoquinones"),
                           as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("UQ", minor_pre_aggregation$species), 
                                        as.character("Ubiquinones"),
                                        as.character(minor_pre_aggregation$lipid_class))


minor_pre_aggregation$lipid_class <- ifelse(grepl("Chl_a|Pheophytin_a|Chlide_a|DivinylChl_a|Pheophytin_a2", minor_pre_aggregation$species), 
                           as.character("Chl_a_etc"),
                           as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Chl_b|Pheophytin_b2", minor_pre_aggregation$species), 
                           as.character("Chl_b_etc"),
                           as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Chl_c1|Chl_c2|Chl_c3|Chl_c2 MGDG 14:0_14:0|Chl_c2 MGDG 18:4_14:0", minor_pre_aggregation$species), 
                           as.character("Chl_c_etc"),
                           as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("19prime_but_fuco|19prime_hex_fuco|Alloxanthin|Crocoxanthin|Ddc_dd|Diatoxanthin|Echinenone|Lutein|Neo_Nosxanthin|Prasinoxanthin|Violaxanthin|Zeaxanthin", minor_pre_aggregation$species), 
                           as.character("Minor Xanthophylls"),
                           as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Fucoxanthin", minor_pre_aggregation$species), 
                           as.character("Fucoxanthin"),
                           as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Carotene", minor_pre_aggregation$species), 
                           as.character("Carotenes"),
                           as.character(minor_pre_aggregation$lipid_class))

minor_aggregation <- minor_pre_aggregation %>% 
  group_by(Variable_Treatment, lipid_class, Treatment, E_hux_Strain) %>% 
  summarize(Summed_Peak_Proportion = sum(Peak_Proportion)) %>% 
  ungroup()

pigment_plot_output <- minor_aggregation # set up a separate df to customize labelling for plot, look at pigment plot output tab for the rest


# separate repletes from treatment samples, 
minor_repletes <- minor_aggregation %>% 
  filter(Treatment == "Replete") %>% 
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    Variable_Treatment, 
    E_hux_Strain, 
    lipid_class
  ) %>% 
  summarise(Replete_Proportions = list(Summed_Peak_Proportion),
            Replete_Average = mean(Summed_Peak_Proportion),
            Replete_St_Dev = sd(Summed_Peak_Proportion)) %>% 
  ungroup() %>% 
  select(-Variable_Treatment)



minor_diff_tests <- minor_aggregation %>% 
  filter(Treatment != "Replete")%>%
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    Variable_Treatment, 
    E_hux_Strain, 
    lipid_class,
    Treatment) %>%
  summarise(Treatment_Proportions = list(Summed_Peak_Proportion),
            Treatment_Average = mean(Summed_Peak_Proportion),
            Treatment_St_Dev = sd(Summed_Peak_Proportion))


minor_diff_stats <- merge(minor_repletes, minor_diff_tests, by = c("E_hux_Strain", "lipid_class")) 
# "FA_total_no_C", "FA_total_no_DB", "degree_oxidation", "lipid_class" "species", 
minor_diff_stats$species <- minor_diff_stats$lipid_class
minor_diff_stats$lipid_class <- rep("Pigment", length(minor_diff_stats$species))

minor_diff_stats <- minor_diff_stats %>% 
  group_by(Variable_Treatment, lipid_class) %>% #species, compound_name
  mutate(Fold_Change = Treatment_Average / Replete_Average,
         Log2_FC = log2(Fold_Change),
         Relative_Difference = Treatment_Average - Replete_Average,
         p_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$p.value,
         t_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$statistic,
         significant = if_else(p_value <= 0.05, true = paste("Significant"), false = paste("Not Significant"))) %>% 
  ungroup()

minor_only_sigs <- minor_diff_stats %>% 
  filter(significant == "Significant")

ggplot(minor_diff_stats, aes(x = species, y = Log2_FC, fill = species))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(rows = vars(E_hux_Strain),
             cols = vars(Treatment))+
  ggtitle("Treatment Proportion Minus Replete Proportion \n Under Nutrient Stress")+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"))+
  ylab("Mean Difference in Peak Proportion")+
  xlab("Lipid Class")+
  geom_abline(slope = 0, intercept = 0)

ggplot(minor_diff_stats, aes(x = Replete_Average, y = Log2_FC, color = species))+
  geom_point()+
  facet_grid(rows = vars(E_hux_Strain),
             cols = vars(Treatment))+
  ggtitle("Fold Changes Under Nutrient Stress")+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"))+
  ylab("Log 2 Fold Change")+
  xlab("Mean Replete Peak Proportion")+
  geom_abline(slope = 0, intercept = 0)

minor_export_n_lim <- minor_diff_stats %>% 
  filter(Treatment == "N_Limited") %>% 
  select(E_hux_Strain, lipid_class, Replete_Average, Replete_St_Dev, Treatment_Average, Treatment_St_Dev, Treatment, Log2_FC, p_value, -Treatment)%>% 
  mutate(Replete_Average = formatC(Replete_Average, format = "e", digits = 2),
         Replete_St_Dev = paste0("±", formatC(Replete_St_Dev, format = "e", digits = 2)),
         
         Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         Treatment_St_Dev = paste0("±", formatC(Treatment_St_Dev, format = "e", digits = 2)),
         Log2_FC = round(Log2_FC, digits = 3),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("Replete Average" = Replete_Average,
         "R_SD" = Replete_St_Dev,
         "N-Limited Average" = Treatment_Average,
         "N-lim SD" = Treatment_St_Dev,
         "N-Lim Log 2 Fold Change" = Log2_FC,
         "Lipid" = lipid_class,
         "Species" = E_hux_Strain, 
         "N-lim p value" = p_value)

minor_export_p_lim <- minor_diff_stats %>% 
  filter(Treatment == "P_Limited") %>% 
  select(E_hux_Strain, lipid_class, Treatment_Average, Treatment_St_Dev, Treatment, Log2_FC, p_value, -Treatment)%>% 
  mutate(Treatment_Average = formatC(Treatment_Average, format = "e", digits = 2),
         Treatment_St_Dev = paste0("±", formatC(Treatment_St_Dev, format = "e", digits = 2)),
         Log2_FC = round(Log2_FC, digits = 3),
         p_value = ifelse(p_value < 0.005, 
                          formatC(p_value, format = "e", digits = 2), 
                          round(p_value, digits = 3)))%>% 
  rename("P-Limited Average" = Treatment_Average,
         "P-lim_SD" = Treatment_St_Dev,
         "P-lim Log 2 Fold Change" = Log2_FC,
         "Lipid" = lipid_class,
         "Species" = E_hux_Strain, 
         "P-lim p value" = p_value)

minor_table_export <- merge(minor_export_n_lim, minor_export_p_lim, by = c("Species", "Lipid"))
write.csv(minor_table_export, file = "minors_table_export.csv")


#######################
## kableExtra Formatting fold change table for export
#######################
library(kableExtra)

fold_change_table <- read.csv("C:/Users/TSQ/Desktop/Daniel Lowenstein/Mayers_Cultures_2017/Final_Plots_Tables/Fold_Change_Table_FDR_Controlled.csv")

knitr::kable(fold_change_table)

######################
## Total peak area per cell across treatments - to see whether total lipids per cell changed under stress
######################

test <- DNPPE_Corrected %>% 
  filter(Time == 0) %>% 
  group_by(Variable_Treatment, E_hux_Strain, Treatment) %>% 
  mutate(Total_Corrected = sum(Corrected_Peak_Area))

ggplot(test, aes(x = Treatment, y = Total_Corrected, color = E_hux_Strain))+
  geom_point(position = position_dodge(width = 0.5))
  


######################
# SQDG/PG and DGCC/PC Ratio
######################

SQDG_PG_Ratio <- No_double_peaks %>% 
  group_by(Variable_Treatment) %>% 
  mutate(LOD = min(Corrected_Per_Cell[Corrected_Per_Cell > 0])/Total_Corrected_Per_Cell) %>% 
  ungroup()  %>% 
  filter(species == "PG" | species == "SQDG" | species == "DGCC" | species == "PC") %>%    
  group_by(Variable_Treatment, 
    Experiment,
    Treatment, 
    E_hux_Strain,
    species) %>% 
  summarize(Agg_Lipid = ifelse(sum(Peak_Proportion) > 0, sum(Peak_Proportion), LOD)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = species, values_from = Agg_Lipid) %>% 
  mutate(Log10_SQDG_PG_Ratio = log10(SQDG/PG), 
         Log10_DGCC_PC_Ratio = log10(DGCC/PC))

SQDG_PG_Ratio$Treatment <- factor(SQDG_PG_Ratio$Treatment, levels = c("Replete", "N_Limited", "P_Limited"))

ggplot(SQDG_PG_Ratio)+
  geom_jitter(width = 0.05, aes(x = Treatment, y = Log10_DGCC_PC_Ratio, color = Treatment))+
  facet_wrap(~ E_hux_Strain)+
  ylab("Log10(DGCC:PC Ratio)")+
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 10))

ggplot(SQDG_PG_Ratio)+
  geom_jitter(width = 0.05, aes(x = Treatment, y = Log10_SQDG_PG_Ratio, color = Treatment))+
  facet_wrap(~ E_hux_Strain)+
  ylab("Log10(SQDG:PG Ratio)")+
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 10))

Treatment_Ratios <- SQDG_PG_Ratio %>% 
  filter(Treatment != "Replete") %>% 
  group_by(E_hux_Strain, Treatment) %>% 
  summarize(Treatment_Ratio_Log10_SQDG_PG = list(Log10_SQDG_PG_Ratio),
            Average_Treatment_Ratio_Log10_SQDG_PG = mean(Log10_SQDG_PG_Ratio),
            Treatment_Ratio_Log10_DGCC_PC = list(Log10_DGCC_PC_Ratio),
            Average_Treatment_Ratio_Log10_DGCC_PC = mean(Log10_DGCC_PC_Ratio))

Replete_Ratios <- SQDG_PG_Ratio %>% 
  filter(Treatment == "Replete") %>% 
  select(-Treatment) %>% 
  group_by(E_hux_Strain) %>% 
  summarize(Replete_Ratio_Log10_SQDG_PG = list(Log10_SQDG_PG_Ratio),
            Average_Replete_Ratio_Log10_SQDG_PG = mean(Log10_SQDG_PG_Ratio),
            Replete_Ratio_Log10_DGCC_PC = list(Log10_DGCC_PC_Ratio),
            Average_Replete_Ratio_Log10_DGCC_PC = mean(Log10_DGCC_PC_Ratio))

Merged_Ratios <- merge(Treatment_Ratios, Replete_Ratios, by = c("E_hux_Strain"))
Ratio_Tests <- Merged_Ratios %>% 
  group_by(E_hux_Strain, Treatment) %>% 
  mutate(SQDG_PG_t_test = t.test(unlist(Replete_Ratio_Log10_SQDG_PG), 
                                 unlist(Treatment_Ratio_Log10_SQDG_PG))$p.value,
         DGCC_PC_t_test = t.test(unlist(Replete_Ratio_Log10_DGCC_PC), 
                                 unlist(Treatment_Ratio_Log10_DGCC_PC))$p.value) %>% 
  select(-Treatment_Ratio_Log10_SQDG_PG, 
         -Treatment_Ratio_Log10_DGCC_PC,
         -Replete_Ratio_Log10_SQDG_PG,
         -Replete_Ratio_Log10_DGCC_PC)

#write.csv(Ratio_Tests, "Nutlim_SQDG_PG_DGCC_PC_Ratios.csv")

Ratio_Test_Stats <- Ratio_Tests %>% 
  select(E_hux_Strain, Treatment, SQDG_PG_t_test, DGCC_PC_t_test) %>% 
  mutate(FDR_SQDG_PG_p_value = p.adjust(SQDG_PG_t_test, method = "bonferroni", n = 8),
         FDR_DGCC_PC_p_value = p.adjust(DGCC_PC_t_test, method = "bonferroni", n = 8))

#######################
## Calculating Average DB Number
#######################













#######################
# DMSP from Erin McParland - but don't we want to find out about DMS, not P?
#######################

DMSP <- read.csv("DMSP Samples Final Kyle_ErinMcParland_2016.csv")
DMSP_csv <- DMSP %>% filter(Time == 0)
DMSP_csv$Treatment <- factor(DMSP_csv$Treatment, levels = c("Replete", "N_Limited", "P_Limited"))

ggplot(DMSP_csv, aes(x = Time, y = Part.DMSP..moles.cell., color = E_hux_Strain))+
  geom_point(position = position_dodge(width = 0.5))+
  facet_grid(cols = vars(Treatment))

#######################
## Mean Double Bonds by Treatment/Species
#######################  
  
Mean_Double_Bonds <-  DNPPE_Corrected %>% 
  filter(species != "DNPPE",
         Experiment == "Nutrient_Limitation",
         Time == 0,
         #Treatment == "Replete",
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         # E_hux_Strain == "I_galbana" | 
         #   E_hux_Strain == "E_huxleyi_374" | 
         #   E_hux_Strain == "P_parvum" | 
         #   E_hux_Strain == "P_gyrans",
         lipid_class == "IP_DAG" | lipid_class == "FFA" | lipid_class == "TAG") %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  #filter(
  # compound_name == "dLCB_GSL_No_FA_OH 37:4 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 38:2 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 38:3 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 38:4 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 38:5 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 39:4 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 39:4 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 39:5 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 39:6 +1O"|
  #        compound_name == "dLCB_GSL_No_FA_OH 41:5 +1O"|
#        compound_name == "dLCB_GSL_No_FA_OH 41:4 +1O",species == "PDPT"
#        )  %>% 
ungroup() %>% 
  group_by(Variable_Treatment, E_hux_Strain, FA_total_no_DB, species, Treatment) %>% 
  summarize(Proportion_by_species = sum(Peak_Proportion)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment, E_hux_Strain, species, Treatment) %>% 
  mutate(Mean_DB_VT = weighted.mean(FA_total_no_DB, Proportion_by_species))


ggplot(Mean_Double_Bonds, aes(x = E_hux_Strain, y = Mean_DB_VT, color = E_hux_Strain))+
  facet_grid(cols = vars(species), rows = vars(Treatment))+
  geom_point(position = position_dodge(width = 0.3))+
  # geom_text()+
  ggtitle("Total DB")+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"))


%>% 
  ungroup() %>% 
  group_by(E_hux_Strain, species, FA_total_no_DB) %>% 
  summarize(Mean_DB = mean(Mean_DB_VT), 
            SD_DB = sd(Mean_DB_VT))
#%>% mutate(DAG = paste(str_extract(string = compound_name, pattern = "(?<= )\\d+:\\d+")))



Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N_Limited", "P_Limited"))

Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("P_carterae", 
                                      "E_huxleyi_1N", 
                                      "E_huxleyi_F", 
                                      "E_huxleyi_370", 
                                      "E_huxleyi_374", 
                                      "E_huxleyi_379", 
                                      "I_galbana", 
                                      "P_parvum", 
                                      "C_ericina", 
                                      "P_globosa", 
                                      "P_globosa_ST", 
                                      "P_antarctica", 
                                      "P_antarctica_ST", 
                                      "P_gyrans"))

agg_GSL <- aggregate(formula = Peak_Proportion ~ E_hux_Strain + species + Treatment + Time + Replicate + Variable_Treatment + Method + FA_total_no_DB, 
                     data = Cultures_GSL,
                     FUN = sum)


ggplot(Cultures_GSL, aes(x = Treatment, y = Mean_DB_VT, color = lipid_class))+
  facet_grid(cols = vars(E_hux_Strain))+
  geom_point(position = position_dodge(width = 0.3))+
  ggtitle("Total DB")+
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"))

+
  geom_text(aes(label = QE_Number))
  

######################
## PDPT figure plots
######################

#replete
Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         # E_hux_Strain == "E_huxleyi_374" |
         #   E_hux_Strain == "I_galbana" |
         #   E_hux_Strain == "P_parvum" |
         #   E_hux_Strain == "P_gyrans" 
         Treatment == "Replete"
  ) %>% 
  group_by(
    Variable_Treatment, 
    species,

    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  filter(species == "PDPT")

Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N_Limited", "P_Limited"))

Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_ericina")] <- paste("H. ericina")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_carterae")] <- paste("P. carterae")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 - 1N")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "I_galbana")] <- paste("I. galbana")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_parvum")] <- paste("P. parvum")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa_ST")] <- paste("P. globosa - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica")] <- paste("P. antarctica")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica_ST")] <- paste("P. antarctica - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_gyrans")] <- paste("P. gyrans")


Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c(
                                      "C. leptoporus",
                                      "E. huxleyi 3268 - 1N", 
                                      "E. huxleyi 3266", 
                                      "E. huxleyi 370", 
                                      "E. huxleyi 374", 
                                      "E. huxleyi 379", 
                                      "H. ericina",
                                      "P. parvum", 
                                      "I. galbana", 
                                      "P. carterae",
                                      "P. globosa", 
                                      "P. globosa - ST", 
                                      "P. antarctica", 
                                      "P. antarctica - ST", 
                                      "P. gyrans"))



a <- ggplot(Cultures_GSL, aes(x = E_hux_Strain, y = Peak_Proportion))+
  #facet_grid(cols = vars(Treatment), scales = "free")+
  geom_boxplot(position = position_dodge(width = 0.5))+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, size = 15, family = "Times"),
        axis.text.y = element_text(size = 15, family = "Times"),
        axis.title.x = element_text(size = 15, family = "Times"),
        axis.title.y = element_text(size = 15, family = "Times"),
        plot.title = element_text(hjust = 0.5, family = "Times"),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 10, family = "Times"))+
  xlab("Species")+
  ylab("Peak Proportion")

# +
#   #ylim(0,0.4)+
#   labs(color = "Species")

# Nut lim experiment

Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         E_hux_Strain == "E_huxleyi_374" |
           E_hux_Strain == "I_galbana" |
           E_hux_Strain == "P_parvum" |
           E_hux_Strain == "P_gyrans") %>% 
  group_by(
    Variable_Treatment, 
    species,
    
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  filter(species == "PDPT")

Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "N_Limited")] <- paste("N-Limited")
Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "P_Limited")] <- paste("P-Limited")

Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N-Limited", "P-Limited"))

Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_ericina")] <- paste("H. ericina")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_carterae")] <- paste("P. carterae")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 - 1N")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "I_galbana")] <- paste("I. galbana")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_parvum")] <- paste("P. parvum")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa_ST")] <- paste("P. globosa - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica")] <- paste("P. antarctica")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica_ST")] <- paste("P. antarctica - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_gyrans")] <- paste("P. gyrans")


Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("P. carterae",
                                      "C. leptoporus",
                                      "E. huxleyi 3268 - 1N", 
                                      "E. huxleyi 3266", 
                                      "E. huxleyi 370", 
                                      "E. huxleyi 374", 
                                      "E. huxleyi 379", 
                                      "I. galbana", 
                                      "P. parvum", 
                                      "C. ericina", 
                                      "P. globosa", 
                                      "P. globosa - ST", 
                                      "P. antarctica", 
                                      "P. antarctica - ST", 
                                      "P. gyrans"))

b <- ggplot(Cultures_GSL, aes(x = E_hux_Strain, y = Peak_Proportion, color = Treatment))+
  #facet_grid(cols = vars(Treatment), scales = "free")+
  geom_boxplot(position = position_dodge(width = 0.85))+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, size = 15, family = "Times"),
        axis.text.y = element_text(size = 15, family = "Times"),
        axis.title.x = element_text(size = 15, family = "Times"),
        axis.title.y = element_text(size = 15, family = "Times"),
        plot.title = element_text(hjust = 0.5, family = "Times"),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 15, family = "Times"),
        legend.title = element_text(size = 15, family = "Times")
        )+
  xlab("Species")+
  ylab("Peak Proportion")+
  #ylim(0,0.4)+
  labs(color = "Treatment")+ 
  scale_color_brewer(palette="Dark2")

library(cowplot)
cowplot::plot_grid(a, b, align = "h")

library(gridExtra)
grid.arrange(a, b, ncol = 2, heights = c(1, 0.8))


###########################
## Growth rates
###########################

growth_rates <- read.csv("Haptophyte_Growth_Rates.csv")
growth_rates$Cell_Count <- as.numeric(as.character(growth_rates$Cell_Count))
growth_rates$Treatment <- factor(growth_rates$Treatment, levels = c("Replete", "N-limited", "P-limited"))

replete_growth_rates <- growth_rates %>% 
  filter(Treatment == "Replete", Day != 0)

ggplot(replete_growth_rates, aes(Day, Cell_Count, color = Replicate, shape = Replicate))+
  geom_path(position = position_dodge(width = 0.5))+
  geom_point(position = position_dodge(width = 0.5))+
  facet_wrap(~Species,
             scales = "free")+
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid", 
                                        color = "gray"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))


nut_lim_growth_rates <- growth_rates %>% 
  filter(Species == "E. huxleyi 374" |
           Species == "P. parvum" |
           Species == "P. gyrans" |
           Species == "I. galbana",
         Treatment != "Replete")

ggplot(nut_lim_growth_rates, aes(Day, Cell_Count, color = Replicate, shape = Replicate))+
  geom_path(position = position_dodge(width = 0.5))+
  geom_point(position = position_dodge(width = 0.5))+
  facet_grid(cols = vars(Treatment), 
             rows = vars(Species),
             scales = "free")+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid", 
                                        color = "gray"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))

#########################
# pigment plot output
#########################


pigment_plot_output <- minor_aggregation
pigment_plot_output$Treatment <- sapply(pigment_plot_output$Treatment, function(x){(str_replace(string = x, pattern = "_", replacement = "-"))})

pigment_plot_output$Treatment <- factor(pigment_plot_output$Treatment, levels = c("Replete", "N-Limited", "P-Limited"))

species_labels <- c("E. huxleyi 374", "I. galbana", "P. gyrans", "P. parvum")
names(species_labels) <- c("E_huxleyi_374", "I_galbana", "P_gyrans", "P_parvum")

pigment_plot_output$lipid_class[pigment_plot_output$lipid_class == "Chl_a_etc"] <- "Chl a pigments" 
pigment_plot_output$lipid_class[pigment_plot_output$lipid_class == "Chl_c_etc"] <- "Chl c pigments" 
pigment_plot_output$lipid_class[pigment_plot_output$lipid_class == "Carotenes"] <- paste0(expression(beta), "-Carotene")


ggplot(pigment_plot_output, aes(x = lipid_class, 
                                y = Summed_Peak_Proportion, 
                                color = Treatment))+
  geom_boxplot(position = position_dodge(width = 0.75))+
  facet_grid(~E_hux_Strain, labeller = labeller(E_hux_Strain = species_labels))+
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid", 
                                        color = "gray"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15))+
  xlab("Species")+
  ylab("Peak Proportion")+
  #ylim(0,0.4)+
  labs(color = "Treatment")+ 
  scale_color_brewer(palette="Dark2")





#########################
## redoing diff stats by species, but aggregating names from the beginning
## so that sums and fold changes will have legitimate st devs and significance tests
#########################


test <- DNPPE_Corrected

test$species <- ifelse(grepl("PQ", test$species), 
                       as.character("Plastoquinones"),
                       as.character(test$species))

test$species <- ifelse(grepl("UQ", test$species), 
                       as.character("Ubiquinones"),
                       as.character(test$species))


test$species <- ifelse(grepl("Chl_a|Pheophytin_a|Chlide_a|DivinylChl_a|Pheophytin_a2", test$species), 
                       as.character("Chl_a_etc"),
                       as.character(test$species))

test$species <- ifelse(grepl("Chl_b|Pheophytin_b2", test$species), 
                       as.character("Chl_b_etc"),
                       as.character(test$species))

test$species <- ifelse(grepl("Chl_c1|Chl_c2|Chl_c3|Chl_c2 MGDG 14:0_14:0|Chl_c2 MGDG 18:4_14:0", test$species), 
                       as.character("Chl_c_etc"),
                       as.character(test$species))

test$species <- ifelse(grepl("19prime_but_fuco|Alloxanthin|Crocoxanthin|Ddc_dd|Diatoxanthin|Echinenone|Lutein|Neo_Nosxanthin|Prasinoxanthin|Violaxanthin|Zeaxanthin", test$species), 
                       as.character("Minor Xanthophylls"),
                       as.character(test$species))

test$species <- ifelse(grepl("Fucoxanthin", test$species), 
                       as.character("Fucoxanthin"),
                       as.character(test$species))

test$species <- ifelse(grepl("Carotene", test$species), 
                       as.character("Carotenes"),
                       as.character(test$species))

test$species <- ifelse(grepl("DGCC", test$species), 
                       as.character("DGCC/LDGCC"),
                       as.character(test$species))

test$species <- ifelse(grepl("MGDG", test$species), 
                       as.character("MGDG/LMGDG"),
                       as.character(test$species))




No_double_peaks <- test %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         G_species == "I_galbana" |
           G_species == "P_parvum" |
           G_species == "P_gyrans" |
           E_hux_Strain == "E_huxleyi_374",
         Time == 0) %>% 
  group_by(
    #compound_name,
    Variable_Treatment, 
    species,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain) %>% 
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup()

# separate repletes from treatment samples, 
Repletes <- No_double_peaks %>% # changed to nodoublepeaksMAJORS for some analyses
  filter(Treatment == "Replete") %>% 
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    #compound_name,
    Variable_Treatment, 
    E_hux_Strain, 
    species,
    #,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation
  ) %>% 
  summarise(Replete_Proportions = list(Peak_Proportion),
            Replete_Average = mean(Peak_Proportion),
            Replete_St_Dev = sd(Peak_Proportion)) %>% 
  ungroup() %>% 
  select(-Variable_Treatment)



diff_tests <- No_double_peaks %>% # changed to nodoublepeaksMAJORS for some analyses
  filter(Time == 0, Treatment != "Replete")%>%
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    #compound_name,
    Variable_Treatment, 
    E_hux_Strain, 
    species,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation,
    # lipid_class,
    Treatment) %>%
  summarise(Treatment_Proportions = list(Peak_Proportion),
            Treatment_Average = mean(Peak_Proportion),
            Treatment_St_Dev = sd(Peak_Proportion))


diff_stats <- merge(Repletes, diff_tests, by = c("E_hux_Strain", "species")) 
# "FA_total_no_C", "FA_total_no_DB", "degree_oxidation", "lipid_class" "species", , "compound_name"


diff_stats <- diff_stats %>% 
  group_by(Variable_Treatment, species) %>% #species, compound_name, lipid_class
  mutate(Fold_Change = Treatment_Average / Replete_Average,
         Log2_FC = log2(Fold_Change),
         Relative_Difference = Treatment_Average - Replete_Average,
         p_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$p.value,
         t_value = t.test(unlist(Replete_Proportions), 
                          unlist(Treatment_Proportions))$statistic,
         significant = if_else(p_value <= 0.05, true = paste("Significant"), false = paste("Not Significant"))) %>% 
  ungroup() 

####

nut_response_all <- diff_stats  %>%
  group_by(species) %>%
  mutate(Max_replete_abundance = max(Replete_Average), Max_treatment_abundance = max(Treatment_Average)) %>%
  filter(Max_replete_abundance >= 0.003 | Max_treatment_abundance >= 0.003) %>%
  select(species, Variable_Treatment, Log2_FC, significant)

nut_response_all$Log2_FC[nut_response_all$significant == "Not Significant"] <- 0
nut_response_all$Log2_FC[nut_response_all$Log2_FC == "NaN"] <- 0

## getting rid of NaNs, which indicates a 0 in both Replete and Treatment
## and -Inf and Inf, which indicates presence/absence -- setting them to -10 and 10, which is an order of magnitude above and below the next largest fold changes
## 
nut_response_all$Log2_FC[nut_response_all$Log2_FC == "Inf"] <- 1.1111111111 # set inf and -inf values to random strings so I can find max and min, then reset them to values above and below the scale
nut_response_all$Log2_FC[nut_response_all$Log2_FC == "-Inf"] <- -1.1111111111 

print(summary(nut_response_all$Log2_FC)) # min = ???, max = ???

nut_response_all$Log2_FC[nut_response_all$Log2_FC == 1.1111111111] <- 7
nut_response_all$Log2_FC[nut_response_all$Log2_FC == -1.1111111111 ] <- -5



nut_response_all <- nut_response_all %>% 
  select(-significant) %>% 
  pivot_wider(names_from = 'Variable_Treatment', values_from = 'Log2_FC')

nut_response_all <- as.data.frame(nut_response_all)
rownames(nut_response_all) <- nut_response_all[,1] 
nut_response_all <- nut_response_all[ ,-1]

nut_response_all_changes <- nut_response_all[-which(rowSums(nut_response_all) == 0),]
nut_response_all_changes[is.na(nut_response_all_changes)] <- 0

library(pheatmap)
library(RColorBrewer)
pheatmap(nut_response_all, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 25, 
         main = "Haptophyte N- and P-Limited Significant Fold Changes \n Log Scale - 'ward D2' clustering - 0.3% peak area abundance cutoff", 
         treeheight_row = 200, treeheight_col = 70,
         legend = TRUE,
         show_rownames = TRUE)

# Testing bootstrap values
# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(nut_response_all, method.hclust="ward.D2",
               method.dist="euclidean")
fit10k <- pvclust(nut_response_all, method.hclust="ward.D2",
                  method.dist="euclidean", nboot = 10000, parallel = TRUE)
#compound_cluster <- fit # ran it first with a flipped cluster, so it calc'ed bootstrap values for the compound cluster--stored in separate variable, compound_cluster--flipping it back to do species/strain
seplot(fit, identify = TRUE)
test <- plot(fit10k)   # dendogram with p values
# add rectangles around groups highly supported by the data
test%>% pvrect(fit, alpha=0.95)

#########################
# Converting Nutrient Response Diff Stats to 
# Cluster Analyses to look at clustering of fold changes
#########################

# setting color gradient for use on both heatmaps
library(RColorBrewer)
colorvector <- seq(-8, 8, by = 0.1)

# make sure current diff stats is by compound (though maybe also check by lipid species)

nut_response_Nlim <- diff_stats %>%
  filter(Treatment == "N_Limited") %>% 
  group_by(compound_name) %>%
  mutate(Max_replete_abundance = max(Replete_Average), Max_treatment_abundance = max(Treatment_Average)) %>%
  filter(Max_replete_abundance >= 0.005 | Max_treatment_abundance >= 0.005) %>%
  select(compound_name, Variable_Treatment, Log2_FC, significant)

nut_response_Nlim$Log2_FC[nut_response_Nlim$significant == "Not Significant"] <- 0
nut_response_Nlim$Log2_FC[nut_response_Nlim$Log2_FC == "NaN"] <- 0

## getting rid of NaNs, which indicates a 0 in both Replete and Treatment
## and -Inf and Inf, which indicates presence/absence -- setting them to -10 and 10, which is an order of magnitude above and below the next largest fold changes
## 
nut_response_Nlim$Log2_FC[nut_response_Nlim$Log2_FC == "Inf"] <- 1.1111111111 # set inf and -inf values to random strings so I can find max and min, then reset them to values above and below the scale
nut_response_Nlim$Log2_FC[nut_response_Nlim$Log2_FC == "-Inf"] <- -1.1111111111 

print(summary(nut_response_Nlim$Log2_FC)) # min = ???, max = ???

nut_response_Nlim$Log2_FC[nut_response_Nlim$Log2_FC == 1.1111111111] <- 8
nut_response_Nlim$Log2_FC[nut_response_Nlim$Log2_FC == -1.1111111111 ] <- -8



nut_response_Nlim <- nut_response_Nlim %>% 
  select(-significant) %>% 
  pivot_wider(names_from = 'Variable_Treatment', values_from = 'Log2_FC')

nut_response_Nlim <- as.data.frame(nut_response_Nlim)
rownames(nut_response_Nlim) <- nut_response_Nlim[,1] 
nut_response_Nlim <- nut_response_Nlim[ ,-1]

nut_response_Nlim_changes <- nut_response_Nlim[-which(rowSums(nut_response_Nlim) == 0),]
nut_response_Nlim_changes[is.na(nut_response_Nlim_changes)] <- 0

a <- pheatmap(nut_response_Nlim, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 30, 
         #main = "Haptophyte N-Limited Fold Change \n Log Scale - 'ward D2' clustering - 0.5% Peak Abundance Cutoff", 
         treeheight_row = 200, treeheight_col = 70,
         legend = FALSE,
         show_rownames = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(colorvector)), 
         breaks = colorvector) 

###
nut_response_Plim <- diff_stats %>%
  filter(Treatment == "P_Limited") %>% 
  group_by(compound_name) %>%
  mutate(Max_replete_abundance = max(Replete_Average), Max_treatment_abundance = max(Treatment_Average)) %>%
  filter(Max_replete_abundance >= 0.005 | Max_treatment_abundance >= 0.005) %>%
  select(compound_name, Variable_Treatment, Log2_FC, significant)

nut_response_Plim$Log2_FC[nut_response_Plim$significant == "Not Significant"] <- 0
nut_response_Plim$Log2_FC[nut_response_Plim$Log2_FC == "NaN"] <- 0

## getting rid of NaNs, which indicates a 0 in both Replete and Treatment
## and -Inf and Inf, which indicates presence/absence -- setting them to -10 and 10, which is an order of magnitude above and below the next largest fold changes
## 
nut_response_Plim$Log2_FC[nut_response_Plim$Log2_FC == "Inf"] <- 1.1111111111 # set inf and -inf values to random strings so I can find max and min, then reset them to values above and below the scale
nut_response_Plim$Log2_FC[nut_response_Plim$Log2_FC == "-Inf"] <- -1.1111111111 

print(summary(nut_response_Plim$Log2_FC)) # min = ???, max = ???

nut_response_Plim$Log2_FC[nut_response_Plim$Log2_FC == 1.1111111111] <- 6
nut_response_Plim$Log2_FC[nut_response_Plim$Log2_FC == -1.1111111111 ] <- -6



nut_response_Plim <- nut_response_Plim %>% 
  select(-significant) %>% 
  pivot_wider(names_from = 'Variable_Treatment', values_from = 'Log2_FC')

nut_response_Plim <- as.data.frame(nut_response_Plim)
rownames(nut_response_Plim) <- nut_response_Plim[,1] 
nut_response_Plim <- nut_response_Plim[ ,-1]

nut_response_Plim_changes <- nut_response_Plim[-which(rowSums(nut_response_Plim) == 0),]
nut_response_Plim_changes[is.na(nut_response_Plim_changes)] <- 0

b <- pheatmap(nut_response_Plim, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 30, 
         #main = "Haptophyte P-Limited Fold Change \n Log Scale - 'ward D2' clustering - 0.5% Peak Abundance Cutoff", 
         treeheight_row = 200, treeheight_col = 70,
         legend = TRUE,
         show_rownames = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(colorvector)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = colorvector) # Sets the breaks of the color scale as in breaksList

library(gridExtra)

grid.arrange(grobs = list(a[[4]], b[[4]]), nrow = 1)

###############################
## N_limited
###############################

N_Lim_Cultures_Cluster <-  DNPPE_Corrected %>% 
  filter(Treatment == "N_Limited",
         Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         G_species != "C_leptoporus",
         Time == 0) %>% # get rid of single samples (i.e. not run in dup or triplicate)
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% # peak proportion per sample
  ungroup() %>%
  group_by(compound_name, Variable_Treatment) %>% # group by compound name so that double peaks are aggregated within each sample
  mutate(Total_by_species = sum(Corrected_Per_Cell), 
         Proportion_by_species = Total_by_species/Total_Corrected_Per_Cell) %>% # combine double peaks in each sample
  ungroup() %>% 
  select(-LOBdbase_mz,
         -lipid_class,
         -species,
         -FA_total_no_C,
         -FA_total_no_DB,
         -degree_oxidation,
         -peak_area,
         -G_species,
         -Treatment,
         -Experiment,
         -Time,
         -Replicate,
         #-Sample_ID,
         -Cells_per_sample,
         -Total_Corrected_Per_Cell,
         -QE_Number,
         -Eluent_Sequence,
         -Max_Blank_Peak,
         -Average_Blank_Peak,
         -Blank_Corrected,
         -Peak_Proportion,
         -Light_Dark, 
         -E_hux_Strain,
         -Total_by_species,
         -Peak_Per_Cell) %>%# pause here and make sure there are matching sample ID's and
  # variable_treatments--there are (2 Samp_ID's per variable treatment- Hummel and Normal
  # now summarize and spread
  group_by(compound_name, Variable_Treatment) %>% 
  summarise(Proportion = mean(Proportion_by_species)) %>% 
  spread("Variable_Treatment", "Proportion")

write.csv(N_Lim_Cultures_Cluster, "N_Lim_Cultures_Cluster.csv")
N_Lim_Cultures_Cluster <- read.csv("N_Lim_Cultures_Cluster.csv")
rownames(N_Lim_Cultures_Cluster) <- N_Lim_Cultures_Cluster[,1] 
N_Lim_Cultures_Cluster <- N_Lim_Cultures_Cluster[ ,-1]

N_Lim_Cultures_Cluster <- N_Lim_Cultures_Cluster[which(rowSums(N_Lim_Cultures_Cluster) > 0),]

log_NLim <- log(N_Lim_Cultures_Cluster)

write.csv(log_NLim, "log_NLim.csv") #write and remove -Inf values....
log_NLim <- read.csv("log_NLim.csv")
log_NLim <- as.data.frame(log_NLim)

rownames(log_NLim) <- log_NLim[,1]
log_NLim <- log_NLim[ ,-1]

pheatmap(log_NLim, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 18, 
         main = "Haptophyte N-Limited Lipid Cluster \n Log Scale - 'ward D2' clustering", 
         treeheight_row = 333, treeheight_col = 100,
         legend = TRUE,
         show_rownames = FALSE)


###############################
## P_limited
###############################

P_Lim_Cultures_Cluster <-  DNPPE_Corrected %>% 
  filter(Treatment == "P_Limited",
         Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         G_species != "C_leptoporus",
         Time == 0) %>% # get rid of single samples (i.e. not run in dup or triplicate)
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% # peak proportion per sample
  ungroup() %>%
  group_by(compound_name, Variable_Treatment) %>% # group by compound name so that double peaks are aggregated within each sample
  mutate(Total_by_species = sum(Corrected_Per_Cell), 
         Proportion_by_species = Total_by_species/Total_Corrected_Per_Cell) %>% # combine double peaks in each sample
  ungroup() %>% 
  select(-LOBdbase_mz,
         -lipid_class,
         -species,
         -FA_total_no_C,
         -FA_total_no_DB,
         -degree_oxidation,
         -peak_area,
         -G_species,
         -Treatment,
         -Experiment,
         -Time,
         -Replicate,
         #-Sample_ID,
         -Cells_per_sample,
         -Total_Corrected_Per_Cell,
         -QE_Number,
         -Eluent_Sequence,
         -Max_Blank_Peak,
         -Average_Blank_Peak,
         -Blank_Corrected,
         -Peak_Proportion,
         -Light_Dark, 
         -E_hux_Strain,
         -Total_by_species,
         -Peak_Per_Cell) %>%# pause here and make sure there are matching sample ID's and
  # variable_treatments--there are (2 Samp_ID's per variable treatment- Hummel and Normal
  # now summarize and spread
  group_by(compound_name, Variable_Treatment) %>% 
  summarise(Proportion = mean(Proportion_by_species)) %>% 
  spread("Variable_Treatment", "Proportion")

#write.csv(P_Lim_Cultures_Cluster, "P_Lim_Cultures_Cluster.csv")
P_Lim_Cultures_Cluster <- read.csv("P_Lim_Cultures_Cluster.csv")
rownames(P_Lim_Cultures_Cluster) <- P_Lim_Cultures_Cluster[,1] 
P_Lim_Cultures_Cluster <- P_Lim_Cultures_Cluster[ ,-1]

P_Lim_Cultures_Cluster <- P_Lim_Cultures_Cluster[which(rowSums(P_Lim_Cultures_Cluster) > 0),]

log_P_Lim <- log(P_Lim_Cultures_Cluster)

#write.csv(log_P_Lim, "log_P_Lim.csv") #write and remove -Inf values....
log_P_Lim <- read.csv("log_P_Lim.csv")
log_P_Lim <- as.data.frame(log_P_Lim)

rownames(log_P_Lim) <- log_P_Lim[,1]
log_P_Lim <- log_P_Lim[ ,-1]

pheatmap(log_P_Lim, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 18, 
         main = "Haptophyte P-Limited Lipid Cluster \n Log Scale - 'ward D2' clustering", 
         treeheight_row = 333, treeheight_col = 100,
         legend = TRUE,
         show_rownames = FALSE)


#####################
## Testing diff stats but indiv treatment samples vs replete average
## to get a log cluster by sample of fold changes
#####################

Treatments <- No_double_peaks %>% 
  filter(Treatment == "N_Limited" | Treatment == "P_Limited")

New_Diff_Stats <- merge(Treatments, Repletes, by = c("E_hux_Strain", "species", "compound_name", "lipid_class"), all = TRUE)
New_Diff_Stats$Replete_Average[is.na(New_Diff_Stats$Replete_Average)] <- 0 # change NAs to 0


New_Diff_Stats <- New_Diff_Stats %>% 
  mutate(Fold_Change = Peak_Proportion / Replete_Average,
                                        Log2_FC = log2(Fold_Change),
                                        Relative_Difference = Peak_Proportion - Replete_Average) %>% 
                                        # ,
                                        # p_value = t.test(unlist(Replete_Proportions), 
                                        #                  unlist(Treatment_Proportions))$p.value,
                                        # t_value = t.test(unlist(Replete_Proportions), 
                                        #                  unlist(Treatment_Proportions))$statistic,
                                        # significant = if_else(p_value <= 0.05, true = paste("Significant"), false = paste("Not Significant"))
  group_by(compound_name) %>% 
  mutate(Average_Species_Abundance = mean(Replete_Average), Max_abundance = max(Replete_Average), Min_abundance = min(Replete_Average)) %>% 
  filter(Max_abundance >= 0.001)


## getting rid of NaNs, which indicates a 0 in both Replete and Treatment
## and -Inf and Inf, which indicates presence/absence -- setting them to -10 and 10, which is an order of magnitude above and below the next largest fold changes
## 
New_Diff_Stats$Log2_FC[New_Diff_Stats$Log2_FC == "NaN"] <- 0
New_Diff_Stats$Log2_FC[New_Diff_Stats$Log2_FC == "Inf"] <- 1.1111111111 # set inf and -inf values to random strings so I can find max and min, then reset them to values above and below the scale
New_Diff_Stats$Log2_FC[New_Diff_Stats$Log2_FC == "-Inf"] <- -1.1111111111 

print(summary(New_Diff_Stats$Log2_FC)) # min = -7.9, max = 6.22

New_Diff_Stats$Log2_FC[New_Diff_Stats$Log2_FC == 1.1111111111] <- 8
New_Diff_Stats$Log2_FC[New_Diff_Stats$Log2_FC == -1.1111111111 ] <- -8



All_New_Diff_Stats <- New_Diff_Stats %>% 
  filter(Treatment == "N_Limited" | Treatment == "P_Limited") %>% 
  select(Variable_Treatment, Log2_FC, compound_name) %>% 
  pivot_wider(names_from = 'Variable_Treatment', values_from = 'Log2_FC')


All_New_Diff_Stats <- as.data.frame(All_New_Diff_Stats)
rownames(All_New_Diff_Stats) <- All_New_Diff_Stats[,1] 
All_New_Diff_Stats <- All_New_Diff_Stats[ ,-1]

All_New_Diff_Stats[is.na(All_New_Diff_Stats)] <- 0 # get rid of NAs (All 0's across a row introduced by spreading the df where there was no compound row for replete or )


# getting rid of compounds that were undetected in every replicate
#All_New_Diff_Stat_changes <- All_New_Diff_Stats[-which(rowSums(All_New_Diff_Stats) == 0),]

pheatmap(All_New_Diff_Stats, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 25, 
         main = "Haptophyte N- and P-Limited Fold Change \n Log Scale - 'ward D2' clustering", 
         treeheight_row = 200, treeheight_col = 70,
         legend = TRUE,
         show_rownames = FALSE)
  
N_Lim_New_Diff_Stats <- New_Diff_Stats %>% 
  filter(Treatment == "N_Limited") %>% 
  select(Variable_Treatment, Log2_FC, species) %>% 
  pivot_wider(names_from = 'Variable_Treatment', values_from = 'Log2_FC')

N_Lim_New_Diff_Stats <- as.data.frame(N_Lim_New_Diff_Stats)
rownames(N_Lim_New_Diff_Stats) <- N_Lim_New_Diff_Stats[,1] 
N_Lim_New_Diff_Stats <- N_Lim_New_Diff_Stats[ ,-1]

N_Lim_New_Diff_Stats[is.na(N_Lim_New_Diff_Stats)] <- 0 # get rid of NAs (All 0's across a row introduced by spreading the df where there was no compound row for replete or )


# getting rid of compounds that were undetected in every replicate
#N_Lim_New_Diff_Stat_changes <- N_Lim_New_Diff_Stats[-which(rowSums(N_Lim_New_Diff_Stats) == 0),]

pheatmap(N_Lim_New_Diff_Stats, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 25, 
         main = "Haptophyte N-Limited Fold Change \n Log Scale - 'ward D2' clustering", 
         treeheight_row = 200, treeheight_col = 70,
         legend = TRUE,
         show_rownames = TRUE)


P_Lim_New_Diff_Stats <- New_Diff_Stats %>% 
  filter(Treatment == "P_Limited") %>% 
  select(Variable_Treatment, Log2_FC, species) %>% 
  pivot_wider(names_from = 'Variable_Treatment', values_from = 'Log2_FC')

P_Lim_New_Diff_Stats <- as.data.frame(P_Lim_New_Diff_Stats)
rownames(P_Lim_New_Diff_Stats) <- P_Lim_New_Diff_Stats[,1] 
P_Lim_New_Diff_Stats <- P_Lim_New_Diff_Stats[ ,-1]

P_Lim_New_Diff_Stats[is.na(P_Lim_New_Diff_Stats)] <- 0 # get rid of NAs (All 0's across a row introduced by spreading the df where there was no compound row for replete or )


# getting rid of compounds that were undetected in every replicate
#P_Lim_New_Diff_Stat_changes <- P_Lim_New_Diff_Stats[-which(rowSums(P_Lim_New_Diff_Stats) == 0),]

pheatmap(P_Lim_New_Diff_Stats, 
         cluster_rows = TRUE, 
         cex=0.6, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2", 
         cex=0.2, 
         fontsize = 25, 
         main = "Haptophyte P-Limited Fold Change \n Log Scale - 'ward D2' clustering", 
         treeheight_row = 200, treeheight_col = 70,
         legend = TRUE,
         show_rownames = TRUE)

##################################
# Unique Lipids per Species and treatment
##################################

No_double_peaks <- DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         G_species == "I_galbana" |
           G_species == "P_parvum" |
           G_species == "P_gyrans" |
           E_hux_Strain == "E_huxleyi_374",
         Time == 0) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class
  ) %>% 
  mutate(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup()

Repletes <- No_double_peaks %>% 
  filter(Treatment == "Replete") %>% 
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    E_hux_Strain, 
    species,
    #,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation
    lipid_class
  ) %>% 
  summarise(Replete_Proportions = list(Peak_Proportion),
            Replete_Average = mean(Peak_Proportion),
            Replete_St_Dev = sd(Peak_Proportion)) %>% 
  ungroup() %>% 
  select(-Variable_Treatment)

Replete_Countsa <- No_double_peaks %>% 
  filter(Corrected_Per_Cell != 0, Treatment == "Replete") %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Peaks_in_Sample = sum(Corrected_Peak_Area)) %>% 
  ungroup%>% 
  count(Variable_Treatment, E_hux_Strain, Total_Peaks_in_Sample) %>% 
  group_by(E_hux_Strain) %>% 
  summarise(Replicates = length(as.list(Variable_Treatment)),
            All_Counts = list(n),
            Average_Detected_Lipids = mean(n),
            SD = sd(n),
            Total_Agg_Peak_Areas = list(Total_Peaks_in_Sample),
            Average_Aggregated_Peak_Area = mean(Total_Peaks_in_Sample),
            SD_Agg_Peaks = sd(Total_Peaks_in_Sample),)



ggplot(Replete_Counts, aes(x = Average_Aggregated_Peak_Area, y = Average_Detected_Lipids))+
  geom_point(size = 2)+
  geom_point(data = Replete_Countsa, aes(x = Average_Aggregated_Peak_Area, y = Average_Detected_Lipids, color = E_hux_Strain), size = 2)+
  geom_text(data = Replete_Countsa, hjust = 0.5, size = 5, aes(label = E_hux_Strain, color = E_hux_Strain))+
  scale_x_continuous(breaks = seq(0, 5e11, 2e10))+
  theme(axis.text.x = element_text(size = 15),
                                                       axis.text.y = element_text(size = 15),
                                                       axis.title.x = element_text(size = 15),
                                                       axis.title.y = element_text(size = 15),
                                                       plot.title = element_text(hjust = 0.5),
                                                       panel.background = element_rect(fill = "white", 
                                                                                       colour = "black", 
                                                                                       size = 0.5, 
                                                                                       linetype = "solid"),
                                                       panel.grid.major = element_line(size = 0.5, 
                                                                                       linetype = "solid", 
                                                                                       color = "gray"),
                                                       legend.text = element_text(size = 15),
                                                       strip.text = element_text(size = 15))

Replete_Counts <- Repletes %>% filter(Replete_Average != 0) %>% # easy first try to get total number of lipids per species in replete cultures, going to do it again to get stdev as well
  count(E_hux_Strain)
Unique_Replete_Count <- No_double_peaks %>% # gets the same value as the above "easy first try" to get replete counts, which really counted the number of values where the average peak value was greater than 0 (i.e. any of the replicates had a positive value)
  filter(Corrected_Per_Cell != 0, Treatment == "Replete") %>% 
  group_by(E_hux_Strain) %>% 
  summarize(total = length(unique(compound_name)))

Unique_Replete_Count_TAG <- No_double_peaks %>% # gets the same value as the above "easy first try" to get replete counts, which really counted the number of values where the average peak value was greater than 0 (i.e. any of the replicates had a positive value)
  filter(Corrected_Per_Cell != 0, Treatment == "Replete", species == "TAG") %>% 
  group_by(E_hux_Strain) %>% 
  summarize(total = length(unique(compound_name)))

###############################
## BLL supp figure
###############################
Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         Treatment == "Replete"
  ) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  filter(compound_name == "BLL 38:5" |
           compound_name == "BLL 38:6" |
           compound_name == "BLL 39:5" |
           compound_name == "BLL 40:5" |
           compound_name == "BLL 40:6" |
           compound_name == "BLL 40:7" |
           compound_name == "BLL 44:12")

Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "N_Limited")] <- paste("N-Limited")
Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "P_Limited")] <- paste("P-Limited")

Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N-Limited", "P-Limited"))

Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_ericina")] <- paste("H. ericina")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_carterae")] <- paste("P. carterae")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 - 1N")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "I_galbana")] <- paste("I. galbana")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_parvum")] <- paste("P. parvum")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa_ST")] <- paste("P. globosa - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica")] <- paste("P. antarctica")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica_ST")] <- paste("P. antarctica - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_gyrans")] <- paste("P. gyrans")


Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("E. huxleyi 3268 - 1N", 
                                      "E. huxleyi 3266", 
                                      "E. huxleyi 370", 
                                      "E. huxleyi 374", 
                                      "E. huxleyi 379", 
                                      "H. ericina", 
                                      "I. galbana", 
                                      "P. carterae",
                                      "P. gyrans",
                                      "P. parvum", 
                                      "P. globosa", 
                                      "P. globosa - ST", 
                                      "P. antarctica", 
                                      "P. antarctica - ST"))

ggplot(Cultures_GSL, aes(x = compound_name, y = Peak_Proportion, color = compound_name))+
  facet_wrap(~E_hux_Strain, nrow = 3)+
  geom_boxplot(position = position_dodge(width = 0.5))+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15))+
  xlab("Species")+
  ylab("Peak Proportion Per Cell")+
  #ylim(0,0.4)+
  labs(color = "Lipid")

###############################
## GADG (GlcADG) supp figure
###############################
Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         Treatment == "Replete"
  ) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  filter(compound_name == "GADG 32:1" | 
           compound_name == "GADG 34:1" | 
           compound_name == "GADG 38:6" | 
           compound_name == "GADG 34:2")

Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "N_Limited")] <- paste("N-Limited")
Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "P_Limited")] <- paste("P-Limited")

Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N-Limited", "P-Limited"))

Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_ericina")] <- paste("H. ericina")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_carterae")] <- paste("P. carterae")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 - 1N")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "I_galbana")] <- paste("I. galbana")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_parvum")] <- paste("P. parvum")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa_ST")] <- paste("P. globosa - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica")] <- paste("P. antarctica")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica_ST")] <- paste("P. antarctica - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_gyrans")] <- paste("P. gyrans")


Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("E. huxleyi 3268 - 1N", 
                                      "E. huxleyi 3266", 
                                      "E. huxleyi 370", 
                                      "E. huxleyi 374", 
                                      "E. huxleyi 379", 
                                      "H. ericina", 
                                      "I. galbana", 
                                      "P. carterae",
                                      "P. gyrans",
                                      "P. parvum", 
                                      "P. globosa", 
                                      "P. globosa - ST", 
                                      "P. antarctica", 
                                      "P. antarctica - ST"))

ggplot(Cultures_GSL, aes(x = compound_name, y = Peak_Proportion, color = compound_name))+
  facet_wrap(~E_hux_Strain, nrow = 3)+
  geom_boxplot(position = position_dodge(width = 0.5))+
  geom_point(position = position_dodge(width = 0.5))+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15))+
  xlab("Species")+
  ylab("Peak Proportion Per Cell")+
  #ylim(0,0.4)+
  labs(color = "Lipid")




###############################
## Replete Table of Aggregated Pigments
###############################
Replete_No_double_peaks <- DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Treatment == "Replete",
         Time == 0,
         E_hux_Strain != "C_leptoporus") %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    # FA_total_no_C, 
    # FA_total_no_DB,
    # degree_oxidation,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class
  ) %>% 
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup()


minor_pre_aggregation <- Replete_No_double_peaks %>% filter(lipid_class == "pigment" |
                                                      lipid_class == "ubiquinone" |
                                                      grepl("plast", lipid_class), 
                                                    Time == 0)

minor_pre_aggregation$lipid_class <- ifelse(grepl("PQ", minor_pre_aggregation$species), 
                                            as.character("Plastoquinones"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("UQ", minor_pre_aggregation$species), 
                                            as.character("Ubiquinones"),
                                            as.character(minor_pre_aggregation$lipid_class))


minor_pre_aggregation$lipid_class <- ifelse(grepl("Chl_a|Pheophytin_a|Chlide_a|DivinylChl_a|Pheophytin_a2", minor_pre_aggregation$species), 
                                            as.character("Chl_a_etc"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Chl_b|Pheophytin_b2", minor_pre_aggregation$species), 
                                            as.character("Chl_b_etc"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Chl_c2|Chl_c3|Chl_c2 MGDG 14:0_14:0|Chl_c2 MGDG 18:4_14:0", minor_pre_aggregation$species), 
                                            as.character("Chl_c_etc"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Alloxanthin|Crocoxanthin|Ddc_dd|Diatoxanthin|Echinenone|Lutein|Neo_Nosxanthin|Prasinoxanthin|Violaxanthin|Zeaxanthin", minor_pre_aggregation$species), 
                                            as.character("Minor Xanthophylls"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Fucoxanthin", minor_pre_aggregation$species), 
                                            as.character("Fucoxanthin"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Carotene", minor_pre_aggregation$species), 
                                            as.character("Carotenes"),
                                            as.character(minor_pre_aggregation$lipid_class))


minor_pre_aggregation$lipid_class <- ifelse(grepl("19prime_but_fuco", minor_pre_aggregation$species), 
                                            as.character("19prime_but_fuco"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("19prime_hex_fuco", minor_pre_aggregation$species), 
                                            as.character("19prime_hex_fuco"),
                                            as.character(minor_pre_aggregation$lipid_class))

minor_pre_aggregation$lipid_class <- ifelse(grepl("Chl_c1", minor_pre_aggregation$species), 
                                            as.character("Chl_c1"),
                                            as.character(minor_pre_aggregation$lipid_class))


Agg_minor_aggregation <- minor_pre_aggregation %>% 
  group_by(Variable_Treatment, lipid_class, Treatment, E_hux_Strain) %>% 
  summarize(Summed_Peak_Proportion = sum(Peak_Proportion)) %>% 
  ungroup()

pigment_plot_output <- Agg_minor_aggregation # set up a separate df to customize labelling for plot, look at pigment plot output tab for the rest



minor_repletes <- pigment_plot_output %>% 
  filter(Treatment == "Replete") %>% 
  mutate(Variable_Treatment = paste0(E_hux_Strain, "_", 
                                     Treatment)) %>% 
  group_by(
    Variable_Treatment, 
    E_hux_Strain, 
    lipid_class
  ) %>% 
  summarise(Replete_Proportions = list(Summed_Peak_Proportion),
            Replete_Average = formatC(mean(Summed_Peak_Proportion), 
                                      format = "e", 
                                      digits = 2),
            Replete_St_Dev = paste0("±", formatC(sd(Summed_Peak_Proportion), 
                                                 format = "e", 
                                                 digits = 2)))  %>% 
  ungroup() %>% 
  select(-Variable_Treatment) %>% 
  rename("Lipid" = lipid_class,
         "G_species" = E_hux_Strain)

 



minor_agg_repletes <- minor_repletes

minor_agg_repletes$Replete_n <- sapply(minor_agg_repletes$Replete_Proportions, 
                                               function(x) {length(x)})

minor_agg_C_ericina <- minor_agg_repletes %>% filter(G_species == "C_ericina") %>% 
  select(-Replete_Proportions) %>% 
  rename("C_ericina_Average" = Replete_Average,
         "C_ericina_SD" = Replete_St_Dev,
         "C_ericina_n" = Replete_n) %>% 
  select(-G_species)
# 
# C_leptoporus <- minor_agg_repletes %>% filter(G_species == "C_leptoporus")%>% 
#   select(-Replete_Proportions)%>% 
#   rename("C_leptoporus_Average" = Replete_Average,
#          "C_leptoporus_SD" = Replete_St_Dev,
#          "C_leptoporus_n" = Replete_n)%>% 
#   select(-G_species)

minor_agg_E_huxleyi_1N <- minor_agg_repletes %>% filter(G_species == "E_huxleyi_1N")%>% 
  select(-Replete_Proportions)%>% 
  rename("E_huxleyi_1N_Average" = Replete_Average,
         "E_huxleyi_1N_SD" = Replete_St_Dev,
         "E_huxleyi_1N_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_E_huxleyi_370 <- minor_agg_repletes %>% filter(G_species == "E_huxleyi_370")%>% 
  select(-Replete_Proportions)%>% 
  rename("E_huxleyi_370_Average" = Replete_Average,
         "E_huxleyi_370_SD" = Replete_St_Dev,
         "E_huxleyi_370_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_E_huxleyi_374    <- minor_agg_repletes %>% filter(G_species == "E_huxleyi_374")%>% 
  select(-Replete_Proportions)%>% 
  rename("E_huxleyi_374_Average" = Replete_Average,
         "E_huxleyi_374_SD" = Replete_St_Dev,
         "E_huxleyi_374_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_E_huxleyi_379    <- minor_agg_repletes %>% filter(G_species == "E_huxleyi_379")%>% 
  select(-Replete_Proportions)%>% 
  rename("E_huxleyi_379_Average" = Replete_Average,
         "E_huxleyi_379_SD" = Replete_St_Dev,
         "E_huxleyi_379_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_E_huxleyi_F      <- minor_agg_repletes %>% filter(G_species == "E_huxleyi_F")%>% 
  select(-Replete_Proportions)%>% 
  rename("E_huxleyi_3266_Average" = Replete_Average,
         "E_huxleyi_3266_SD" = Replete_St_Dev,
         "E_huxleyi_3266_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_I_galbana       <- minor_agg_repletes %>% filter(G_species == "I_galbana")%>% 
  select(-Replete_Proportions)%>% 
  rename("I_galbana_Average" = Replete_Average,
         "I_galbana_SD" = Replete_St_Dev,
         "I_galbana_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_P_antarctica     <- minor_agg_repletes %>% filter(G_species == "P_antarctica")%>% 
  select(-Replete_Proportions)%>% 
  rename("P_antarctica_Average" = Replete_Average,
         "P_antarctica_SD" = Replete_St_Dev,
         "P_antarctica_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_P_antarctica_ST  <- minor_agg_repletes %>% filter(G_species == "P_antarctica_ST")%>% 
  select(-Replete_Proportions)%>% 
  rename("P_antarctica_ST_Average" = Replete_Average,
         "P_antarctica_ST_SD" = Replete_St_Dev,
         "P_antarctica_ST_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_P_carterae       <- minor_agg_repletes %>% filter(G_species == "P_carterae")%>% 
  select(-Replete_Proportions)%>% 
  rename("P_carterae_Average" = Replete_Average,
         "P_carterae_SD" = Replete_St_Dev,
         "P_carterae_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_P_globosa       <- minor_agg_repletes %>% filter(G_species == "P_globosa")%>% 
  select(-Replete_Proportions)%>% 
  rename("P_globosa_Average" = Replete_Average,
         "P_globosa_SD" = Replete_St_Dev,
         "P_globosa_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_P_globosa_ST     <- minor_agg_repletes %>% filter(G_species == "P_globosa_ST")%>% 
  select(-Replete_Proportions)%>% 
  rename("P_globosa_ST_Average" = Replete_Average,
         "P_globosa_ST_SD" = Replete_St_Dev,
         "P_globosa_ST_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_P_gyrans         <- minor_agg_repletes %>% filter(G_species == "P_gyrans")%>% 
  select(-Replete_Proportions)%>% 
  rename("P_gyrans_Average" = Replete_Average,
         "P_gyrans_SD" = Replete_St_Dev,
         "P_gyrans_n" = Replete_n)%>% 
  select(-G_species)

minor_agg_P_parvum        <- minor_agg_repletes %>% filter(G_species == "P_parvum")%>% 
  select(-Replete_Proportions)%>% 
  rename("P_parvum_Average" = Replete_Average,
         "P_parvum_SD" = Replete_St_Dev,
         "P_parvum_n" = Replete_n)%>% 
  select(-G_species)

#Clade_Table <- merge(C_ericina, C_leptoporus, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_C_ericina, minor_agg_E_huxleyi_1N, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_E_huxleyi_370, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_E_huxleyi_374, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_E_huxleyi_379, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_E_huxleyi_F, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_I_galbana, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_P_antarctica, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_P_antarctica_ST, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_P_carterae, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_P_globosa, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_P_globosa_ST, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_P_gyrans, all = TRUE, by = "Lipid")
minor_agg_Clade_Table <- merge(minor_agg_Clade_Table, minor_agg_P_parvum, all = TRUE, by = "Lipid")

write.csv(minor_agg_Clade_Table, "Haptophyte_Agg_Pigments_Table.csv")


#######################
## Plastoquinones
#######################

PQs <- Replete_No_double_peaks %>% 
  filter(grepl("PQ", compound_name))

ggplot(PQs, aes(x = compound_name, y = Peak_Proportion))+
  geom_boxplot()+
  facet_grid(~E_hux_Strain)







########################
## DGCC
########################

Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         Treatment == "Replete"
  ) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  filter(compound_name == "DGCC 32:5" |
           compound_name == "DGCC 34:4" |
           compound_name == "DGCC 34:5" |
           compound_name == "DGCC 36:5" |
           compound_name == "DGCC 36:6" |
           compound_name == "DGCC 40:7" |
           compound_name == "DGCC 44:12"|
           compound_name == "DGCC 38:6" |
           compound_name == "DGCC 38:10" |
           compound_name == "DGCC 40:10" |
           compound_name == "DGCC 40:11" |
           compound_name == "DGCC 42:11")

Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "N_Limited")] <- paste("N-Limited")
Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "P_Limited")] <- paste("P-Limited")

Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N-Limited", "P-Limited"))

Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_ericina")] <- paste("H. ericina")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_carterae")] <- paste("P. carterae")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 - 1N")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "I_galbana")] <- paste("I. galbana")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_parvum")] <- paste("P. parvum")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa_ST")] <- paste("P. globosa - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica")] <- paste("P. antarctica")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica_ST")] <- paste("P. antarctica - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_gyrans")] <- paste("P. gyrans")


Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("E. huxleyi 3268 - 1N", 
                                      "E. huxleyi 3266", 
                                      "E. huxleyi 370", 
                                      "E. huxleyi 374", 
                                      "E. huxleyi 379", 
                                      "H. ericina", 
                                      "I. galbana", 
                                      "P. carterae",
                                      "P. gyrans",
                                      "P. parvum", 
                                      "P. globosa", 
                                      "P. globosa - ST", 
                                      "P. antarctica", 
                                      "P. antarctica - ST"))

ggplot(Cultures_GSL, aes(x = compound_name, y = Peak_Proportion, color = compound_name))+
  facet_wrap(~E_hux_Strain, nrow = 3)+
  geom_boxplot(position = position_dodge(width = 0.5))+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15))+
  xlab("Species")+
  ylab("Peak Proportion Per Cell")+
  #ylim(0,0.4)+
  labs(color = "Lipid")

##############################
# DGTS/DGTA
##############################

Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         Treatment == "Replete"
  ) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  filter(compound_name == "DGTS_DGTA 32:1" |
           compound_name == "DGTS_DGTA 36:5" |
           compound_name == "DGTS_DGTA 30:0" |
           compound_name == "DGTS_DGTA 34:2"
           )

Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "N_Limited")] <- paste("N-Limited")
Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "P_Limited")] <- paste("P-Limited")

Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N-Limited", "P-Limited"))

Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_ericina")] <- paste("H. ericina")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_carterae")] <- paste("P. carterae")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 - 1N")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "I_galbana")] <- paste("I. galbana")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_parvum")] <- paste("P. parvum")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa_ST")] <- paste("P. globosa - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica")] <- paste("P. antarctica")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica_ST")] <- paste("P. antarctica - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_gyrans")] <- paste("P. gyrans")


Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("E. huxleyi 3268 - 1N", 
                                      "E. huxleyi 3266", 
                                      "E. huxleyi 370", 
                                      "E. huxleyi 374", 
                                      "E. huxleyi 379", 
                                      "H. ericina", 
                                      "I. galbana", 
                                      "P. carterae",
                                      "P. gyrans",
                                      "P. parvum", 
                                      "P. globosa", 
                                      "P. globosa - ST", 
                                      "P. antarctica", 
                                      "P. antarctica - ST"))

ggplot(Cultures_GSL, aes(x = compound_name, y = Peak_Proportion, color = compound_name))+
  facet_wrap(~E_hux_Strain, nrow = 3)+
  geom_boxplot(position = position_dodge(width = 0.5))+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15))+
  xlab("Species")+
  ylab("Peak Proportion Per Cell")+
  #ylim(0,0.4)+
  labs(color = "Lipid")



############################
## PDPT Supp Fig
############################
Cultures_GSL <-   DNPPE_Corrected %>% 
  filter(Experiment == "Nutrient_Limitation", 
         species != "DNPPE", 
         Replicate != "Box4_A",
         Replicate != "Box4_B",
         Replicate != "A",
         Time == 0,
         G_species != "C_leptoporus",
         Treatment == "Replete"
  ) %>% 
  group_by(
    compound_name,
    Variable_Treatment, 
    species,
    Sample_ID,
    Time,
    Experiment,
    Treatment, 
    E_hux_Strain, lipid_class,
    Replicate
  ) %>% 
  
  summarise(Corrected_Per_Cell = sum(Corrected_Per_Cell)) %>% 
  ungroup() %>% 
  group_by(Variable_Treatment) %>% 
  mutate(Total_Corrected_Per_Cell = sum(Corrected_Per_Cell), 
         Peak_Proportion = Corrected_Per_Cell / Total_Corrected_Per_Cell) %>% 
  ungroup() %>% 
  filter(compound_name == "PDPT 34:5" |
           compound_name == "PDPT 38:6" |
           compound_name == "PDPT 36:6" |
           compound_name == "PDPT 32:5" |
           compound_name == "PDPT 40:11" |
           compound_name == "PDPT 40:7" |
           compound_name == "PDPT 44:12")

Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "N_Limited")] <- paste("N-Limited")
Cultures_GSL$Treatment[which(Cultures_GSL$Treatment == "P_Limited")] <- paste("P-Limited")

Cultures_GSL$Treatment <- factor(Cultures_GSL$Treatment, levels = c("Replete", "N-Limited", "P-Limited"))

Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_ericina")] <- paste("H. ericina")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_carterae")] <- paste("P. carterae")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "C_leptoporus")] <- paste("C. leptoporus")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_1N")] <- paste("E. huxleyi 3268 - 1N")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_F")] <- paste("E. huxleyi 3266")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_370")] <- paste("E. huxleyi 370")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_374")] <- paste("E. huxleyi 374")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "E_huxleyi_379")] <- paste("E. huxleyi 379")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "I_galbana")] <- paste("I. galbana")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_parvum")] <- paste("P. parvum")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa_ST")] <- paste("P. globosa - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica")] <- paste("P. antarctica")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_antarctica_ST")] <- paste("P. antarctica - ST")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_globosa")] <- paste("P. globosa")
Cultures_GSL$E_hux_Strain[which(Cultures_GSL$E_hux_Strain == "P_gyrans")] <- paste("P. gyrans")


Cultures_GSL$E_hux_Strain <- factor(Cultures_GSL$E_hux_Strain, 
                                    c("E. huxleyi 3268 - 1N", 
                                      "E. huxleyi 3266", 
                                      "E. huxleyi 370", 
                                      "E. huxleyi 374", 
                                      "E. huxleyi 379", 
                                      "H. ericina", 
                                      "I. galbana", 
                                      "P. carterae",
                                      "P. gyrans",
                                      "P. parvum", 
                                      "P. globosa", 
                                      "P. globosa - ST", 
                                      "P. antarctica", 
                                      "P. antarctica - ST"))

ggplot(Cultures_GSL, aes(x = compound_name, y = Peak_Proportion, color = compound_name))+
  facet_wrap(~E_hux_Strain, nrow = 3)+
  geom_boxplot(position = position_dodge(width = 0.5))+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "gray"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15))+
  xlab("Species")+
  ylab("Peak Proportion Per Cell")+
  #ylim(0,0.4)+
  labs(color = "Lipid")

