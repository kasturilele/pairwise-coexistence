#new script that consolidates all the model plots + rejection sampling in one script

setwd("~/Documents/Research_local/community/model/new-pairwise") #changed the directory for the new data
library(ggplot2)
library(reshape2)
library(mctoolsr)
library(vegan)
library(data.table)
library(dplyr)

#open the file 

new_data <- as.data.frame(read.csv("~/Documents/Research_local/community/model/new-pairwise/Outnew_1000.csv", header = T, row.names = 1))
new_data_plot <- as.data.table(read.csv("~/Documents/Research_local/community/model/new-pairwise/Outnew_1000.csv", header = T))
new_data <- round(new_data)
new_data_plot[,c(2:10)] <- round(new_data_plot[,c(2:10)]) #to get rid of the one spurious negative value that somehow snuck into the data

colour_names <- c("F.sanfranciscensis" = "#d9d9d9",
                  "L.brevis" = "#882255",
                  "L.plantarum" = "#DDCC77",
                  "A.malorum" = "#CC6677",
                  "C.paralimentarius" = "#aa4499",
                  "S.cerevisiae" = "#117733",
                  "W.anomalus" = "#88ccee",
                  "K.humilis" = "#332288",
                  "K.servazzii" = "#44aa99")

#-------- bray-curtis and nmds ------

#community composition
#distance matrix - bray-curtis dissimilarity
distance_model <- as.matrix(vegdist(new_data, "bray")) #distance matrix

#do nmds to visualize this above stuff
nmds_model <- metaMDS(distance_model, k = 2, maxit = 999, try = 500)
plot(nmds_model)

#-------- make a dendrogram, get the order and use it to order the data frame before plotting -------

#transpose data frame
trans_data <- transpose(new_data)
#trans_data <- transpose(subsampled)

#redefine row and column names
rownames(trans_data) <- colnames(new_data)
colnames(trans_data) <- rownames(new_data)

#create clusters
dm_model_data<-calc_dm(trans_data, method = "bray")
hc_model <- hclust(dm_model_data, method = "ward.D2")
plot(hc_model)

dd_model <- as.dendrogram(hc_model)
plot(dd_model)

#make a data frame with the cluster order
model_clustorder<-hc_model$labels[c(hc_model$order)] 
model_clustorder<-as.data.frame(model_clustorder)
model_clustorder<-setDT(model_clustorder, keep.rownames = TRUE)[]

model_clusters <- model_clustorder %>% rename(order = rn, Population = model_clustorder)

#-------plot -------

# change data from wide format to long format
new_data_long <- melt(data = new_data_plot, id.vars = c('Population'), variable.name = 'Species', value.name = "Abundance")

new_data_long <- merge(new_data_long, model_clusters, by = 'Population')

new_data_long$order<-as.numeric(new_data_long$order)
new_data_long<-arrange(new_data_long, order)

p_ordered <- ggplot(data = new_data_long, mapping = aes(x = reorder(Population, order) , y = Abundance, fill = Species)) +
  theme_minimal() +
  geom_bar(stat="identity", position="fill", width = 1.0) + #set width to 1.0 to remove the white space between bars
  scale_fill_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 16),
        legend.position = 'none')

p_ordered

ggsave("~/Documents/Research_local/Plots/model_data/new_model_plots/new_model_subsampled_dendrogram.pdf", width=12.42, height=2.61, units = "in")

#-------- compare experimental communities to the new model community --------

#load the old files and try to manipulate them
compare_meta <- as.data.frame(read.csv("~/Documents/Research_local/community/Rotation2_metadata.csv", header = T)) #meta
compare_meta <- compare_meta[3:26,]
compare_meta <- compare_meta[c(1,2)]

compare_df <- as.data.frame(read.csv("~/Documents/Research_local/community/Rotation2_CFU.csv", header = T, row.names = 1))
compare_df <- compare_df[3:26,]

Sample_ID <- new_data_plot$Population
Source <- rep('M', times = 1000)

model_meta <- cbind(Sample_ID,Source)
model_meta <- as.data.frame(model_meta)

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum",
                  "C.paralimentarius","S.cerevisiae","W.anomalus","K.humilis","K.servazzii")

both_df <- rbind(new_data,compare_df)
both_meta <- rbind(model_meta,compare_meta)

#---------- run all the analyses from before on this new data set ----------

#1 - diversity indices - alpha diversity
both_meta$Shannon <- diversity(both_df,"shannon")

#plot data wrt treatment (l/s)
both_shannon <- ggplot(both_meta, aes(x = Source, y = Shannon)) +
  theme_bw() +
  geom_boxplot(aes(color = Source)) +
  geom_jitter(aes(color = Source),width = 0.2, height = 0)
both_shannon

#2 - community composition
#distance matrix - bray-curtis dissimilarity

both_distance <- as.matrix(vegdist(both_df, "bray")) #distance matrix
#do a permanova
both_adonis <- adonis2(both_distance~Source, data = both_meta, permutations = 1000)

#do nmds to visualize this above stuff
nmds_both <- metaMDS(both_distance, k = 2, maxit = 999, try = 500)
plot(nmds_both)
#save the coordinates
coords_both <- as.data.frame(nmds_both$points)

#3 - plot these
#add columns to data frame 
coords_both$Source = both_meta$Source

p2 = ggplot(coords_both, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 21, size = 2.0, stroke = 1.0, aes(colour = Source)) +
  scale_color_manual(values = c("#DED9B3","#3D2C7D")) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
#legend.position = 'none')


p2

#1 - making a subset of the distance matrix to work with

#open the file with information on replicates used for each run
new_data_reps <- as.data.frame(read.csv("testOut_1000_reps.csv", header = T, row.names = 1))

#for the new model with both single and pairwise parameters being permuted - 
data_reps_S <- as.data.frame(read.csv("Outnew_1000_reps_single.csv", header = T, row.names = 1))
data_reps_P <- as.data.frame(read.csv("Outnew_1000_reps_paired.csv", header = F, row.names = 1))

#open the files with the parameters
paired_params <- read.csv("~/Documents/Research_local/gauseR/new_outputs/paired_parameters_1.csv", header = TRUE, sep = ",")
single_params <- read.csv("~/Documents/Research_local/gauseR/single_parameters.csv", header = TRUE, sep = ",")

both_distance_sub <- as.data.frame(both_distance[1001:1024,1:1000])
ModelID <- colnames(both_distance_sub)

mindist <- rep(0, times = 1000)

i = 1
for(i in 1:1000){
  vector_temp <- both_distance_sub[i]
  minval <- min(vector_temp)
  mindist[i] <- minval
  
}
ordered_long <- as.data.frame(cbind(ModelID,mindist))
#make a data frame that has info about which parameters were picked for the chosen model communities
x <- 1
all_modelIDs_S <- tibble() #for the new model
all_modelIDs_P <- tibble()

for(x in 1:length(ordered_short)){
  y <- ordered_short[x]
  temp_S <- data_reps_S[y,] #for the new model
  temp_P <- data_reps_P[y,]
  all_modelIDs_S <- rbind(all_modelIDs_S, temp_S)
  all_modelIDs_P <- rbind(all_modelIDs_P, temp_P)
}

#for writing the parameters chosen in a csv file
write.table(all_modelIDs_S, file = "parameters_subset_S10.csv", append = TRUE,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)
write.table(all_modelIDs_P, file = "parameters_subset_P10.csv", append = TRUE,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)

#-------- redoing the NMDS plot, but marking the chosen communities in a different colour ----------

chosen_data <- new_data_plot[ordered_short,]$Population
unchosen_data <- new_data_plot[others,]$Population

Source <- rep('SM', times = length(ordered_short))
Source_other <- rep('M', times = length(others))

chosen_data_meta <- as.data.frame(cbind(chosen_data,Source))
chosen_data_meta <- chosen_data_meta %>% rename(Sample_ID = chosen_data)

others_meta <- as.data.frame(cbind(unchosen_data,Source_other))
others_meta <- others_meta %>% rename(Sample_ID = unchosen_data, Source = Source_other)

chosen_both_meta <- rbind(chosen_data_meta,others_meta,compare_meta)
both_ordered <- rbind(new_data[ordered_short,],new_data[others,],compare_df)

#2 - community composition
#distance matrix - bray-curtis dissimilarity

chosen_distance <- as.matrix(vegdist(both_ordered, "bray")) #distance matrix
#do a permanova
chosen_adonis <- adonis2(chosen_distance~Source, data = chosen_both_meta, permutations = 1000)

#do nmds to visualize this above stuff
nmds_chosen <- metaMDS(chosen_distance, k = 2, maxit = 999, try = 500)
plot(nmds_chosen)

#trying a 3d nmds to see if it looks different from the 2d nmds
nmds_chosen_3d <- metaMDS(chosen_distance, k = 3, maxit = 999, try = 500)
plot(nmds_chosen_3d)
#save the coordinates
coords_chosen <- as.data.frame(nmds_chosen$points)

#3 - plot these
#add columns to data frame 
coords_chosen$Source = chosen_both_meta$Source

p3 = ggplot(coords_chosen, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 21, size = 2.0, stroke = 1.0, aes(colour = Source)) +
  scale_color_manual(values = c("#DED9B3", "#3D2C7D","#FF7206")) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
#legend.position = 'none')


p3

#plot the chosen model compositions
chosen_plot <- new_data_plot[ordered_short,]

# change data from wide format to long format
chosen_long <- melt(data = chosen_plot, id.vars = c('Population',), variable.name = 'Species', value.name = "Abundance")

p_chosen <- ggplot(data = chosen_long, mapping = aes(x = factor(Population, levels = chosen_data) , y = Abundance, fill = Species)) +
  theme_bw() +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = 'none')

p_chosen

#--------- nmds for the old model+rejection sampling model -----------

#1 - assembling the data and metadata
subsampled <- as.data.frame(read.csv("testOut_100_subsampled_10.csv", header = T, row.names = 1))
subsampled_plot <- as.data.table(read.csv("testOut_100_subsampled_10.csv", header = T))

Sample_ID <- new_data_plot$Population
Source <- rep('M', times = 1000)

subsampled_ID <- subsampled_plot$Population
subSource <- rep('SM', times = 100)

sample_meta <- as.data.frame(cbind(Sample_ID,Source))
subsample_meta <- as.data.frame(cbind(subsampled_ID,subSource))
subsample_meta <- subsample_meta %>% rename(Sample_ID = subsampled_ID, Source = subSource)

sample_data <- rbind(new_data,subsampled,compare_df)
sample_meta <- rbind(sample_meta,subsample_meta,compare_meta)

#2 - community composition and nmds
#distance matrix - bray-curtis dissimilarity

sample_distance <- as.matrix(vegdist(sample_data, "bray")) #distance matrix
#do a permanova
sample_adonis <- adonis2(sample_distance~Source, data = sample_meta, permutations = 1000)

#do nmds to visualize this above stuff
nmds_sample <- metaMDS(sample_distance, k = 2, maxit = 999, try = 500)
plot(nmds_sample)

#save the coordinates
coords_sample <- as.data.frame(nmds_sample$points)

#3 - plot these
#add columns to data frame 
coords_sample$Source = sample_meta$Source

p3 = ggplot(coords_sample, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 21, size = 2.0, stroke = 1.0, aes(colour = Source)) +
  scale_color_manual(values = c("#DED9B3","#3D2C7D","#FF7206")) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
#legend.position = 'none')

p3




