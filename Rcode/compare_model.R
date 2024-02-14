#rewrite model comparison script with new model predictions
setwd("~/Documents/Research_local/new_analysis")
library(ggplot2)
library(vegan)
library(tidyr)

#open the file 
new_data_main <- as.data.frame(read.csv("~/Documents/Research_local/new_analysis/model/endpoints_new1.csv", header = T, row.names = 1))
new_data_plot <- as.data.table(read.csv("~/Documents/Research_local/new_analysis/model/endpoints_new1.csv", header = T))
new_data_main <- round(new_data_main)
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

#processing the data so that it is relative abundance instead of absolute abundance
new_data <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(new_data) <- colnames(new_data_main)
for(i in 1:1000){
  temp_data <- new_data_main[i,]
  temp_data <- temp_data/sum(temp_data)
  new_data <- rbind(new_data,temp_data)
}

#-------- bray-curtis and nmds ------

#community composition
#distance matrix - bray-curtis dissimilarity
distance_model <- vegdist(new_data, "bray") #distance matrix

#can we use the vegan distance model for hierarchical clustering? yes
hc_model <- hclust(distance_model, method = "ward.D2")
plot(hc_model)

dd_model <- as.dendrogram(hc_model)
plot(dd_model) #1242X521

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

p_ordered #1242X261

#-------- compare experimental communities to the new model community --------

#new community assembly data - 
compare_df <- as.data.frame(read.csv("~/Documents/Research_local/community/new_CFU_relative.csv", header = T, row.names = 1))
compare_long <- as.data.frame(read.csv("~/Documents/Research_local/community/new_CFU.csv", header = T))

compare_df <- compare_df[1:10,]
compare_long <- compare_long[1:10,c(1,3:11)]
  
Sample_ID <- c(0:9)
Source <- compare_df$Origin

compare_meta <- as.data.frame(cbind(Sample_ID, Source))
compare_df <- compare_df[,2:10]

Sample_ID <- new_data_plot$Population
Source <- rep('M', times = 1000)

model_meta <- as.data.frame(cbind(Sample_ID,Source))

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
  scale_colour_manual(values = c("#3D2C7D","#DED9B3")) +
  geom_jitter(aes(color = Source),width = 0.2, height = 0) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
both_shannon

t.test(Shannon ~ Source, both_meta) #t-test for finding out differences between shannon indices

#2 - community composition
#distance matrix - bray-curtis dissimilarity

both_distance <- as.matrix(vegdist(both_df, "bray")) #distance matrix
#do a permanova
both_adonis <- adonis2(both_distance~Source, data = both_meta, permutations = 1000)
both_adonis

#do nmds to visualize this above stuff
nmds_both <- metaMDS(both_distance, k = 2, maxit = 999, try = 500)
plot(nmds_both)
#save the coordinates
coords_both <- as.data.frame(nmds_both$points)

#3 - plot these
#add columns to data frame 
coords_both$Source = both_meta$Source

p_nmds = ggplot(coords_both, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 21, size = 2.0, stroke = 1.0, aes(colour = Source)) +
  scale_color_manual(values = c("#3D2C7D","#DED9B3")) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
p_nmds

#-------- additional analysis - count #of species/species rank by cluster--------
#1 count number of bacteria/yeast in the community

count_b <- data.frame(new_data) %>%
  rowwise() %>% 
  do(data.frame(sum(which(. > 0.001) <= 5)))

colnames(count_b) <- "num_bac"

count_y = data.frame(new_data) %>%
  rowwise() %>% 
  do(data.frame(sum(which(. > 0.001) > 5)))

colnames(count_y) <- "num_yeast"

count_b_y <- cbind(count_b,count_y)
rownames(count_b_y) <- rownames(new_data)

#counting median, max, and min no. of bacteria and yeasts
median(as.numeric(count_b$num_bac))
max(as.numeric(count_b$num_bac))
min(as.numeric(count_b$num_bac))

median(as.numeric(count_y$num_yeast))
max(as.numeric(count_y$num_yeast))
min(as.numeric(count_y$num_yeast))

hist(count_b$num_bac)
hist(count_y$num_yeast)

#2 rank each species by abundance in model
ranked_data <- data.frame(new_data) %>%
  rowwise() %>% 
  do(data.frame(t(rank(-unlist(.)))))

rownames(ranked_data) <- rownames(new_data)

#sum how many communities for aeach species have which rank
sum_ranked_data <- data.frame(Species=character(),
                              rank=double(),
                              sum_rank=integer())

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","C.paralimentarius","A.malorum",
                  "S.cerevisiae","W.anomalus","K.humilis","K.servazzii")
counter <- 1
i <- 1
for(i in 1:9){
  ranked_col <- ranked_data[[i]]
  ranks <- unique(ranked_col)
  L <- length(ranks)
  j <- 1
  for(j in 1:L){
    rank_num <- which(ranked_col == ranks[j])
    sum_ranked_data[counter,1] <- strain_names[i]
    sum_ranked_data[counter,2] <- ranks[j]
    sum_ranked_data[counter,3] <- length(rank_num)
    counter <- counter + 1
  }
}

#convert the data into wide format
sum_ranked_wide <- pivot_wider(sum_ranked_data, rank, sum_rank) #used the tidyr method instead of the reshape2 method as reshape2 has been deprecated
sum_ranked_wide <- sum_ranked_wide %>% replace(is.na(.), 0)

#plot a graph of the above data to visualise which species is most abundant in most communities
p1 <- ggplot(data = sum_ranked_data, mapping = aes(x = rank, y = sum_rank, colour = Species))+
  theme_bw() +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  scale_colour_manual(values = colour_names) +
  theme(panel.grid.minor.x = element_blank())
p1

p1 <- ggplot(data = ranked_data, mapping = aes(x = rank, y = sum_rank, colour = Species))+
  theme_bw() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  scale_colour_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
p1

#--------trial plot 1- nmds with colour gradient based on hierarchical clustering---------
coords_nmds_model <- coords_both[1:1000,]
coords_nmds_obs <- coords_both[1001:1010,]
coords_nmds_model$Population <- rownames(coords_nmds_model)
nmds_plot_data <- merge(coords_nmds_model, model_clusters, by = 'Population')
nmds_plot_data$order <- as.numeric(nmds_plot_data$order)/1000

#plot nmds coloured by gradient
p_grad = ggplot(nmds_plot_data, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 21, size = 2.0, stroke = 1.0, alpha = 0.75, aes(colour = order)) +
  geom_point(data = coords_nmds_obs, mapping = aes(x = MDS1, y = MDS2), shape = 16, size = 3.0, colour = "#000000") +
  scale_color_viridis_c(option = "B") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
p_grad

#subset of communities that are closest to the model communities
both_distance_sub <- as.data.frame(both_distance[1001:1010,1:1000])
ModelID <- colnames(both_distance_sub)

mindist <- rep(0, times = 1000)

i = 1
for(i in 1:1000){
  vector_temp <- both_distance_sub[i]
  minval <- min(vector_temp)
  mindist[i] <- minval
  
}
ordered_long <- as.data.frame(cbind(ModelID,mindist))
ordered_short <- order(ordered_long$mindist)[1:20] #pick the top 10 model iterations

#making the chosen model IDs into a new data frame
points_chosen <- nmds_plot_data[ordered_short,]

p_grad2 = ggplot(nmds_plot_data, aes(x = MDS1, y = MDS2)) +
  theme_bw() +
  geom_point(shape = 21, size = 2.0, stroke = 1.0, alpha = 0.75, aes(colour = order)) +
  geom_point(data = points_chosen, mapping = aes(x = MDS1, y = MDS2, colour = order),shape = 21, size = 2.5, stroke = 2.0) +
  geom_point(data = coords_nmds_obs, mapping = aes(x = MDS1, y = MDS2), shape = 16, size = 3.0, colour = "#000000") +
  scale_color_viridis_c(option = "B") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none')
p_grad2

#clustered plot of relative abundance
chosen_data <- rbind(new_data[ordered_short,],compare_df)
chosen_data_plot <- rbind(new_data_plot[ordered_short,],compare_long)

distance_chosen <- vegdist(chosen_data, "bray") #distance matrix

#can we use the vegan distance model for hierarchical clustering? yes
hc_chosen <- hclust(distance_chosen, method = "ward.D2")
dd_chosen <- as.dendrogram(hc_chosen)
plot(dd_chosen) #1242X521

#make a data frame with the cluster order
chosen_clustorder<-hc_chosen$labels[c(hc_chosen$order)] 
chosen_clustorder<-as.data.frame(chosen_clustorder)
chosen_clustorder<-setDT(chosen_clustorder, keep.rownames = TRUE)[]

chosen_clusters <- chosen_clustorder %>% rename(order = rn, Population = chosen_clustorder)

#-------plot -------

# change data from wide format to long format
chosen_long <- melt(data = chosen_data_plot, id.vars = c('Population'), variable.name = 'Species', value.name = "Abundance")

chosen_long <- merge(chosen_long, chosen_clusters, by = 'Population')

chosen_long$order<-as.numeric(chosen_long$order)
chosen_long<-arrange(chosen_long, order)

p_ordered2 <- ggplot(data = chosen_long, mapping = aes(x = reorder(Population, order) , y = Abundance, fill = Species)) +
  theme_minimal() +
  geom_bar(stat="identity", position="fill", width = 1.0) + #set width to 1.0 to remove the white space between bars
  scale_fill_manual(values = colour_names) +
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 16),
        legend.position = 'none')

p_ordered2 #1242X261

