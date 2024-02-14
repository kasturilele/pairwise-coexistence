#script to look at the new parameter estimates + plot them on existing growth curves
setwd("~/Documents/Research_local/new_analysis/new_python_outputs")

#using ggplot to plot nicer? looking graphs
library(ggplot2)
library(dplyr)
library(readxl)
library(cowplot)

#--------- strain to strain names function ---------
#the role of this function is to return the strain names when strain numbers are an input
#this is to give row names to the parameter outputs
number_to_name <- function(strain_input){
  strain_numbers <- c("17B2","0092a","232","460","550","253","163","228","177")
  strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum",
                    "C.paralimentarius","S.cerevisiae","W.anomalus","K.humilis","K.servazzii")
  a <- strain_input
  b <- strain_names[match(strain_input, strain_numbers)]
  return(b)
}

#--------- plots ------------
strain_order <- c("17B2","0092a","232","460","550","253","163","228","177")
strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","C.paralimentarius","A.malorum",
                  "S.cerevisiae","W.anomalus","K.humilis","K.servazzii")

#strain-specific colours
colour_names1 <- c("17B2" = "#d9d9d9",
                  "0092a" = "#882255",
                  "232" = "#DDCC77",
                  "460" = "#CC6677",
                  "550" = "#aa4499",
                  "253" = "#117733",
                  "163" = "#88ccee",
                  "228" = "#332288",
                  "177" = "#44aa99")

colour_names <- c("F.sanfranciscensis" = "#d9d9d9",
                  "L.brevis" = "#882255",
                  "L.plantarum" = "#DDCC77",
                  "A.malorum" = "#CC6677",
                  "C.paralimentarius" = "#aa4499",
                  "S.cerevisiae" = "#117733",
                  "W.anomalus" = "#88ccee",
                  "K.humilis" = "#332288",
                  "K.servazzii" = "#44aa99",
                  "Median data point" = "#000000")

#opening files
single_params <- read.table(file = "est_single.txt",
                            sep = ",",
                            header = TRUE)
Species <- number_to_name(single_params$Strains)
single_params <- cbind(single_params,Species)
single_params_old <- read.csv("~/Documents/Research_local/gauseR/new_outputs/single_parameters_new.csv")

paired_params <- read.table(file = "est_pairs_new.txt",
                            sep = ",",
                            header = TRUE)
paired_params_old <- read.csv("~/Documents/Research_local/gauseR/new_outputs/paired_parameters_new.csv")
Species_1 <- number_to_name(paired_params$Strain_1)
Species_2 <- number_to_name(paired_params$Strain_2)
paired_params <- cbind(paired_params,Species_1,Species_2)

#checking if the single params are the same in the pairwise estimates
temp_r1 <- paired_params[1:50,]$r2
temp_r2 <- paired_params[51:100,]$r2

temp_a1 <- paired_params[1:50,]$a22
temp_a2 <- paired_params[51:100,]$a22

#plot 1 - distribution of single parameters
r1_graph <- ggplot() + 
  theme_bw() +
  geom_point(data = single_params, mapping = aes(x = factor(Species, level = strain_names), y = r, colour = Species), shape = 16, size = 5.0, stroke = 2.0, alpha = 0.025) +
  geom_point(data = single_params_old, mapping = aes(x = factor(Strains, level = strain_names), y = r, colour = Strains), size = 2.0) +
  scale_color_manual(values = colour_names) +
  labs(x = "Species", y = "calculated r value") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") 
r1_graph  
#1048X640

a11_graph <- ggplot() + 
  theme_bw() +
  geom_point(data = single_params, mapping = aes(x = factor(Species, level = strain_names), y = a11, colour = Species), shape = 16, size = 5.0, stroke = 2.0, alpha = 0.025) +
  geom_point(data = single_params_old, mapping = aes(x = factor(Strains, level = strain_names), y = a11, colour = Strains), size = 2.0) +
  scale_color_manual(values = colour_names) +
  labs(x = "Species", y = "calculated a11 value") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") 
a11_graph 

#making a new data frame with all pairwise data
paired_temp1 <- paired_params[c("Species_1","Species_2","a12")]
colnames(paired_temp1) <- c("Species_i","Species_j","aij")
paired_temp2 <- paired_params[c("Species_2","Species_1","a21")]
colnames(paired_temp2) <- c("Species_i","Species_j","aij")
paired_params_combi <- rbind(paired_temp1,paired_temp2)
data_ID <- rep("aij", times = 3600)
paired_params_combi <- cbind(paired_params_combi, data_ID)
  
paired_temp1 <- paired_params_old[c("Strain_1","Strain_2","a12")]
colnames(paired_temp1) <- c("Species_i","Species_j","aij")
paired_temp2 <- paired_params_old[c("Strain_2","Strain_1","a21")]
colnames(paired_temp2) <- c("Species_i","Species_j","aij")
paired_params_old_combi <- rbind(paired_temp1,paired_temp2)
data_ID <- rep("aij", times = 210)
paired_params_old_combi <- cbind(paired_params_old_combi, data_ID)

aij_graph <- ggplot() + 
  theme_bw() +
  geom_point(data = paired_params_combi, mapping = aes(x = data_ID, y = aij, colour = Species_j), shape = 16, size = 5.0, stroke = 2.0, alpha = 0.25) +
  #geom_point(data = paired_params_old_combi, mapping = aes(x = data_ID, y = aij, colour = Species_j), size = 2.0) +
  scale_color_manual(values = colour_names) +
  labs(x = "Species", y = "calculated r value") + 
  facet_grid(factor(Species_j, level = strain_names) ~factor(Species_i, level = strain_names), scale = "free_y") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") 

aij_graph #1048X960

#plot new paired params on top of old ones
paired_params_new <- read.table(file = "est_pairs_new2.txt",
                                sep = ",",
                                header = TRUE)

Species_1 <- number_to_name(paired_params_new$Strain_1)
Species_2 <- number_to_name(paired_params_new$Strain_2)
paired_params_new <- cbind(paired_params_new,Species_1,Species_2)

#making a new data frame with all pairwise data
paired_temp1 <- paired_params[c("Species_1","Species_2","a12")]
colnames(paired_temp1) <- c("Species_i","Species_j","aij")
paired_temp2 <- paired_params[c("Species_2","Species_1","a21")]
colnames(paired_temp2) <- c("Species_i","Species_j","aij")
paired_params_combi <- rbind(paired_temp1,paired_temp2)
data_ID <- rep("old", times = 3600)
paired_params_combi <- cbind(paired_params_combi, data_ID)

paired_temp1 <- paired_params_new[c("Species_1","Species_2","a12")]
colnames(paired_temp1) <- c("Species_i","Species_j","aij")
paired_temp2 <- paired_params_new[c("Species_2","Species_1","a21")]
colnames(paired_temp2) <- c("Species_i","Species_j","aij")
paired_params_new_combi <- rbind(paired_temp1,paired_temp2)
data_ID <- rep("new", times = 720)
paired_params_new_combi <- cbind(paired_params_new_combi, data_ID)
paired_params_new_combi <- rbind(paired_params_combi, paired_params_new_combi)

aij_graph2 <- ggplot() + 
  theme_bw() +
  #geom_point(data = paired_params_combi, mapping = aes(x = data_ID, y = aij, colour = Species_j), shape = 21, size = 5.0, stroke = 2.0, alpha = 0.25) +
  geom_boxplot(data = paired_params_new_combi, mapping = aes(x = data_ID, y = aij, colour = Species_j)) +
  scale_color_manual(values = colour_names) +
  scale_fill_manual(values = colour_names) +
  labs(x = "Species", y = "calculated r value") + 
  facet_grid(factor(Species_j, level = strain_names) ~factor(Species_i, level = strain_names), scale = "free_y") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none") 
aij_graph2
