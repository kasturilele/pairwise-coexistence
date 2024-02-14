#the package for reading data out of excel files
library(readxl)
library(ggplot2)
library(cowplot)

Strain_order <- c("17B2","0092a","232","460","550","253","163","228","177")
#make a box plot using ggplot2

colour_names_2 <- c("163" = "#88CCEE",
                    "228" = "#332288",
                    "253" = "#117733",
                    "177" = "#44AA99",
                    "17B2" = "#D9D9D9",
                    "0092a" = "#882255",
                    "232" = "#DDCC77",
                    "460" = "#CC6677",
                    "550" = "#AA4499")

#------- plot 1: single species growth curves -----------

GRdata_init <- read_excel("Growth_curve_CFU.xlsx", sheet = "single_all")
GRdata_plot <- GRdata_init
GRdata_plot$Strain <- factor(GRdata_plot$Strain, levels=Strain_order)
GRdata_model_single <- read.csv("~/Documents/Research_local/new_analysis/Growth_curves_single.csv", header = TRUE)
colnames(GRdata_model_single) <- c("Time","Strain","Replicate","CFU_strain_1")
single_list <- list()

for(i in 1:9){
  GRdata_subset <- subset(GRdata_plot, Strain==Strain_order[i])
  model_subset <- subset(GRdata_model_single, Strain==Strain_order[i])
  single_plot <- ggplot() + 
    theme_bw()+
    geom_point(data = GRdata_subset, mapping = aes(x = Time,y = TOTAL_CFUS_well, group = Replicate, colour = Strain)) + #plots the points for strain 1
    geom_line(data = model_subset, mapping = aes(x = Time, y = CFU_strain_1, group = Replicate, colour = Strain), alpha = 0.25) +
    scale_color_manual(values = colour_names_2) +
    labs(title = strain_names[i], x = "time (hrs)") + 
    theme(plot.title = element_text(hjust = 0.5), #centers the title
          legend.position = "none")
  single_list[[length(single_list)+1]] <- single_plot
}

#plot all the plots in a grid
single_all <- plot_grid(
  plotlist = single_list,
  ncol = 3,
  byrow = TRUE
)
single_all #1048X720

#--------plot 2: pairwise growth curves---------
GRdatapw <- read_excel("Growth_curve_CFU.xlsx", sheet = "pairwise_all")
GRdatamodel <- read.csv("~/Documents/Research_local/new_analysis/Growth_curves_pairs_newpars2.csv", header = TRUE)

GRdatapw$Strain_1 <- factor(GRdatapw$Strain_1, levels=Strain_order)
GRdatapw$Strain_2 <- factor(GRdatapw$Strain_2, levels=Strain_order)

GRdatamodel$Strain_1 <- factor(GRdatamodel$Strain_1, levels=Strain_order)
GRdatamodel$Strain_2 <- factor(GRdatamodel$Strain_2, levels=Strain_order)

grplotlist <- list()

ggplot(data = GRdatapw, mapping = aes(x = Time, group = Replicate)) +
  theme_bw()+
  geom_point(aes(y = CFU_Strain_1, colour = Strain_1)) + #plots the points for strain 1
  geom_line(data = GRdatamodel, mapping = aes(x = Time, y = CFU_strain_1, group = Replicate, colour = Strain_1), alpha = 0.25) +
  geom_point(aes(y = CFU_Strain_2, colour = Strain_2)) + #plots strain 2
  geom_line(data = GRdatamodel, mapping = aes(x = Time, y = CFU_strain_2, group = Replicate, colour = Strain_2), alpha = 0.25) +
  scale_color_manual(values = colour_names_2) +
  #scale_y_log10() + #sets the y scale to a log scale instead of a linear scale
  facet_grid(Strain_2~Strain_1) +
  #, scale = "free_y") + #plots all the strains separately
  labs(x = "time (hrs)", y = "total CFU/well") + 
  theme(plot.title = element_text(hjust = 0.5), #centers the title
        legend.position = "none")

ggsave("~/Documents/Research_local/Plots/new_params/pairwise_comparison_5.pdf", width=13.44, height=10.11, units = "in")


