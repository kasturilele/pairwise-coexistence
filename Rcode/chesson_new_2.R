#script for calculating pairwise coexistence from chesson's framework
#new : moving the bits from the old script that are to be retained over here
#change to new: changed the function that calculates preds of coexistence to be organised by species id

setwd("~/Documents/Research_local/Plots/new_params")

library(ggplot2)
library(dplyr)
library(cowplot)
#reading data from the pairwise parameters sheet - new data
paired_data <- read.table(file = "~/Documents/Research_local/new_analysis/new_python_outputs/est_pairs_new.txt",
                          sep = ",",
                          header = TRUE)

strain_names <- c("F.sanfranciscensis","L.brevis","L.plantarum","A.malorum","C.paralimentarius","S.cerevisiae", "W.anomalus","K.humilis","K.servazzii")
strain_order <- c("17B2","0092a","232","460","550","253","163","228","177")

#-----------predictions of coexistence-------------
coexist <- data.frame(Strain_1=character(),
                      Strain_2=character(),
                      rho=double(),
                      f2_f1=double())
pred_coexist <- data.frame(Strain_1=character(),
                           Strain_2=character(),
                           coexist = integer(),
                           S1_win = integer(),
                           S2_win = integer(),
                           priority = integer())
c <- 1
for(i in 1:8){
  S1 <- strain_order[i]
  k <- i+1
  for(j in k:9){
    S2 <- strain_order[j]
    flip <- FALSE
    temp_pars <- subset(paired_data, Strain_1 == S1 & Strain_2 == S2)
    if(nrow(temp_pars) < 1){
      temp_pars <- subset(paired_data, Strain_1 == S2 & Strain_2 == S1)
      flip <- TRUE #S1 and S2 are flipped, this affects calculation of f2/f1
    }
    #create new blank data frame
    temp_coexist <- data.frame(Strain_1=character(),
                               Strain_2=character(),
                               rho=double(),
                               f2_f1=double())
    for(r in 1:10){
      #calculate rho and f2/f1 for all pairs
      temp_rho <- (temp_pars[r,"a12"]/temp_pars[r,"a11"])*(temp_pars[r,"a21"]/temp_pars[r,"a22"])
      temp1_f2f1 <- (temp_pars[r,"a11"]/temp_pars[r,"a22"])*(temp_pars[r,"a12"]/temp_pars[r,"a21"])
      temp_f2f1 <- sqrt(abs(temp1_f2f1))*(temp_pars[r,"r2"]/temp_pars[r,"r1"])
      
      
      if(flip == TRUE){
        temp_coexist[r,1] <- temp_pars[r, "Strain_2"]
        temp_coexist[r,2] <- temp_pars[r, "Strain_1"]
        temp_coexist[r,3] <- sqrt(abs(temp_rho))
        temp_coexist[r,4] <- (1/temp_f2f1)
      }else{
        temp_coexist[r,1] <- temp_pars[r, "Strain_1"]
        temp_coexist[r,2] <- temp_pars[r, "Strain_2"]
        temp_coexist[r,3] <- sqrt(abs(temp_rho))
        temp_coexist[r,4] <- (temp_f2f1)
      }
    }
    
    #from the filled in data frame, calculate predictions of coexistence
    rho <- temp_coexist$rho
    f2_f1 <- temp_coexist$f2_f1
    var_coexist <- f2_f1 > rho & f2_f1 < 1/rho
    var_sp1 <- f2_f1 < rho & f2_f1 < 1/rho
    var_sp2 <- f2_f1 > rho & f2_f1 > 1/rho
    var_prior <- f2_f1 < rho & f2_f1 > 1/rho
    
    pred_coexist[c,1] <- S1
    pred_coexist[c,2] <- S2
    pred_coexist[c,3] <- sum(var_coexist)
    pred_coexist[c,4] <- sum(var_sp1)
    pred_coexist[c,5] <- sum(var_sp2)
    pred_coexist[c,6] <- sum(var_prior)
    
    c <- c+1
    
    coexist <- rbind(coexist, temp_coexist)
  }
}

#------------plots -------------
#plot 1 - separate bacteria-bacteria and yeast-yeast competitions, superimposed on OG chesson plot
fun_default <- function(x) {x}
fun_inverse <- function(x) {1/x}
fun_vert <- function(x,y) {y <- 1}

colour_b_y <- c("17B2" = "#DDCC77",
                "0092a" = "#DDCC77",
                "232" = "#DDCC77",
                "460" = "#DDCC77",
                "550" = "#DDCC77",
                "253" = "#332288",
                "163" = "#332288",
                "228" = "#332288",
                "177" = "#332288")

plot_chesson <- ggplot(data = coexist, mapping = aes(x = rho, y = f2_f1)) +
  theme_bw() +
  geom_point(colour = '#212121', shape = 21, size = 1.6, stroke = 1.0, alpha = 0.5) +
  #geom_point(aes(colour = Strain_1, fill = Strain_2), shape = 21, size = 1.6, stroke = 1.0, alpha = 1) +
  scale_fill_manual(values = colour_b_y) +
  scale_color_manual(values = colour_b_y) +
  xlim(0, 3) +
  #ylim(0, 10) +
  scale_y_log10() +
  stat_function(fun = fun_inverse, inherit.aes = F) +
  stat_function(fun = fun_default, inherit.aes = F) +
  theme(axis.text = element_text(size = 16),
        legend.position = "none",
        axis.title = element_blank()) #removing the facet labels - they are not necessary

#make pie charts for the proportion of species that coexist (or not?)
category <- colnames(pred_coexist[,3:6])
count <- 1
plot_list <- list()

for(i in 1:8){
  S1 = strain_names[i]
  for(j in 2:9){
    if(j <= i){
      plot_list[[length(plot_list)+1]] <- pNew
      plot_list[length(plot_list)] <- list(NULL)
      print(length(plot_list))
    }
    else{
      S2 = strain_names[j]
      cat(S1, ", ", S2, "\n")
      print(count)
      amount <- c(t(pred_coexist[count,3:6]))
      
      temp_data <- tibble(category = category,
                          amount = amount)
      
      pNew <- ggplot(temp_data, aes(x="", y=amount, fill=category)) +
        theme_classic() +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        geom_text(aes(label = amount), position = position_stack(vjust=0.5), size = 2.5) +
        labs(x = NULL, y = NULL) +
        scale_fill_manual(values = c("#ED6B3F","#648FFF","#DC267F","#FFB000")) +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position="none")
      plot_list[[length(plot_list)+1]] <- pNew
      count = count + 1
    }
  }
}

#plot all the bar plots in a grid
p_all <- plot_grid(
  plotlist = plot_list,
  ncol = 8,
  byrow = TRUE
)
p_all #1048X1048

