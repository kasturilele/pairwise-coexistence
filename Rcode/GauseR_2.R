#writing a script to get parameters from data in a loop
#so that I can process all my data at once

setwd("~/Documents/Research_local/gauseR")

library(readxl)
library(gauseR)

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

#------------single strain ------------

# load data
GRdata_solo_1 <- read_excel("Growth_curve_CFU.xlsx", sheet = "single_all")
strain_1 <- c("17B2","0092a","232","460","550","253","163","228","177")

timepoints_all <- GRdata_solo_1$Time

#strain_temp <- c("228","460")

i <- 1
counter <- 1 #to add data to a new row of the single_params data frame
single_params <- data.frame(Strains=character(),
                            Replicate=integer(),
                            N0=double(),
                            r=double(),
                            a11=double(),
                            goodness_of_fit=double()
                            #N0_lower=double(),
                            #r_lower=double(),
                            #a11_lower=double(),
                            #N0_upper=double(),
                            #r_upper=double(),
                            #a11_upper=double()
                            )
for(i in 1:9){
  GRdata_temp1 <- subset(GRdata_solo_1, Strain==strain_1[i])
  j <- 1
  for(j in 1:6){
    tryCatch({
      species_CFU <- subset(GRdata_temp1, Replicate==j)$TOTAL_CFUS_well
      timepoints <- subset(GRdata_temp1, Replicate==j)$Time #this gives a vector for timepoints that can be used as an input for GauseR
      names(species_CFU) = number_to_name(strain_1[i])
      pps <- c(0,0,0)
      
      #this is the wrapper function that fits the trajectory and optimizes it over time
      gause_out<-gause_wrapper(time=timepoints, species=species_CFU)
      
      single_params[counter,1] <- number_to_name(strain_1[i])
      single_params[counter,2] <- j
      single_params[counter,3:5] <- gause_out[["parameter_intervals"]][["mu"]]
      #single_params[i,6:8] <- gause_out[["parameter_intervals"]][["lower_sd"]]
      #single_params[i,9:11] <- gause_out[["parameter_intervals"]][["upper_sd"]]
      pps <- gause_out[["parameter_intervals"]][["mu"]]
      #testing the goodness of fit
      pred_short <- get_logistic(time = timepoints, N0 = pps[1], r = pps[2], K = -(pps[2]/pps[3]))
      
      GoF <- test_goodness_of_fit(observed = species_CFU, predicted = pred_short)
      print(GoF)
      single_params[counter,6] <- GoF
      
      counter <- counter + 1
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  } 
}

write.table(single_params, file = "~/Documents/Research_local/gauseR/new_outputs/single_parameters_new.csv", append = TRUE,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)

#------------ single strains while considering median of data --------------
i <- 1
counter <- 1 #to add data to a new row of the single_params data frame
single_params <- data.frame(Strains=character(),
                            N0=double(),
                            r=double(),
                            a11=double()
                            #N0_lower=double(),
                            #r_lower=double(),
                            #a11_lower=double(),
                            #N0_upper=double(),
                            #r_upper=double(),
                            #a11_upper=double()
)
for(i in 1:9){
  GRdata_temp1 <- subset(GRdata_solo_1, Strain==strain_1[i])
  j <- 1
  for(j in 1:9){
    GRdata_median1 <-subset(GRdata_temp1, Time==timepoints[j])$TOTAL_CFUS_well
    species_CFU[j] <- median(GRdata_median1)
  }
  print(species_CFU)
  names(species_CFU) = number_to_name(strain_1[i])
  #this is the wrapper function that fits the trajectory and optimizes it over time
  gause_out<-gause_wrapper(time=timepoints, species=species_CFU)
  
  single_params[counter,1] <- number_to_name(strain_1[i])
  single_params[counter,2:4] <- gause_out[["parameter_intervals"]][["mu"]]
  #single_params[i,6:8] <- gause_out[["parameter_intervals"]][["lower_sd"]]
  #single_params[i,9:11] <- gause_out[["parameter_intervals"]][["upper_sd"]]
  counter <- counter + 1
  
}

#------------ two strains (median) -------------

# load data
GRdata_paired <- read_excel("Growth_curve_CFU.xlsx", sheet = "Pairwise2")
strain_1 <- unique(GRdata_paired$Strain_1)
strain_2 <- unique(GRdata_paired$Strain_2)

timepoints_2 <- unique(GRdata_paired$Time) #this gives a vector for timepoints that can be used as an input for GauseR

k <- 1
counter <- 1
paired_params <- data.frame(Strain_1=character(),
                            Strain_2=character(),
                            N10=double(),
                            N20=double(),
                            r1=double(),
                            r2=double(),
                            a11=double(),
                            a12=double(),
                            a21=double(),
                            a22=double())
for(k in 1:length(strain_1)){
  GRdata_temp3 <- subset(GRdata_paired, Strain_1==strain_1[k])
  i <- 1
  for(i in 1:length(strain_2)){
    GRdata_temp1 <- subset(GRdata_temp3, Strain_2==strain_2[i])
    if(nrow(GRdata_temp1)!=0){
      species_CFU <- data.frame(Strain_1=double(),Strain_2=double())
      j <- 1
      for(j in 1:9){
        GRdata_median1 <-subset(GRdata_temp1, Time==timepoints_2[j])$CFU_Strain_1
        GRdata_median2 <-subset(GRdata_temp1, Time==timepoints_2[j])$CFU_Strain_2
        species_CFU[j,] <- c(median(GRdata_median1),median(GRdata_median2))
      }
      colnames(species_CFU) <- c(number_to_name(strain_1[k]),number_to_name(strain_2[i]))
      gause_out<-gause_wrapper(time=timepoints_2, species=species_CFU)
      paired_params[counter,1] <- number_to_name(strain_1[k])
      paired_params[counter,2] <- number_to_name(strain_2[i])
      paired_params[counter,3:10] <- gause_out[["parameter_intervals"]][["mu"]]
      counter <- counter + 1
    }
  }
}

write.table(paired_params, file = "~/Documents/Research_local/paired_parameters.csv", append = TRUE,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)

#------------ two strains replicates separate -------------

# load data
GRdata_paired <- read_excel("Growth_curve_CFU.xlsx", sheet = "Pairwise")
strain_1 <- unique(GRdata_paired$Strain_1)
strain_2 <- unique(GRdata_paired$Strain_2)

timepoints_2 <- unique(GRdata_paired$Time) #this gives a vector for timepoints that can be used as an input for GauseR

k <- 1
counter <- 1
paired_params <- data.frame(Strain_1=character(),
                            Strain_2=character(),
                            Replicate=integer(),
                            N10=double(),
                            N20=double(),
                            r1=double(),
                            r2=double(),
                            a11=double(),
                            a12=double(),
                            a21=double(),
                            a22=double())
for(k in 1:length(strain_1)){
  GRdata_temp3 <- subset(GRdata_paired, Strain_1==strain_1[k])
  i <- 1
  for(i in 1:length(strain_2)){
    GRdata_temp1 <- subset(GRdata_temp3, Strain_2==strain_2[i])
    if(nrow(GRdata_temp1)!=0){
      j <- 1
      for(j in 1:3){
        GRdata_strain1 <-subset(GRdata_temp1, replicate==j)$CFU_Strain_1
        GRdata_strain2 <-subset(GRdata_temp1, replicate==j)$CFU_Strain_2
        species_CFU <- cbind(GRdata_strain1,GRdata_strain2)
        colnames(species_CFU) <- c(number_to_name(strain_1[k]),number_to_name(strain_2[i]))
        gause_out<-gause_wrapper(time=timepoints_2, species=species_CFU)
        paired_params[counter,1] <- number_to_name(strain_1[k])
        paired_params[counter,2] <- number_to_name(strain_2[i])
        paired_params[counter,3] <- j
        paired_params[counter,4:11] <- gause_out[["parameter_intervals"]][["mu"]]
        counter <- counter + 1
      }
    }
  }
}

write.table(paired_params, file = "~/Documents/Research_local/paired_parameters_2.csv", append = TRUE,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)

#for writing the gause_out data in csv file
write.table(gause_out[["out"]], file = "~/Documents/Research_local/test_out.csv", append = TRUE,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)

#--------- graphs of single params ----------
library(ggplot2)

df_N0<-data.frame(Strains=single_params$Strains,
               means=single_params$N0,
               lower=single_params$N0_lower,
               upper=single_params$N0_upper)

ggplot(df_N0, aes(x=Strains, y=means, fill=Strains)) +
  geom_bar(position=position_dodge(), stat="identity",colour='black') +
  scale_y_log10() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  labs(title = "N0 values from model fitted to data", x = "Strain names", y = "log(N0)") + 
  theme(plot.title = element_text(hjust = 0.5)) #centers the title

df_r<-data.frame(Strains=single_params$Strains,
                  means=single_params$r,
                  lower=single_params$r_lower,
                  upper=single_params$r_upper)

ggplot(df_r, aes(x=Strains, y=means, fill=Strains)) +
  geom_bar(position=position_dodge(), stat="identity",colour='black') +
  #scale_y_log10() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  labs(title = "r values from model fitted to data", x = "Strain names", y = "r") + 
  theme(plot.title = element_text(hjust = 0.5)) #centers the title

df_a11<-data.frame(Strains=single_params$Strains,
                 means=single_params$a11,
                 lower=single_params$a11_lower,
                 upper=single_params$a11_upper)

ggplot(df_a11, aes(x=Strains, y=means, fill=Strains)) +
  geom_bar(position=position_dodge(), stat="identity",colour='black') +
  scale_y_reverse() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  labs(title = "self-interaction coeffcients from model fitted to data", x = "Strain names", y = "a11 (y axis is reversed)") + 
  theme(plot.title = element_text(hjust = 0.5)) #centers the title

