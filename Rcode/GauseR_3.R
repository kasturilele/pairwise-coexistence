# new script based on gauser, to play around with the parameter extractions

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

# ------ two strains gauseR ----------

# load data
GRdata_single <- read.csv("single_parameters.csv", header = TRUE, sep = ",")
GRdata_paired <- read_excel("Growth_curve_CFU.xlsx", sheet = "pairwise_all")
strain_1 <- unique(GRdata_paired$Strain_1)
strain_2 <- unique(GRdata_paired$Strain_2)

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
                            a22=double(),
                            goodness_of_fit=double())
for(k in 1:length(strain_1)){
  GRdata_temp3 <- subset(GRdata_paired, Strain_1==strain_1[k])
  GRsingle1 <- subset(GRdata_single, Strains==number_to_name(strain_1[k]))
  i <- 1
  for(i in 1:length(strain_2)){
    GRdata_temp1 <- subset(GRdata_temp3, Strain_2==strain_2[i])
    GRsingle2 <- subset(GRdata_single, Strains==number_to_name(strain_2[i]))
    if(nrow(GRdata_temp1)!=0){
      j <- 1
      for(j in 1:3){
        tryCatch({
          timepoints <- unique(GRdata_temp1$Time) #this gives a vector for timepoints that can be used as an input for GauseR
          GRdata_strain1 <-subset(GRdata_temp1, replicate==j)$CFU_Strain_1
          GRdata_strain2 <-subset(GRdata_temp1, replicate==j)$CFU_Strain_2
          species_CFU <- cbind(GRdata_strain1,GRdata_strain2) #creating a dataframe that maintains the proper format for GauseR
          colnames(species_CFU) <- c(number_to_name(strain_1[k]),number_to_name(strain_2[i])) #adding column names for species identity
          
          # gause_out<-gause_wrapper(time=timepoints, species=species_CFU,
          #                          parm_signs = c(1,1,-1,-1,-1,-1),
          #                          N_starting = c(GRsingle1$N0[j],GRsingle2$N0[j]),
          #                          r_starting = c(GRsingle1$r[j],GRsingle2$r[j]),
          #                          A_starting = c(GRsingle1$a11[j],NA, NA, GRsingle2$a11[j]),
          #                          keeptimes = TRUE)
          
          gause_out<-gause_wrapper(time=timepoints, species=species_CFU,
                                  parm_signs = c(1,1,-1,-1,-1,-1),
                                  keeptimes = TRUE)
          
          paired_params[counter,1] <- number_to_name(strain_1[k])
          paired_params[counter,2] <- number_to_name(strain_2[i])
          paired_params[counter,3] <- j
          paired_params[counter,4:11] <- gause_out[["parameter_intervals"]][["mu"]]
          
          #calculate goodness of fit
          GoF <- test_goodness_of_fit(observed = species_CFU, predicted = as.data.frame(gause_out["out"])[c(2,3)])
          print(GoF)
          
          paired_params[counter,12] <- GoF
          
          counter <- counter + 1
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
}

#for writing the outputs of gauseR in a csv file
write.table(paired_params, file = "~/Documents/Research_local/gauseR/new_outputs/paired_parameters_new.csv", append = TRUE,
            quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)
