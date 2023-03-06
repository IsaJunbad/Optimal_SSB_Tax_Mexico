# Main Code for "Creating an optimal tax on sugar-sweetened beverages to 
#                maximize societal benefits in Mexico."

# Authors: Juan Carlos Salgado Hernandez, PhD
#          Ana Basto-Abreu, PhD
#          Luis Alberto Moreno Aguilar, PhD
#          Isabel Junquera Badilla, BSc
#          M. Arantxa Colchero, PhD
#          Tonatiuh Barrientos-Gutierrez, PhD

# Last modified the 14 July 2022

# SETUP ########################################################################

#Clean up the environment 
rm(list=ls())      

#Set working directory to the location of the file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load required libraries 
if (!require(devtools)){install.packages("devtools")}
if (!require(bw)){devtools::install_github("INSP-RH/bw", build_vignettes = TRUE)}
if (!require(survey)){install.packages("survey")}
if (!require(tidyverse)){install.packages("tidyverse")}
if (!require(xlsx)){install.packages("xlsx")}

# DATA #########################################################################

#Read data
DATA <- read_csv("DATA_OptimalTax_SSB-Clean.csv") 

#Add variables needed for analisis 
years <- 10
days <- (0:years)*365
WeightCutoff = 30

# TAX SCHEMES ##################################################################

Schemes <- c("UnitedKingdom", "SouthAfrica","Portugal", "Chile",  "Peru", "Ireland", 
             
             "Ecuador", "Thailand", "India", "Kiribati", "Bahrain",
             
             "Philippines", "BerkeleyEEUU", "PhiladelphiaEEUU", "BoulderEEUU",
             
             "Mexico")

Tax_Change <- c(0.069836348, #UnitedKingdom
                0.074335871, #SouthAfrica
                0.114175398, #Portugal
                0.1457658  , #Chile
                0.212508725, #Peru
                0.265940562, #Ireland
                
                0.093740064, #Ecuador
                0.13123609 , #Thailand
                0.262472179, #India
                0.374960256, #Kiribati
                0.46870032 , #Bahrain
                
                0.121862083, #Philippines
                0.221226551, #Berkeley(EEUU)
                0.331839827, #Philadelphia(EEUU)
                0.442453102, #Boulder(EEUU)
                
                0.062599359) #Mexico

Reformulation_Change <- c( 0.0818, #UnitedKingdom
                           0.1496, #SouthAfrica
                           0.0818, #Portugal
                           0.0945, #Chile
                           0.1244, #Peru
                           0.0818, #Ireland
                           
                           0, #Ecuador
                           0, #Thailand
                           0, #India
                           0, #Kiribati
                           0, #Bahrain
                           
                           0, #Philippines
                           0, #Berkeley(EEUU)
                           0, #Philadelphia(EEUU)
                           0, #Boulder(EEUU)
                           
                           0) #Mexico

# ANALISIS VARIABLES ###########################################################

for(j in 1:length(Tax_Change)){
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Change_consumption_Kcal_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Change_consumption_Ml_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_Kcal_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Final_weight_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_weight_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Final_BMI_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_BMI_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Obesity_final_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_of_obesity_", Schemes[j])
  
  DATA[, ncol(DATA) + 1] <-  double(nrow(DATA))
  names(DATA)[ncol(DATA)] <- paste0("Reduction_relative_obesity_", Schemes[j])
}

# SIMULATIONS ##################################################################

for (j in 1:length(Tax_Change)) {
  
  DATA[, paste0("Reduction_Kcal_", Schemes[j])] <- -1*DATA$SSB_cal*(1-Reformulation_Change[j])*Tax_Change[j]*DATA$Elasticity - DATA$SSB_cal*Reformulation_Change[j]
  
  #We define the matrix of daily changes in energy intake
  mat_reduc    <- matrix(rep(unlist(DATA[, paste0("Reduction_Kcal_", Schemes[j])]), years + 1), ncol = years  +1)
  EIchange_SSB <- energy_build(mat_reduc, days, "Stepwise_R")
  
  rm(mat_reduc) #Eliminate mat_reduc to save memory
  
  ### We run the Hall model 
  
  SSB_model <- adult_weight(bw = DATA$Weight, ht = (DATA$Height)/100,
                                       age = DATA$Age, sex = DATA$Sex,
                                       EIchange = EIchange_SSB, days = max(days))
  rm(EIchange_SSB) #Eliminate EIchange_SSB to save memory
  
  ### We save antropometric data 
  
  DATA[, paste0("Change_consumption_Ml_", Schemes[j])] <- -1*Tax_Change[j]*DATA$Elasticity*DATA$SSB_ml
  
  DATA[, paste0("Change_consumption_Kcal_", Schemes[j])] <- -1*Tax_Change[j]*DATA$Elasticity
  
  DATA[, paste0("Final_weight_", Schemes[j])] <- SSB_model$Body_Weight[, 365*years ] #Body weight at day max(days) +`1 (last day)`
  
  DATA[, paste0("Reduction_weight_", Schemes[j])] <- DATA[, paste0("Final_weight_", Schemes[j])] - DATA$Weight 
  
  DATA[, paste0("Final_BMI_", Schemes[j])]  <- SSB_model$Body_Mass_Index[, 365*years ] #Body weight at day max(days) +`1 (last day)`
  
  DATA[, paste0("Reduction_BMI_", Schemes[j])]  <- DATA[, paste0("Final_BMI_", Schemes[j])] - DATA$BMI
  
  DATA[, paste0("Obesity_final_", Schemes[j])] <- ifelse(DATA[, paste0("Final_BMI_", Schemes[j])] >= WeightCutoff, 1, 0)
  
  DATA[, paste0("Reduction_of_obesity_", Schemes[j])] <- DATA[, paste0("Obesity_final_", Schemes[j])] - DATA$Obesity
  
  DATA[, paste0("Reduction_relative_obesity_", Schemes[j])] <- DATA[, paste0("Reduction_of_obesity_", Schemes[j])]
  
  rm(SSB_model) #We eliminate the antropometric simulations, since we saved the ones we need 
  
  print(Schemes[j]) #We print the scheme we are in
  
}

# SURVEY DESIGN ################################################################

#Survey design 
Design_SSB <- svydesign(id= ~ID, strata= ~Strata, weights= ~Weight_svy, PSU= ~PSU, 
                        data= DATA, nest = TRUE)
options(survey.lonely.psu = "adjust")


#We calculate the mean of the obesity prevalence and use it to make the relative change in obesity prevalence 
mean_Obesity <- svymean(~Obesity, Design_SSB, na.rm = TRUE)
for (j in 1:length(Tax_Change)) {
  
  DATA[, paste0("Reduction_relative_obesity_", Schemes[j])] <- DATA[, paste0("Reduction_relative_obesity_", Schemes[j])]/mean_Obesity
  
}

#We run the survey design with the relative change in obesity prevalence
Design_SSB <- svydesign(id= ~ID, strata= ~Strata, weights= ~Weight_svy, PSU= ~PSU, 
                        data= DATA, nest = TRUE)
options(survey.lonely.psu = "adjust")

# RESULTS ######################################################################

#Create new matrix to save population results
SimReduc <- (matrix(nrow = length(Tax_Change), ncol = 21))

#Fill in the matrix with the population results
for (j in 1:length(Tax_Change)) {
  SimReduc[j, 1] <- paste0(round(Tax_Change[j]*100, digits = 1), "%")

  SimReduc[j, 2] <- paste0(Reformulation_Change[j]*100, "%")
  
  ### Consumption 
  Change_consumption_Kcal_ <- paste0("svymean(~Change_consumption_Kcal_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 3] <- paste0(round(eval(parse(text = Change_consumption_Kcal_))[1]*100, digits = 1), "%")
  
  ### Calories
  Reduction_Kcal <- paste0("svymean(~Reduction_Kcal_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 4] <- round(eval(parse(text = Reduction_Kcal))[1], digits = 1)
  
  CIKCAL <- paste0("confint(", Reduction_Kcal, ")")
  SimReduc[j, 5] <- paste0("( ", round(eval(parse(text = CIKCAL))[1], digits = 1), 
                           " , ", round(eval(parse(text = CIKCAL))[2], digits = 1), " )")
  
  ### Weight
  weightfinal <- paste0("svymean(~Reduction_weight_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 6] <- round(eval(parse(text = weightfinal))[1], digits = 1)
  
  CIW <- paste0("confint(", weightfinal, ")")
  SimReduc[j, 7] <- paste0("( ", round(eval(parse(text = CIW))[1], digits = 1),
                           " , ", round(eval(parse(text = CIW))[2], digits = 1), " )")
  
  ### BMI
  BMIfinal <- paste0("svymean(~Reduction_BMI_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 8] <- round(eval(parse(text = BMIfinal))[1], digits = 1)
  
  CIBMI <- paste0("confint(", BMIfinal, ")")
  SimReduc[j, 9] <-  paste0("( ", round(eval(parse(text = CIBMI))[1], digits = 1),
                            " , ", round(eval(parse(text = CIBMI))[2], digits = 1), " )")
  
  ### Obesity
  obesity_final <- paste0("svymean(~Obesity_final_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 10] <- round(eval(parse(text = obesity_final))[1]*100, digits = 1)
  
  CIO <- paste0("confint(", obesity_final, ")")
  SimReduc[j, 11] <- paste0("( ", round(eval(parse(text = CIO))[1]*100, digits = 1),
                            " , ", round(eval(parse(text = CIO))[2]*100, digits = 1), " )")
  
  ### Reduction obesity
  obes_reducfinal <- paste0("svymean(~Reduction_of_obesity_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 12] <- round(eval(parse(text = obes_reducfinal))[1]*100, digits = 1)
  
  CIOR <- paste0("confint(", obes_reducfinal, ")")
  SimReduc[j, 13] <-  paste0("( ",  round(eval(parse(text = CIOR))[1]*100, digits = 1),
                             " , ", round(eval(parse(text = CIOR))[2]*100, digits = 1), " )")
  
  ### Relative Reduction obesity
  obes_reducRelfinal <- paste0("svymean(~Reduction_relative_obesity_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 14] <- round(eval(parse(text = obes_reducRelfinal))[1]*100, digits = 1)
  
  CIORr <- paste0("confint(", obes_reducRelfinal, ")")
  SimReduc[j, 15] <- paste0("( ", round(eval(parse(text = CIORr))[1]*100, digits = 1),
                            " , ", round(eval(parse(text = CIORr))[2]*100, digits = 1), " )")
  
  
  #We calculate the mean and confidence interval of initial consumption 
  SimReduc[j, 16] <- round(svymean(~SSB_cal, Design_SSB, na.rm = TRUE), digits = 1)
  SimReduc[j, 17] <- paste0("( ", round(confint(svymean(~SSB_cal, Design_SSB, na.rm = TRUE))[1], digits = 1),  
                            " , ", round(confint(svymean(~SSB_cal, Design_SSB, na.rm = TRUE))[2], digits = 1), " )")
  
  #We calculate the mean and confidence interval of initial obesity prevalence
  SimReduc[j, 18] <- round(svymean(~Obesity, Design_SSB, na.rm = TRUE)*100 , digits = 1)
  SimReduc[j, 19] <- paste0("( ", round(confint(svymean(~Obesity, Design_SSB, na.rm = TRUE))[1]*100, digits = 1),  
                            " , ", round(confint(svymean(~Obesity, Design_SSB, na.rm = TRUE))[2]*100, digits = 1), " )")

  ### Consumption ML
  Cambio_ConsumptionMl_ <- paste0("svymean(~Change_consumption_Ml_", Schemes[j], ", Design_SSB, na.rm = TRUE)")
  SimReduc[j, 20] <- round(eval(parse(text = Cambio_ConsumptionMl_))[1], digits = 1)
  
  CIMl <- paste0("confint(", Cambio_ConsumptionMl_, ")")
  SimReduc[j, 21] <-  paste0("( ", round(eval(parse(text = CIMl))[1], digits = 1),
                            " , ", round(eval(parse(text = CIMl))[2], digits = 1), " )")
  
  print(Schemes[j]) #We print the scheme we are in
  }

#Turn the SimReduc matrix into a data-set
SimReduc <- data.frame(SimReduc)

#Add the schemes 
SimReduc <- cbind(Schemes,SimReduc) 

#Add the col names 
names(SimReduc) <- c("Country","Price increase", "Reformulation", "Decreased consumption (%)",
                     
                     "Reduction of Kcal", "Confidence interval reduction Kcal",
                     "Reduction of weight", "Confidence interval reduction weight",
                     "Reduction of BMI" , "Confidence interval reduction BMI",
                     "Obesity post-tax scheme", "Confidence interval obesity post-tax scheme",
                     "Reduction of obesity", "Confidence interval reduction of obesity",
                     "Relative reduction of obesity", "Confidence interval relative reduction of obesity",
                     "Consumption pre-tax", "Confidence interval consumption pre-tax",
                     "Obesity pre-tax", "Confidence interval obesity pre-tax",
                     "Decreased consumption in ml", "Confidence decreased consumption in ml")



# SAVE #########################################################################

write.xlsx(SimReduc, 
           row.names=FALSE,
           file = paste0("SSB 2018_",Sys.Date(),"_","MAIN.xlsx"),
           append = FALSE)
