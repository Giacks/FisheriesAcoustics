##########################################################################################
#####            CREATING THE AGE-LENGTH KEY AND DERIVING N AT AGE FOR 2024           ####
##########################################################################################

## Copyright 2024 Technical University of Denmark
## Author: Giacomo Gardella <giaga@aqua.dtu.dk>


##########################################################################################
## The age length key is estimated from survey data from 2018-2023 and used in          ##
## in combination with the length distribution estimated from Christoffer's model to    ##
## split the cod densities estimated from echo-interation and echo-counting among the   ##
## relevant ages.                                                                       ##
##########################################################################################

# Load relevant libraries
library(readxl)
library(mgcv)
library(ggplot2)
library(tidyr)
library(FSA)

##########################################################################################

# Set working directory
setwd("C:/Giacomo Gardella/Comparisons LSSS")

# Load survey data - now can load ALK directly
dfNuukIndex <- as.data.frame(read_excel("C:/Giacomo Gardella/Comparisons LSSS/CODICOD_ALK_2018_2023.xlsx", sheet = 1))
names(dfNuukIndex) = c("Length_cm", "Age1", "Age2", "Age3", "Age4", "Age5", "Age6", "Age7","Age8","Age9","Age10")
head(dfNuukIndex)
dim(dfNuukIndex)

# Add 0.5 cm to all lengths, since cod are measured to the closest cm
dfNuukIndex$Length_cm = dfNuukIndex$Length_cm+0.5

# Rename
dfNuukIndex_Prop = dfNuukIndex


##############         Interpolation

# Find indices where row total is ==0. Length is index + 1, since the first row corresponds to length = 0
dfNuukIndex_Prop$LengthTot <- dfNuukIndex_Prop$Age1+dfNuukIndex_Prop$Age2+dfNuukIndex_Prop$Age3+dfNuukIndex_Prop$Age4+dfNuukIndex_Prop$Age5+
  dfNuukIndex_Prop$Age6+dfNuukIndex_Prop$Age7+dfNuukIndex_Prop$Age8+dfNuukIndex_Prop$Age9+dfNuukIndex_Prop$Age10
indices <- which(dfNuukIndex_Prop$LengthTot==0)

# Interpolation criterium
# Up to 14 cm - these are 0-1 year olds, put in age class 1
# Above 83 cm - put age 10 for now, better interpolation needed

# Iterate over these indices
for (i in indices) {
  if (i < 13) {
    dfNuukIndex_Prop$Age1[i] = 1
  } 
  else if (i > 83) {
  dfNuukIndex_Prop$Age10[i] = 1
  # This threshold will need to be smoothened out
    
  }
}

# Remove useless rows
dfNuukIndex_Prop = dfNuukIndex_Prop[,c(1:11)]


##########     GAM - TEST

# Transform df in long format
dfNuukIndex_Prop_Long <- dfNuukIndex_Prop %>%
  pivot_longer(cols = starts_with("Age"),
               names_to = "Age",
               names_prefix = "Age",
               values_to = "Prop")
dfNuukIndex_Prop_Long = as.data.frame(dfNuukIndex_Prop_Long)
head(dfNuukIndex_Prop_Long)

# Make ordered factor column
dfNuukIndex_Prop_Long$fAge =  as.factor(dfNuukIndex_Prop_Long$Age)

# Test GAM Plot
ggplot(dfNuukIndex_Prop_Long, aes(x=Length_cm, y=Prop, colour = fAge))+
  geom_point()+
  geom_smooth(method="gam",
              method.args=list(family="binomial"),
              formula=y~s(x),
              se=FALSE)+
  labs(x = "Length (cm)",
    y = "Proportion",
    colour = "Age (years)") +
  theme_classic()
#ggsave("alk plot gam.png", height=4, width=5.5)


####################################################################################
######                        CREATE GAM TO SMOOTHEN ALK                  ##########
####################################################################################


#Create GAM
gam_alk = gam(Prop~s(Length_cm)*fAge, data=dfNuukIndex_Prop_Long, family = "quasibinomial", se=FALSE)

dfPred_ALK <- data.frame(Length_cm=0:150,
                        pred_age1=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(1,151))), type="response"),
                        pred_age2=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(2,151))), type="response"),
                        pred_age3=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(3,151))), type="response"),
                        pred_age4=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(4,151))), type="response"),
                        pred_age5=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(5,151))), type="response"),
                        pred_age6=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(6,151))), type="response"),
                        pred_age7=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(7,151))), type="response"),
                        pred_age8=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(8,151))), type="response"),
                        pred_age9=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(9,151))), type="response"),
                        pred_age10=predict(gam_alk, data.frame(Length_cm=0:150, fAge=as.factor(rep(10,151))), type="response"))

dfPred_ALK_Norm <- dfPred_ALK
pred_sums <- rowSums(dfPred_ALK_Norm[, -1]) # Sum predictions across all ages for each length

for (i in 2:ncol(dfPred_ALK_Norm)) { # Start at column 2 since column 1 is 'Length'
  dfPred_ALK_Norm[, i] <- dfPred_ALK_Norm[, i] / pred_sums
}

# Area plot
dfPred_ALK_Norm_Mat <- as.matrix(dfPred_ALK_Norm[,2:11])
dfPred_ALK_Norm_Mat_Num = dfPred_ALK_Norm_Mat
colnames(dfPred_ALK_Norm_Mat_Num) <- 1:ncol(dfPred_ALK_Norm_Mat_Num)
alkPlot(dfPred_ALK_Norm_Mat_Num,"area")

# Transform to long format
dfPred_ALK_Norm_Long <- dfPred_ALK_Norm %>%
  pivot_longer(cols = starts_with("pred_age"), 
               names_to = "Age", 
               values_to = "Proportion")

# Plot with length distribution
ggplot(data = NULL) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age1",], 
            aes(x = Length_cm, y = Proportion, color = "Age 1"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age2",], 
            aes(x = Length_cm, y = Proportion, color = "Age 2"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age3",], 
            aes(x = Length_cm, y = Proportion, color = "Age 3"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age4",], 
            aes(x = Length_cm, y = Proportion, color = "Age 4"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age5",], 
            aes(x = Length_cm, y = Proportion, color = "Age 5"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age6",], 
            aes(x = Length_cm, y = Proportion, color = "Age 6"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age7",], 
            aes(x = Length_cm, y = Proportion, color = "Age 7"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age8",], 
            aes(x = Length_cm, y = Proportion, color = "Age 8"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age9",], 
            aes(x = Length_cm, y = Proportion, color = "Age 9"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age10",], 
            aes(x = Length_cm, y = Proportion, color = "Age 10"), lwd = 1.2) +
  scale_x_continuous(limits = c(0,80),breaks=seq(0,80, by=5))+
  labs(x = "Length (cm)", y = "Proportion", color = "") +
  theme_classic(base_size = 14) +
  scale_color_manual(
    values = c(
      "Age 1" = "darkred",
      "Age 2" = "#F8766D", 
      "Age 3" = "orange", 
      "Age 4" = "#B79F00", 
      "Age 5" = "#00BA38", 
      "Age 6" = "darkolivegreen", 
      "Age 7" = "lightblue", 
      "Age 8" = "#619CFF", 
      "Age 9" = "purple", 
      "Age 10" = "purple3"
    ),
    labels = c(
      "Age 1", "Age 10", "Age 2", "Age 3", "Age 4", "Age 5", "Age 6", 
      "Age 7", "Age 8", "Age 9")
  )
#ggsave(filename="ALK corrected no LD.png", bg="white", width=8, height=4, dpi=700)




# Plot with length distribution
ggplot(data = NULL) +
  geom_line(data = dfLdist_long, 
            aes(x = Length_cm, y = Proportion, group = Cell_ID, color = "Length distribution used for each cell"), 
            lwd = 1.1) +
  geom_line(data = dfMeanLdist_long, 
            aes(x = Length_cm, y = Proportion, color = "Mean length distribution"), 
            lwd = 1.5) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age1",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 1"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age2",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 2"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age3",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 3"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age4",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 4"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age5",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 5"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age6",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 6"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age7",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 7"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age8",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 8"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age9",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 9"), lwd = 1.2) +
  geom_line(data = dfPred_ALK_Norm_Long[dfPred_ALK_Norm_Long$Age == "pred_age10",], 
            aes(x = Length_cm, y = Proportion / 4.1, color = "Age 10"), lwd = 1.2) +
  scale_x_continuous(limits = c(0,80),breaks=seq(0,80, by=5))+
  scale_y_continuous(limits = c(0, 0.25)) +
  labs(x = "Length (cm)", y = "Proportion", color = "") +
  theme_classic(base_size = 14) +
  scale_color_manual(
    values = c(
      "Length distribution used for each cell" = "grey", 
      "Mean length distribution" = "black", 
      "Age 1" = "darkred",
      "Age 2" = "#F8766D", 
      "Age 3" = "orange", 
      "Age 4" = "#B79F00", 
      "Age 5" = "#00BA38", 
      "Age 6" = "darkolivegreen", 
      "Age 7" = "lightblue", 
      "Age 8" = "#619CFF", 
      "Age 9" = "purple", 
      "Age 10" = "purple3"
    ),
    labels = c(
      "Age 1", "Age 10", "Age 2", "Age 3", "Age 4", "Age 5", "Age 6", 
      "Age 7", "Age 8", "Age 9", 
      "Length distribution\nused for each cell", 
      "Mean length\ndistribution"
    )
  )
#ggsave(filename="ALK corrected overlayed to LD.png", bg="white", width=8, height=4, dpi=700)




###########################################################################################
######                                                                              #######
######                      Loading length distribution files                       #######
######                                                                              #######
###########################################################################################

# Run EchoIntegration scrip and check dfAcousticData_Integration_int 

dim(dfAcousticData_Integration_int)
head(dfAcousticData_Integration_int)

dfAcousticData_int_Age <- dfAcousticData_Integration_int

# Remember to rename files in directories so that 0025- is 0025 and 0385+ is 0385

# Create veraible that has unique value for each cell
dfAcousticData_int_Age$Cell_ID <- paste(
  "N", sprintf("%05.2f", dfAcousticData_int_Age$Rect_Lat),
  "W", sprintf("%05.2f", abs(dfAcousticData_int_Age$Rect_Lon)),
  "D", sprintf("%04d", dfAcousticData_int_Age$DepthMean_m), sep="")

dfAcousticData_int_Age$LD_ID = paste(
  "N", sprintf("%05.2f", dfAcousticData_int_Age$Rect_Lat),
  "W", sprintf("%05.2f", abs(dfAcousticData_int_Age$Rect_Lon)),
  "D", ifelse(dfAcousticData_int_Age$DepthMean_m <= 25, "0025",
              ifelse(dfAcousticData_int_Age$DepthMean_m >= 385, "0385", sprintf("%04d", dfAcousticData_int_Age$DepthMean_m))),
  sep=""
)

dfAcousticData_int_Age$LD_Area_ID = paste(
  dfAcousticData_int_Age$Area,
  "D", ifelse(dfAcousticData_int_Age$DepthMean_m <= 25, "0025",
              ifelse(dfAcousticData_int_Age$DepthMean_m >= 385, "0385", sprintf("%04d", dfAcousticData_int_Age$DepthMean_m))),
  sep=""
)

length(unique(dfAcousticData_int_Age$Cell_ID))
length(unique(dfAcousticData_int_Age$LD_ID))
length(unique(dfAcousticData_int_Age$LD_Area_ID))

# Add empty variables for ages
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age1 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age2 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age3 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age4 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age5 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age6 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age7 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age8 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age9 <- NA
dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age10 <- NA

# Set working directory
setwd("C:/Giacomo Gardella/Length distribution - Chris/New files/By area/results_A1")
# List all .xlsx files in the directory
file_list_LD_Area <- list.files(pattern = "*.xlsx")

# Set working directory
setwd("C:/Giacomo Gardella/Length distribution - Chris/New files/By cell/results_A2")
# List all .xlsx files in the directory
file_list_LD_Cell <- list.files(pattern = "*.xlsx")

#List of LD_IDs
LD_ID_list = substr(file_list_LD_Cell,38,54)


####### For loop to calculate N at age for each cell based on length distribution and age-length key

for (i in 1:length(dfAcousticData_int_Age$Cell_ID)) {
  
  # Define Cell ID
  cell_id_processed <- dfAcousticData_int_Age$Cell_ID[i]
  
  # Define LD ID
  LD_ID_processed <- dfAcousticData_int_Age$LD_ID[dfAcousticData_int_Age$Cell_ID == cell_id_processed]
  
  if (LD_ID_processed %in% LD_ID_list) { # Use LD output from Cell model
    
    # Find the file that ends with LD_ID_processed
    file_name <- file_list_LD_Cell[grepl(paste(LD_ID_processed, "\\.xlsx$", sep = ""), file_list_LD_Cell)]
    
    # Read the .xlsx file and correct the name
    file_data <- as.data.frame(read_excel(paste("C:/Giacomo Gardella/Length distribution - Chris/New files/By cell/results_A2/", file_name, sep = ""), sheet = 1))
    names(file_data) <- c("Length_cm", "prop")
    file_data$Length_cm <- as.numeric(file_data$Length_cm)
    
    # Add N of observations
    n_obs_data <- as.data.frame(read_excel(paste("C:/Giacomo Gardella/Length distribution - Chris/New files/By cell/results_A2/", file_name, sep = ""), sheet = 4))
    file_data$n_obs <- sum(n_obs_data$NumberOfObservations)/2 #Divide by 2 because the last row of the observations is the total number of observations
    
    # Add LD_ID
    file_data$LD_ID <- substr(file_name, 38, 54)
    
    if (file_data$n_obs[1] > 20) {
      # Create vector of densities
      length_n_vec <- dfAcousticData_int_Age$Density_N_nmi2x10m_cod[dfAcousticData_int_Age$Cell_ID == cell_id_processed] * file_data$prop[file_data$Length_cm <= 149]
      n_150_group <- sum(dfAcousticData_int_Age$Density_N_nmi2x10m_cod[dfAcousticData_int_Age$Cell_ID == cell_id_processed] * file_data$prop[file_data$Length_cm > 149])
      length_n_vec <- c(length_n_vec, n_150_group)
      
      # Multiply by ALK matrix
      n_at_age <- length_n_vec %*% dfPred_ALK_Norm_Mat
      
      # Overwrite results in dfAcousticData_int_Age
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age1[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[1]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age2[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[2]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age3[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[3]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age4[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[4]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age5[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[5]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age6[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[6]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age7[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[7]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age8[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[8]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age9[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[9]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age10[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age[10]
      
    } else {
      # Use output from Area model when n_obs <= 20
      Area_LD_ID_processed <- dfAcousticData_int_Age$LD_Area_ID[dfAcousticData_int_Age$Cell_ID == cell_id_processed]
      
      # Find the file that ends with Area_LD_ID_processed
      file_name_Area <- file_list_LD_Area[grepl(paste(Area_LD_ID_processed, "\\.xlsx$", sep = ""), file_list_LD_Area)]
      
      # Read the .xlsx file and correct the name
      file_data_Area <- as.data.frame(read_excel(paste("C:/Giacomo Gardella/Length distribution - Chris/New files/By area/results_A1/", file_name_Area, sep = ""), sheet = 1))
      names(file_data_Area) <- c("Length_cm", "prop")
      file_data_Area$Length_cm <- as.numeric(file_data_Area$Length_cm)
      
      # Create vector of densities
      length_n_vec_Area <- dfAcousticData_int_Age$Density_N_nmi2x10m_cod[dfAcousticData_int_Age$Cell_ID == cell_id_processed] * file_data_Area$prop[file_data_Area$Length_cm <= 149]
      n_150_group_Area <- sum(dfAcousticData_int_Age$Density_N_nmi2x10m_cod[dfAcousticData_int_Age$Cell_ID == cell_id_processed] * file_data_Area$prop[file_data_Area$Length_cm > 149])
      length_n_vec_Area <- c(length_n_vec_Area, n_150_group_Area)
      
      # Multiply by ALK matrix
      n_at_age_Area <- length_n_vec_Area %*% dfPred_ALK_Norm_Mat
      
      # Overwrite results in dfAcousticData_int_Age
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age1[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[1]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age2[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[2]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age3[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[3]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age4[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[4]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age5[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[5]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age6[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[6]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age7[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[7]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age8[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[8]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age9[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[9]
      dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age10[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[10]
    }
    
  } else {
    # Use output from Area model if LD_ID_processed is not in LD_ID_list
    Area_LD_ID_processed <- dfAcousticData_int_Age$LD_Area_ID[dfAcousticData_int_Age$Cell_ID == cell_id_processed]
    
    # Find the file that ends with Area_LD_ID_processed
    file_name_Area <- file_list_LD_Area[grepl(paste(Area_LD_ID_processed, "\\.xlsx$", sep = ""), file_list_LD_Area)]
    
    # Read the .xlsx file and correct the name
    file_data_Area <- as.data.frame(read_excel(paste("C:/Giacomo Gardella/Length distribution - Chris/New files/By area/results_A1/", file_name_Area, sep = ""), sheet = 1))
    names(file_data_Area) <- c("Length_cm", "prop")
    file_data_Area$Length_cm <- as.numeric(file_data_Area$Length_cm)
    
    # Create vector of densities
    length_n_vec_Area <- dfAcousticData_int_Age$Density_N_nmi2x10m_cod[dfAcousticData_int_Age$Cell_ID == cell_id_processed] * file_data_Area$prop[file_data_Area$Length_cm <= 149]
    n_150_group_Area <- sum(dfAcousticData_int_Age$Density_N_nmi2x10m_cod[dfAcousticData_int_Age$Cell_ID == cell_id_processed] * file_data_Area$prop[file_data_Area$Length_cm > 149])
    length_n_vec_Area <- c(length_n_vec_Area, n_150_group_Area)
    
    # Multiply by ALK matrix
    n_at_age_Area <- length_n_vec_Area %*% dfPred_ALK_Norm_Mat
    
    # Overwrite results in dfAcousticData_int_Age
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age1[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[1]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age2[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[2]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age3[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[3]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age4[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[4]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age5[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[5]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age6[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[6]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age7[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[7]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age8[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[8]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age9[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[9]
    dfAcousticData_int_Age$Density_N_nmi2x10m_cod_Age10[dfAcousticData_int_Age$Cell_ID == cell_id_processed] <- n_at_age_Area[10]
    
  }
}

head(dfAcousticData_int_Age)

# Replace NA with 0
names_age_col = c("Density_N_nmi2x10m_cod_Age1","Density_N_nmi2x10m_cod_Age2","Density_N_nmi2x10m_cod_Age3",
                  "Density_N_nmi2x10m_cod_Age4","Density_N_nmi2x10m_cod_Age5", "Density_N_nmi2x10m_cod_Age6","Density_N_nmi2x10m_cod_Age7",
                  "Density_N_nmi2x10m_cod_Age8","Density_N_nmi2x10m_cod_Age9","Density_N_nmi2x10m_cod_Age10")

dfAcousticData_int_Age[, names_age_col] <- lapply(dfAcousticData_int_Age[, names_age_col], function(x) replace(x, is.na(x), 0))


# Count cod in the fjord by age

#Aggregate by rectangle
dfAcousticData_AreaDensity_Age <- aggregate(dfAcousticData_int_Age[, c("Density_N_nmi2x10m_cod_Age1",
                                                                       "Density_N_nmi2x10m_cod_Age2",
                                                                       "Density_N_nmi2x10m_cod_Age3",
                                                                       "Density_N_nmi2x10m_cod_Age4",
                                                                       "Density_N_nmi2x10m_cod_Age5",
                                                                       "Density_N_nmi2x10m_cod_Age6",
                                                                       "Density_N_nmi2x10m_cod_Age7",
                                                                       "Density_N_nmi2x10m_cod_Age8",
                                                                       "Density_N_nmi2x10m_cod_Age9",
                                                                       "Density_N_nmi2x10m_cod_Age10")],
                                            by = list(dfAcousticData_int_Age$Rect_Lat, dfAcousticData_int_Age$Rect_Lon),FUN = sum)

# Quality assurance step: make map by rectangle and make sure its the same as fromthe EchoIntegrationScript

mean_AreaDensity_fishnmi2_age <- sapply(dfAcousticData_AreaDensity_Age[, c("Density_N_nmi2x10m_cod_Age1",
                                                                           "Density_N_nmi2x10m_cod_Age2",
                                                                           "Density_N_nmi2x10m_cod_Age3",
                                                                           "Density_N_nmi2x10m_cod_Age4",
                                                                           "Density_N_nmi2x10m_cod_Age5",
                                                                           "Density_N_nmi2x10m_cod_Age6",
                                                                           "Density_N_nmi2x10m_cod_Age7",
                                                                           "Density_N_nmi2x10m_cod_Age8",
                                                                           "Density_N_nmi2x10m_cod_Age9",
                                                                           "Density_N_nmi2x10m_cod_Age10")], mean)
mean_AreaDensity_fishnmi2_age

#Scale up to entire fjord
setwd("C:/Giacomo Gardella/Comparisons lsss")
dfBathymetryData = read.csv("Godthåbsfjord_Ameralik_strata_20200414.csv", sep=";")
numArea_nm2 = (sum(dfBathymetryData$areakm2)*(1000^2))/(1852^2)

Cod_total <- mean_AreaDensity_fishnmi2_age*numArea_nm2
sum(Cod_total) #13995992

WISC_total <- Cod_total*0.55
WISC_total

# Before cod stock split
Cod_total_beforeSplit <- Cod_total
names(Cod_total_beforeSplit) = c("N_Age1","N_Age2","N_Age3","N_Age4","N_Age5","N_Age6","N_Age7","N_Age8","N_Age9","N_Age10")
Cod_total_beforeSplit
sum(Cod_total_beforeSplit)


# Vertical distribution
head(dfAcousticData_int_Age)

dfAcousticData_int_Age10_Layer = aggregate(Density_N_nmi2x10m_cod_Age10~DepthMean_m, dfAcousticData_int_Age, mean)
dfAcousticData_int_Age10_Layer$Density_N_per_m3 = dfAcousticData_int_Age10_Layer$Density_N_nmi2x10m_cod_Age10/((1852^2)*10)

smoother_int_Age10 <-sreg(I(rev(dfAcousticData_int_Age10_Layer$DepthMean_m)),dfAcousticData_int_Age10_Layer$Density_N_per_m3, lambda = 0.2)
predVD_int_Age10 <- predict(smoother_int_Age10)


# Vertical profile
plot(dfAcousticData_int_Age10_Layer$Density_N_per_m3,dfAcousticData_int_Age10_Layer$DepthMean_m, col="white",
     ylim=rev(range(dfAcousticData_int_Age10_Layer$DepthMean_m)),
     xlab=expression("Volume density of Age 10 cod (cod"~m^3~")"), ylab="Depth (m)")

lines(predVD_int_Age10,dfAcousticData_int_Age10_Layer$DepthMean_m, lwd=6, col="purple")


# Maps for ages

# Age 1
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age1)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 21000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 1", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=1434782", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age1.png", bg="white", width=5.3, height=4, dpi=500)

# Age 2
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age2)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 27000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 2", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=2622436", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age2.png", bg="white", width=5.3, height=4, dpi=500)

# Age 3
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age3)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 24000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 3", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=1800033", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age3.png", bg="white", width=5.3, height=4, dpi=500)

# Age 4
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age4)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 24000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 4", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=1520468", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age4.png", bg="white", width=5.3, height=4, dpi=500)

# Age 5
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age5)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 24000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 5", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=1081682", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age5.png", bg="white", width=5.3, height=4, dpi=500)

# Age 6
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age6)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 24000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 6", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=877388", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age6.png", bg="white", width=5.3, height=4, dpi=500)

# Age 7
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age7)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 21000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 7", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=826533", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age7.png", bg="white", width=5.3, height=4, dpi=500)

# Age 8
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age8)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 21000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 8", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=941139", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age8.png", bg="white", width=5.3, height=4, dpi=500)

# Age 9
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age9)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 41000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 9", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=826533", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age9.png", bg="white", width=5.3, height=4, dpi=500)


# Age 10
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age10)) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 39000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 10", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=1324067", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age10.png", bg="white", width=5.3, height=4, dpi=500)

# Age 10 - proportion
ggplot() +
  geom_tile(data = expand.grid(Rect_Lon = unique(dfAcousticData_AreaDensity_Age$Group.2),
                               Rect_Lat = unique(dfAcousticData_AreaDensity_Age$Group.1)),
            aes(x = Rect_Lon, y = Rect_Lat),
            fill = "lightblue1") +
  geom_tile(data = dfAcousticData_AreaDensity_Age, aes(x = Group.2, y = Group.1, fill = Density_N_nmi2x10m_cod_Age10/) +
  scale_fill_gradient(low = "yellow", high = "red",
                      limits = c(0, 39000),  
                      breaks = seq(0, 21000, by = 5000)) +
  geom_sf(data = nuuk_countries, color = "black") +
  coord_sf(xlim = c(nuuk_bbox$xmin, nuuk_bbox$xmax), ylim = c(nuuk_bbox$ymin, nuuk_bbox$ymax), expand = FALSE) +
  labs(title="Age 10", x = "Longitude", y = "Latitude", fill=paste("Areal cod density
( cod",expression(nmi^2),")",sep=" "))+
  scale_x_continuous(breaks = seq(nuuk_bbox$xmin, nuuk_bbox$xmax, by = 0.3))+
  annotate("text", x = -50.85, y = 64.05, label = "Tot=1324067", size = 4, hjust = 0) + 
  theme_classic()
#ggsave(filename="density maps integration Age10.png", bg="white", width=5.3, height=4, dpi=500)
