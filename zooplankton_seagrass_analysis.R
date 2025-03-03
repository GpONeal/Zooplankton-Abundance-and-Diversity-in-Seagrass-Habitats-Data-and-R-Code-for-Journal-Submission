###ZOOPLANKTON TUBE TRAP ANALYSIS



# This R file contains the data manipulation and analysis for zooplankton abundances, diversity, and seasonality between seagrass species.

# Couple steps before proceeding:
# Step (1) Enable Soft wrapping - It will save you the headache
# Step (2) Ctrl + F and replace ?.. (Question mark ..) with ï.. (i with diaeresis) - importing the data may cause it to change some column names to start with ï.. and for some reason saving MAY change it to ?...

#Lets get started!
library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
#Import the file tube_trap_data

tube_trap_data <- read.csv(" ")

#Make a copy of the data 

tube_data <- tube_trap_data

#Make a summary table for taxon seen in each sample by seagrass#
trap_stats <- tube_data %>%
  group_by(Seagrass) %>%
  summarise(across(where(is.numeric), list(
    mean = ~mean(. , na.rm = TRUE),
    sd = ~sd(. , na.rm = TRUE),
    min = ~min(. , na.rm = TRUE),
    max = ~max(. , na.rm = TRUE))))

# View the summary statistics

print(trap_stats, n = Inf, width = Inf)

#Calculate Shannon Diversity
calculate_shannon <- function(row) {
  # Remove zeros
  row_non_zero <- row[row > 0]
  # Calculate proportions
  proportions <- row_non_zero / sum(row_non_zero)
  # Calculate Shannon Diversity Index
  shannon_index <- -sum(proportions * log(proportions))
  return(shannon_index)
}

# Apply the Shannon Diversity Index function to each row and create a new column
trap_shannon <- tube_data %>%
  rowwise() %>%
  mutate(Shannon_Index = calculate_shannon(c_across(8:22)))

#Summarize Shhannon indexies
shannon_stats <- trap_shannon %>%
  group_by(Seagrass,) %>%
  summarise(
    Mean = mean(Shannon_Index, na.rm = TRUE),
    SD = sd(Shannon_Index, na.rm = TRUE),
    Min = min(Shannon_Index, na.rm = TRUE),
    Max = max(Shannon_Index, na.rm = TRUE),
  )

print(shannon_stats)

#ANOVA: Seagrass

trap_anova_result <- aov(Shannon_Index ~ Seagrass, data = trap_shannon)

summary(trap_anova_result)

summary(trap_anova_result)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Seagrass     2  0.276  0.1382   1.117  0.334
#Residuals   58  7.172  0.1236

#ANOVA: Site

trap_anova_result <- aov(Shannon_Index ~ Site, data = trap_shannon)

summary(trap_anova_result)

#            Df Sum Sq Mean Sq F value Pr(>F)  
#ï..Site      4  1.180  0.2950   2.636 0.0434 *
#Residuals   56  6.268  0.1119 

#PostHoc shows only Brewers is significant. This makes sense because it only has Halophila stipulacea within it. Run two-way anova with seagrass type.


#Two-way ANOVA: Seagrass & Site


trap_2anova_result <- aov(Shannon_Index ~ Seagrass * Site, data = trap_shannon)

summary(trap_2anova_result)
#                  Df Sum Sq Mean Sq F value Pr(>F) 
#Seagrass          2  0.276  0.1382   1.297 0.2829  
#ï..Site           4  1.047  0.2617   2.456 0.0582 
#Seagrass:ï..Site  6  1.010  0.1683   1.579 0.1738  
#Residuals        48  5.115  0.1066          

#Not Significant

#Plot the Shannon indexies
ggplot(trap_shannon, aes(x = Seagrass, y = Shannon_Index)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", color = "black") +  # Create the boxplot
  labs(title = "", 
       x = "Seagrass Species", 
       y = "Shannon Diversity Index") +
  
  scale_x_discrete(labels = c("HS" = "H. stipulacea", 
                              "SF" = "S. filiforme", 
                              "TT" = "T. testudinum")) +
  theme_minimal()

#EDITED COLORS

ggplot(trap_shannon, aes(x = Seagrass, y = Shannon_Index, fill = Seagrass)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Create the boxplot
  labs(title = "", 
       x = "Seagrass Species", 
       y = "Shannon Diversity Index") +
  
  scale_x_discrete(labels = c("HS" = expression(italic("H. stipulacea")), 
                              "SF" =expression(italic ("S. filiforme")), 
                              "TT" =expression(italic("T. testudinum")))) +
  
  scale_fill_manual(values = c("HS" = "orange", 
                               "SF" = "lightgreen", 
                               "TT" = "deepskyblue")) +  
  
  theme_minimal() 
  

  #Calculate average abundance of taxa per sample
average_abundances_by_seagrass <- tube_data %>%
  group_by(Seagrass) %>%  
  summarise(across(Harpacticoid:Fish_Larvae, mean, na.rm = TRUE))

# Pivot the data
long_data <- average_abundances_by_seagrass %>%
  pivot_longer(cols = Harpacticoid:Fish_Larvae, 
               names_to = "Taxon", 
               values_to = "Average_Abundance")

# Calculate mean and standard error (SE)
abundance_stats <- tube_data %>%
  group_by(Seagrass) %>%
  summarise(across(Harpacticoid:Fish_Larvae, 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
                   .names = "{col}_{fn}"))

# Pivot the data for mean and standard error
long_data <- abundance_stats %>%
  pivot_longer(cols = -Seagrass, 
               names_to = c("Taxon", ".value"), 
               names_sep = "_") %>%
  rename(Average_Abundance = mean, Std_Error = se)

# The data is now in a tidy format with mean and standard error for each taxon and seagrass




# Calculate mean and standard error (SE)
abundance_stats <- tube_data %>%
  group_by(Seagrass) %>%
  summarise(across(Harpacticoid:Fish_Larvae, 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
                   .names = "{col}_{fn}"))

# Pivot the data for mean and standard error
long_data <- abundance_stats %>%
  pivot_longer(cols = -Seagrass, 
               names_to = c("Taxon", ".value"), 
               names_sep = "_") %>%
  rename(Average_Abundance = mean, Std_Error = se)

# Fix taxon names
long_data <- long_data %>%
  mutate(Taxon = gsub("[_\\.]", " ", Taxon)) %>%
  mutate(Taxon = gsub("Fish eggs", "Fish Eggs", Taxon)) %>%
  mutate(Taxon = gsub("Ascidicea Larva", "Ascidian Larvae", Taxon))

# Combine specific taxa into "Other" category
long_data <- long_data %>%
  mutate(Taxon = case_when(
    Taxon %in% c("Cumacea", "Ascidian Larvae", "Isopoda", "Fish Larvae") ~ "Other",
    TRUE ~ Taxon
  ))

# Calculate overall mean abundance per taxon for ordering
taxon_order <- long_data %>%
  group_by(Taxon) %>%
  summarise(Overall_Mean_Abundance = mean(Average_Abundance, na.rm = TRUE)) %>%
  arrange(desc(Overall_Mean_Abundance)) %>%
  pull(Taxon)

# Reorder taxa, placing "Other" at the end
long_data <- long_data %>%
  mutate(Taxon = fct_relevel(factor(Taxon, levels = taxon_order), "Other", after = Inf))

long_data <- long_data %>%
  filter(!is.na(Average_Abundance) & !is.na(Std_Error))

ggplot(long_data, aes(x = Taxon, y = Average_Abundance, fill = Seagrass)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.80, color = "black") +
  geom_errorbar(aes(ymin = Average_Abundance - Std_Error, ymax = Average_Abundance + Std_Error), 
                position = position_dodge(width = 0.9), width = 0.25, color = "black") +
  scale_fill_manual(values = c("HS" = "orange", "SF" = "lightgreen", "TT" = "deepskyblue")) +  
  labs(x = "Taxon", y = "Abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)
  )
















# Edit Names
long_data <- long_data %>%
  mutate(Taxon = gsub("[_\\.]", " ", Taxon)) %>% 
  mutate(Taxon = gsub("Fish eggs", "Fish Eggs", Taxon)) %>% 
  mutate(Taxon = gsub("Ascidicea Larva", "Ascidian Larvae", Taxon))
  

long_data <- long_data %>%
  mutate(Taxon = case_when(
    Taxon %in% c("Cumacea", "Ascidian Larvae", "Isopoda", "Fish Larvae") ~ "Other",
    TRUE ~ Taxon
  ))
long_data <- long_data %>%
  mutate(Taxon = fct_reorder(Taxon, Average_Abundance, .desc = TRUE))


#Plot Average Abundance per Seagrass
long_data <- long_data %>%
  arrange(desc(Average_Abundance)) %>%
  mutate(Taxon = factor(Taxon, levels = unique(Taxon)))



long_data <- long_data %>%
  mutate(Taxon = case_when(
    Taxon %in% c("Cumacea", "Ascidian Larvae", "Isopoda", "Fish Larvae") ~ "Other",
    TRUE ~ Taxon
  )) %>%
  group_by(Taxon) %>%
  summarise(Average_Abundance = mean(Average_Abundance, na.rm = TRUE),
            Std_Error = mean(se, na.rm = TRUE)) %>%  # Aggregate standard error
  ungroup() %>%
  mutate(Taxon = fct_reorder(Taxon, mean, .desc = TRUE))

# Pivot the data for mean and standard deviation
long_data <- abundance_stats %>%
  pivot_longer(cols = -Seagrass, 
               names_to = c("Taxon", ".value"), 
               names_sep = "_")

# Edit Names
long_data <- long_data %>%
  mutate(Taxon = gsub("[_\\.]", " ", Taxon)) %>% 
  mutate(Taxon = gsub("Fish eggs", "Fish Eggs", Taxon)) %>% 
  mutate(Taxon = gsub("Ascidicea Larva", "Ascidian Larvae", Taxon))

# Group certain taxa as 'Other' and reorder factors
long_data <- long_data %>%
  mutate(Taxon = case_when(
    Taxon %in% c("Cumacea", "Ascidian Larvae", "Isopoda", "Fish Larvae") ~ "Other",
    TRUE ~ Taxon
  ))  
  ungroup() %>%
  mutate(Taxon = fct_reorder(Taxon, Average_Abundance, .desc = TRUE))

# Now you can plot using the standard deviation in long_data
  ggplot(long_data, aes(x = Taxon, y = Average_Abundance, fill = Seagrass)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.80, color = "black") +
    geom_errorbar(aes(ymin = Average_Abundance - Std_Error, ymax = Average_Abundance + Std_Error), 
                  position = position_dodge(width = 0.9), width = 0.25, color = "black") +
    scale_fill_manual(values = c("HS" = "orange", "SF" = "lightgreen", "TT" = "deepskyblue")) +  
    labs(x = "Taxon", y = "Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  




#graph it
ggplot(long_data, aes(x = Taxon, y = mean, fill = Seagrass)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.80) +  # Adjust position and width for spacing
  scale_fill_manual(values = c("HS" = "orange", 
                               "SF" = "lightgreen", 
                               "TT" = "deepskyblue")) +  
  labs( 
       x = "Taxon", 
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))  



ggplot(long_data, aes(x = Taxon, y = Average_Abundance, fill = Seagrass)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.80, color = "black") +  # Adjust width and add border around bars
  geom_errorbar(aes(ymin = Average_Abundance - Std_Dev, ymax = Average_Abundance + Std_Dev), 
                position = position_dodge(width = 0.9), width = 0.25, color = "black") +  # Add error bars for std
  scale_fill_manual(values = c("HS" = "orange", 
                               "SF" = "lightgreen", 
                               "TT" = "deepskyblue")) +  
  labs( 
    x = "Taxon", 
    y = "Abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase x-axis label size
    panel.border = element_rect(color = "black", fill = NA, size = 1.5))



amp_anova_result <- aov(Amphipod ~ Seagrass, data = tube_data)
dec_anova_result <- aov(Decapod ~ Seagrass, data = tube_data)
iso_anova_result <- aov(Isopoda ~ Seagrass, data = tube_data)
har_anova_result <- aov(Harpacticoid ~ Seagrass, data = tube_data)

cal_anova_result <- aov(Calanoida ~ Seagrass, data = tube_data)
cyc_anova_result <- aov(Cyclopoida ~ Seagrass, data = tube_data)
mys_anova_result <- aov(Mysida ~ Seagrass, data = tube_data)
cum_anova_result <- aov(Cumacea ~ Seagrass, data = tube_data)
tan_anova_result <- aov(Tanaidacea ~ Seagrass, data = tube_data)

app_anova_result <- aov(Appendicularians ~ Seagrass, data = tube_data)
nem_anova_result <- aov(Nematoda ~ Seagrass, data = tube_data)
egg_anova_result <- aov(Fish.eggs ~ Seagrass, data = tube_data)
lar_anova_result <- aov(Fish_Larvae ~ Seagrass, data = tube_data)
ost_anova_result <- aov(Ostracod ~ Seagrass, data = tube_data)
asc_anova_result <- aov(Ascidicea_Larva ~ Seagrass, data = tube_data)

summary(amp_anova_result)
summary(iso_anova_result) # Isopods Significant
#F = 3.25 df= 2, p-value=0.0459
summary(har_anova_result) #Harpacticoid Significant
#F=4.877, df=2, p-value=0.011
summary(cal_anova_result)
summary(cyc_anova_result)
summary(mys_anova_result)
summary(tan_anova_result)
summary(app_anova_result) #Appendicularians Significant
#F=3.688, df = 2, p-value=0.0311
summary(nem_anova_result)
summary(egg_anova_result)
summary(lar_anova_result)
summary(ost_anova_result)
summary(asc_anova_result)#Ascideans Significant
#F=3.322, df =2, p-value=0.043

iso_tukey <- TukeyHSD(iso_anova_result)
print(iso_tukey)

#SF-HS p-value=0.0422061

har_tukey <- TukeyHSD(har_anova_result)
print(har_tukey)

#SF-HS p-value=0.0149888

app_tukey <- TukeyHSD(app_anova_result)
print(app_tukey)
#SF-HS p-value = 0.0332832
asc_tukey <- TukeyHSD(asc_anova_result)
print(asc_tukey)
#SF-HS p-value=0.0332832


#Seasonal Abundance

samp_anova_result <- aov(Amphipod ~ Season, data = tube_data)
sdec_anova_result <- aov(Decapod ~ Season, data = tube_data)
siso_anova_result <- aov(Isopoda ~ Season, data = tube_data)
shar_anova_result <- aov(Harpacticoid ~ Season, data = tube_data)

scal_anova_result <- aov(Calanoida ~ Season, data = tube_data)
scyc_anova_result <- aov(Cyclopoida ~ Season, data = tube_data)
smys_anova_result <- aov(Mysida ~ Season, data = tube_data)
scum_anova_result <- aov(Cumacea ~ Season, data = tube_data)
stan_anova_result <- aov(Tanaidacea ~ Season, data = tube_data)

sapp_anova_result <- aov(Appendicularians ~ Season, data = tube_data)
snem_anova_result <- aov(Nematoda ~ Season, data = tube_data)
segg_anova_result <- aov(Fish.eggs ~ Season, data = tube_data)
slar_anova_result <- aov(Fish_Larvae ~ Season, data = tube_data)
sost_anova_result <- aov(Ostracod ~ Season, data = tube_data)
sacs_anova_result <- aov(Ascidicea_Larva ~ Season, data = tube_data)


summary(samp_anova_result)
summary(siso_anova_result)



summary(shar_anova_result) 
summary(scal_anova_result)
summary(cycs_anova_result)
summary(scum_anova_result)
summary(stan_anova_result)
summary(sapp_anova_result) 
summary(snem_anova_result) 
summary(segg_anova_result)
summary(slar_anova_result)


summary(sost_anova_result)
summary(sacs_anova_result)

### NMDS Zooplankton in Seagrass Species

library(ggplot2)
library(vegan)
library(tidyverse)

#Import the data set zooplankton_nmds

zooplankton_nmds <- read.csv(" ")

#make Data only numeric, remove non-numeric data

plank_nmd <- zooplankton_nmds

plank_nmd$Naupaulis.stage.III = NULL
plank_nmd$Stomatopod = NULL
plank_nmd$Naupaulis.stage.I. = NULL
plank_nmd$Chaetognath = NULL
plank_nmd$Polychaete = NULL
plank_nmd$Bivalve = NULL
plank_nmd$Chiton =NULL
plank_nmd$Snail.Larvae = NULL

print(names(plank_nmd))

# Change column names
names(plank_nmd) <- c("ï..Sample.ID", "Date", "Site", "Seagrass", "Amphipoda", "Cumacea", "Tanaidacea","Mysida", "Isopoda","Decapoda","Calanoida", "Cyclopoida","Harpacticoida", "Ostracoda", "Appendicularia", "Ascidian Larvae", "Fish Larvae", "Fish Eggs", "Nematoda")



non_num_plank <- plank_nmd

non_num_plank$Date = NULL
non_num_plank$ï..Sample.ID = NULL
non_num_plank$Site = NULL
non_num_plank$Seagrass = NULL


view(non_num_plank)
grass_NMDS <- metaMDS(non_num_plank, k = 2, trymax = 100) 



grass_NMDS

#Stress: 0.2533287

stressplot(grass_NMDS)

#Non-metric fit R^2 = 0.936
#Linear Fit = 0.686

#Make a plotfram to pull the non numeric data
plotframe<-plank_nmd%>%
  select(ï..Sample.ID, Site, Seagrass, Date)

##pulling the species locations out of the NMDS (x and y)
fam <- as.data.frame(grass_NMDS$species)


# this line is pulling the species names and linking them to their location
fam$sp <- rownames(grass_NMDS$species)


#this is pulling the points for each individual transect out of the nmds
seagrass_type <- as.data.frame(grass_NMDS$points)

#we are pulling the seagrass categories out of our plotframe that we made earlier so we can make the 95% confidence intervals for our categories depth/percent halophila
seagrass_type$Seagrass <- plotframe[, 3]

#need to keep this code because you cant pull the points in the next line unless if this ordiplot exists
ordiplot(grass_NMDS, type = "none") 

#this is storing the info to make the 95% confidence intervals as an object from the ordiplot
finallord<-ordiellipse(grass_NMDS, seagrass_type$Seagrass, conf = 0.95, label = TRUE)

#hidden vegan function that calculates ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#convert zone to a factor
seagrass_type$Seagrass <- factor(seagrass_type$Seagrass)

#create a dataframe of points that will draw the ellipses
final_df_ell <- data.frame()

for(g in levels(seagrass_type$Seagrass)){
  final_df_ell <- rbind(final_df_ell, cbind(as.data.frame(with(seagrass_type[seagrass_type$Seagrass==g,],
                                                               veganCovEllipse(finallord[[g]]$cov,finallord[[g]]$center,finallord[[g]]$scale)))
                                            ,zone=g))
}

#plot
finalnMDSplot<-ggplot(seagrass_type) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(Seagrass)), size = 3) + geom_path(data = final_df_ell, aes(x = NMDS1, y = NMDS2, group = zone, color=zone), linewidth = 1, linetype = 1)+
  theme_classic() +
  labs(color = "Seagrass", fill = "Seagrass") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(panel.border = element_rect(fill = NA),
        legend.position = "bottom")+ geom_text(data = fam, aes(x = MDS1, y = MDS2, label = sp), alpha = 0.3, size - 6)+scale_fill_manual(name = "Category", labels = c("Halophila stipulacea", "Thalassia testudinum","Syrigodium filiforme"),values = c("orange","deepskyblue","lightgreen"))+scale_color_manual(name = "Category", labels = c("Halophila stipulacea", "Thalassia testudinum","Syrigodium filiforme"), values = c("orange","deepskyblue","lightgreen"))+annotate(geom="text", x=1.5, y=1.2, label="stress=0.234")

finalnMDSplot <- ggplot(seagrass_type) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(Seagrass)), size = 6) + 
  geom_path(data = final_df_ell, aes(x = NMDS1, y = NMDS2, group = zone, color=zone), linewidth = 1, linetype = 1) +
  theme_classic() +
  labs(color = "Seagrass", fill = "Seagrass") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(panel.border = element_rect(fill = NA),
        legend.position = "bottom") +
  geom_text(data = fam, aes(x = MDS1, y = MDS2, label = sp), alpha = 0.3, size = 6) +  # Increase the size here
  scale_fill_manual(name = "Category", labels = c("Halophila stipulacea", "Thalassia testudinum", "Syrigodium filiforme"), values = c("orange","deepskyblue","lightgreen")) +
  scale_color_manual(name = "Category", labels = c("Halophila stipulacea", "Thalassia testudinum", "Syrigodium filiforme"), values = c("orange","deepskyblue","lightgreen")) +
  annotate(geom="text", x=1.5, y=1.2, label="stress=0.253")

finalnMDSplot <- ggplot(seagrass_type) +
  geom_point(aes(x = MDS1, y = MDS2, color = as.factor(Seagrass)), size = 2) + 
  geom_path(data = final_df_ell, aes(x = NMDS1, y = NMDS2, group = zone, color = zone), linewidth = 1, linetype = 1) +
  theme_classic() +
  labs(color = "Seagrass", fill = "Seagrass") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "bottom",
    axis.text = element_text(size = 15),  # Adjust text size for axis
    legend.text = element_text(size = 15),  # Adjust text size for legend
    axis.title = element_text(size = 17)  # Adjust title size for axis
  ) +
  geom_text(data = fam, aes(x = MDS1, y = MDS2, label = sp), color = "black", alpha = 1, size = 3) +  # Set zooplankton taxon labels to black and fully opaque
  scale_fill_manual(
    name = "Category",
    labels = c(expression(italic("Halophila stipulacea")), 
               expression(italic("Thalassia testudinum")), 
               expression(italic("Syringodium filiforme"))),
    values = c("orange", "deepskyblue", "lightgreen")
  ) +
  scale_color_manual(
    name = "Category",
    labels = c(expression(italic("Halophila stipulacea")), 
               expression(italic("Thalassia testudinum")), 
               expression(italic("Syringodium filiforme"))),
    values = c("orange", "deepskyblue", "lightgreen")
  ) +
  annotate(geom = "text", x = 1.11, y = 1.2, label = "stress = 0.2533")

finalnMDSplot
