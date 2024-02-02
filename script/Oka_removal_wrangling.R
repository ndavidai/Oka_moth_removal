### Data cleaning ###

#### egg data ####
#loading data 
removal_counts <- read.csv("input/Oka_removal_counts.csv")

# quick visualizations
summary(removal_counts)
str(removal_counts)

library(janitor) #janitor cleans up column names.It removes all unique characters and replaces spaces with _.
#piping through `dplyr`
removal_counts_1 <- removal_counts %>%
  clean_names() #Cleans names of an object (usually a data.frame)

## there are 2 columns not needed (Notes, Entered_by)- let's remove them ##
library(dplyr)
removal_counts_1 <-removal_counts_1 %>% select(1:8)

# if any column names need replacing
colnames(removal_counts_1)[colnames(removal_counts_1)=="dbh_diameter_cm"] <- "dbh_cm"
colnames(removal_counts_1)[colnames(removal_counts_1)=="species_genus"] <- "tree_species"

# quick visualizations
summary(removal_counts_1)
str(removal_counts_1)

# change species, removal height, quant removed to 'factor'; change DBH, dist from site to 'numerical'
removal_counts_1$dbh_cm <- as.numeric(removal_counts_1$dbh_cm)
removal_counts_1$distance_from_site_m <- as.numeric(removal_counts_1$distance_from_site_m)
removal_counts_1$tree_species <- as.factor(removal_counts_1$tree_species)
removal_counts_1$removal_height_cm <- as.factor(removal_counts_1$removal_height_cm)
removal_counts_1$quantity_removed <- as.factor(removal_counts_1$quantity_removed)

str(removal_counts_1)

unique(removal_counts_1$dbh_cm)

# replace [?] with 'NA'
removal_counts_1 <- replace(removal_counts_1, removal_counts_1=='?', 'NA')

unique(removal_counts_1$tree_species)

unique(removal_counts_1$removal_height_cm)
unique(removal_counts_1$quantity_removed)
# replace 'undetected' with 'none'
removal_counts_1 <- replace(removal_counts_1, removal_counts_1=='undetected', 'none')

unique(removal_counts_1$distance_from_site_m)


unique(removal_counts_1$dbh_cm)
unique(removal_counts_1$tree_species)
unique(removal_counts_1$removal_height_cm)
unique(removal_counts_1$quantity_removed)
unique(removal_counts_1$distance_from_site_m)


#visualize data
#library(esquisse)
#esquisser(egg_mass_counts_1)

library(esquisse)
esquisser(removal_counts_1)

library(ggplot2)


write.csv(removal_counts_1, "Oka_removal_counts_1.csv")






### Extra code from Riikka's original script ###
egg_counts_1 <- replace(egg_counts_1, egg_counts_1==',', NA)
# replace "no entry" with NA
egg_counts_1 <- replace(egg_counts_1, egg_counts_1=='no entry', NA)
# replace "none" with NA ? maybe not as this seems categorical
egg_counts_1 <- replace(egg_counts_1, egg_counts_1=='none', 0)
unique(egg_counts_1$Egg.Count)
unique(egg_counts_1$Caterpillar.Count)

# if any column names need replacing
colnames(egg_counts_1)[colnames(egg_counts_1)=="Site.Name"] <- "Site_Name"

#### caterpillar data ####
#loading data 
caterp_counts <- read.csv("MSB Caterpillar Counts.csv")

# quick visualizations
summary(caterp_counts)
str(caterp_counts)
unique(caterp_counts$Egg.Count)
unique(caterp_counts$Site.Name)

#library(janitor) #janitor cleans up column names.It removes all unique characters and replaces spaces with _.
#piping through `dplyr`
caterp_counts <- caterp_counts %>%
  clean_names() #Cleans names of an object (usually a data.frame) 

## let's remove notes column ##
library(dplyr)
caterp_counts_1 <-caterp_counts %>% select(-13)
# let's check for inconsistencies in names etc
unique(caterp_counts_1$dbh_cm)
unique(caterp_counts_1$ )

## there are some rows that are completely empty- we want to replace those with NA##
# replace empty cells with NA
caterp_counts_1 <- replace(caterp_counts_1, caterp_counts_1=='', NA)
# replace "no entry" with NA (or with something else?)
caterp_counts_1 <- replace(caterp_counts_1, caterp_counts_1=='no entry', NA)
# replace "none" with NA ? maybe not as this seems categorical
caterp_counts_1 <- replace(caterp_counts_1, caterp_counts_1=='none', 0)

unique(caterp_counts_1$dbh_cm) #see "no entry "
# replace "no entry " with NA (or with something else?)
caterp_counts_1 <- replace(caterp_counts_1, caterp_counts_1=='no entry ', NA)
unique(caterp_counts_1$tree_id) #you may want to change these to separate common name and scientific name columns

## replace space between words with underscore
library(tidyverse)
caterp_counts_1 <- caterp_counts_1 %>%
  mutate(tree_id = str_replace(tree_id, " ", "_"))

# replace "none " with "none"
caterp_counts_1 <- replace(caterp_counts_1, caterp_counts_1=='none ', "none")

summary(caterp_counts_1)
str(caterp_counts_1)

## remove % and replace with nothing
caterp_counts_1 <- caterp_counts_1 %>%
  mutate(humidity = str_replace(humidity, "%", ""))

# if any column names need replacing
colnames(caterp_counts_1)[colnames(caterp_counts_1)=="humidity"] <- "humidity_percent"

## some numeric columns are character here- let's change them
caterp_counts_1$dbh_cm <- as.numeric(caterp_counts_1$dbh_cm)
caterp_counts_1$humidity_percent <- as.numeric(caterp_counts_1$humidity_percent)

summary(caterp_counts_1)
str(caterp_counts_1)

egg_mass_counts_1 %>%
  filter(!is.na(temp_c)) %>%
  filter(!is.na(humidity_percent)) %>%
  filter(!is.na(in_center)) %>%
  filter(!is.na(count_height)) %>%
  filter(!is.na(egg_count)) %>%
  filter(!is.na(caterpillar_count)) %>%
  filter(!is.na(entered_by)) %>%
  ggplot() +
  aes(x = dbh_cm, y = egg_count) +
  geom_boxplot(fill = "#2FD6BE") +
  scale_x_continuous(trans = "log") +
  coord_flip() +
  theme_minimal()
