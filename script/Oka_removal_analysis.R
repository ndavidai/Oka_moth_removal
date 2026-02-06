
#loading data 
removal_counts <- read.csv("input/Oka_removal_counts_clean.csv")
hatch_rates <- read.csv("input/Oka_hatch_rates_clean.csv")


library(janitor)                       
library(dplyr)
library(ggplot2)
library(esquisse)
library(ordinal)
library(tidyr)
library(ordinal)
library(DHARMa)
library(Matrix)
library(lme4)
library(glmmTMB)
library(forcats)

#install.packages("TMB")
#install.packages("glmmTMB", type = "source")



### removal counts analyses -------------------------------------------------

summary(removal_counts)
str(removal_counts)

# replace 'NA' with 'zero'
# add "zero" as a valid factor level
levels(removal_counts$quantity_removed) <-
  c(levels(removal_counts$quantity_removed), "zero")

# replace NA with "zero"
removal_counts$quantity_removed[
  is.na(removal_counts$quantity_removed)] <- "zero"


##visual explorations
ggplot(removal_counts, aes(distance_from_site_m, quantity_removed)) +
  geom_boxplot() +
  facet_wrap(~ tree_species)

ggplot(removal_counts, aes(removal_height_cm, quantity_removed)) +
  geom_boxplot()


##separate site ID into a separate variable that does not include tree #

removal_counts$Site_ID <- factor(
  sub("-.*$", "", removal_counts$id)
)

##check data structure
table(removal_counts$quantity_removed, removal_counts$tree_species)
table(removal_counts$quantity_removed, removal_counts$removal_height_cm)
table(removal_counts$quantity_removed, removal_counts$dbh_cm)
table(removal_counts$quantity_removed, removal_counts$distance_from_site_m)
table(removal_counts$quantity_removed, removal_counts$Site_ID)


### Ordinal mixed-effects model ---------------------------------------------

removal_counts$quantity_removed <- factor(
  removal_counts$quantity_removed,
  levels = c("zero", "low", "medium", "high"),
  ordered = TRUE
)


#convert 'quantity removed' to 'factor'
str(removal_counts)
removal_counts$quantity_removed <- as.factor(removal_counts$quantity_removed)

##with site as a random effect 
m0 <- clmm(
  quantity_removed ~ tree_species +
    removal_height_cm +
    distance_from_site_m +
    dbh_cm + 
    (1 | Site_ID),
  data = removal_counts)
summary(m0)

##adding interaction terms, one at a time

##Species × Height - Are egg masses more common higher on certain tree species?
m1 <- clmm(
  quantity_removed ~ tree_species * removal_height_cm +
    distance_from_site_m +
    dbh_cm +
    (1 | Site_ID),
  data = removal_counts)
summary(m1)

##Distance × Height - more eggs near human disturbance?
m2 <- clmm(
  quantity_removed ~ tree_species +
    removal_height_cm * distance_from_site_m +
    dbh_cm +
    (1 | Site_ID),
  data = removal_counts)
summary(m2)

anova(m0, m1, m2)

#getting a Hessian warning in m1, with the tree species and removal height
#interaction term, because some of the combinations
#are empty.  So will 

##collapse tree species into genera levels, which will
##preserve biological meaning and statistical stability, but improve structure

removal_counts <- removal_counts %>%
  mutate(
    tree_genus = case_when(
      tree_species %in% c("Oak", "Rea Oak", "White Oak") ~ "Oak",
      tree_species %in% c("Pine", "White Pine") ~ "Pine",
      tree_species %in% c("Maple", "Red Maple", "Silver Maple") ~ "Maple",
      tree_species %in% c("Beech", "American Beech") ~ "Beech",
      tree_species %in% c("Birch", "White Birch") ~ "Birch",
      tree_species %in% c("Ash", "Dead tree", "Hackberry") ~ "Other"
    ),
    tree_genus = factor(tree_genus, levels = 
          c("Beech", "Birch", "Maple", "Oak", "Pine", "Other"), ordered = TRUE)
  )

#check data structure
table(removal_counts$quantity_removed, removal_counts$tree_genus)

#convert 'tree genus' to 'character'
str(removal_counts)
removal_counts$tree_genus <- as.character(removal_counts$tree_genus)
str(removal_counts)

##repeat all 3 models from above, with tree genera instead of species
m3 <- clmm(
  quantity_removed ~ tree_genus +
    removal_height_cm +
    distance_from_site_m +
    dbh_cm + 
    (1 | Site_ID),
  data = removal_counts)
summary(m3)


##Species × Height 
m4 <- clmm(
  quantity_removed ~ tree_genus * removal_height_cm +
    distance_from_site_m +
    dbh_cm +
    (1 | Site_ID),
  data = removal_counts)
summary(m4)

##Distance × Height 
m5 <- clmm(
  quantity_removed ~ tree_genus +
    removal_height_cm * distance_from_site_m +
    dbh_cm +
    (1 | Site_ID),
  data = removal_counts)
summary(m5)

anova(m3, m4, m5)

##with tree ID as a random effect

m6 <- clmm(
  quantity_removed ~ tree_genus +
    removal_height_cm +
    distance_from_site_m +
    dbh_cm + 
    (1 | id),
  data = removal_counts)
summary(m6)

##Species × Height 
m7 <- clmm(
  quantity_removed ~ tree_genus * removal_height_cm +
    distance_from_site_m +
    dbh_cm +
    (1 | id),
  data = removal_counts)
summary(m7)

##Distance × Height 
m8 <- clmm(
  quantity_removed ~ tree_genus +
    removal_height_cm * distance_from_site_m +
    dbh_cm +
    (1 | id),
  data = removal_counts)
summary(m8)

anova(m3, m4, m5, m6, m7, m8)

#other interactions that could be considered...
#Species × Height - Are egg masses more common higher on certain tree species?
#Distance × Tree Species - Are moths selecting different host trees closer to camps?\

### Prediction probability plot ---------------------------------------------

##Create prediction grid
#str(removal_counts$tree_genus)
#str(removal_counts$removal_height_cm)

#levels(removal_counts$tree_genus)
#levels(removal_counts$removal_height_cm)

#removal_counts$tree_genus <- factor(removal_counts$tree_genus)
#removal_counts$removal_height_cm <- factor(removal_counts$removal_height_cm)

#levels(removal_counts$tree_genus)
#levels(removal_counts$removal_height_cm)

#newdat <- expand.grid(
 # tree_genus = levels(removal_counts$tree_genus),
  #removal_height_cm = levels(removal_counts$removal_height_cm),
  #distance_from_site_m = mean(removal_counts$distance_from_site_m, na.rm = TRUE),
 # dbh_cm = mean(removal_counts$dbh_cm, na.rm = TRUE))

#newdat$tree_genus <- factor(
 # newdat$tree_genus,
#  levels = levels(removal_counts$tree_genus))

#newdat$removal_height_cm <- factor(
 # newdat$removal_height_cm,
#  levels = levels(removal_counts$removal_height_cm))

#stopifnot(
 # nrow(newdat) > 0,
#  all(levels(newdat$tree_genus) == levels(removal_counts$tree_genus)),
#  all(levels(newdat$removal_height_cm) == levels(removal_counts$removal_height_cm)))


#nrow(newdat)
# is 18

#names(beta)

#X <- model.matrix(
 # delete.response(terms(m4)),
#  newdat)

#X <- X[, names(beta), drop = FALSE]

#stopifnot(
 # nrow(X) == nrow(newdat),
#  ncol(X) == length(beta))

#theta <- m4$Theta
#eta <- as.vector(X %*% beta)

#levels(removal_counts$quantity_removed)
# "zero" < "low" < "medium" < "high"

#With 4 ordered categories, clmm estimates 3 thresholds
#theta[1] : zero | low
#theta[2] : low | medium
#theta[3] : medium | high

#The cumulative probabilities are:𝑃(𝑌≤𝑘)=logit−1(𝜃𝑘 −𝜂)

#Convert linear predictor → probabilities

#logit <- function(x) 1 / (1 + exp(-x))

#P_le_high   <- logit(theta[1] - eta)
#P_le_low    <- logit(theta[2] - eta)
#P_le_medium <- logit(theta[3] - eta)

#p_high   <- P_le_high
#p_low    <- P_le_low    - P_le_high
#p_medium <- P_le_medium - P_le_low
#p_zero   <- 1 - P_le_medium

#check - should be very close to 1
#range(p_high + p_low + p_medium + p_zero)
# every row's probably sums to 1

#bind predictions to grid
#pred_df <- cbind(
 # newdat,
#  high   = p_high,
 # low    = p_low,
#  medium = p_medium,
 # zero   = p_zero)

#Convert to long format for ggplot
#pred_long <- pred_df |>
#  pivot_longer(
 #   cols = c("high", "medium", "low", "zero"),
#    names_to = "removal_class",
 #   values_to = "probability")

#pred_long$removal_class <- factor(
 # pred_long$removal_class,
#  levels = c("high", "medium", "low", "zero"))

#Make the predicted probability plot
#ggplot(pred_long,
 #      aes(x = removal_height_cm,
 #          y = probability,
#  fill = removal_class)) +
 # geom_col(color = "black", width = 0.7) +
#  facet_wrap(~ tree_genus) +
 # scale_y_continuous(labels = scales::percent) +
#  labs(
 #   x = "Removal height",
#    y = "Predicted probability",
 #   fill = "Egg mass removal") +
#  theme_bw() +
 # theme(
#    axis.text.x = element_text(angle = 45, hjust = 1),
 #   strip.text = element_text(size = 8))

### Hatch rate analysis -----------------------------------------------------

#order hatch for an ordinal model

hatch_rates <- hatch_rates %>%
  mutate(
    tree_genus = factor(tree_genus),
    removal_height_cm = factor(removal_height_cm),
    site_ID = factor(site_ID),
    hatch = factor(hatch,
                   levels = c("Low", "Medium", "High", "Very high"),
                   ordered = TRUE)
  )


### Egg mass removed --------------------------------------------------------

m_mass <- lmer(
  log(mass_g + 0.01) ~ tree_genus + removal_height_cm + (1 | site_ID),
  data = hatch_rates
)

summary(m_mass)
anova(m_mass)


### Hatch rates -------------------------------------------------------------

str(hatch_rates)
# replace 'NA' with 'zero'
# add "zero" as a valid factor level
#levels(hatch_rates$hatch) <-
 # c(levels(hatch_rates$hatch), "zero")

# replace NA with "zero"
#hatch_rates$hatch[
#  is.na(hatch_rates$hatch)] <- "zero"

# change genus, removal height, quant removed to 'factor' 
hatch_rates$removal_height_cm <- as.factor(hatch_rates$removal_height_cm)
hatch_rates$hatch <- as.factor(hatch_rates$hatch)
hatch_rates$tree_genus <- as.factor(hatch_rates$tree_genus)

str(hatch_rates)


#m_hatch <- clmm(
 # hatch ~ tree_genus + removal_height_cm + (1 | site_ID),
  #data = hatch_rates)

#summary(m_hatch)


#getting a Hessian warning because some genus x height x site combinations
#are empty.  So will 

##collapse hatch into 2 ordered levels, which will
##preserve biological meaning and statistical stability, but improve structure

hatch_rates <- hatch_rates %>%
 mutate(
    hatch2 = case_when(
      hatch %in% c("zero", "Low", "Medium") ~ "Low",
      hatch %in% c("High", "Very high") ~ "High" ),
    hatch2 = factor(hatch2, levels = c("Low", "High"), ordered = TRUE))

#check data structure
#table(hatch_rates$hatch2, hatch_rates$tree_genus)
#table(hatch_rates$hatch2, hatch_rates$removal_height_cm)

#refit the model with the collapsed data
#m_hatch2 <- clmm(
 # hatch2 ~ tree_genus + removal_height_cm + (1 | site_ID),
  #data = hatch_rates)
#summary(m_hatch2)

#still getting warnings, so will 

##collapse rare genera into 1 group

hatch_rates$tree_group <- fct_collapse(
  hatch_rates$tree_genus,
  Other = c("Birch", "Beech", "Dead tree", "Hackberry", "inanimate", "Ash")
)

table(hatch_rates$tree_group, hatch_rates$hatch2)

m_hatch2 <- clmm(
  hatch2 ~ tree_group + removal_height_cm + (1 | site_ID),
  data = hatch_rates
)
summary(m_hatch2)


### Hatch counts ------------------------------------------------------------

##check for overdispersion

mean(hatch_rates$hatch_count)
var(hatch_rates$hatch_count)

m_count <- glmmTMB(
  hatch_count ~ tree_genus + removal_height_cm + (1 | site_ID),
  family = nbinom2,
  data = hatch_rates
)
summary(m_count)

##Overdispered (=0.559)

#residual diagnostics

sim_res <- simulateResiduals(m_count)
plot(sim_res)

#same model, but with the genera groupings from ‘hatch levels’ 
m_count2 <- glmmTMB(
  hatch_count ~ tree_group + removal_height_cm + (1 | site_ID),
  family = nbinom2,
  data = hatch_rates
)
summary(m_count2)

sim_res2 <- simulateResiduals(m_count2)
plot(sim_res2)


## Hatch SUCCESS -----------------------------------------------------------

#calculate rate of eggs hatched per gram removed
hatch_rates <- hatch_rates %>%
  mutate(
    prop_hatched = hatch_count / mass_g
  )

#summarize by tree genus and by tree group
summary_by_genus <- hatch_rates %>%
  group_by(tree_genus) %>%
  summarise(
    total_mass = sum(mass_g),
    total_hatched = sum(hatch_count),
    prop_hatched = total_hatched / total_mass,
    n = n()
  ) %>%
  arrange(desc(prop_hatched))

summary_by_group <- hatch_rates %>%
  group_by(tree_group) %>%
  summarise(
    total_mass = sum(mass_g),
    total_hatched = sum(hatch_count),
    prop_hatched = total_hatched / total_mass,
    n = n()
  ) %>%
  arrange(desc(prop_hatched))

ggplot(hatch_rates, aes(x = tree_genus, y = prop_hatched)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ylab("Proportion of eggs hatched") +
  xlab("Tree genus") +
  theme_minimal()

ggplot(hatch_rates, aes(x = tree_group, y = prop_hatched)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ylab("Proportion of eggs hatched") +
  xlab("Tree group") +
  theme_minimal()



m_rate <- glmmTMB(
  hatch_count ~ tree_group + offset(log(mass_g)) + (1 | site_ID),
  family = nbinom2, data = hatch_rates)

summary(m_rate)



