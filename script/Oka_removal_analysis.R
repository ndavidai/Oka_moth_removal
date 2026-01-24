
#loading data 
removal_counts <- read.csv("input/Oka_removal_counts_clean.csv")
hatch_rates <- read.csv("input/Oka_hatch_rates_clean.csv")

library(janitor)                       
library(dplyr)
library(ggplot2)
library(esquisse)
library(ordinal)
library(lme4)
library(tidyr)


# removal counts analyses -------------------------------------------------

# replace 'NA' with 'zero'
# add "zero" as a valid factor level
levels(removal_counts$quantity_removed) <-
  c(levels(removal_counts$quantity_removed), "zero")

# replace NA with "zero"
removal_counts$quantity_removed[
  is.na(removal_counts$quantity_removed)] <- "zero"


#visual explorations
ggplot(removal_counts, aes(distance_from_site_m, quantity_removed)) +
  geom_boxplot() +
  facet_wrap(~ tree_species)

ggplot(removal_counts, aes(removal_height_cm, quantity_removed)) +
  geom_boxplot()


# Ordinal mixed-effects model ---------------------------------------------

removal_counts$quantity_removed <- factor(
  removal_counts$quantity_removed,
  levels = c("zero", "low", "medium", "high"),
  ordered = TRUE
)


#convert 'quantity removed' to 'factor'
str(removal_counts)
removal_counts$quantity_removed <- as.factor(removal_counts$quantity_removed)

m0 <- clmm(
  quantity_removed ~ tree_species +
    removal_height_cm +
    distance_from_site_m +
    dbh_cm + 
    (1 | id),
  data = removal_counts)

##adding interaction terms, one at a time

#Species × Height - Are egg masses more common higher on certain tree species?
m1 <- clmm(
  quantity_removed ~ tree_species * removal_height_cm +
    distance_from_site_m +
    dbh_cm +
    (1 | id),
  data = removal_counts)
summary(m1)

#Distance × Height - more eggs near human disturbance?
m2 <- clmm(
  quantity_removed ~ tree_species +
    removal_height_cm * distance_from_site_m +
    dbh_cm +
    (1 | id),
  data = removal_counts)
summary(m2)

anova(m0, m1, m2)

#other interactions that could be considered...
#Species × Height - Are egg masses more common higher on certain tree species?
#Distance × Tree Species - Are moths selecting different host trees closer to camps?\

# Prediction probability plot ---------------------------------------------

#Create prediction grid
newdat <- removal_counts %>%
  distinct(tree_species, removal_height_cm) %>%
  mutate(
    distance_from_site_m = mean(removal_counts$distance_from_site_m, na.rm = TRUE),
    dbh_cm = mean(removal_counts$dbh_cm, na.rm = TRUE)
  )

nrow(newdat)
# should be ~40

names(beta)

X <- model.matrix(
  delete.response(terms(m1)),
  newdat
)

X <- X[, colnames(X) %in% names(beta), drop = FALSE]

stopifnot(length(beta) == ncol(X))

eta <- as.vector(X %*% beta)
theta <- m1$Theta

levels(removal_counts$quantity_removed)
# "zero" < "low" < "medium" < "high"

#With 4 ordered categories, clmm estimates 3 thresholds
#theta[1] : zero | low
#theta[2] : low | medium
#theta[3] : medium | high

#The cumulative probabilities are:𝑃(𝑌≤𝑘)=logit−1(𝜃𝑘 −𝜂)

#Convert linear predictor → probabilities
logit <- function(x) 1 / (1 + exp(-x))

P_le_high   <- logit(theta[1] - eta)
P_le_low    <- logit(theta[2] - eta)
P_le_medium <- logit(theta[3] - eta)

p_high   <- P_le_high
p_low    <- P_le_low    - P_le_high
p_medium <- P_le_medium - P_le_low
p_zero   <- 1 - P_le_medium

#check - should be very close to 1
range(p_high + p_low + p_medium + p_zero)
# every row's probably sums to 1

#bind predictions to grid
pred_df <- cbind(
  newdat,
  high   = p_high,
  low    = p_low,
  medium = p_medium,
  zero   = p_zero
)

#Convert to long format for ggplot
pred_long <- pred_df |>
  pivot_longer(
    cols = c("high", "medium", "low", "zero"),
    names_to = "removal_class",
    values_to = "probability"
  )

pred_long$removal_class <- factor(
  pred_long$removal_class,
  levels = c("high", "medium", "low", "zero")
)


#Make the predicted probability plot
ggplot(pred_long,
       aes(x = removal_height_cm,
           y = probability,
           fill = removal_class)) +
  geom_col(color = "black", width = 0.7) +
  facet_wrap(~ tree_species) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Removal height",
    y = "Predicted probability",
    fill = "Egg mass removal"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)
  )



