# does latitudes decreases with climate change? ---------------------------
# latitude x climate change -----------------------------------------------

library(lmerTest)

m <- lmerTest::lmer(latitude ~ scenario * spGroup + (1 | species), data = centroids)
anova(m)

m2 <- lmerTest::lmer(latitude ~ scenario + spGroup + (1 | species), data = centroids)
anova(m, m2) # significant, do not remove interaction

summary(m)

# does suitable area decreases with climate change? -----------------------

areas$scenario <- factor(areas$scenario, # area in R_data line 33
                        levels = c(
                                 "present",
                                 "future_11",
                                 "future_12",
                                 "future_21",
                                 "future_22")
                        ) 

library(lmerTest)

areas$species <- factor(areas$species, levels = c("Bothrops fonsecai", "Bothrops cotiara", "Bothrops itapetiningae", "Bothrops alternatus"))
areas$scenario <- factor(areas$scenario, levels = c("present","future_11", "future_12", "future_21", "future_22"))
areas$spGroups <- factor(areas$spGroup, levels = c("Forest & Mountain","Open Areas & Plateau"))
areas$year <- factor(areas$year, levels = c("Present","2050", "2100"))
areas$futureScenario <- factor(areas$futureScenario, levels = c("SSP 370","SSP 585"))

mod_misto <- lmerTest::lmer(log(area_km2) ~ scenario * spGroup + (1 | species), data = areas) # with interaction
anova(mod_misto) # only scenario is significant

mod_misto2 <- lmerTest::lmer(log(area_km2) ~ scenario + spGroup + (1 | species), data = areas) # without interaction
anova(mod_misto, mod_misto2) # keep without interaction model 2

mod_misto3 <- lmerTest::lmer(log(area_km2) ~ scenario + (1 | species), data = areas) # without species group
anova(mod_misto2, mod_misto3) # keep model with species groups model 2

summary(mod_misto2)

# extracting predictions --------------------------------------------------

newdat <- expand.grid(
  scenario = unique(areas$scenario),
  species = unique(areas$species)
)

newdat$pred_log <- predict(mod_misto, newdata = newdat, re.form = ~(1|species))
newdat$pred_area <- exp(newdat$pred_log)

  
library(ggplot2)

ggplot(
  newdat,
  aes(
    x = year,
    y = pred_area,
    color = species,
    group = interaction(species, futureScenario),
    linetype = futureScenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Predicted Area (km²)",
    color = "Species",
    linetype = "Scenario"
  )
