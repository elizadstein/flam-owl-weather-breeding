library(R2jags)
library(tidyverse)
library(lubridate)
library(mcp)
library(tidybayes)

### READ IN DATA ###

mass_climate <- read.csv("../data/mass_data.csv")

# Split into all adult (a), adult male (m), adult female (f), and nestling (n)
mass_climate_a <- filter(mass_climate, age == "AD" | age == "U")
mass_climate_m <- filter(mass_climate, age == "AD" & sex == "M")
mass_climate_f <- filter(mass_climate, age == "AD" & sex == "F")
mass_climate_n <- filter(mass_climate, age == "N" | age == "F")

owlet.wet <- mass_climate_n %>%
  filter(wet_dry == "Wet")

owlet.dry <- mass_climate_n %>%
  filter(wet_dry == "Dry")

owlet.warm <- mass_climate_n %>%
  filter(warm_cold == "Warm")

owlet.cold <- mass_climate_n %>%
  filter(warm_cold == "Cold")



### PRECIP ###

## Global model for wet/dry years:
m.mass = list(
  mass ~  1 + (1|band) + day, # slope of seg_1
  1 + (1|band) ~ 0 + day     # slope of seg_2
)



# Run model for wet years
m.mass.wet = mcp(m.mass, data = owlet.wet, 
                 iter = 4000, adapt = 2000)     #20,000 iter, 4,000 burn-in

# Run model for dry years
m.mass.dry = mcp(m.mass, data = owlet.dry, 
                 iter = 8000, adapt = 2000)     #20,000 iter, 4,000 burn-in



### TEMP ###


# Run model for warm years
m.mass.warm = mcp(m.mass, data = owlet.warm, 
                  iter = 8000, adapt = 800)     #20,000 iter, 4,000 burn-in

# Run model for cold years
m.mass.cold = mcp(m.mass, data = owlet.cold, 
                  iter = 8000, adapt = 800)     #20,000 iter, 4,000 burn-in




### EVALUATE OVERLAP ###

### Wet/Dry

# Print summary
print(m.mass.wet)
print(m.mass.dry)

# Store as df
m.mass.wet.df <- data.frame("wet_dry" = "Wet", unclass(summary(m.mass.wet)), check.names = FALSE, stringsAsFactors = FALSE)
m.mass.wet.df$name <- c("Changepoint", "Slope 1", "Slope 2", "Intercept 1", "Intercept 2", "Sigma")


m.mass.dry.df <- data.frame("wet_dry" = "Dry", unclass(summary(m.mass.dry)), check.names = FALSE, stringsAsFactors = FALSE)
m.mass.dry.df$name <- c("Changepoint", "Slope 1", "Slope 2", "Intercept 1", "Intercept 2", "Sigma")

model.output.owlet.precip <- rbind(m.mass.wet.df, m.mass.dry.df)
model.output.owlet.precip <- rename(model.output.owlet.precip, climate = wet_dry)


### Warm/Cold

print(m.mass.warm)
print(m.mass.cold)

# Store as df
m.mass.warm.df <- data.frame("warm_cold" = "Warm", unclass(summary(m.mass.warm)), check.names = FALSE, stringsAsFactors = FALSE)
m.mass.warm.df$name <- c("Changepoint", "Slope 1", "Slope 2", "Intercept 1", "Intercept 2", "Sigma")


m.mass.cold.df <- data.frame("warm_cold" = "Cold", unclass(summary(m.mass.cold)), check.names = FALSE, stringsAsFactors = FALSE)
m.mass.cold.df$name <- c("Changepoint", "Slope 1", "Slope 2", "Intercept 1", "Intercept 2", "Sigma")

model.output.owlet.temp <- rbind(m.mass.warm.df, m.mass.cold.df)
model.output.owlet.temp <- rename(model.output.owlet.temp, climate = warm_cold)


# extract 3 chains from m.mass.wet
chain.1.mass.wet <- as.data.frame(m.mass.wet$mcmc_post[[1]])
chain.1.mass.wet$chain <- as.factor(1)
names(chain.1.mass.wet) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.2.mass.wet <- as.data.frame(m.mass.wet$mcmc_post[[2]])
chain.2.mass.wet$chain <- as.factor(2)
names(chain.2.mass.wet) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.3.mass.wet <- as.data.frame(m.mass.wet$mcmc_post[[3]])
chain.3.mass.wet$chain <- as.factor(3)
names(chain.3.mass.wet) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

# combine into one wet dataframe
df.mass.wet <- rbind(chain.1.mass.wet, chain.2.mass.wet, chain.3.mass.wet)
df.mass.wet$wet_dry <- "Wet"


# extract 3 chains from m.mass.dry
chain.1.mass.dry <- as.data.frame(m.mass.dry$mcmc_post[[1]])
chain.1.mass.dry$chain <- as.factor(1)
names(chain.1.mass.dry) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.2.mass.dry <- as.data.frame(m.mass.dry$mcmc_post[[2]])
chain.2.mass.dry$chain <- as.factor(2)
names(chain.2.mass.dry) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.3.mass.dry <- as.data.frame(m.mass.dry$mcmc_post[[3]])
chain.3.mass.dry$chain <- as.factor(3)
names(chain.3.mass.dry) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

# combine into one dry dataframe
df.mass.dry <- rbind(chain.1.mass.dry, chain.2.mass.dry, chain.3.mass.dry)
df.mass.dry$wet_dry <- "Dry"


# combine to create new df for all nest age posteriors
posteriors.mass.precip <- rbind(df.mass.wet, df.mass.dry)
names(posteriors.mass.precip) <- c("Changepoint", "Segment_1", "Segment_2", "Intercept_1", "Intercept_2", "Sigma", "Chain", "Wet_Dry")

# extract 3 chains from m.mass.wet
chain.1.mass.warm <- as.data.frame(m.mass.warm$mcmc_post[[1]])
chain.1.mass.warm$chain <- as.factor(1)
names(chain.1.mass.warm) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.2.mass.warm <- as.data.frame(m.mass.warm$mcmc_post[[2]])
chain.2.mass.warm$chain <- as.factor(2)
names(chain.2.mass.warm) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.3.mass.warm <- as.data.frame(m.mass.warm$mcmc_post[[3]])
chain.3.mass.warm$chain <- as.factor(3)
names(chain.3.mass.warm) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

# combine into one wet dataframe
df.mass.warm <- rbind(chain.1.mass.warm, chain.2.mass.warm, chain.3.mass.warm)
df.mass.warm$warm_cold <- "Warm"


# extract 3 chains from m.mass.dry
chain.1.mass.cold <- as.data.frame(m.mass.cold$mcmc_post[[1]])
chain.1.mass.cold$chain <- as.factor(1)
names(chain.1.mass.cold) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.2.mass.cold <- as.data.frame(m.mass.cold$mcmc_post[[2]])
chain.2.mass.cold$chain <- as.factor(2)
names(chain.2.mass.cold) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

chain.3.mass.cold <- as.data.frame(m.mass.cold$mcmc_post[[3]])
chain.3.mass.cold$chain <- as.factor(3)
names(chain.3.mass.cold) <- c("cp", "seg1", "seg2", "int1", "int2", "sigma", "chain")

# combine into one dry dataframe
df.mass.cold <- rbind(chain.1.mass.cold, chain.2.mass.cold, chain.3.mass.cold)
df.mass.cold$warm_cold <- "Cold"


# combine to create new df for all nest age posteriors
posteriors.mass.temp <- rbind(df.mass.warm, df.mass.cold)
names(posteriors.mass.temp) <- c("Changepoint", "Segment_1", "Segment_2", "Intercept_1", "Intercept_2", "Sigma", "Chain", "Warm_Cold")


### PREDICT ###

# Precipitation
f.mass.wet <- fitted(m.mass.wet) # get fitted values
f.mass.wet$wet_dry <- "Wet"

f.mass.dry <- fitted(m.mass.dry) # get fitted values
f.mass.dry$wet_dry <- "Dry"

f.mass.precip <- rbind(f.mass.wet, f.mass.dry) # Merge into one df

# Temperature
f.mass.warm <- fitted(m.mass.warm) # get fitted values
f.mass.warm$warm_cold <- "Warm"

f.mass.cold <- fitted(m.mass.cold) # get fitted values
f.mass.cold$warm_cold <- "Cold"

f.mass.temp <- rbind(f.mass.warm, f.mass.cold) # Merge into one df


### QUANTILES ###

## Custom function to get quantiles (prob) from a given vector (vec):
getQuant <- function(vec, prob){
  quantile(vec, probs = prob, names = FALSE)
}

### Apply getQuant for 80% and 95% and output data.frame

## Precipitation
df.quant.precip <- cbind(t(data.frame(
  lower_80 = apply(X = df.mass.wet[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.1),
  upper_80 = apply(X = df.mass.wet[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.9),
  lower_95 = apply(X = df.mass.wet[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.025),
  upper_95 = apply(X = df.mass.wet[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.975), row.names = c("wet_cp", "wet_seg1", "wet_seg2"))),
  t(data.frame(
    lower_80 = apply(X = df.mass.dry[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.1),
    upper_80 = apply(X = df.mass.dry[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.9),
    lower_95 = apply(X = df.mass.dry[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.025),
    upper_95 = apply(X = df.mass.dry[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.975), row.names = c("dry_cp", "dry_seg1", "dry_seg2")))
)

# Reorder columns for ease of viewing
df.quant.precip <- df.quant.precip[, c("wet_cp", "dry_cp", "wet_seg1", "dry_seg1", "wet_seg2", "dry_seg2")]



## Temperature

df.quant.temp <- cbind(t(data.frame(
  lower_80 = apply(X = df.mass.warm[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.1),
  upper_80 = apply(X = df.mass.warm[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.9),
  lower_95 = apply(X = df.mass.warm[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.025),
  upper_95 = apply(X = df.mass.warm[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.975), row.names = c("warm_cp", "warm_seg1", "warm_seg2"))),
  t(data.frame(
    lower_80 = apply(X = df.mass.cold[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.1),
    upper_80 = apply(X = df.mass.cold[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.9),
    lower_95 = apply(X = df.mass.cold[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.025),
    upper_95 = apply(X = df.mass.cold[,1:3], MARGIN = 2, FUN = getQuant, prob = 0.975), row.names = c("cold_cp", "cold_seg1", "cold_seg2")))
)

# Reorder columns for ease of viewing
df.quant.temp <- df.quant.temp[, c("warm_cp", "cold_cp", "warm_seg1", "cold_seg1", "warm_seg2", "cold_seg2")]

df.quant.nestling <- t(cbind(df.quant.precip, df.quant.temp))


### WRITE MODEL OUTPUTS ###

# saveRDS(m.mass.wet, "m_mass_wet.rda")
# saveRDS(m.mass.dry, "m_mass_dry.rda")
# saveRDS(m.mass.warm, "m_mass_warm.rda")
# saveRDS(m.mass.cold, "m_mass_cold.rda")
# 
# write.csv(f.mass.wet, "f_mass_wet.csv")
# write.csv(f.mass.dry, "f_mass_dry.csv")
# write.csv(f.mass.warm, "f_mass_warm.csv")
# write.csv(f.mass.cold, "f_mass_cold.csv")
# 
# write.csv(df.quant.temp, "df_quant_temp.csv")
# write.csv(df.quant.precip, "df_quant_precip.csv")