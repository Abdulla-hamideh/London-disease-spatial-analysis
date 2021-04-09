# loading the necessary library 
library(sp)
library(CARBayes)
library(coda)
library(R2OpenBUGS)
library(spdep)
library(rgeos)

# loading the dataset 
load("london.RData")
head(london@data)

# 1.A) creating a subset of the assigned data 
london_s <- st_as_sf(london)
london_sub <- london_s[london_s$year == 2005,]
london_2005 <- as_Spatial(london_sub)

# 1.B) calculating Sir  
sir = london_2005$Observed/london_2005$Expected
# adding sir to the dataset
london_2005$Sir = sir

#1.C) plotting the sir
spplot(london_2005,"Sir", colorkye = T, main=list(label="SIR of admissions for respiratory diseases in London 2005"))

# 1.D) plotting the relationship between the sir and the other 3 variables in one image 
par(mfrow = c(1,3))
# first variable 
plot(london_2005$Sir~london_2005$PM25, data = as.data.frame(london_2005), col = "blue", main = "Association between SIR and PM25")
# fitting a line 
abline(lm(london_2005$Sir~london_2005$PM25), col = "blue",lty=2, lwd=3)
# second variable 
plot(london_2005$Sir~london_2005$JSA, data = as.data.frame(london_2005), col = "blue", main = "Association between SIR and JSA")
# fitting a line 
abline(lm(london_2005$Sir~london_2005$JSA), col = "blue",lty=2, lwd=3)
# third variable 
plot(london_2005$Sir~london_2005$Price, data = as.data.frame(london_2005), col = "blue", main = "Association between SIR and Price")
# fitting a line 
abline(lm(london_2005$Sir~london_2005$Price), col = "blue",lty=2, lwd=3)


# 2.A) Assign data.
model = list(Y = london_2005$Observed, E = london_2005$Expected, PM25 = as.numeric(scale(london_2005$PM25)), JSA = as.numeric(scale(london_2005$JSA)), Price= as.numeric(scale( london_2005$Price)), N = length(london_2005$Expected))

Bet = list(list(beta0 = 0, beta1 =0, beta2 = 0, beta3 =0))

## run MCMC sample using bugs function.
model.bugs = bugs(data = model,
                 model.file = "Model.txt",
                 parameters.to.save = c("beta0", "beta1", "beta2", "beta3"),
                 inits = NULL,
                 n.iter = 1000,
                 n.burnin = 500,
                 n.thin = 1,
                 n.chains = 4,
                 bugs.directory = "C:/Users/abdul/Documents/WinBUGS14")
