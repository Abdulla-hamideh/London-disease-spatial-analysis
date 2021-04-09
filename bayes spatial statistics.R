# loading the necessary library 
library(sp)
library(CARBayes)
library(coda)
library(R2WinBUGS)
library(spdep)
library(rgeos)
library(devtools)
remotes::install_github("davidaarmstrong/asmcjr")
library(remotes)


# loading the dataset 
load("london.RData")
head(london@data)

# 1.A) creating a subset of the assigned data 
london@data <-london@data[london@data$year == 2005,] 
head(london@data)
# 1.B) calculating Sir  
sir = london$Observed/london$Expected
# adding sir to the dataset
london$Sir = sir

#1.C) plotting the sir
spplot(london,"Sir", colorkye = T, main=list(label="SIR of admissions for respiratory diseases in London 2005"))

# 1.D) plotting the relationship between the sir and the other 3 variables in one image 
par(mfrow = c(1,3))
# first variable 
plot(london$Sir~london$PM25, data = as.data.frame(london), col = "blue", main = "Association between SIR and PM25")
# fitting a line 
abline(lm(london$Sir~london$PM25), col = "blue",lty=2, lwd=3)
# second variable 
plot(london$Sir~london$JSA, data = as.data.frame(london), col = "blue", main = "Association between SIR and JSA")
# fitting a line 
abline(lm(london$Sir~london$JSA), col = "blue",lty=2, lwd=3)
# third variable 
plot(london$Sir~london$Price, data = as.data.frame(london), col = "blue", main = "Association between SIR and Price")
# fitting a line 
abline(lm(london$Sir~london$Price), col = "blue",lty=2, lwd=3)


# 2.A) Assign data.
data = list(Y = london$Observed, E = london$Expected, PM25 = as.numeric(scale(london$PM25)), JSA = as.numeric(scale(london$JSA)), Price= as.numeric(scale( london$Price)), N = length(london$Expected))
Bet = list(list(beta0 = 0, beta1 =0, beta2 = 0, beta3 =0))

## run MCMC sample using bugs function.
model.fit = R2WinBUGS::bugs(data = data,
                 model.file = "Model.txt",
                 parameters.to.save =c("beta0", "beta1", "beta2", "beta3"),
                 inits = rep(Bet,4),
                 n.iter = 5000,
                 n.burnin = 2500,
                 n.thin = 2,
                 n.chains = 4,
                 bugs.directory = "C:/Users/abdul/Documents/WinBUGS14")

# plotting the trace plot for beta 0-3
par(mfrow = c(2,2))
plot(model.fit$sims.list$beta0, type ="l", ylab = "beta0")
plot(model.fit$sims.list$beta1, type ="l", ylab = "beta1")
plot(model.fit$sims.list$beta2, type ="l", ylab = "beta2")
plot(model.fit$sims.list$beta3, type ="l", ylab = "beta3")
mtext("Traceplots", side = 3, line = -3, outer = TRUE)
## Posterior median and credible intervals.
print(model.fit, digits.summary = 3)

###2b: Use Gelman-Rubin diagostic (Rhat) and geweke diagnostic plots to assess convergence.

#### Gelman-rubin, checking RHAT
par(mfrow = c(1,1))
model.fit$summary[,8]
hist(model.fit$summary[,8], main = " Distribution of Rhat", ylab = " Frequency", xlab = "Rhat")
### it is between 1 -1.2 --> converge

asmcjr::geweke.ggplot(model.fit)


### geweke convergence
par(mfrow = c(2,2))
geweke.plot(model)
geweke.plot(as.mcmc(as.data.frame(model.fit$sims.list$beta0)))
geweke.plot(as.mcmc(as.data.frame(model.fit$sims.list$beta1)))
geweke.plot(as.mcmc(as.data.frame(model.fit$sims.list$beta2)))
geweke.plot(as.mcmc(as.data.frame(model.fit$sims.list$beta3)))
mtext("Traceplots", side = 3, line = -3, outer = TRUE)
###2c: Calculate pearson residuals and check 2 model assumptions
# Overdispersion

model.fit$summary[,1]
mu =london$Expected*exp(-0.25679050 + 0.15936450*as.numeric(scale(london$PM25)) +0.13717505*as.numeric(scale(london$JSA)) -0.13717505*as.numeric(scale( london$Price)))
pearson = (london$Observed - mu)/sqrt(mu)

## check overdispersion - plot fitted values and residual

par(mfrow = c(1,1))
plot(pearson ~ mu, xlab = "mu", ylab = " res", main = " residuals vs mu", col = "mediumorchid4")
var(pearson)

## check independence - moran's i statistics
check <- nb2listw(poly2nb(london))
moran =moran.test(pearson, check,alternative = "two.sided")
moran
moran$estimate[1]
moran$p.value


# 3.a  
#proudce the adjacency matrix, and various pieces of associated info Bugs needs
W = nb2mat(poly2nb(london), style = "B")
inds = lapply(1:nrow(W), function(i) which (W[i, ] ==1))
Adj = Reduce("c", inds)
Num.Adj = rowSums(W)
SumNumNeigh = sum(Num.Adj)


# combine all of the data, and constatns into a single list
car_data= list(observed = london$Observed, expected = london$Expected, N = nrow(london), PM25 = as.numeric(scale(london$PM25)), JSA = as.numeric(scale(london$JSA)), Price= as.numeric(scale( london$Price)), Adj = Adj, Num = Num.Adj, SumNumNeigh = SumNumNeigh)
car_bet = list(list(beta0 =0,beta1 =0, beta2 = 0, beta3 =0))

car.bugs = R2WinBUGS::bugs(data=car_data, 
                 model.file = "model2.txt", 
                 parameters.to.save =  c("beta0", "beta1", "beta2", "beta3", "phi"), 
                 inits = rep(car_bet,5), 
                 n.iter = 10000, 
                 n.burnin = 5000, 
                 n.thin = 1, 
                 n.chains = 5,
                 bugs.directory = "C:/Users/abdul/Documents/WinBUGS14")

par(mfrow = c(2,3))
plot(car.bugs$sims.list$beta0, type ="l", ylab = " beta0")
plot(car.bugs$sims.list$beta1, type ="l", ylab = " beta1")
plot(car.bugs$sims.list$beta2, type ="l", ylab = " beta2")
plot(car.bugs$sims.list$beta3, type ="l", ylab = " beta3")
plot(car.bugs$sims.list$phi[,1], type ="l", ylab = " phi.1")
plot(car.bugs$sims.list$phi[,50], type ="l", ylab = " phi.50")
plot(car.bugs$sims.list$phi[,100], type ="l", ylab = " phi.100")
plot(car.bugs$sims.list$phi[,150], type ="l", ylab = " phi.150")
plot(car.bugs$sims.list$phi[,300], type ="l", ylab = " phi.300")
plot(car.bugs$sims.list$phi[,600], type ="l", ylab = " phi")
mtext("Traceplots", side = 3, line = -3, outer = TRUE)
print(car.bugs, digits.summary = 3)
car.bugs$summary[2,1:7]

## 3.b - Use Gelman-Rubin diagostic (Rhat) and geweke diagnostic plots to assess convergence.

#### Gelman-rubin, checking RHAT
par(mfrow = c(1,1))
hist(car.bugs$summary[,8],main = " Distribution of Rhat", ylab = " Frequency", xlab = "Rhat",col= c("slategray1", "turquoise1","lightblue", "lavender", "lightcyan1", "cyan"))
### it is between 1 -1.2 --> converge.

### geweke convergence
asmcjr::geweke.ggplot(car.bugs)
par(mfrow = c(2,2))
geweke.plot(as.mcmc(as.data.frame(car.bugs$sims.list$beta0)))
geweke.plot(as.mcmc(as.data.frame(car.bugs$sims.list$beta1)))
geweke.plot(as.mcmc(as.data.frame(car.bugs$sims.list$beta2)))
geweke.plot(as.mcmc(as.data.frame(car.bugs$sims.list$beta3)))
geweke.plot(as.mcmc(as.data.frame(car.bugs$sims.list$phi[,c(1,50,150,300)])))


###3c: Calculate pearson residuals and check 2 model assumption 
car.bugs$summary[,1]
phi.mean = apply(car.bugs$sims.list$phi, 2, mean)
phi.mean

#########
mu.car_model = london$Expected*exp(-2.657904e-01 + 1.847929e-01*as.numeric(scale(london$PM25)) +1.133079e-01*as.numeric(scale(london$JSA))-4.579765e-02*as.numeric(scale( london$Price))+phi.mean)
pearson.car_model = (london$Observed - mu.car_model)/sqrt(mu.car_model)
var(pearson.car_model)

## check overdispersion - plot fitted values and residual
par(mfrow = c(1,1))
plot(pearson.car_model ~ mu.car_model, xlab = "Mu", ylab = " residual", main = "Pearson residuals vs mu", col = "hotpink3")
var(pearson.car_model)

## check independence - moran's i statistics
Wnb.car <- nb2listw(poly2nb(london))
moran.car = moran.test(pearson.car_model, Wnb.car)
moran.car
moran.car$estimate[1]
# summary table 
names(model.fit)
names(car.bugs)


Model = c("Poisson Regression", "Poisson Car")
DIC = c(model.fit$DIC,car.bugs$DIC) 
MoranI = c(moran$estimate[1],moran.car$estimate[1])
p_value = c(moran$p.value,moran.car$p.value)
PS_Var = c(var(pearson), var(pearson.car_model))



PM25_mean_exp = c(exp(model.fit$summary[2,1]),exp(car.bugs$summary[2,1]))
PM25_95Low_exp = c(exp(model.fit$summary[2,3]),exp(car.bugs$summary[2,3]))
PM25_95_HI_exp = c(exp(model.fit$summary[2,7]),exp(car.bugs$summary[2,7]))
PM25_mean = c(model.fit$summary[2,1],car.bugs$summary[2,1])
PM25_95Low = c(model.fit$summary[2,3],car.bugs$summary[2,3])
PM25_95_HI = c(model.fit$summary[2,7],car.bugs$summary[2,7])

PM25_mean_exp
PM25_95Low_exp
PM25_95_HI_exp
exp(car.bugs$summary[1:4,c(1,3,7)])


summary.table.converted = data.frame(Model, DIC, MoranI, p_value, PS_Var, PM25_mean_exp, PM25_95Low_exp, PM25_95_HI_exp)
summary.table.converted

summary.table = data.frame(Model, DIC, MoranI, p_value, PS_Var, PM25_mean, PM25_95Low, PM25_95_HI)
summary.table

PM25_mean
PM25_95_HI

exp(0.1593506)
log(1.172749)
