
# Final Project
# Forest Fire area burned Bayesian analysis

library(tidyverse)
library(mvtnorm)
library(coda)
library(xtable)

# read in the forest fire data
forestfire <- read.csv(header=TRUE,"forestfires.csv")

# EDA
par(mfrow=c(1,2))
hist(forestfire$area[forestfire$area>0],breaks=100,probability = TRUE,main="Histogram of Areas",xlab="area")
x <- seq(0.01,400,by=.01)
lines(x,dlnorm(x,1.5,1),col="red")
hist(log(forestfire$area[forestfire$area>0]),breaks=100,probability = TRUE,main="Histogram of Log Areas",xlab="log(area)")
par(mfrow=c(1,1))

# Look at data greater than 0
datnot0 <- forestfire[forestfire$area>0,]
plot(~log(area)+FFMC+DMC+ISI+temp+RH+wind+rain,data=datnot0)
plot(log(area)~rain,data=datnot0)
hist(log(forestfire$area+1),breaks=50)

# Just a quick lm :)
summary(lm(log(area)~. - X -Y,data=datnot0))

# Group months by season and days by weekends
# Maybe change the groupings
forestfire$Weekend <- as.factor(as.numeric(forestfire$day %in% c("sat","sun")))
fire_season <- c("jul","aug","sep")
fall_winter <- c("oct","nov","dec","jan","feb")
spring_summer <- c("mar","apr","may","jun")
forestfire$Season <- numeric(nrow(forestfire))
forestfire$Season[forestfire$month %in% fire_season] <- 1
forestfire$Season[forestfire$month %in% fall_winter] <- 2
forestfire$Season[forestfire$month %in% spring_summer] <- 3
forestfire$Season <- as.factor(forestfire$Season)

forestfire %>%
  group_by(Season) %>%
  summarise(n=n())


# The model used for my final analysis is model 9. Results from this model are found staring on lines 822 and 1145
# Backward DIC selection is used to determine which model out of the 10 to use.

# Multiple models used to do backward DIC selection
# Model 1: #####
# area_i ~ (1-pi) * delta_0(area_i) + pi * Lognormal(mu_i, lambda)
# pi ~ Beta(2,2)
# mu_i = Beta0 + Beta1-Beta2 * season + Beta4 * weekday + Beta5 * FFMC + Beta6 * DMC +
#   Beta7 * CD + Beta8 * ISI + Beta9 * temp + Beta10 * RH + Beta11 * wind + Beta12 * rain
# lambda ~ gamma(2,2)
# Beta_i ~iid Normal(0,.001)

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

# Functions for the likelihood and posterior
log_likelihood <- function(betas,lambda,pi) {
  log(1-pi)*n_with_0 + log(pi)*n_without_0 + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_pi <- function(pi) {
  dbeta(pi,2,2)
}
log_posterior <- function(betas,lambda,pi) {
  log_likelihood(betas,lambda,pi) + log_prior_lambda(lambda) + log_prior_betas(betas) + log_prior_pi(pi)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(0.5,.01,numeric(12))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
model1_draws <- matrix(nrow=Ndraws,numeric(Ndraws*14))
model1_draws[1,] <- center

accept <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for pi using Gibbs sampling
  if (i %% 3==1) {
    new_pi <- rbeta(1,n_without_0 + 2, n_with_0 + 2)
    model1_draws[i,] <- model1_draws[i-1,]
    model1_draws[i,1] <- new_pi
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model1_draws[i-1,3:14])^2))
    model1_draws[i,] <- model1_draws[i-1,]
    model1_draws[i,2] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model1_draws[i-1,3:14] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior(propose_new_vals,model1_draws[i-1,2],model1_draws[i-1,1])-log_posterior(model1_draws[i-1,3:14],model1_draws[i-1,2],model1_draws[i-1,1])
    model1_draws[i,] <- model1_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model1_draws[i,3:14] <- propose_new_vals
      accept <- accept + 1
    }
  }
}

accept/(Ndraws/3)
# Trace plots for each parameter
par(mfrow=c(2,2))
plot(1:Ndraws,model1_draws[,1],type="l")
plot(1:Ndraws,model1_draws[,2],type="l")
plot(1:Ndraws,model1_draws[,3],type="l")
plot(1:Ndraws,model1_draws[,4],type="l")
plot(1:Ndraws,model1_draws[,5],type="l")
plot(1:Ndraws,model1_draws[,6],type="l")
plot(1:Ndraws,model1_draws[,7],type="l")
plot(1:Ndraws,model1_draws[,8],type="l")
plot(1:Ndraws,model1_draws[,9],type="l")
plot(1:Ndraws,model1_draws[,10],type="l")
plot(1:Ndraws,model1_draws[,11],type="l")
plot(1:Ndraws,model1_draws[,12],type="l")
plot(1:Ndraws,model1_draws[,13],type="l")
plot(1:Ndraws,model1_draws[,14],type="l")
par(mfrow=c(1,1))
effectiveSize(model1_draws)

# It looks like I need to burn the first 1000
burned_draws_1 <- model1_draws[1001:Ndraws,]

# Check to see if the values are similar to using an lm
lm_fit <- lm(log_areas~Xmat -1)
lm_fit$coefficients
apply(burned_draws_1,2,mean)
cbind(lm_fit$coefficients,apply(burned_draws_1,2,mean)[3:14])

par(mfrow=c(2,2))
hist(burned_draws[,1])
hist(burned_draws[,2])
hist(burned_draws[,3])
hist(burned_draws[,4])
hist(burned_draws[,5])
hist(burned_draws[,6])
hist(burned_draws[,7])
hist(burned_draws[,8])
hist(burned_draws[,9])
hist(burned_draws[,10])
hist(burned_draws[,11])
hist(burned_draws[,12])
hist(burned_draws[,13])
hist(burned_draws[,14])
par(mfrow=c(1,1))
#####

# Model 2: #####
# area_i ~ (1-p_i) * delta_0(area_i) + p_i * Lognormal(mu_i, lambda)
# Logistic regression for p_i
# p_i = 1 / (1+exp(-X_i*alpha))
# alpha_i ~iid Normal(0,.001)
# mu_i = Beta0 + Beta1-Beta2 * season + Beta4 * weekday + Beta5 * FFMC + Beta6 * DMC +
#   Beta7 * CD + Beta8 * ISI + Beta9 * temp + Beta10 * RH + Beta11 * wind + Beta12 * rain
# lambda ~ gamma(2,2)
# Beta_i ~iid Normal(0,.001)

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

# Functions for the likelihood and posterior
log_likelihood_m2 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0%*%alphas) / (1+exp(-Xmat_with_0%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m2 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m2 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m2 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m2 <- function(betas,lambda,alphas) {
  log_likelihood_m2(betas,lambda,alphas) + log_prior_lambda_m2(lambda) + log_prior_betas_m2(betas) + log_prior_alphas_m2(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(24))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
model2_draws <- matrix(nrow=Ndraws,numeric(Ndraws*25))
model2_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model2_draws[i-1,14:25] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m2(model2_draws[i-1,2:13],model2_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m2(model2_draws[i-1,2:13],model2_draws[i-1,1],model2_draws[i-1,14:25])
    model2_draws[i,] <- model2_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model2_draws[i,14:25] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model2_draws[i-1,2:13])^2))
    model2_draws[i,] <- model2_draws[i-1,]
    model2_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model2_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m2(propose_new_vals,model2_draws[i-1,1],model2_draws[i-1,14:25])-log_posterior_m2(model2_draws[i-1,2:13],model2_draws[i-1,1],model2_draws[i-1,14:25])
    model2_draws[i,] <- model2_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model2_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}

accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model2_draws[,15],type="l")

burned_draws_2 <- model2_draws[1001:Ndraws,]
round(apply(model2_draws,2,mean),4)

effectiveSize(burned_draws_2)

# Look at coefficients from a logistic regression model to see if the values are the same for the alphas
greater_than_0 <- as.numeric(forestfire$area>0)
logistic_fit <- glm(greater_than_0~Xmat_pre-1, family="binomial")
summary(logistic_fit)

cbind(round(apply(model2_draws,2,mean),4),round(c(0.4,lm_fit$coefficients,logistic_fit$coefficients),4 ))
hist(model2_draws[,21])

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_2[,14:25],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Use backward selection for the logistic regression coefficients

# Model 5: 
# Eliminate DC #####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]

# Functions for the likelihood and posterior
log_likelihood_m5 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m5 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m5 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m5 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m5 <- function(betas,lambda,alphas) {
  log_likelihood_m5(betas,lambda,alphas) + log_prior_lambda_m5(lambda) + log_prior_betas_m5(betas) + log_prior_alphas_m5(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(23))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))
model5_draws <- matrix(nrow=Ndraws,numeric(Ndraws*24))
model5_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model5_draws[i-1,14:24] + B_logistic%*%rnorm(11)
    log_mh_ratio <- log_posterior_m5(model5_draws[i-1,2:13],model5_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m5(model5_draws[i-1,2:13],model5_draws[i-1,1],model5_draws[i-1,14:24])
    model5_draws[i,] <- model5_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model5_draws[i,14:24] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model5_draws[i-1,2:13])^2))
    model5_draws[i,] <- model5_draws[i-1,]
    model5_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model5_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m5(propose_new_vals,model5_draws[i-1,1],model5_draws[i-1,14:24])-log_posterior_m5(model5_draws[i-1,2:13],model5_draws[i-1,1],model5_draws[i-1,14:24])
    model5_draws[i,] <- model5_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model5_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}
accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model5_draws[,15],type="l")

burned_draws_5 <- model5_draws[1001:Ndraws,]

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_5[,14:24],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Model 6: 
# Eliminate DC, and RH #####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]

# Functions for the likelihood and posterior
log_likelihood_m6 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m6 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m6 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m6 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m6 <- function(betas,lambda,alphas) {
  log_likelihood_m6(betas,lambda,alphas) + log_prior_lambda_m6(lambda) + log_prior_betas_m6(betas) + log_prior_alphas_m6(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(22))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))
model6_draws <- matrix(nrow=Ndraws,numeric(Ndraws*23))
model6_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model6_draws[i-1,14:23] + B_logistic%*%rnorm(10)
    log_mh_ratio <- log_posterior_m6(model6_draws[i-1,2:13],model6_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m6(model6_draws[i-1,2:13],model6_draws[i-1,1],model6_draws[i-1,14:23])
    model6_draws[i,] <- model6_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model6_draws[i,14:23] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model6_draws[i-1,2:13])^2))
    model6_draws[i,] <- model6_draws[i-1,]
    model6_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model6_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m6(propose_new_vals,model6_draws[i-1,1],model6_draws[i-1,14:23])-log_posterior_m6(model6_draws[i-1,2:13],model6_draws[i-1,1],model6_draws[i-1,14:23])
    model6_draws[i,] <- model6_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model6_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}
accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model6_draws[,15],type="l")

burned_draws_6 <- model6_draws[1001:Ndraws,]

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_6[,14:23],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Model 3: 
# Eliminate DC, DMC, and RH #####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - DMC - RH,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]

# Functions for the likelihood and posterior
log_likelihood_m3 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m3 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m3 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m3 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m3 <- function(betas,lambda,alphas) {
  log_likelihood_m3(betas,lambda,alphas) + log_prior_lambda_m3(lambda) + log_prior_betas_m3(betas) + log_prior_alphas_m3(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(21))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))
model3_draws <- matrix(nrow=Ndraws,numeric(Ndraws*22))
model3_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model3_draws[i-1,14:22] + B_logistic%*%rnorm(9)
    log_mh_ratio <- log_posterior_m3(model3_draws[i-1,2:13],model3_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m3(model3_draws[i-1,2:13],model3_draws[i-1,1],model3_draws[i-1,14:22])
    model3_draws[i,] <- model3_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model3_draws[i,14:22] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model3_draws[i-1,2:13])^2))
    model3_draws[i,] <- model3_draws[i-1,]
    model3_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model3_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m3(propose_new_vals,model3_draws[i-1,1],model3_draws[i-1,14:22])-log_posterior_m3(model3_draws[i-1,2:13],model3_draws[i-1,1],model3_draws[i-1,14:22])
    model3_draws[i,] <- model3_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model3_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}
accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model3_draws[,15],type="l")

burned_draws_3 <- model3_draws[1001:Ndraws,]

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_3[,14:22],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Model 7: 
# Eliminate DC, DMC, RH, and rain #####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - DMC - RH - rain,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]

# Functions for the likelihood and posterior
log_likelihood_m7 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m7 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m7 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m7 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m7 <- function(betas,lambda,alphas) {
  log_likelihood_m7(betas,lambda,alphas) + log_prior_lambda_m7(lambda) + log_prior_betas_m7(betas) + log_prior_alphas_m7(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(20))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))
model7_draws <- matrix(nrow=Ndraws,numeric(Ndraws*21))
model7_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model7_draws[i-1,14:21] + B_logistic%*%rnorm(8)
    log_mh_ratio <- log_posterior_m7(model7_draws[i-1,2:13],model7_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m7(model7_draws[i-1,2:13],model7_draws[i-1,1],model7_draws[i-1,14:21])
    model7_draws[i,] <- model7_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model7_draws[i,14:21] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model7_draws[i-1,2:13])^2))
    model7_draws[i,] <- model7_draws[i-1,]
    model7_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model7_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m7(propose_new_vals,model7_draws[i-1,1],model7_draws[i-1,14:21])-log_posterior_m7(model7_draws[i-1,2:13],model7_draws[i-1,1],model7_draws[i-1,14:21])
    model7_draws[i,] <- model7_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model7_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}
accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model7_draws[,15],type="l")

burned_draws_7 <- model7_draws[1001:Ndraws,]

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_7[,14:21],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Model 4:
# Eliminate DC, DMC, RH, rain, and weekend now
#####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - DMC - RH - Weekend - rain,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]

# Functions for the likelihood and posterior
log_likelihood_m4 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m4 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m4 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m4 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m4 <- function(betas,lambda,alphas) {
  log_likelihood_m4(betas,lambda,alphas) + log_prior_lambda_m4(lambda) + log_prior_betas_m4(betas) + log_prior_alphas_m4(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(19))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))
model4_draws <- matrix(nrow=Ndraws,numeric(Ndraws*20))
model4_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model4_draws[i-1,14:20] + B_logistic%*%rnorm(7)
    log_mh_ratio <- log_posterior_m4(model4_draws[i-1,2:13],model4_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m4(model4_draws[i-1,2:13],model4_draws[i-1,1],model4_draws[i-1,14:20])
    model4_draws[i,] <- model4_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model4_draws[i,14:20] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model4_draws[i-1,2:13])^2))
    model4_draws[i,] <- model4_draws[i-1,]
    model4_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model4_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m4(propose_new_vals,model4_draws[i-1,1],model4_draws[i-1,14:20])-log_posterior_m4(model4_draws[i-1,2:13],model4_draws[i-1,1],model4_draws[i-1,14:20])
    model4_draws[i,] <- model4_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model4_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}
accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model4_draws[,15],type="l")

burned_draws_4 <- model4_draws[1001:Ndraws,]

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_4[,14:20],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Model 8:
# Eliminate DC, DMC, RH, rain, weekend, and ISI now
#####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - DMC - RH - ISI - rain - Weekend,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]

# Functions for the likelihood and posterior
log_likelihood_m8 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m8 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m8 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m8 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m8 <- function(betas,lambda,alphas) {
  log_likelihood_m8(betas,lambda,alphas) + log_prior_lambda_m8(lambda) + log_prior_betas_m8(betas) + log_prior_alphas_m8(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(18))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))
model8_draws <- matrix(nrow=Ndraws,numeric(Ndraws*19))
model8_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model8_draws[i-1,14:19] + B_logistic%*%rnorm(6)
    log_mh_ratio <- log_posterior_m8(model8_draws[i-1,2:13],model8_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m8(model8_draws[i-1,2:13],model8_draws[i-1,1],model8_draws[i-1,14:19])
    model8_draws[i,] <- model8_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model8_draws[i,14:19] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model8_draws[i-1,2:13])^2))
    model8_draws[i,] <- model8_draws[i-1,]
    model8_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model8_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m8(propose_new_vals,model8_draws[i-1,1],model8_draws[i-1,14:19])-log_posterior_m8(model8_draws[i-1,2:13],model8_draws[i-1,1],model8_draws[i-1,14:19])
    model8_draws[i,] <- model8_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model8_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}
accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model8_draws[,15],type="l")

burned_draws_8 <- model8_draws[1001:Ndraws,]

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_8[,14:19],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Model 9:
# Eliminate DC, DMC, RH, rain, weekend, ISI, and temp now
#####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - DMC - RH - ISI - rain - Weekend - temp,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
# I change these hyperparameters to check different priors
lambda_hyper <- 2
beta_alpha_hyper <- 0.00001


# Functions for the likelihood and posterior
log_likelihood_m9 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m9 <- function(lambda) {
  dgamma(lambda,lambda_hyper,rate=lambda_hyper,log=TRUE)
}
log_prior_betas_m9 <- function(betas) {
  sum(dnorm(betas,0,1/beta_alpha_hyper,log=TRUE))
}
log_prior_alphas_m9 <- function(alphas) {
  sum(dnorm(alphas,0,1/beta_alpha_hyper,log=TRUE))
}
log_posterior_m9 <- function(betas,lambda,alphas) {
  log_likelihood_m9(betas,lambda,alphas) + log_prior_lambda_m9(lambda) + log_prior_betas_m9(betas) + log_prior_alphas_m9(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(17))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))

set.seed(12)
doModel9 <- function(center_val) {
  model9_draws <- matrix(nrow=Ndraws,numeric(Ndraws*18))
  model9_draws[1,] <- center_val
  
  accept_alpha <- 0
  accept_beta <- 0
  for (i in 2:Ndraws) {
    # Update for alphas using MCMC since it is not conjugate
    if (i %% 3==1) {
      propose_new_alpha_vals <- model9_draws[i-1,14:18] + B_logistic%*%rnorm(5)
      log_mh_ratio <- log_posterior_m9(model9_draws[i-1,2:13],model9_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m9(model9_draws[i-1,2:13],model9_draws[i-1,1],model9_draws[i-1,14:18])
      model9_draws[i,] <- model9_draws[i-1,]
      if (log(runif(1)) < log_mh_ratio) {
        model9_draws[i,14:18] <- propose_new_alpha_vals
        accept_alpha <- accept_alpha + 1
      }
    } else if(i %% 3==2) {
      # Update for lambda using Gibbs sampling
      new_lambda <- rgamma(1,n_without_0/2+lambda_hyper,lambda_hyper + 0.5*sum((log_areas-Xmat%*%model9_draws[i-1,2:13])^2))
      model9_draws[i,] <- model9_draws[i-1,]
      model9_draws[i,1] <- new_lambda
    } else {
      # Metropolis Hastings for the Betas using multivariate draws
      propose_new_vals <- model9_draws[i-1,2:13] + B%*%rnorm(12)
      log_mh_ratio <- log_posterior_m9(propose_new_vals,model9_draws[i-1,1],model9_draws[i-1,14:18])-log_posterior_m9(model9_draws[i-1,2:13],model9_draws[i-1,1],model9_draws[i-1,14:18])
      model9_draws[i,] <- model9_draws[i-1,]
      if (log(runif(1)) < log_mh_ratio) {
        model9_draws[i,2:13] <- propose_new_vals
        accept_beta <- accept_beta + 1
      }
    }
  }
  alpha_accept_rate <- accept_alpha/(Ndraws/3)
  beta_accept_rate <- accept_beta/(Ndraws/3)
  
  burned_draws_9 <- model9_draws[1001:Ndraws,]
  list(alpha_accept_rate,beta_accept_rate,burned_draws_9,model9_draws)
}
set.seed(12)
m9_moderate_center <- doModel9(center)
m9_high_center <- doModel9(c(.001,rep(0.1,12),rep(0,5)))
m9_low_center <- doModel9(c(3,rep(-0.1,12),rep(0,5)))
# Different priors 
# lambda 0.1, beta 1
prior2 <- doModel9(center)
# expected values 
prior2_means <- apply(prior2[[3]],2,mean)
prior3 <- doModel9(center)
prior3_means <- apply(prior3[[3]],2,mean)

# trace plots
par(mfrow=c(2,2))
plot(1:Ndraws,m9_low_center[[4]][,1],type="l",main="Trace Plot for Lambda",ylab="Posterior Draw")
lines(1:Ndraws,m9_high_center[[4]][,1],col="red")
lines(1:Ndraws,m9_moderate_center[[4]][,1],col="blue")
plot(1:Ndraws,m9_low_center[[4]][,3],type="l",main="Trace Plot for Beta1",ylab="Posterior Draw")
lines(1:Ndraws,m9_high_center[[4]][,3],col="red")
lines(1:Ndraws,m9_moderate_center[[4]][,3],col="blue")
plot(1:Ndraws,m9_low_center[[4]][,6],type="l",main="Trace Plot for Beta4",ylab="Posterior Draw")
lines(1:Ndraws,m9_high_center[[4]][,6],col="red")
lines(1:Ndraws,m9_moderate_center[[4]][,6],col="blue")
plot(1:Ndraws,m9_low_center[[4]][,15],type="l",main="Trace Plot for Alpha1",ylab="Posterior Draw")
lines(1:Ndraws,m9_high_center[[4]][,15],col="red")
lines(1:Ndraws,m9_moderate_center[[4]][,15],col="blue")
par(mfrow=c(1,1))

# posterior distribution trace plot
posterior_values <- apply(m9_moderate_center[[4]],1, function(x) log_posterior_m9(x[2:13],x[1],x[14:18]))
par(mfrow=c(1,2))
plot(1:Ndraws,posterior_values,type="l",main="Without Burning",ylab="log posterior of draws")
plot(1:49000,posterior_values[1001:50000],type="l",main="Burning",ylab="log posterior of draws")
par(mfrow=c(1,1))
# Check to see if I'm getting similar results
par(mfrow=c(1,2))
acf(m9_moderate_center[[3]][,1],main="ACF for Lambda Draws")
acf(m9_moderate_center[[3]][,3],main="ACF for Beta1 Draws")
plot(1:49000,m9_moderate_center[[3]][,1],type="l",main="Trace Plot for Lambda",ylab="draws")
plot(1:49000,m9_moderate_center[[3]][,3],type="l",main="Trace Plot for Beta1",ylab="draws")
par(mfrow=c(1,1))
burned_draws_9 <- m9_moderate_center[[3]]
effect_sizes_m9 <- t(effectiveSize(burned_draws_9))
names_of_col <- c("Lambda","Beta0","Beta1","Beta2","Beta3","Beta4","Beta5","Beta6","Beta7","Beta8","Beta9","Beta10","Beta11","Alpha0","Alpha1","Alpha2","Alpha3","Alpha4")
colnames(effect_sizes_m9) <- names_of_col
rownames(effect_sizes_m9) <- "Effective Sample Size"
m9_moderate_center[[1]]
m9_moderate_center[[2]]
xtable(effect_sizes_m9)

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_9[,14:18],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# Model 10:
# Eliminate DC, DMC, RH, rain, weekend, ISI, temp, and Season now
#####

# Prepare the data for the analysis
Xmat_pre <- model.matrix(~. -X -Y -area -month - day,data=forestfire)
rows_with_0 <- which(forestfire$area==0)
Xmat_without_0 <- Xmat_pre[-rows_with_0,]
Xmat <- Xmat_without_0
Xmat_with_0 <- Xmat_pre[rows_with_0,]
n_with_0 <- length(rows_with_0)
n_without_0 <- nrow(forestfire) - n_with_0
areas_pre <- forestfire$area
areas <- areas_pre[-rows_with_0]
log_areas <- log(areas)

Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - DMC - RH - ISI - rain - Weekend - temp - Season,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]

# Functions for the likelihood and posterior
log_likelihood_m10 <- function(betas,lambda,alphas) {
  pi_not_0 <- 1 / (1+exp(-Xmat_without_0_logistic%*%alphas))
  pi_is_0 <- exp(-Xmat_with_0_logistic%*%alphas) / (1+exp(-Xmat_with_0_logistic%*%alphas))
  sum(log(pi_is_0)) + sum(log(pi_not_0)) + sum(dlnorm(areas,Xmat%*%betas,1/sqrt(lambda),log=TRUE))
}
log_prior_lambda_m10 <- function(lambda) {
  dgamma(lambda,2,rate=2,log=TRUE)
}
log_prior_betas_m10 <- function(betas) {
  sum(dnorm(betas,0,1/.001,log=TRUE))
}
log_prior_alphas_m10 <- function(alphas) {
  sum(dnorm(alphas,0,1/.001,log=TRUE))
}
log_posterior_m10 <- function(betas,lambda,alphas) {
  log_likelihood_m10(betas,lambda,alphas) + log_prior_lambda_m10(lambda) + log_prior_betas_m10(betas) + log_prior_alphas_m10(alphas)
}

# MCMC
# Do multivariate draws for the betas
# Gibbs sampler for pi
Ndraws <- 50000
center <- c(.01,numeric(15))
# Covariance matrix for multivariate beta draws
cov_draws <- 1.5 * solve(t(Xmat)%*%Xmat)
B <- t(chol(cov_draws))
cov_draws_logistic <- 1.5 * solve(t(Xmat_without_0_logistic)%*%Xmat_without_0_logistic)
B_logistic <- t(chol(cov_draws_logistic))
model10_draws <- matrix(nrow=Ndraws,numeric(Ndraws*16))
model10_draws[1,] <- center

accept_alpha <- 0
accept_beta <- 0
set.seed(12)
for (i in 2:Ndraws) {
  # Update for alphas using MCMC since it is not conjugate
  if (i %% 3==1) {
    propose_new_alpha_vals <- model10_draws[i-1,14:16] + B_logistic%*%rnorm(3)
    log_mh_ratio <- log_posterior_m10(model10_draws[i-1,2:13],model10_draws[i-1,1],propose_new_alpha_vals)-log_posterior_m10(model10_draws[i-1,2:13],model10_draws[i-1,1],model10_draws[i-1,14:16])
    model10_draws[i,] <- model10_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model10_draws[i,14:16] <- propose_new_alpha_vals
      accept_alpha <- accept_alpha + 1
    }
  } else if(i %% 3==2) {
    # Update for lambda using Gibbs sampling
    new_lambda <- rgamma(1,n_without_0/2+2,2+0.5*sum((log_areas-Xmat%*%model10_draws[i-1,2:13])^2))
    model10_draws[i,] <- model10_draws[i-1,]
    model10_draws[i,1] <- new_lambda
  } else {
    # Metropolis Hastings for the Betas using multivariate draws
    propose_new_vals <- model10_draws[i-1,2:13] + B%*%rnorm(12)
    log_mh_ratio <- log_posterior_m10(propose_new_vals,model10_draws[i-1,1],model10_draws[i-1,14:16])-log_posterior_m10(model10_draws[i-1,2:13],model10_draws[i-1,1],model10_draws[i-1,14:16])
    model10_draws[i,] <- model10_draws[i-1,]
    if (log(runif(1)) < log_mh_ratio) {
      model10_draws[i,2:13] <- propose_new_vals
      accept_beta <- accept_beta + 1
    }
  }
}
accept_alpha/(Ndraws/3)
accept_beta/(Ndraws/3)
# Trace plot
plot(1:Ndraws,model10_draws[,15],type="l")

burned_draws_10 <- model10_draws[1001:Ndraws,]

# Find logistic variables that have the lowest "p-value"
apply(burned_draws_10[,14:16],2,function(x) min(mean(x>0),1-mean(x>0)))
#####

# DIC for the models #####
# Model 1
posterior_mean1 <- apply(burned_draws_1,2,mean)
pdic1 <- 2 * (log_likelihood(posterior_mean1[3:14],posterior_mean1[2],posterior_mean1[1]) - mean( apply(burned_draws_1,1,function(x) log_likelihood(x[3:14],x[2],x[1]))))
dic1 <- -2*log_likelihood(posterior_mean1[3:14],posterior_mean1[2],posterior_mean1[1]) + 2*pdic1
# Model 2
posterior_mean2 <- apply(burned_draws_2,2,mean)
pdic2 <- 2 * (log_likelihood_m2(posterior_mean2[2:13],posterior_mean2[1],posterior_mean2[14:25]) - mean( apply(burned_draws_2,1,function(x) log_likelihood_m2(x[2:13],x[1],x[14:25]))))
dic2 <- -2*log_likelihood_m2(posterior_mean2[2:13],posterior_mean2[1],posterior_mean2[14:25]) + 2*pdic2
# Model 5
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean5 <- apply(burned_draws_5,2,mean)
pdic5 <- 2 * (log_likelihood_m5(posterior_mean5[2:13],posterior_mean5[1],posterior_mean5[14:24]) - mean( apply(burned_draws_5,1,function(x) log_likelihood_m5(x[2:13],x[1],x[14:24]))))
dic5 <- -2*log_likelihood_m5(posterior_mean5[2:13],posterior_mean5[1],posterior_mean5[14:24]) + 2*pdic5
# Model 6
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean6 <- apply(burned_draws_6,2,mean)
pdic6 <- 2 * (log_likelihood_m6(posterior_mean6[2:13],posterior_mean6[1],posterior_mean6[14:23]) - mean( apply(burned_draws_6,1,function(x) log_likelihood_m6(x[2:13],x[1],x[14:23]))))
dic6 <- -2*log_likelihood_m6(posterior_mean6[2:13],posterior_mean6[1],posterior_mean6[14:23]) + 2*pdic6
# Model 3
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH - DMC,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean3 <- apply(burned_draws_3,2,mean)
pdic3 <- 2 * (log_likelihood_m3(posterior_mean3[2:13],posterior_mean3[1],posterior_mean3[14:22]) - mean( apply(burned_draws_3,1,function(x) log_likelihood_m3(x[2:13],x[1],x[14:22]))))
dic3 <- -2*log_likelihood_m3(posterior_mean3[2:13],posterior_mean3[1],posterior_mean3[14:22]) + 2*pdic3
# Model 7
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH - DMC - rain,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean7 <- apply(burned_draws_7,2,mean)
pdic7 <- 2 * (log_likelihood_m7(posterior_mean7[2:13],posterior_mean7[1],posterior_mean7[14:21]) - mean( apply(burned_draws_7,1,function(x) log_likelihood_m7(x[2:13],x[1],x[14:21]))))
dic7 <- -2*log_likelihood_m7(posterior_mean7[2:13],posterior_mean7[1],posterior_mean7[14:21]) + 2*pdic7
# Model 4
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH - DMC - rain - Weekend,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean4 <- apply(burned_draws_4,2,mean)
pdic4 <- 2 * (log_likelihood_m4(posterior_mean4[2:13],posterior_mean4[1],posterior_mean4[14:20]) - mean( apply(burned_draws_4,1,function(x) log_likelihood_m4(x[2:13],x[1],x[14:20]))))
dic4 <- -2*log_likelihood_m4(posterior_mean4[2:13],posterior_mean4[1],posterior_mean4[14:20]) + 2*pdic4
# Model 8
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH - DMC - rain - Weekend - ISI,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean8 <- apply(burned_draws_8,2,mean)
pdic8 <- 2 * (log_likelihood_m8(posterior_mean8[2:13],posterior_mean8[1],posterior_mean8[14:19]) - mean( apply(burned_draws_8,1,function(x) log_likelihood_m8(x[2:13],x[1],x[14:19]))))
dic8 <- -2*log_likelihood_m8(posterior_mean8[2:13],posterior_mean8[1],posterior_mean8[14:19]) + 2*pdic8
# Model 9
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH - DMC - rain - Weekend - ISI - temp,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean9 <- apply(burned_draws_9,2,mean)
pdic9 <- 2 * (log_likelihood_m9(posterior_mean9[2:13],posterior_mean9[1],posterior_mean9[14:18]) - mean( apply(burned_draws_9,1,function(x) log_likelihood_m9(x[2:13],x[1],x[14:18]))))
dic9 <- -2*log_likelihood_m9(posterior_mean9[2:13],posterior_mean9[1],posterior_mean9[14:18]) + 2*pdic9
# Model 10
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH - DMC - rain - Weekend - ISI - temp - Season,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
posterior_mean10 <- apply(burned_draws_10,2,mean)
pdic10 <- 2 * (log_likelihood_m10(posterior_mean10[2:13],posterior_mean10[1],posterior_mean10[14:16]) - mean( apply(burned_draws_10,1,function(x) log_likelihood_m10(x[2:13],x[1],x[14:16]))))
dic10 <- -2*log_likelihood_m10(posterior_mean10[2:13],posterior_mean10[1],posterior_mean10[14:16]) + 2*pdic10

# Table of DICs
dic_tab <- t(t(c("All Variables"=dic2,
             "Eliminate DC"=dic5,
             "RH"=dic6,
             "DMC"=dic3,
             "rain"=dic7,
             "weekend"=dic4,
             "ISI"=dic8,
             "temp"=dic9,
             "season"=dic10,
             "none"=dic1)))
colnames(dic_tab) <- "DIC"

xtable(dic_tab)
dic_tab
#####


# I will use model 9 because it has the lowest DIC out of all the models

## Test general goodness-of-fit in 'residuals'
test1 <- function(observed, expected) {  # Test general goodness-of-fit in 'residuals'
  sum( (observed - expected)^2 / expected )
}

# Model 9
Xmat_pre_logistic <- model.matrix(~. -X -Y -area -month - day - DC - RH - DMC - rain - Weekend - ISI - temp,data=forestfire)
Xmat_without_0_logistic <- Xmat_pre_logistic[-rows_with_0,]
Xmat_with_0_logistic <- Xmat_pre_logistic[rows_with_0,]
set.seed(15)
post_pred_dat_9 <- lapply(1:nrow(burned_draws_9), function(i) cbind(rbinom(517,1,1 / (1 + exp(-Xmat_pre_logistic%*%burned_draws_9[i,14:18]))) * rlnorm(517,Xmat_pre %*% burned_draws_9[i,2:13],1/(sqrt(burned_draws_9[i,1]))), Xmat_pre))
## Test general goodness-of-fit in 'residuals'

expected_ys <- lapply(1:nrow(burned_draws_9), function(i)  1 / (1 + exp(-Xmat_pre_logistic%*%burned_draws_9[i,14:18])) * exp(Xmat_pre%*%burned_draws_9[i,2:13] + 0.5/burned_draws_9[i,1]))

statisticFromPredictive <- sapply(1:nrow(burned_draws_9), function(i) test1(post_pred_dat_9[[i]][,1], expected_ys[[i]] ))
statisticFromObserved   <- sapply(1:nrow(burned_draws_9), function(i) test1(forestfire$area, expected_ys[[i]] ))
p.value1 <- mean(statisticFromPredictive >= statisticFromObserved)

# residual plot for only the data that wasn't 0
estimated_parms <- apply(burned_draws_9,2,mean)
estimated_parms_quant <- t(apply(burned_draws_9,2,quantile,prob=c(0.025,0.975)))

standardized_residuals <- log(areas)-Xmat%*%estimated_parms[2:13] / (1/sqrt(estimated_parms[1]))
plot(Xmat%*%estimated_parms[2:13],standardized_residuals,pch=19,xlab="Fitted Values",ylab="Standardized Residuals",main="Plot of the Residuals")
abline(h=0)

# Plot to show linearity
plot(~., data = datnot0)

estimated_parms_mat <- t(t(estimated_parms))
estimated_parms_mat <- cbind(estimated_parms_mat,estimated_parms_quant)
rownames(estimated_parms_mat) <- names_of_col
colnames(estimated_parms_mat) <- c("Expected Value","2.5%","97.5%")

# plot for all the variables
par(mfrow=c(2,3))
hist(burned_draws_9[,1])
hist(burned_draws_9[,2])
hist(burned_draws_9[,3])
abline(v=0,col="red")
hist(burned_draws_9[,4])
abline(v=0,col="red")
hist(burned_draws_9[,5])
abline(v=0,col="red")
hist(burned_draws_9[,6])
abline(v=0,col="red")
hist(burned_draws_9[,7])
abline(v=0,col="red")
hist(burned_draws_9[,8])
abline(v=0,col="red")
hist(burned_draws_9[,9])
abline(v=0,col="red")
hist(burned_draws_9[,10])
abline(v=0,col="red")
hist(burned_draws_9[,11])
abline(v=0,col="red")
hist(burned_draws_9[,12])
abline(v=0,col="red")
hist(burned_draws_9[,13])
abline(v=0,col="red")
hist(burned_draws_9[,14])
hist(burned_draws_9[,15])
abline(v=0,col="red")
hist(burned_draws_9[,16])
abline(v=0,col="red")
hist(burned_draws_9[,17])
abline(v=0,col="red")
hist(burned_draws_9[,18])
abline(v=0,col="red")
par(mfrow=c(1,1))

cbind(1:18,apply(burned_draws_9,2,function(x) min(mean(x>0),1-mean(x>0))))

par(mfrow=c(2,3))
hist(burned_draws_9[,4],main=NULL,xlab="Beta DMC")
abline(v=0,col="red")
hist(burned_draws_9[,6],main=NULL,xlab="Beta ISI")
abline(v=0,col="red")
hist(burned_draws_9[,11],main=NULL,xlab="Beta Weekend1")
abline(v=0,col="red")
hist(burned_draws_9[,12],main=NULL,xlab="Beta Season2")
abline(v=0,col="red")
hist(burned_draws_9[,16],main=NULL,xlab="Alpha Wind")
abline(v=0,col="red")
hist(burned_draws_9[,18],main=NULL,xlab="Alpha Season3")
abline(v=0,col="red")
par(mfrow=c(1,1))




# Frequentist analysis
lm_fit <- lm(log_areas~Xmat -1)
summary(lm_fit)
# Logistic regression
ind_greater_than_0 <- as.numeric(forestfire$area>0)
glm_fit <- glm(ind_greater_than_0~Xmat_pre_logistic-1,family="binomial")
summary(glm_fit)
freq_parm <-c(1/1.513^2,lm_fit$coefficients,glm_fit$coefficients)

# Table for the results
the_print <- cbind(estimated_parms_mat,freq_parm,prior2_means,prior3_means)

xtable(the_print,digits=4)
