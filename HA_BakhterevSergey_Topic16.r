#Simulating Geometric Brownian motion with mean-reversion
gbm_mr = function(x_0, alpha, mu, sigma, len, n){
  x = rep(0, n) #x will contain values of stochastic process
  # x[t] - value of stochastic process at time (t-1)*len/n
  x[1] = x_0 #initializing process at time zero
  dt = len/n #difference in time
  for (i in (2:n)){
    dx = (alpha*(mu - x[i-1]) - 1/2*sigma^2)*dt + sigma*sqrt(dt)*rnorm(1, 0, 1)
    x[i] = x[i-1] + dx
  }
  return = x
}

#Plot of the simulation
t = seq(0, 0.25, 0.25/365)
length(t)
plot(t, gbm_mr(100, 100, 100, 0.3, 0.25, 366), type="l", lty=1, main="Mean-reversion Paths", ylab="x_t")


#Simulating Geometric brownian motion with mean-reversion and stochastic volatility
gbm_sv = function(x_0, mu, alpha,len, n, a, avg_v, ksi, rho){
  dt = len/n
  x = rep(0, n)
  v = rep(0, n)
  x[1] = x_0
  v[0] = avg_v
  eps_1 = rnorm(n, 0, 1)
  eps_2 = rnorm(n, 0, 1)
  for (i in (2:n)){
    dv = a * (avg_v - v[i-1])*dt + ksi*sqrt(v[i-1]*dt)*(rho*eps_1[i] + sqrt(1 - rho^2)*eps_2[i])
    v[i] = v[i-1] + dv
    dx = (alpha*(mu - x[i-1]) - 1/2*v[i]**2)*dt + (v[i])*sqrt(dt)*rnorm(1, 0, 1)
    x[i] = x[i-1] + dx
  }
  return = x
}


#Plot of the simulation
t = seq(0, 0.25, 0.25/365)
plot(t, gbm_sv(100, 100, 100, 0.25, 366, 10, 0.09, 0.5, 0.0), type="l",
     lty=1, main="Stochastic Volatility Paths", ylab="x_t")




#Simulating Geometric brownian motion with jumps
gbm_jumps = function(x_0, r, phi, kappa, sigma, gamma, len, n){
  dt = len/n
  x = rep(0, n)
  x[1] = x_0
  eps_1 = rnorm(n, 0, 1)
  eps_2 = rnorm(n, 0, 1)
  u = rnorm(n, 0, 1)
  for (i in (2:n)){
    dx = (r - phi*kappa - 1/2*sigma**2)*dt + sigma*sqrt(dt)*eps_1[i] +
      (kappa + gamma*eps_2[i])*as.integer(u[i] > phi*dt)
    x[i] = x[i-1] + dx
  }
  return = x
}


#Plot of the simulation
t = seq(0, 0.25, 0.25/(365))
plot(t, gbm_jumps(100, 0.03, 100, 0, 0.3, 0.02, 0.25, 366), type="l",
     lty=1, main="GBM with jumps", ylab="x_t")


#Calibration Geometric Brownian motion with mean-reversion
y = gbm_mr(100, 100, 100, 0.3, 0.25, 366)

plot(t, y,  type="l",
     lty=1, main="GBM with mean-reversion", ylab="x_t")
n = length(y)

len = 0.25
calibrate_gbm_mr = function(y, len, n){
  dt = len/n
  y_0 = y[1]
  y_i = y[1:(length(y)-1)]
  y_i1 = y[2:(length(y))]
  b = lm(y_i1 ~ y_i)$coefficient[1]
  a = lm(y_i1 ~ y_i)$coefficient[2]
  alpha = - log(a)/dt
  avg_y = b/(1-a)
  sigma = sqrt(anova(lm(y_i1 ~ y_i))[2 , 3])*sqrt(-2*log(a)/(dt*(1-a^2)))
  
  return = c(y_0, alpha, avg_y, sigma)
}

v4 = calibrate_gbm_mr(y, len, n)
  
lines(t, gbm_mr(v4[1], v4[2], v4[3], v4[4], len, 366), col = "red") 

#Calibration GBM with stochastic volatility

u = gbm_sv(100, 100, 100, 0.25, 366, 10, 0.09, 0.5, 0.0)

t2 = seq(0, 0.25, 0.25/365)
n = length(u)
len = 0.25
gbm_sv_garch = function(x_0, alpha, mu, v, len, n){
  x = rep(0, n)
  x[1] = x_0
  dt = len/n
  for (i in (2:n)){
    dx = (alpha*(mu - x[i-1]) - 1/2*v[i]**2)*dt + v[i]*sqrt(dt)*rnorm(1, 0, 1)
    x[i] = x[i-1] + dx
  }
  return = x
}
calibrate_gbm_sv = function(u, len, n){
  dt = len/n
  u_0 = u[1]
  u.garchFit=garchFit(~ garch(1,1), data = u, trace=F)
  v1 = sqrt(volatility(u.garchFit))
  u_i = u[1:(length(u)-1)]
  u_i1 = u[2:(length(u))]
  b = lm(u_i1 ~ u_i)$coefficient[1]
  a = lm(u_i1 ~ u_i)$coefficient[2]
  alpha = - log(a)/dt
  mu = b/(1-a)
  return = list(c(u_0, alpha, mu), v1)
}

v3 = calibrate_gbm_sv(u, len, n)

plot(t2, u, type = "l",
     lty=1, main="GBM with stochastic volatility", ylab="x_t")

u_est = gbm_sv_garch(v3[[1]][1], v3[[1]][2], v3[[1]][3], v3[[2]], 0.25, 366)
lines(t2, u_est, col = "red")


#Calibration of GBM with jumps
z = gbm_jumps(100, 0.05, 100, 0, 0.3, 0.02, 0.25, 366)
z_0 = z[1]
n = length(z)

t = seq(0, 0.25, length = 366)
plot(t, z,  type="l",
     lty=1, main="GBM with jumps", ylab="x_t")



calibrate_gbm_jumps = function(price, len, n){
  price_0 = price[1]
  dprice = rep(0, n)
  dt = len/n
  for (i in 2:n){
    dprice[i] = price[i] - price[i-1]
  }
  
  diff = 1
  count = 0
  count1 = 0
  count2 = 0
  price_new = price
  jumps2 = rep(0, n)
  while(diff != 0){
    jump_detect = sd(dprice, na.rm = T)*3
    count2 = count2 + count
    count1 = count
    count = 0
    jumps = rep(0,n)
    for (j in (1 : length(dprice))){
      if (abs(dprice[j]) > jump_detect){
        count = count + 1
        jumps[j] = dprice[j]
        price_new = price_new[- price_new[j]]
      }
      
    }
    jump = jumps[jumps != 0]
    jumps2 = append(jumps2, jump)
    dprice = dprice[! dprice %in% jump]
    diff = count1 - count
  }
  jumps2 = jumps2[jumps2 != 0]
  phi = count2/len
  gamma = sd(jumps2, na.rm = T)
  
  
  price_i = price_new[1:(length(price_new)-1)]
  price_i1 = price_new[2:(length(price_new))]
  a = lm(price_i1 ~ price_i)$coefficient[2]
  dt = len/n
  alpha = - log(a)/dt
  sigma = sqrt(anova(lm(price_i1 ~ price_i))[2 , 3])*sqrt(-2*log(a)/(len*(1-a**2)))
  const = rep(1, length(dprice))
  r = 1/2*sigma**2 + lm(dprice ~ const)$coefficient[1]/dt
  return = c(price_0, r, phi, gamma, sigma)
}

v = calibrate_gbm_jumps(z,len,n)

plot(t, z,  type="l",
     lty=1, main="GBM with jumps", ylab="x_t")
lines(t, gbm_jumps(v[1], v[2], v[3], 0, v[5], v[4], 0.25, 366), col = "red")


#gbm with jumps and seasonality
gbm_jumps_sv = function(x_0, r, phi, kappa, sigma, gamma, len, n){
  dt = len/n
  x = rep(0, n)
  x[1] = x_0
  eps_1 = rnorm(n, 0, 1)
  eps_2 = rnorm(n, 0, 1)
  u = rnorm(n, 0, 1)
  for (i in (2:n)){
    dx = (r - phi*kappa - 1/2*sigma[i]**2)*dt + sigma[i]*sqrt(dt)*eps_1[i] +
      (kappa + gamma*eps_2[i])*as.integer(u[i] > phi*dt)
    x[i] = x[i-1] + dx
  }
  return = x
}

calibrate_gbm_jumps_sv = function(price, len, n){
  price_0 = price[1]
  dprice = rep(0, n)
  dt = len/n
  for (i in 2:n){
    dprice[i] = price[i] - price[i-1]
  }
  
  diff = 1
  count = 0
  count1 = 0
  count2 = 0
  price_new = price
  jumps2 = rep(0, n)
  while(diff != 0){
    jump_detect = sd(dprice, na.rm = T)*3
    count2 = count2 + count
    count1 = count
    count = 0
    jumps = rep(0,n)
    for (j in (1 : length(dprice))){
      if (abs(dprice[j]) > jump_detect){
        count = count + 1
        jumps[j] = dprice[j]
        price_new = price_new[- price_new[j]]
      }
      
    }
    jump = jumps[jumps != 0]
    jumps2 = append(jumps2, jump)
    dprice = dprice[! dprice %in% jump]
    diff = count1 - count
  }
  jumps2 = jumps2[jumps2 != 0]
  phi = count2/len
  gamma = sd(jumps2, na.rm = T)
  
  
  price_i = price_new[1:(length(price_new)-1)]
  price_i1 = price_new[2:(length(price_new))]
  a = lm(price_i1 ~ price_i)$coefficient[2]
  dt = len/n
  price.garchFit=garchFit(~ garch(1,1), data = price_new, trace=F)
  sigma = sqrt(volatility(price.garchFit))
  sigma_avg = mean(sigma)
  const = rep(1, length(dprice))
  r = 1/2*sigma_avg**2 + lm(dprice ~ const)$coefficient[1]/dt
  return = list(c(price_0, r, phi, gamma), sigma)
}


sigma = rnorm(366, 0.126232, 0.01)

season = gbm_jumps_sv(100, 0.05, 100, 0, sigma, 0.02, 0.25, 366)
t = seq(0, 0.25, length = 366)

plot(t, season,  type="l",
     lty=1, main="GBM with jumps and seasonality ", ylab="x_t")

v2 = calibrate_gbm_jumps_sv(season, 0.25, 366)
lines(t, gbm_jumps_sv(v2[[1]][1], v2[[1]][2], v2[[1]][3], 0, v2[[2]], v2[[1]][4], 0.25, 366), col = "red")



#Real data
#Spot prices for electricity prices in USA 2019
#SP15 EZ Gen DA LMP Peak

install.packages("readxl", dep = T)
library(readxl)
data = read_excel("C:/Users/defaultuser0/R/ice_electric-2019final.xlsx")
log_price = log(data$`Wtd avg price $/MWh`[1055:1300])
n = length(log_price)
t1 = seq(0, 1, length = n)
plot(t1, log_price, type = "l")


coef = calibrate_gbm_mr(log_price, 1, n)
#Simulating by GBM with mean-reversion
plot(t1, log_price, type = "l")
lines(t1, gbm_mr(coef[1], coef[2], coef[3], coef[4],1, n), col = "red")

t_append = seq(0, 1.25, length = as.integer(5*n/4))
#forecasting for half of the year by GBM with mean-reversion
forecast = append(log_price, gbm_mr(log_price[n], coef[2], coef[3], coef[4], 0.25, n/4))
plot(t_append, forecast, type = "l")


#coefficients for other models
coef2 = calibrate_gbm_sv(log_price, 1, n)
coef3 = calibrate_gbm_jumps(log_price, 1, n)
coef4 = calibrate_gbm_jumps_sv(log_price, 1, n)
