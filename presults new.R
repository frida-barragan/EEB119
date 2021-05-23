library(deSolve)

IC <- c(
  N = 3490,
  S = 1000,
  E = 1000,
  I1 = 100,
  I2 = 50,
  I3 = 340,
  R = 0 #recovered must be 0 initially
) #Fill in the initial values

times = seq(0,20, 1) #Create a sequence of times

pars = c(
  a = 0.095,
  a1 = 0.065,
  a2 = 0.018,
  p = 0.4,
  mu = 0.05,
  z = 0.35,
  x = 0.25,
  e1 = 0.1,
  e2 = 0.25,
  k1 = 0.5,
  k2 = 0.75,
  #i directly changed 1-k to 0.5 in the equations above
  gamma = 0.6,
  alpha = 0.28,
  beta2 = 0.05,
  beta3 = 0.03,
  D = 0.005,
  D0 = 0.01,
  rho = 0.01,
  rho1 = 0.005
)

SIRODE <- function(t, vars, pars) {
  with(as.list(c(vars,pars)), { #This let's you use S, I, beta, and gamma directly 
    N = S + E+ I1 + I2 + I3 + R
    
    dSdt = (1 - p)*mu*N + (1-z)*mu*N + (1-x)*mu*N- D*S - ((e1 + k1)*I1*S)/N - ((e2 + k2)*I2*S)/N - alpha*S - gamma*S
    
    dEdt = p*mu*N + z*mu*N + x*mu*N - D*E - a*E + ((e1 + k1)*I1*S)/N +  ((e2 + k2)*I2*S)/N + (0.5)*R + alpha*S + gamma*S
    
    dI1dt = -D*I1 +  a*E - a1*I1 - (rho*I1*S)/N + beta2*I2
    
    dI2dt = -D0*I2 + a1*I1 - a2*I2 - beta2*I2 + beta3*I3
    
    dI3dt = -(D0*I3) + a2*I2 - beta3*I3
    
    dRdt = -(D*R) + (rho1*I1*S)/N - (0.5)*R
    list(c(N,dSdt,dEdt,dI1dt,dI2dt,dI3dt,dRdt))
  }) }
output <- lsoda(IC, times, SIRODE, pars) 
plot(output[,1],output[,2],type="l", xlab = "time", ylab="population")
lines(output[,1],output[,3],col="blue") #Infected
lines(output[,1],output[,4],col="red")
lines(output[,1],output[,5],col="green")
lines(output[,1],output[,6],col="yellow")
lines(output[,1],output[,7],col="purple")
lines(output[,1],output[,8],col="black")

legend("left", c("total population","susceptibles", "exposed", "overweight", "obese", "extremely obese", "recovered"),
       col = c("blue", "red", "green", "yellow", "purple", "black", "pink"), lty = 1, bty = "n")
