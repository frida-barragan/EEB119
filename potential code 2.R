library(deSolve)

IC <- c(
  N = 3650,
  S = 1000,
  E = 1000,
  I1 = 500,
  I2 = 100,
  I3 = 1050,
  R = 0 #recovered must be 0 initially
) #Fill in the initial values

times = seq(0,10, 1) #Create a sequence of times

pars = c(
  a = 0.01,
  a1 = 0.005,
  a2 = 0.018,
  p = 0.4,
  mu = 0.05,
  z = 0.35,
  x = 0.25,
  e1 = 0.001,
  e2 = 0.25,
  k1 = 0.05,
  k2 = 0.0075,
  #i directly changed 1-k to 0.5 in the equations below
  gamma = 0.6,
  alpha = 0.0028,
  beta2 = 0.0005,
  beta3 = 0.0003,
  D = 0.05,
  D0 = 0.56,
  rho = 0.001,
  rho1 = 0.005
)

SIRODE <- function(t, vars, pars) {
  with(as.list(c(vars,pars)), { #This let's you use S, I, beta, and gamma directly 
    dNdt = (dSdt + dEdt+ dI1dt + dI2dt + dI3dt + dRdt)
    
    dSdt = ((1 - p)*mu*dNdt + (1-z)*mu*dNdt+ (1-x)*mu*dNdt- D*S - ((e1 * k1)*I1*S)/dNdt - ((e2 * k2)*I2*S)/dNdt - alpha*S - gamma*S)
    
    dEdt = (p*mu*dNdt + z*mu*dNdt + x*mu*dNdt - D*E - a*E + ((e1* k1)*I1*S)/dNdt +  ((e2 * k2)*I2*S)/dNdt + (0.95)*R + alpha*S + gamma*S)
    
    dI1dt = (-D*I1 +  a*E - a1*I1 - (rho*I1*S)/dNdt + beta2*I2)
    
    dI2dt = (-D0*I2 + a1*I1 - a2*I2 - beta2*I2 + beta3*I3)
    
    dI3dt = (-(D0*I3) + a2*I2 - beta3*I3)
    
    dRdt = (-(D*R) + (0.95)*I1*S/dNdt - (rho1)*R)
    list(c(dNdt,dSdt,dEdt,dI1dt,dI2dt,dI3dt,dRdt))
  }) }
output <- lsoda(IC, times, SIRODE, pars) 
output
plot(output[,1],output[,2],type="l", xlab = "time", ylab="population", col = "pink")
lines(output[,1],output[,3],col="blue") #Infected
lines(output[,1],output[,4],col="red")
lines(output[,1],output[,5],col="green")
lines(output[,1],output[,6],col="yellow")
lines(output[,1],output[,7],col="purple")
lines(output[,1],output[,8],col="black")

legend("left", c("total population","susceptibles", "exposed", "overweight", "obese", "extremely obese", "recovered"),
       col = c("pink","blue", "red", "green", "yellow", "purple", "black"), lty = 1, bty = "n")