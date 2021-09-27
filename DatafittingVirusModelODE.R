library(simecol)
library(dplyr)
library(sfsmisc) # for displaying better scale labels

"
save your csv file in the following format:
Patient |     V      | Day
A       |    10^8    |  7
A       |    10^7    |  10
A       |    10^5    |  16
B       |    10^8    |  8
B       |    10^2    |  10
B       |    10^3    |  15
etc.
V stands for viral load, might also be A for antibodies
First line is the header
Fitting can also be performed for Antibody compartment - just replace all $V with $A e.g. yobs=data.frame(A=c(patient$A))
"

corona = new("odeModel",
               main=function(time, init, parms){
                 with(as.list(c(init, parms)),{
                   if(V < 1){b1 = 0} else {b1 = b}
                   if(V < 1){p1 = 0} else {p1 = p}
                   dV = p*I-y*V*A-b*U*V
                   dU = r*U*(1-(U/Umax))-b*U*V
                   dL = b*U*V-k*L
                   dI = k*L-b*I*C
                   dC = w*I-f*C
                   dA = x*C-q*A
                   list(c(dV, dU, dL, dI, dC, dA))})},
             parms=c(Umax = 4e8,
                     p = 9,
                     r = 0.5,
                     k = 0.8,
                     b = 4e-8,
                     y = 7e-2,
                     z = 5e-8,
                     w = 0.09,
                     x = 2e-5,
                     q = 0.01,
                     f = 0.09),
               times=c(seq(0, 30, by=0.1)),
               init=c(V=10, U=4e8, L=0, I=0, C=0, A=0), solver="lsoda")

wolfel = read.csv2("/Users/YOUR_NAME/Desktop/XXXX.csv") # loading data from Excel (example path)
patient = wolfel %>% filter(Patient == "A") # subsetting for single patient data

obstime=patient$Day # setting x coordinate for data points
yobs=data.frame(V=c(patient$V)) # setting y coordinate for data points
parms(corona)[c("p", "y")] = c(10, 1e-2) # setting starting values for parameters of interest (needed for 'starting orientation' for the algorithm)

fitModel = fitOdeModel(corona, method="L-BFGS-B", lower=c(p=0, y=0), whichpar=c("p", "y"), obstime, yobs) # telling the algorithm which parameters to fit and making sure the parameters don't go below zero

parms(corona)[c("p", "y")] = fitModel$par # giving new parameters into a 2nd simulation
result = out(sim(corona)) # saving the simulation result in a new variable

plot(result$time, result$V, type="l", main="Patient #NUM", col="red", ylim=c(1, 2e8), yaxt="n", xlab="Days after infection", ylab="Virus copies / mL")
points(obstime, yobs$V, pch=16, col="red2")
eaxis(2) # Plotting the fitted curve with the data points, to which it is fitted

ssq = ssqOdeModel(p=c(p=fitModel$par[1], y=fitModel$par[2]), corona, obstime, yobs) # calculation of the sum of squares
fitModel$par # outputs result of fitting in the console
ssq # outputs result of sum-of-squares calculation in the console