library(simecol) # solves differential equations
library(ggplot2) # for creating visually pleasing plots
require(scales) # for logarithmic axes
library(tidyr) # for transforming wide date format into long data format
library(dplyr)

corona = new("odeModel",
             main=function(time, init, parms){
               with(as.list(c(init, parms)),{
                 if(V < 1){b1 = 0} else {b1 = b}
                 if(V < 1){p1 = 0} else {p1 = p} # Heaviside functions to stop unwanted reactivations of "dormant" virus
                 dV = p1*I-y*V*A-b*U*V # Here, the differential equations must be given
                 dU = r*U*(1-(U/Umax))-b1*U*V
                 dL = b1*U*V-k*L
                 dI = k*L-z*I*C
                 dC = w*I-f*C
                 dA = x*C-q*A
                 list(c(dV, dU, dL, dI, dC, dA))})},
             parms=c(Umax = 4e8, # Here, all parameter values must be entered
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
             times=c(seq(0, 60, by=0.1)),
             # Here, the duration and intervals of simulation must be set
             init=c(V=100, U=4e8, L=0, I=0, C=0, A=0), solver="lsoda")
outc = out(sim(corona)) # The output table is stored under a variable

outc = outc %>% rename( # Here, all column titles are renamed for the plot legend
  Virus = V, Target = U, Latent = L, Infected = I, Cytotoxic = C, Antibody = A)

long_outc = gather(outc, Compartment, Value, Virus:Antibody)
# Here, wide data is turned into long data, which is necessary when using ggplot2
long_outc$Value = pmax(long_outc$Value, 1) # Here, all values under 1 are left out in order to get a tidy logarithmic y-axis formating

model = ggplot(long_outc) + # Here, the plot is stored under a variable
  geom_line(aes(x = time, y = Value, group = Compartment, colour = Compartment)) +
  labs(x = "Days after infection", y = "Compartment value")

model + # Here, the plot is drawn with specific modifications
  theme_bw(base_size = 25) +
  geom_hline(yintercept=1e2, linetype="dashed", color="grey") +
  geom_vline(xintercept=5, linetype="dashed", color="grey") +
  geom_vline(xintercept=55, linetype="dashed", color="grey") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))