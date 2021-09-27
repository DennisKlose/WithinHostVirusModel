library(simecol) # solves differential equations
library(ggplot2) # for creating visually pleasing plots
require(scales) # for logarithmic axes
library(tidyr) # for transforming wide date format into long data format
library(dplyr) # for renaming column titles

m = c(0, 0.001, 0.01, 0.1, 1) #Plasma transfer no/yes

corona = new("odeModel",
             main=function(time, init, parms){
               with(as.list(c(init, parms)),{
                 if(V > 1){b1 = b} else {b1 = 0}
                 if(time < 5){a1 = 0} else {a1 = a} # Heaviside function in order to stop asymptotic curve
                 dV = p*I-y*V*A-b*U*V # Here, the differential equations must be given
                 dU = r*U*(1-(U/Umax))-b1*U*V
                 dL = b1*U*V-k*L
                 dI = k*L-z*I*C
                 dC = w*I-f*C
                 dA = x*C-q*A+a1*P
                 dP = -a1*P
                 list(c(dV, dU, dL, dI, dC, dA, dP))})},
             parms=c(Umax = 4e8, # Here, all parameter values must be set
                     p = 9,
                     r = 0.5,
                     k = 0.8,
                     b = 4e-8,
                     y = 7e-2,
                     z = 5e-8,
                     w = 0.09,
                     x = 2e-5,
                     q = 0.01,
                     f = 0.09,
                     a = m[1]),
             times=c(seq(0, 20, by=0.1)),
             # Here, the duration and intervals of simulation must be set
             init=c(V=10, U=4e8, L=0, I=0, C=0, A=0, P=100000), solver="lsoda")
# Here, the start values are set
outc_m_list = list()

for (i in 1:length(m)){ # This is needed for creating multiple graphs with different 'a'
  parms(corona)["a"] = m[i]
  outc_m = out(sim(corona))
  outc_m$Plasma = m[i]
  outc_m_list[[i]] = outc_m}

big_outc_m = dplyr::bind_rows(outc_m_list) # merging the data
big_outc_m[,9] = as.factor(big_outc_m[,9]) # for correct legend display in ggplot2
big_outc_m$V = pmax(big_outc_m$V, 1)

param = ggplot(big_outc_m) +
  geom_line(aes(x = time, y = V, group = Plasma, colour = Plasma)) +
  labs(x = "Days after infection", y = "Virus copies/ mL")

my_color = brewer.pal(n = 5, "Greens")[2:6] 

param + theme_bw(base_size = 20) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values=my_color)