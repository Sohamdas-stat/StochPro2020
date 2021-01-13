########### verify absorption probabilities ############
L = 20 #how many proportions of A considered
M = 50 #number of simulations to get proportion of times of absorption at A
N = 500 #Population size
p = seq(0.01, 0.99, length = L)

library(GillespieSSA)
transitions = cbind(c(-1,1), c(1,-1))


proportion.fixation = NULL

for(i in 1:L){
  initial.A = floor(N*p[i])
  initial.a = N - initial.A
  n.abs.A = 0
  for(j in 1:M){
    #set.seed(j)
    r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
            parms = NULL, tf = Inf, method = ssa.btl())
    tau_disc = length(r$data[,1])
    if(r$data[tau_disc, 2] == N){n.abs.A = n.abs.A + 1}
    #rm(simu)
  }
  proportion.fixation = c(proportion.fixation, (n.abs.A/M))
  
}





plot(p, proportion.fixation, pch = 20, cex = 0.8, col = "darkgreen", 
     xlab = "Proportion of allele A initially", ylab = "Proportion of times absorbed at A")

abline(0, 1, lwd = 2, lty = "dotted")


