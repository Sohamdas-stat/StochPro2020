install.packages("GillespieSSA")

### Showing that fixation probability = initial proportion ###

########### verify absorption probabilities ############
L = 50 #how many proportions of A considered
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
  }
  proportion.fixation = c(proportion.fixation, (n.abs.A/M))
  
}





plot(p, proportion.fixation, pch = 20, cex = 0.8, col = "blue", 
     xlab = "Proportion of allele A initially", ylab = "Proportion of times absorbed at A")

abline(0, 1, lwd = 1.5, lty = "dotted")

###### simulating and plotting Moran process ##############

library(GillespieSSA)

transitions = cbind(c(-1,1), c(1,-1))

N = 1000
p = 0.8


initial.A = floor(N*p)
initial.a = N - initial.A

r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
        parms = NULL, tf = 105, method = ssa.otl())

plot(r$data[,1], r$data[,2]/N, cex = 0.4, pch = 20, type = 'l', col = "darkmagenta",
     xlab = "time", ylab = "X(t)/2N", ylim = c(0,1), xlim = c(0,100))


for(i in 1:29){
  initial.A = floor(N*p)
  initial.a = N - initial.A
  
  r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
          parms = NULL, tf = 105, method = ssa.otl())
  
  lines(r$data[,1], r$data[,2]/N, cex = 0.4, pch = 20, col = "darkmagenta",
        xlab = "time", ylab = "X(t)/2N", ylim = c(0,1), xlim = c(0,100))
}

########## Drawing histogram of tau ################

### draw histograms of rows of the matrix tau_mat

############# stopping time tau ###########

N = 500 #population size
p = c(0.01, 0.1, 0.5, 0.8)

M = 1000


library(GillespieSSA)
transitions = cbind(c(-1,1), c(1,-1))
tau_mat = matrix(nrow = length(p), ncol = M)

myfunc = function(x){
  r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
          parms = NULL, tf = Inf, method = ssa.otl())
  tau_disc = length(r$data[,1])
  tau = r$data[tau_disc, 1]
  return(tau)
}

for(i in 1:length(p)){
  initial.A = floor(N*p[i])
  initial.a = N - initial.A
  vec = rep(1, M)
  tau_mat[i, ] = sapply(vec, myfunc)
}

########### Plotting expected tau against starting proportion of A #####

############# stopping time tau ###########

N = 500 #population size
p = seq(0.01, 0.99, length = 20)
M = 200


library(GillespieSSA)
transitions = cbind(c(-1,1), c(1,-1))
tau_mat2 = matrix(nrow = length(p), ncol = M)

myfunc = function(x){
  r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
          parms = NULL, tf = Inf, method = ssa.otl())
  tau_disc = length(r$data[,1])
  tau = r$data[tau_disc, 1]
  return(tau)
}

for(i in 1:length(p)){
  initial.A = floor(N*p[i])
  initial.a = N - initial.A
  vec = rep(1, M)
  tau_mat2[i, ] = sapply(vec, myfunc)
}


meanv = NULL
for(i in 1:length(p)){
  meanv = c(meanv, mean(tau_mat2[i, ]))
}

plot(p, meanv, pch = 20, cex = 2, col = "darkmagenta", xlab = "Starting proportion of allele A", 
     ylab = "Expected stopping time")

theo_meanv = NULL
for(i in 1:length(p)){
  x = -2*N * (p[i]*log(p[i]) + (1 - p[i])*log(1-p[i]))
  theo_meanv = c(theo_meanv, x)
}

points(p, theo_meanv, pch = 20, col = "red")
