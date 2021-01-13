############# stopping time tau ###########
start_time = Sys.time()
N = 1000 #population size
p = seq(0.01, 0.99, length = 25)
M = 100


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
end_time = Sys.time()

end_time - start_time

meanv = NULL
for(i in 1:length(p)){
  meanv = c(meanv, mean(tau_mat[i, ]))
}

plot(p, meanv, pch = 20, cex = 2, col = "darkmagenta", xlab = "Starting proportion of allele A", 
     ylab = "Expected stopping time")

theo_meanv = NULL
for(i in 1:length(p)){
  x = -2*N * (p[i]*log(p[i]) + (1 - p[i])*log(1-p[i]))
  theo_meanv = c(theo_meanv, x)
}

points(p, theo_meanv, pch = 20, col = "red")

