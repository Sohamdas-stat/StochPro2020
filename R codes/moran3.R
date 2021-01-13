############# stopping time tau ###########

N = 500 #population size
p = seq(0.1, 0.9, by = 0.5)

M = 100


library(GillespieSSA)
transitions = cbind(c(-1,1), c(1,-1))

tau_mat = matrix(nrow = length(p), ncol = M)

for(i in 1:length(p)){
  initial.A = floor(N*p[i])
  initial.a = N - initial.A
  for(j in 1:M){
    r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
            parms = NULL, tf = Inf, method = ssa.otl())
    tau_disc = length(r$data[,1])
    tau = r$data[tau_disc, 1]
    tau_mat[i,j] = tau
  }
}

par(mfrow = c(1,2))
for(i in 1:length(p)){
  hist(tau_mat[i, ], breaks = 20)
}

mean(tau_mat[2,])
