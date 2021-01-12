install.packages("adaptivetau")
library(adaptivetau)
transitions = cbind(c(-1,1), c(1,-1))

driftRateF = function(x, params, t){
  rep(x[1]*x[2]/sum(x), 2)
}

#seedr = 123
#set.seed(seedr)
initial.A = 500
initial.a = 1500
r = ssa.adaptivetau(c(initial.A, initial.a), transitions, 
                    driftRateF, params = NULL, tf = Inf)

tau_disc = length(r[,1])
tau = r[tau_disc, 1]
hetero = c()
for(i in 1:tau_disc){
  h = 2*r[i,2]*r[i,3]/((initial.a+initial.A)^2)
  hetero = c(hetero, h)
}

plot(r[,1], r[,2], cex = 0.4, pch = 20, type = 'p')
plot(r[,1], hetero, cex = 0.4, pch = 20, type = 'p')
