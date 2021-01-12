library(GillespieSSA)

transitions = cbind(c(-1,1), c(1,-1))

driftRateF = function(x, params, t){
  rep(x[1]*x[2]/sum(x), 2)
}

N = 2000
initial.A = 500
initial.a = 1500
r = ssa(c(Y1 = 500, Y2 = 1500),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
                    parms = NULL, tf = Inf, method = ssa.etl())
plot(r$data[,1], r$data[,2], cex = 0.4, pch = 20, type = 'p')
