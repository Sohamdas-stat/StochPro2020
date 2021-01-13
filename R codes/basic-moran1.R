library(GillespieSSA)

transitions = cbind(c(-1,1), c(1,-1))



N = 2000
initial.A = 500
initial.a = 1500
set.seed(1)
r1 = ssa(c(Y1 = 500, Y2 = 1500),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
        parms = NULL, tf = Inf, method = ssa.d())
tau1 = r1$data[ length(r1$data[,1]), 1]

set.seed(1)
r2 = ssa(c(Y1 = 500, Y2 = 1500),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
        parms = NULL, tf = Inf, method = ssa.etl())
tau2 = r2$data[ length(r2$data[,1]), 1]

plot(r$data[,1], r$data[,2], cex = 0.4, pch = 20, type = 'p')
