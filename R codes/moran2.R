library(GillespieSSA)

transitions = cbind(c(-1,1), c(1,-1))

draw_A = function(N, p ){
  initial.A = floor(N*p)
  initial.a = N - initial.A
  
  r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
          parms = NULL, tf = Inf, method = ssa.otl())
  
  plot(r$data[,1], r$data[,2], cex = 0.4, pch = 20, type = 'b', xlab = "time", ylab = "A(t)")
}

draw_A(2000, 0.5)

draw_hetero = function(N, p ){
  initial.A = floor(N*p)
  initial.a = N - initial.A
  
  r = ssa(c(Y1 = initial.A, Y2 = initial.a),  a = c("{Y1}*{Y2}/N", "{Y1}*{Y2}/N") ,transitions, 
          parms = NULL, tf = Inf, method = ssa.etl())
  
  h = NULL
  for(i in 1:length(r$data[,1])){
    het = 2*r$data[i,2]*r$data[i,3]/(N^2)
    h = c(h, het)
  }
  plot(r$data[,1], h, cex = 0.4, pch = 20, type = 'p', xlab = "time", ylab = "Heterozygosity(t)")
}

draw_hetero(2000, 0.5)
