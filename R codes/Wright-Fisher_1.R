library(tictoc)

pop.size = c(50,100,1000,5000)
start.p = c(0.01,0.1,0.5,0.8)
n.gen = 100
n.rep = 50

wright.fisher <- function(){
  wf.tab = data.frame()
  for(N in pop.size){
    print(N)
    for(p_0 in start.p){
      for(k in 1:n.rep){
        p = p_0
        for(j in 1:n.gen){
          X = rbinom(1,2*N,p)
          p = X/(2*N)
          h = 2*p*(1-p)*(2*N)/(2*N-1)
          wf.tab = rbind(wf.tab,c(k,N,j,p_0,p,h))
        }
      }
    }
  }
  names(wf.tab) = c("replicate","N","gen","start.p","p","hetzy")
  print(dim(wf.tab))
  return(wf.tab)
}

tic()
table = wright.fisher()
toc()

plot.p <- function(table){
 quartz(width = 10, height = 7)
 layout(matrix(1:16, nrow = 4, byrow = T))	
  for(pop in pop.size){
    for(p in start.p){
      count_0 = 0
      count_1 = 0
      plot(NULL,xlim=c(0,100),ylim=c(0,1),xlab="gen",ylab="p")
      for(k in c(1:n.rep)){
        n.row = which(table$replicate == k & table$N == pop & table$start.p == p)
        tab = table[n.row,]
        print(dim(tab))
        lines(c(1:100),tab$p,lwd=1,lty=1,col=rgb(0,0.5,1))
        count_0 = count_0 + isTRUE(tab$p[50] == 0)
        count_1 = count_1 + isTRUE(tab$p[50] == 1)
      }
      tau_0 = count_0/50
      tau_1 = count_1/50
      #			legend("topleft",legend=c(paste("Pr_0=",tau_0,", Pr_1=",tau_1)))
    }
  }
}

plot.p(table)