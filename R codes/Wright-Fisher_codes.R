library(tictoc)

##	Plotting p_t vs t for different starting frequency and different population sizes

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
 par(mfrow=c(4,4))
  for(pop in pop.size){
    for(p in start.p){
      count_0 = 0
      count_1 = 0
      plot(NULL,xlim=c(0,100),ylim=c(0,1),xlab="gen(t)",ylab="p_t")
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


##	Plotting probability of fixation vs starting probability

pop.size = c(50,100,500,1000)
start.p = seq(0.01,0.99,by=0.02)
n.rep = 200

wright.fisher <- function(){
	tau.tab = data.frame()
	for(N in pop.size){
		print(N)
		for(p_0 in start.p){
#			row = NULL
			for(k in 1:n.rep){
				p = p_0
				count = 0
				while(p*(1-p) != 0){
					X = rbinom(1,2*N,p)
					p = X/(2*N)
					count = count + 1
#					row = c(row,count)
				}
				tau.tab = rbind(tau.tab,c(N,p_0,X,count))
			}
		}
	}
	names(tau.tab) = c("pop.size","start.p","X_fix","Tau")
	return(tau.tab)
}

tic()
table = wright.fisher()
toc()

plot.fixprob <- function(table){
	quartz(width = 10, height = 7)
	layout(matrix(1:4, nrow = 2, byrow = T))
	for(pop in pop.size){
	plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="starting probability",ylab="prob. of fixation to 2N")
	lines(x=c(0,1),y=c(0,1),col=1,lty=1,lwd=1)	
		for(p in start.p){
			n.row = which(table$pop.size == pop & table$start.p == p)
			tab = table[n.row,]
			p.fix = length(which(tab$X_fix == 0))/n.rep
			points(p,1-p.fix,pch=20,col=rgb(0.4,0.2,0.8))
		}
	}
}

plot.fixprob(table)

##	Plotting fixation time vs starting probability

plot.tau <- function(table){
	quartz(width = 10, height = 7)
	layout(matrix(1:1, nrow = 1, byrow = T))
	plot(NULL,xlim=c(0,1),ylim=c(0,3000),xlab="starting probability",ylab="expected fixation time")
	for(i in 1:4){
		pop = pop.size[i]
		for(p in start.p){
			n.row = which(table$pop.size == pop & table$start.p == p)
			tab = table[n.row,]
			points(p,mean(tab$Tau),pch=20,col=c("red","blue","green","black")[i])
		}
	}
}

plot.tau(table)


##	Plotiing histogram of fixation time

pop.size = c(50,200,500)
start.p = c(0.01,0.1,0.5,0.8)
n.rep = 1000

wright.fisher <- function(){
	tau.tab = data.frame()
	for(N in pop.size){
		print(N)
		for(p_0 in start.p){
			for(k in 1:n.rep){
				p = p_0
				count = 0
				while(p*(1-p) != 0){
					X = rbinom(1,2*N,p)
					p = X/(2*N)
					count = count + 1
				}
				tau.tab = rbind(tau.tab,c(count,X))
			}
			
		}
	}
	return(tau.tab)
}

tic()
table = wright.fisher()
toc()
names(table) = c("Tau","Fix")

hist.tau <- function(table){
	quartz(width = 10, height = 7)
	layout(matrix(1:12, nrow = 3, byrow = T))	
	for(i in 1:3){
		for(j in 1:4){
			rows = seq(4000*(i-1)+1000*(j-1)+1,4000*(i-1)+1000*j,by=1)
			tab = table[rows,]
			p.fix = length(which(tab$Fix == 0))/1000
			hist(as.numeric(tab$Tau), breaks=seq(0,max(table)+200,200), prob=T, ylim=c(0,0.002), main="\n", xlab="Tau", col=rgb(0.8,0.8,0))
			legend(x=c(2500,8000),y=c(0.0008,0.002),legend=c(paste("Mean = ",mean(tab$Tau),"\nFix_prob: 1->",round(1-p.fix,3),"\n                0->",round(p.fix,3))))
		}
	}
}

hist.tau(table)
