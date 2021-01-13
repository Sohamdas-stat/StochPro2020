library(tictoc)

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

xyz <- function(table){
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

xyz(table)
