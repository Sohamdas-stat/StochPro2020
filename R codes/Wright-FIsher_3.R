library(tictoc)

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

plot.tau <- function(table){
#	quartz(width = 10, height = 7)
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

plot.tau(table)
