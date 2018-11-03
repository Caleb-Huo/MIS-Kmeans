S <- 2
K <- 3
G <- 1000
g1 <- 50
g2 <- 50
n0 <- 20
n <- K*n0
labels <- cut(1:n,breaks=K,labels=FALSE)

set.seed(32611)
S1 <- matrix(rnorm(G*n), nrow=G, ncol=n)
S2 <- matrix(rnorm(G*n), nrow=G, ncol=n)

S1[1:g1, labels==1] <- S1[1:g1, labels==1] + 2
S1[1:g1, labels==3] <- S1[1:g1, labels==3] - 2
S1[g1 + 1:g2, labels==1] <- S1[g1 + 1:g2, labels==1] - 2
S1[g1 + 1:g2, labels==2] <- S1[g1 + 1:g2, labels==2] + 2


S2[1:g1, labels==2] <- S2[1:g1, labels==2] + 2
S2[1:g1, labels==1] <- S2[1:g1, labels==1] - 2
S2[g1 + 1:g2, labels==2] <- S2[g1 + 1:g2, labels==2] - 2
S2[g1 + 1:g2, labels==3] <- S2[g1 + 1:g2, labels==3] + 2


S = list(t(S1),t(S2))
groups <- Map('c',1:g1,g1 + 1:g2)

res <- MISKmeans(d = S, K = 3, gamma = 0.4, group = groups)
inferLambda(res)
inferAlpha(res)

if(F){
	d = S; K = 3; gamma = 0.4; lambda = 0.5; alpha = 0.5; group = groups; nstart = 20; 
	    wsPre = NULL; iniWbound = 20; penaltyInfo = NULL; silent = FALSE; maxiter = 20; sampleSizeAdjust = FALSE
		

d, K, groupInfo=groupInfoIni, Cs, ws=wsPre, tss.x, lambda, sampleSizeAdjust = sampleSizeAdjust, silent=silent

}

names(res)

order_S1 <- order(res$Cs[[1]])
order_S2 <- order(res$Cs[[2]])

S1_sub <- S1[res$ws!=0, order_S1]
S2_sub <- S2[res$ws!=0, order_S2]

col_S1 <- palette()[res$Cs[[1]][order_S1]]
col_S2 <- palette()[res$Cs[[2]][order_S2]]

heatmap(S1_sub, Rowv=NA, Colv=NA, ColSideColors=col_S1)
heatmap(S2_sub, Rowv=NA, Colv=NA, ColSideColors=col_S2)


