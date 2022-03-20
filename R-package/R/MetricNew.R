library("mvtnorm")

set.seed(1)

#################indicator function
indic = function(X, u, v){
  r = sqrt(sum((u-v)^2))
  return(1*(sum((X-u)^2)<=r^2))
}


##################Empirical metric dist. fun.
#input: data: a matrix with row num being the dimension, and col num being the sample size
Fn = function(u, v, data){
  size = length(data[1,])
  a = rep(0, size)
  a = apply(data, 2, indic, u, v)
  return(sum(a)/size)
}

###########################################We first estimate the global quantile by
###########################################generating a large sample of size N 
#################generating r.v.s
N = 2000
K = 2

#####K-dimensional normal dist.
XX = t(rmvnorm(N, mean = rep(0, K), sigma = diag(K)))

#mixture of 3 K-dimensional normal dists
mu1 = c(-3, 0)
mu2 = c(3, 0)
mu3 = c(0, -2.5)
sigma1 = matrix(c(5, -4, -4, 5), 2, 2)
sigma2 = matrix(c(5, 4, 4, 5), 2, 2)
sigma3 = matrix(c(4, 0, 0, 1), 2, 2)
XX1 = rmvnorm(N, mean = mu1, sigma = sigma1)
XX2 = rmvnorm(N, mean = mu2, sigma = sigma2)
XX3 = rmvnorm(N, mean = mu3, sigma = sigma3)
u = runif(N)
XX = (u<3/8)*XX1 + ((u>=3/8)&(u<3/4))*XX2 + (u>=3/4)*XX3
#XX = (u<1/2)*XX1 + (u>=1/2)*XX2
XX = t(XX)

plot(XX[1, ], XX[2, ], asp = 1, cex = 1, xlab = ' ', ylab = ' ')

u = c(0, 0)
v = c(5, 5)

Fn(u, v, XX)

##################compute 0th global quantile function
#########n^{-1}\sum Fn(X, u): function for performing the maximization
average_Fn = function(u, data){
  size = length(data[1,])
  a = rep(0, size)
  a = apply(data, 2, Fn, u, data)
  return(sum(a)/size)
}

average_Fn(c(10,10), XX)


#########n^{-1}\sum_{i=1}^n Fn(X_i, X_j) for j = 1, ..., n
aver_all = apply(XX, 2, average_Fn, XX)

Rank = rank(aver_all)

tau = c(0.2, 0.5, 0.8)
epsilon = 0.02

q1_index = which((Rank >= tau[1]*N)&(Rank < (tau[1]+epsilon)*N))
q1 = XX[,q1_index]

q2_index = which((Rank >= tau[2]*N)&(Rank < (tau[2]+epsilon)*N))
q2 = XX[,q2_index]

q3_index = which((Rank >= tau[3]*N)&(Rank < (tau[3]+epsilon)*N))
q3 = XX[,q3_index]


############################Empirical quantiles
n = 500


#####K-dimensional normal dist.
X = t(rmvnorm(n, mean = rep(0, K), sigma = diag(K)))

#mixture of 2 K-dimensional normal dists
X1 = rmvnorm(n, mean = mu1, sigma = sigma1)
X2 = rmvnorm(n, mean = mu2, sigma = sigma2)
X3 = rmvnorm(n, mean = mu3, sigma = sigma3)
u = runif(n)
X = (u<3/8)*X1 + ((u>=3/8)&(u<3/4))*X2 + (u>=3/4)*X3
X = t(X)

aver_sample = apply(X, 2, average_Fn, X)

Rank_WW = rank(aver_sample)

qn1_index = which((Rank_WW >= tau[1]*n)&(Rank_WW < tau[1]*n+1))
qn1 = X[, qn1_index]

qn2_index = which((Rank_WW >= tau[2]*n)&(Rank_WW < tau[2]*n+1))
qn2 = X[, qn2_index]

qn3_index = which((Rank_WW >= tau[3]*n)&(Rank_WW < tau[3]*n+1))
qn3 = X[, qn3_index]

par(mfrow = c(1,1), mar = c(2, 2, 1,1))
plot(X[1, ], X[2, ], asp = 1, cex = 1, xlab = ' ', ylab = ' ')

#tau[1] sample quantile
points(qn1[1], qn1[2], col = 'blue', pch = 2, lwd=1.5, cex = 1.5, asp = 1)
#estimation of tau[1] global quantile
polygon(q1[1, order(atan2(q1[1,]-mean(q1[1,]),q1[2,]-mean(q1[2,])))], q1[2,order(atan2(q1[1,]-mean(q1[1,]),q1[2,]-mean(q1[2,])))], border = "blue", lwd = 1.5, lty = 1)
#points(q1[1,], q1[2,], col = 'green', pch = 2, lwd=2, cex = 1, asp = 1)

#tau[2] sample quantile
points(qn2[1], qn2[2], col = 'purple', pch = 8, lwd=1.5, cex = 1.5, asp = 1)
#estimation of tau[2] global quantile
polygon(q2[1, order(atan2(q2[1,]-mean(q2[1,]),q2[2,]-mean(q2[2,])))], q2[2,order(atan2(q2[1,]-mean(q2[1,]),q2[2,]-mean(q2[2,])))], border = "purple", lwd = 1.5, lty = 1)
#points(q2[1,], q2[2,], col = 'yellow', pch = 2, lwd=2, cex = 1, asp = 1)


#tau[3] sample quantile
points(qn3[1], qn3[2], col = 'red', pch = 0, lwd=1.5, cex = 1.5, asp = 1)
#estimation of tau[3] global quantile
polygon(q3[1, order(atan2(q3[1,]-mean(q3[1,]),q3[2,]-mean(q3[2,])))], q3[2,order(atan2(q3[1,]-mean(q3[1,]),q3[2,]-mean(q3[2,])))], border = "red", lwd = 1.5, lty = 1)
#points(q3[1,], q3[2,], col = 'black', pch = 2, lwd=2, cex = 1, asp = 1)



legend('bottomleft', plot= T, 
       legend = c("0.2 empirical quantile", "0.2 quantile", 
                  "0.5 empirical quantile", "0.5 quantile", 
                  "0.8 empirical quantile", "0.8 quantile"), 
       pch = c(2, NA, 8, NA, 0, NA), 
       col = c('blue', 'blue', 'purple', 'purple', 'red', 'red'),
       lty = c(NA, 1, NA, 1, NA, 1),
       cex = 1)

##############################plot of estimation of global quantile contour
par(mfrow = c(1,1), mar = c(2, 2, 1,1), cex = 1)
plot(XX[1, ], XX[2, ], asp = 1, cex = 0.5, xlab = ' ', ylab = ' ')

#estimation of tau[1] global quantile
polygon(q1[1, order(atan2(q1[1,]-mean(q1[1,]),q1[2,]-mean(q1[2,])))], q1[2,order(atan2(q1[1,]-mean(q1[1,]),q1[2,]-mean(q1[2,])))], border = "blue", lwd = 1.5, lty = 1)
#estimation of tau[2] global quantile
polygon(q2[1, order(atan2(q2[1,]-mean(q2[1,]),q2[2,]-mean(q2[2,])))], q2[2,order(atan2(q2[1,]-mean(q2[1,]),q2[2,]-mean(q2[2,])))], border = "purple", lwd = 1.5, lty = 1)
#estimation of tau[3] global quantile
polygon(q3[1, order(atan2(q3[1,]-mean(q3[1,]),q3[2,]-mean(q3[2,])))], q3[2,order(atan2(q3[1,]-mean(q3[1,]),q3[2,]-mean(q3[2,])))], border = "red", lwd = 1.5, lty = 1)

points(q1[1,], q1[2,], col = 'blue', pch = 0, lwd=1.5, cex = 1, asp = 1)
points(q2[1,], q2[2,], col = 'purple', pch = 0, lwd=1.5, cex = 1, asp = 1)
points(q3[1,], q3[2,], col = 'red', pch = 0, lwd=1.5, cex = 1, asp = 1)

legend('bottomleft', plot= T, 
       legend = c("0.2 metric quantile", 
                  "0.5 metric quantile", 
                  "0.8 metric quantile"), 
       col = c('blue', 'purple', 'red'),
       lty = c(1, 1, 1),
       cex = 1)

