#simulate truncated normal
trunc_norm <- function(mu=0, sd=1, l, u, size){
  x <- runif(size)
  left <- pnorm(l,mean=mu,sd=sd)
  right <- pnorm(u, mean=mu, sd=sd)
  quan <- left + x*(right-left)
  return(qnorm(quan, mean=mu, sd=sd))
}

#Gibbs sampler
Gibbs <- function(X, Y, D, Pi, init_beta, lam=1, a=0.1, b=0.01, niter=5000, BurnIn=3000){
  m = ncol(X)
  n = nrow(X)
  Xobs = X[D]
  Yobs = Y[D]
  beta = init_beta
  mat_var = lam*diag(m) + crossprod(X, diag(1/Pi)%*%X)
  r <- eigen(mat_var, symmetric=TRUE)
  V <- r$vectors %*% diag(1/sqrt(r$values))
  
  InvPiX = diag(1/Pi)%*%X
  yPiy = crossprod(Y, diag(1/Pi)%*%Y)
  y_pred = numeric(length=n)
  mean_post = numeric(length=n)
  sig_post = numeric(length=1)
  overall_mean = numeric(length=niter-BurnIn)
  for(i in 1:niter){
    post_a = a + (sum(D/Pi)+m)/2
    post_b = (crossprod(beta, mat_var%*%beta) - 2*crossprod(Y,InvPiX%*%beta)+yPiy+2*b)/2
    sig = 1/rgamma(n=1, shape=post_a, rate=post_b)
    post_mu = solve(mat_var)%*%crossprod(InvPiX, Y)
    post_sd = sqrt(sig)*V
    beta = post_sd%*%rnorm(m) + post_mu
    if(i>BurnIn){
      temp = rnorm(n, mean=X%*%beta, sd=sqrt(sig)) 
      y_pred = y_pred + temp 
      mean_post = mean_post + X%*%beta
      sig_post = sig_post + sig
      overall_mean[i-BurnIn] = mean(temp)
    }
  }
  y_pred = y_pred/(niter-BurnIn)
  mse1 = mean((y_pred-Y)^2)
  mean_post = mean_post/(niter-BurnIn)
  mse2 = mean((mean_post-Y)^2)
  sig_post = sig_post/(niter-BurnIn)
  l = quantile(overall_mean, 0.025)
  u = quantile(overall_mean, 0.975)
  return(c(mean(mean_post), sig_post, mse1, mse2, l, u))
}

#true model case
res_mean = numeric(1000)
res_sig = numeric(1000)
mse1 = numeric(1000)
mse2 = numeric(1000)
l = numeric(1000)
u = numeric(1000)
set.seed(1)
for(i in 1:1000){
  beta0 <- c(2.5,1.5)
  beta1 <- c(-2,1.6)
  X1 <- trunc_norm(l=-2,u=1,size=100)
  X2 <- trunc_norm(l=1,u=3, size=100)
  X <- cbind(X1,X2)
  y <- X%*%beta1 + rnorm(sd=1.5,100)
  p <- exp(X%*%beta0)/(1+exp(X%*%beta0))
  index <- rbinom(n=100,size=1,prob=p)
  fit <- glm(index~X, family=binomial(link="logit"))
  Pi <- 1/(1+exp(-predict(fit)))
  temp <- Gibbs(X, y, as.logical(index), Pi, c(0,0), lam=0.01)
  res_mean[i] = temp[1]
  res_sig[i] = temp[2]
  mse1[i] = temp[3]
  mse2[i] = temp[4]
  l[i] = temp[5]
  u[i] = temp[6]
}

tar = (dnorm(1)-dnorm(3))/(pnorm(3)-pnorm(1))*1.6+(dnorm(-2)-dnorm(1))/(pnorm(1)-pnorm(-2))*(-2)
sum(tar>l&tar<u)
mean(mse1)
mean(mse2)

#wrong model case
res_mean2 = numeric(1000)
res_sig2 = numeric(1000)
mse12 = numeric(1000)
mse22 = numeric(1000)
l2 = numeric(1000)
u2 = numeric(1000)
set.seed(1)
for(i in 1:1000){
  beta0 <- c(2.5,1.5)
  X1 <- trunc_norm(l=-2,u=1,size=100)
  X2 <- trunc_norm(l=1,u=3, size=100)
  X <- cbind(X1,X2)
  y <- 1.5*(2*X1-1)^3 - (X2+1.7)^3 + rnorm(100)
  p <- exp(X%*%beta0)/(1+exp(X%*%beta0))
  index <- rbinom(n=100,size=1,prob=p)
  fit <- glm(index~X, family=binomial(link="logit"))
  Pi <- 1/(1+exp(-predict(fit)))
  temp <- Gibbs(X, y, as.logical(index), Pi, c(0,0), lam=0.01)
  res_mean2[i] = temp[1]
  res_sig2[i] = temp[2]
  mse12[i] = temp[3]
  mse22[i] = temp[4]
  l2[i] = temp[5]
  u2[i] = temp[6]
}


x1 <- trunc_norm(l=-2,u=1,size=10000)
x2 <- trunc_norm(l=1,u=3, size=10000)

tar = mean(1.5*(2*x1-1)^3 - (x2+1.7)^3)
sum(tar>l2&tar<u2)
mean(mse12)
mean(mse22)
