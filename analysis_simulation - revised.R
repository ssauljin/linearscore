par(mfrow=c(1,2), mar=c(4.1, 4.1, 3.1, 1.1), 
    oma=c(0,0,0,0), mgp=c(2,0.8,0))

library(pROC)
library(nleqslv)
library(cvTools)
library(MASS)

R    <- 10
set.seed(1000)

N <- 10000
n <- 500

U  <- rbinom(n, 1, 0.1)
X1 <- U * mvrnorm(n = n, mu=c(4.95 ,2.2), Sigma=matrix((c(0.5,0,0,0.3)),nrow=2)) +
   (1-U)* mvrnorm(n = n, mu=c(2.5  ,2.6), Sigma=matrix((c(0.2,0,0,0.1)),nrow=2))

X0 <-     mvrnorm(n = N, mu=c(4.95 ,1.5), Sigma=matrix((c(1  ,0,0,0.3)),nrow=2))

#plot( X0, pch='.', col="blue")
#lines(X1, pch='.', col="red", type="p")

dt     <- cbind(          X0, rep(0, nrow(X0)))
dt     <- rbind(dt, cbind(X1, rep(1, nrow(X1))))
dt     <- as.data.frame(dt)
colnames(dt) <- c("Xo","Xt","Y")
head(dt)
folds <- cvFolds(NROW(dt), K=R)

rels_se_logit_test <- rep(0, 100)
rels_se_exp9_test  <- rep(0, 100)
rels_se_exp5_test  <- rep(0, 100)
rels_se_exp1_test  <- rep(0, 100)
rels_se_log9_test  <- rep(0, 100)
rels_se_log5_test  <- rep(0, 100)
rels_se_log1_test  <- rep(0, 100)

rels_sp_logit_test <- rep(0, 100)
rels_sp_exp9_test  <- rep(0, 100)
rels_sp_exp5_test  <- rep(0, 100)
rels_sp_exp1_test  <- rep(0, 100)
rels_sp_log9_test  <- rep(0, 100)
rels_sp_log5_test  <- rep(0, 100)
rels_sp_log1_test  <- rep(0, 100)

auc_logit_test <- 0
auc_exp9_test  <- 0
auc_exp5_test  <- 0
auc_exp1_test  <- 0
auc_log9_test  <- 0
auc_log5_test  <- 0
auc_log1_test  <- 0

log9_minus <- function(z) {
  gamma = 0.9
  result <- -0.5/gamma*log(1/(1+gamma*exp(z)))-
    0.5*log((1-gamma)*exp(-z)/(1+exp(-z)*(1-gamma)))
  return(result) }

log5_minus <- function(z) {
  gamma = 0.5
  result <- -0.5/gamma*log(1/(1+gamma*exp(z)))-
    0.5*log((1-gamma)*exp(-z)/(1+exp(-z)*(1-gamma)))
  return(result) }

log1_minus <- function(z) {
  gamma = 0.1
  result <- -0.5/gamma*log(1/(1+gamma*exp(z)))-
    0.5*log((1-gamma)*exp(-z)/(1+exp(-z)*(1-gamma)))
  return(result) }


log9_plus <- function(z) {
  gamma = 0.9
  result <- -0.5*log(gamma*exp(z)/(1+exp(z)*gamma))-
    0.5/(1-gamma)*log(1/(1+(1-gamma)*exp(-z)))
  return(result) }

log5_plus <- function(z) {
  gamma = 0.5
  result <- -0.5*log(gamma*exp(z)/(1+exp(z)*gamma))-
    0.5/(1-gamma)*log(1/(1+(1-gamma)*exp(-z)))
  return(result) }

log1_plus <- function(z) {
  gamma = 0.1
  result <- -0.5*log(gamma*exp(z)/(1+exp(z)*gamma))-
    0.5/(1-gamma)*log(1/(1+(1-gamma)*exp(-z)))
  return(result) }


for (i in 1:R) {

  data <- dt[folds$subsets[folds$which != i],]
  test <- dt[folds$subsets[folds$which == i],]
  
  fit_logit_loss <- glm(Y ~ ., family=binomial(link = "logit"), data=data)
  tX       <- model.matrix(fit_logit_loss)
  
  fun_exp9_loss <- function(beta) {
    lambda = 0.9
    result <- sum(data$Y *(1-lambda)*exp(  -lambda *(tX %*% beta))+
              (1-data$Y)*   lambda *exp((1-lambda)*(tX %*% beta)))
    return(result) }

  fun_exp9_eqn <- function(beta) {
    lambda = 0.9
    result <- t(tX) %*% (-data$Y *exp(  -lambda *(tX %*% beta))+
                       (1-data$Y)*exp((1-lambda)*(tX %*% beta)))*
      lambda*(1-lambda)
    return(result) }
  
  # par9 <- matrix(NA, ncol=100, nrow=ncol(tX))
  # val9 <- rep(NA, 100)
  # for (j in 1:100) {
  #   set.seed(j)
  #   fit_exp9_loss <- optim(coef(fit_logit_loss)+runif(ncol(tX))/500-0.001, fun_exp9_loss)
  #   par9[,j] <- fit_exp9_loss$par
  #   val9[j]  <- fit_exp9_loss$val }
  #zz <- fit_exp9_loss$par
  #summary(fit_exp9_loss$par - zz)
  fit_exp9_loss <- optim(coef(fit_logit_loss), fun_exp9_loss)
  
  fit_exp9_eqn <- nleqslv(fit_exp9_loss$par, fun_exp9_eqn) 
  
  fun_log9_loss <- function(beta) {
    lambda = 0.9
    z <- tX %*% beta
    result <- sum(data$Y *log9_plus( z)+
               (1-data$Y)*log9_minus(z))
    return(result) }

  fit_log9_loss <- optim(coef(fit_logit_loss), fun_log9_loss)
  
  fun_exp5_loss <- function(beta) {
    lambda = 0.5
    result <- sum(data$Y *(1-lambda)*exp(  -lambda *(tX %*% beta) )+
                    (1-data$Y)*   lambda *exp((1-lambda)*(tX %*% beta)))
    return(result) }
  
  fun_exp5_eqn <- function(beta) {
    lambda = 0.5
    result <- t(tX) %*% (-data$Y *exp(  -lambda *(tX %*% beta))+
                           (1-data$Y)*exp((1-lambda)*(tX %*% beta)))*
      lambda*(1-lambda)
    return(result) }
  
  fit_exp5_loss <- optim(coef(fit_logit_loss), fun_exp5_loss)
  
  fit_exp5_eqn <- nleqslv(fit_exp5_loss$par, fun_exp5_eqn) 
  
  # par5 <- matrix(NA, ncol=100, nrow=ncol(tX))
  # val5 <- rep(NA, 100)
  # for (j in 1:100) {
  #   set.seed(j)
  #   fit_exp5_loss <- optim(coef(fit_logit_loss)+runif(ncol(tX))/500-0.001, fun_exp5_loss)
  #   par5[,j] <- fit_exp5_loss$par
  #   val5[j]  <- fit_exp5_loss$val }

  fun_log5_loss <- function(beta) {
    lambda = 0.5
    z <- tX %*% beta
    result <- sum(data$Y *log5_plus( z)+
               (1-data$Y)*log5_minus(z))
    return(result) }
  
  fit_log5_loss <- optim(coef(fit_logit_loss), fun_log5_loss)
  
    
  fun_exp1_loss <- function(beta) {
    lambda = 0.1
    result <- sum(data$Y *(1-lambda)*exp(  -lambda *(tX %*% beta) )+
                    (1-data$Y)*   lambda *exp((1-lambda)*(tX %*% beta)))
    return(result) }
  
  fun_exp1_eqn <- function(beta) {
    lambda = 0.1
    result <- t(tX) %*% (-data$Y *exp(  -lambda *(tX %*% beta))+
                           (1-data$Y)*exp((1-lambda)*(tX %*% beta)))*
      lambda*(1-lambda)
    return(result) }
  
  fit_exp1_loss <- optim(coef(fit_logit_loss), fun_exp1_loss)
  
  fit_exp1_eqn <- nleqslv(fit_exp1_loss$par, fun_exp1_eqn) 
  
  # par1 <- matrix(NA, ncol=100, nrow=ncol(tX))
  # val1 <- rep(NA, 100)
  # for (j in 1:100) {
  #   set.seed(j)
  #   fit_exp1_loss <- optim(coef(fit_logit_loss)+runif(ncol(tX))/500-0.001, fun_exp1_loss)
  #   par1[,j] <- fit_exp1_loss$par
  #   val1[j]  <- fit_exp1_loss$val }
  
  fun_log1_loss <- function(beta) {
    lambda = 0.1
    z <- tX %*% beta
    result <- sum(data$Y *log1_plus( z)+
               (1-data$Y)*log1_minus(z))
    return(result) }
  
  fit_log1_loss <- optim(coef(fit_logit_loss), fun_log1_loss)
  
  fit_logit_test <- glm(Y ~ ., family=binomial(link = "logit"), data=test)
  tS       <- model.matrix(fit_logit_test)
  
  roc_logit_test <- roc(test$Y ~ predict(fit_logit_loss, test))
  roc_exp9_test  <- roc(test$Y ~ as.vector(tS %*% fit_exp9_eqn$x))
  roc_exp5_test  <- roc(test$Y ~ as.vector(tS %*% fit_exp5_eqn$x))
  roc_exp1_test  <- roc(test$Y ~ as.vector(tS %*% fit_exp1_eqn$x))
  roc_log9_test  <- roc(test$Y ~ as.vector(tS %*% fit_log9_loss$par))
  roc_log5_test  <- roc(test$Y ~ as.vector(tS %*% fit_log5_loss$par))
  roc_log1_test  <- roc(test$Y ~ as.vector(tS %*% fit_log1_loss$par))

    
  auc_logit_test <- auc_logit_test + auc(roc_logit_test)/R
  auc_exp9_test  <- auc_exp9_test  + auc(roc_exp9_test )/R
  auc_exp5_test  <- auc_exp5_test  + auc(roc_exp5_test )/R
  auc_exp1_test  <- auc_exp1_test  + auc(roc_exp1_test )/R
  auc_log9_test  <- auc_log9_test  + auc(roc_log9_test )/R
  auc_log5_test  <- auc_log5_test  + auc(roc_log5_test )/R
  auc_log1_test  <- auc_log1_test  + auc(roc_log1_test )/R
  
    
  rels_sp_logit_test <- rels_sp_logit_test +
    as.vector(coords(roc=roc_logit_test, x=(900:999)/1000,
                      input="sensitivity", ret="specificity")[,1])/R
  rels_se_logit_test <- rels_se_logit_test +
    as.vector(coords(roc=roc_logit_test, x=(900:999)/1000,
                      ret="sensitivity", input="specificity")[,1])/R
  rels_sp_exp9_test <- rels_sp_exp9_test +
    as.vector(coords(roc=roc_exp9_test, x=(900:999)/1000,
                      input="sensitivity", ret="specificity")[,1])/R
  rels_se_exp9_test <- rels_se_exp9_test +
    as.vector(coords(roc=roc_exp9_test, x=(900:999)/1000,
                      ret="sensitivity", input="specificity")[,1])/R
  
  rels_sp_log9_test <- rels_sp_log9_test +
    as.vector(coords(roc=roc_log9_test, x=(900:999)/1000,
                     input="sensitivity", ret="specificity")[,1])/R
  rels_se_log9_test <- rels_se_log9_test +
    as.vector(coords(roc=roc_log9_test, x=(900:999)/1000,
                     ret="sensitivity", input="specificity")[,1])/R
  
  rels_sp_exp5_test <- rels_sp_exp5_test +
    as.vector(coords(roc=roc_exp5_test, x=(900:999)/1000,
                      input="sensitivity", ret="specificity")[,1])/R
  rels_se_exp5_test <- rels_se_exp5_test +
    as.vector(coords(roc=roc_exp5_test, x=(900:999)/1000,
                      ret="sensitivity", input="specificity")[,1])/R
  
  rels_sp_log5_test <- rels_sp_log5_test +
    as.vector(coords(roc=roc_log5_test, x=(900:999)/1000,
                     input="sensitivity", ret="specificity")[,1])/R
  rels_se_log5_test <- rels_se_log5_test +
    as.vector(coords(roc=roc_log5_test, x=(900:999)/1000,
                     ret="sensitivity", input="specificity")[,1])/R
  
  rels_sp_exp1_test <- rels_sp_exp1_test +
    as.vector(coords(roc=roc_exp1_test, x=(900:999)/1000,
                      input="sensitivity", ret="specificity")[,1])/R
  rels_se_exp1_test <- rels_se_exp1_test +
    as.vector(coords(roc=roc_exp1_test, x=(900:999)/1000,
                      ret="sensitivity", input="specificity")[,1])/R

  rels_sp_log1_test <- rels_sp_log1_test +
    as.vector(coords(roc=roc_log1_test, x=(900:999)/1000,
                     input="sensitivity", ret="specificity")[,1])/R
  rels_se_log1_test <- rels_se_log1_test +
    as.vector(coords(roc=roc_log1_test, x=(900:999)/1000,
                     ret="sensitivity", input="specificity")[,1])/R
}


rels_se_exp9s  <- rels_se_exp9_test*100
rels_se_exp5s  <- rels_se_exp5_test*100
rels_se_exp1s  <- rels_se_exp1_test*100
rels_se_log9s  <- rels_se_log9_test*100
rels_se_log5s  <- rels_se_log5_test*100
rels_se_log1s  <- rels_se_log1_test*100
rels_se_logits <- rels_se_logit_test*100


rels_sp_exp9s  <- rels_sp_exp9_test*100
rels_sp_exp5s  <- rels_sp_exp5_test*100
rels_sp_exp1s  <- rels_sp_exp1_test*100
rels_sp_log9s  <- rels_sp_log9_test*100
rels_sp_log5s  <- rels_sp_log5_test*100
rels_sp_log1s  <- rels_sp_log1_test*100
rels_sp_logits <- rels_sp_logit_test*100

ggcol <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }


plot( 1-100:1/1000, 1-rels_sp_exp9s/100  , col=ggcol(4)[1],
      type="l", xlab="True Positive Rate", ylab="False Positive Rate",
      main="Classification with high-sensitivity", ylim=c(0,0.6), cex.main=0.8)
lines(1-100:1/1000, 1-rels_sp_exp5s/100  , ylim = c(78,95), col=ggcol(4)[2])
lines(1-100:1/1000, 1-rels_sp_exp1s/100  , ylim = c(78,95), col=ggcol(4)[3])
lines(1-100:1/1000, 1-rels_sp_logits/100,                   col=ggcol(4)[4])
mtext("(Example 1)",cex=0.7)
legend("topleft", legend=c(expression(paste(lambda, " = ", 0.9)),
                           expression(paste(lambda, " = ", 0.5)), 
                           expression(paste(lambda, " = ", 0.1)), "logistic"),
       col=ggcol(4), lty=1, cex=0.8)

plot( 1-100:1/1000, 1-rels_se_exp9s/100  ,  col=ggcol(4)[1],
      type="l", xlab="True Negative Rate", ylab="False Negative Rate",
      main="Classification with high-specificity", cex.main=0.8)

lines(1-100:1/1000, 1-rels_se_exp5s/100  , ylim = c(64,88), col=ggcol(4)[2])
lines(1-100:1/1000, 1-rels_se_exp1s/100  , ylim = c(64,88), col=ggcol(4)[3])
lines(1-100:1/1000, 1-rels_se_logits/100,                   col=ggcol(4)[4])
mtext("(Example 1)",cex=0.7)
legend("topleft", legend=c(expression(paste(lambda, " = ", 0.9)),
                           expression(paste(lambda, " = ", 0.5)), 
                           expression(paste(lambda, " = ", 0.1)), "logistic"),
       col=ggcol(4), lty=1, cex=0.8)



par(mfrow=c(1,2), mar=c(4.1, 4.1, 3.1, 1.1))

ggcol <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }


plot( 1-100:1/1000, 1-rels_sp_log9s/100  , col=ggcol(4)[1],
      type="l", xlab="True Positive Rate", ylab="False Positive Rate",
      main="Classification with high-sensitivity", ylim=c(0,0.6), cex.main=0.8)
lines(1-100:1/1000, 1-rels_sp_log5s/100  , ylim = c(78,95), col=ggcol(4)[2])
lines(1-100:1/1000, 1-rels_sp_log1s/100  , ylim = c(78,95), col=ggcol(4)[3])
lines(1-100:1/1000, 1-rels_sp_logits/100,                   col=ggcol(4)[4])
mtext("(Example 2)",cex=0.7)
legend("topleft", legend=c(expression(paste(gamma, " = ", 0.9)),
                           expression(paste(gamma, " = ", 0.5)), 
                           expression(paste(gamma, " = ", 0.1)), "logistic"),
       col=ggcol(4), lty=1, cex=0.8)

plot( 1-100:1/1000, 1-rels_se_log9s/100  ,  col=ggcol(4)[1],
      type="l", xlab="True Negative Rate", ylab="False Negative Rate",
      main="Classification with high-specificity", cex.main=0.8)

lines(1-100:1/1000, 1-rels_se_log5s/100  , ylim = c(64,88), col=ggcol(4)[2])
lines(1-100:1/1000, 1-rels_se_log1s/100  , ylim = c(64,88), col=ggcol(4)[3])
lines(1-100:1/1000, 1-rels_se_logits/100,                   col=ggcol(4)[4])
mtext("(Example 2)",cex=0.7)
legend("topleft", legend=c(expression(paste(gamma, " = ", 0.9)),
                           expression(paste(gamma, " = ", 0.5)), 
                           expression(paste(gamma, " = ", 0.1)), "logistic"),
       col=ggcol(4), lty=1, cex=0.8)


auc_logit_test
auc_exp9_test
auc_exp5_test
auc_exp1_test
auc_log9_test
auc_log5_test
auc_log1_test

par(mfrow=c(1,1), mar=c(4.1, 4.1, 3.1, 1.1), 
    oma=c(0,0,0,0), mgp=c(2,0.8,0))

plot(1, 1, col = "white",
     xlab="False Positive Rate",
     ylab="True Positive Rate",
     xlim=c(0,1), ylim=c(0,1), main="ROC curves")
axis(side=1, at=(0:10)/10)
axis(side=2, at=(0:10)/10)
polygon(x = c(0.0, 1.0, 1.0, 0.0),                           # X-Coordinates of polygon
        y = c(0.9, 0.9, 1.0, 1.0),                             # Y-Coordinates of polygon
        col = "#cdcdcd", border=NA)  
lines(y=  roc_logit_test$sensitivities,
     x=1-roc_logit_test$specificities,
     xlab="False Positive Rate",
     ylab="True Positive Rate",
     type='l', col=ggcol(4)[4])
                                   # Color of polygon
abline(0,1, lty=2)
#abline(h=0.9, lty=4)
lines(y=  roc_exp9_test$sensitivities,
      x=1-roc_exp9_test$specificities,
      xlab="False Positive Rate",
      ylab="True Positive Rate",
      type='l' , col=ggcol(4)[1])
legend("bottomright", legend=c(expression(paste(lambda, " = ", 0.9)),
                           "logistic"),
       col=ggcol(4)[c(1,4)], lty=1, cex=0.8)

ind <- (1:10)*10-9

sptable <- rbind(
  cbind(rels_sp_exp9s[ind], rels_sp_exp5s[ ind],
        rels_sp_exp1s[ind], rels_sp_logits[ind]),
 c(mean(rels_sp_exp9s), mean(rels_sp_exp5s ),
   mean(rels_sp_exp1s), mean(rels_sp_logits)))


colnames(sptable) <- c("$\\lambda=0.9$", "$\\lambda=0.5$",
                       "$\\lambda=0.1$", "Logistic")
rownames(sptable) <- c("$90\\%$", "$91\\%$", "$92\\%$",
                       "$93\\%$", "$94\\%$", "$95\\%$",
                       "$96\\%$", "$97\\%$", "$98\\%$",
                       "$99\\%$", "Mean")

kable(sptable,digits=2,booktabs = T,
      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
      bottomrule="\\hhline{=====}",escape = FALSE) %>% kable_styling(latex_options = c("hold_position", "scale_down")) 


setable <- rbind(
  cbind(rels_se_exp9s[ind], rels_se_exp5s[ind],
        rels_se_exp1s[ind], rels_se_logits[ind]),
  c(mean(rels_se_exp9s), mean(rels_se_exp5s ),
    mean(rels_se_exp1s), mean(rels_se_logits)))

colnames(setable) <- c("$\\lambda=0.9$", "$\\lambda=0.5$",
                       "$\\lambda=0.1$", "Logistic")
rownames(setable) <- c("$90\\%$", "$91\\%$", "$92\\%$",
                       "$93\\%$", "$94\\%$", "$95\\%$",
                       "$96\\%$", "$97\\%$", "$98\\%$",
                       "$99\\%$", "Mean")

kable(setable,digits=2,booktabs = T,
      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
      bottomrule="\\hhline{=====}",escape = FALSE) %>% kable_styling(latex_options = c("hold_position", "scale_down")) 