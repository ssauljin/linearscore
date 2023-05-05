coeftables <- NULL

par(mfrow=c(1,2), mar=c(4.1, 4.1, 3.1, 1.1), 
    oma=c(0,0,0,0), mgp=c(2,0.8,0))

#### LGPIF data ####
library(pROC)
library(MASS)
library(nleqslv)
library(numDeriv)

load("data.RData")



#par(mfrow=c(1,2), mar=c(4.1, 4.1, 3.1, 1.1))

data$FreqBC[data$FreqBC>10] <- 10
FreqBCcount <- table(data$FreqBC)
rownames(FreqBCcount)[11] <- "10+"

FreqBCprop <- FreqBCcount/sum(FreqBCcount)*100

#barplot(FreqBCcount, main="Distribution of FreqBC",
#        xlab="Count")

R <- 1

data <- data[,c(9:13, 16:18, 20)]
data$Y <- (data$FreqBC>4)*1
data$FreqBC <- NULL

load("dataout.RData")

test <- dataout[,c(9:13, 16:18, 20)]
test$Y <- (test$FreqBC>4)*1
test$FreqBC <- NULL
rm(dataout)

rels_se_logit_test <- rep(0, 100)
rels_se_exp9_test  <- rep(0, 100)
rels_se_exp5_test  <- rep(0, 100)
rels_se_exp1_test  <- rep(0, 100)

rels_sp_logit_test <- rep(0, 100)
rels_sp_exp9_test  <- rep(0, 100)
rels_sp_exp5_test  <- rep(0, 100)
rels_sp_exp1_test  <- rep(0, 100)

auc_logit_test <- 0
auc_exp9_test  <- 0
auc_exp5_test  <- 0
auc_exp1_test  <- 0

for (i in 1:R) {
  
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
  
  fit_logit_test <- glm(Y ~ ., family=binomial(link = "logit"), data=test)
  tS       <- model.matrix(fit_logit_test)
  
  roc_logit_test <- roc(test$Y ~ predict(fit_logit_loss, test))
  roc_exp9_test  <- roc(test$Y ~ as.vector(tS %*% fit_exp9_eqn$x))
  roc_exp5_test  <- roc(test$Y ~ as.vector(tS %*% fit_exp5_eqn$x))
  roc_exp1_test  <- roc(test$Y ~ as.vector(tS %*% fit_exp1_eqn$x))

  auc_logit_test <- auc_logit_test + auc(roc_logit_test)/R
  auc_exp9_test  <- auc_exp9_test  + auc(roc_exp9_test )/R
  auc_exp5_test  <- auc_exp5_test  + auc(roc_exp5_test )/R
  auc_exp1_test  <- auc_exp1_test  + auc(roc_exp1_test )/R

  
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
  
  rels_sp_exp5_test <- rels_sp_exp5_test +
    as.vector(coords(roc=roc_exp5_test, x=(900:999)/1000,
                     input="sensitivity", ret="specificity")[,1])/R
  rels_se_exp5_test <- rels_se_exp5_test +
    as.vector(coords(roc=roc_exp5_test, x=(900:999)/1000,
                     ret="sensitivity", input="specificity")[,1])/R
  
  rels_sp_exp1_test <- rels_sp_exp1_test +
    as.vector(coords(roc=roc_exp1_test, x=(900:999)/1000,
                     input="sensitivity", ret="specificity")[,1])/R
  rels_se_exp1_test <- rels_se_exp1_test +
    as.vector(coords(roc=roc_exp1_test, x=(900:999)/1000,
                     ret="sensitivity", input="specificity")[,1])/R
  
}


rels_se_exp9s  <- rels_se_exp9_test*100
rels_se_exp5s  <- rels_se_exp5_test*100
rels_se_exp1s  <- rels_se_exp1_test*100
rels_se_logits <- rels_se_logit_test*100

rels_sp_exp9s  <- rels_sp_exp9_test*100
rels_sp_exp5s  <- rels_sp_exp5_test*100
rels_sp_exp1s  <- rels_sp_exp1_test*100
rels_sp_logits <- rels_sp_logit_test*100


ggcol <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }

plot( 1-100:1/1000, 1-rels_sp_exp9s/100  , ylim = c(0.1,1), col=ggcol(4)[1],
      type="l", xlab="True Positive Rate", ylab="False Positive Rate",
      main="Classification with high-sensitivity")
lines(1-100:1/1000, 1-rels_sp_exp5s/100  , ylim = c(64,88), col=ggcol(4)[2])
lines(1-100:1/1000, 1-rels_sp_exp1s/100  , ylim = c(64,88), col=ggcol(4)[3])
lines(1-100:1/1000, 1-rels_sp_logits/100,                   col=ggcol(4)[4])
mtext("(h=4)",cex=0.9)
legend("topleft", legend=c(expression(paste(lambda, " = ", 0.9)),
                          expression(paste(lambda, " = ", 0.5)),
                          expression(paste(lambda, " = ", 0.1)), "logistic"),
      col=ggcol(4), lty=1, cex=0.8)

# plot( 1-100:1/1000, 1-rels_se_exp9s/100  , ylim = c(0.1,1), col=ggcol(4)[1],
#       type="l", xlab="True Negative Rate", ylab="False Positive Rate",
#       main="Classification with high-specificity")
# lines(1-100:1/1000, 1-rels_se_exp5s/100  , ylim = c(78,95), col=ggcol(4)[2])
# lines(1-100:1/1000, 1-rels_se_exp1s/100  , ylim = c(78,95), col=ggcol(4)[3])
# lines(1-100:1/1000, 1-rels_se_logits/100,                   col=ggcol(4)[4])
# legend("topleft", legend=c(expression(paste(lambda, " = ", 0.9)),
#                           expression(paste(lambda, " = ", 0.5)),
#                           expression(paste(lambda, " = ", 0.1)), "logistic"),
#       col=ggcol(4), lty=1, cex=0.8)


library(knitr)
library(kableExtra)
options(knitr.kable.NA = '')
options(knitr.table.format = "latex")

coeftable <- cbind(
  fit_exp9_eqn$x, fit_exp5_eqn$x,
  fit_exp1_eqn$x, fit_logit_loss$coefficients)


colnames(coeftable) <- c("$\\lambda=0.9$", "$\\lambda=0.5$",
                         "$\\lambda=0.1$", "Logistic")

kable(coeftable,digits=4,booktabs = T,
      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
      bottomrule="\\hhline{=====}",escape = FALSE) %>% kable_styling(latex_options = c("hold_position", "scale_down")) 

coeftables <- cbind(coeftables, coeftable)



round(c(auc_exp9_test, auc_exp5_test,
        auc_exp1_test, auc_logit_test),4)


kable(coeftables,digits=4,booktabs = T,
      linesep = c("", "", "","",  "","","", "","", "", "","\\hline"),
      bottomrule="\\hhline{=================}",escape = FALSE) %>% 
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  add_header_above(c(" "=1, "\\$d=1$" = 4,"$d=2$" = 4,"$d=3$" = 4,"$d=4$" = 4))