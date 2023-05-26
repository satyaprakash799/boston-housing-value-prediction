library(MASS)
data(Boston)
help(Boston)
pairs(Boston)
Boston_log <- Boston
Boston_log$log_medv <- log(Boston$medv)
pairs(Boston_log[c("log_medv", "age", "rm", "zn", "dis", "rad", "crim", 
                   "ptratio", "chas", "indus", "nox", "tax", "black", "lstat")])
cor(Boston_log)
dim(Boston_log)
sum(is.na.data.frame(Boston_log))
sapply(Boston_log,class)
data_len = nrow(Boston_log)
pred_len = ncol(Boston_log)-1

library(leaps)
set.seed(799)
bss.fit = regsubsets(log_medv~.-medv, data = Boston_log, nvmax = 15)
bss.fit.summary = summary(bss.fit) ; bss.fit.summary

par(mfrow=c(2,2))
plot(bss.fit.summary$rss, xlab="Number of Variables",ylab="RSS", type = "l")
which.min(bss.fit.summary$rss)
points(13, bss.fit.summary$rss[13], col ="red",cex =2, pch =20)

plot(bss.fit.summary$adjr2, xlab="Number of Variables",ylab="Adj-R2", type = "l")
which.max(bss.fit.summary$adjr2)
points(12, bss.fit.summary$adjr2[12], col ="red",cex =2, pch =20)

plot(bss.fit.summary$cp, xlab="Number of Variables",ylab="cp", type = "l")
which.min(bss.fit.summary$cp)
points(11, bss.fit.summary$cp[11], col ="red",cex =2, pch =20)

plot(bss.fit.summary$bic, xlab="Number of Variables",ylab="bic", type = "l")
which.min(bss.fit.summary$bic)
points(10, bss.fit.summary$bic[10], col ="red",cex =2, pch =20)

coef(bss.fit, 11)
coef(bss.fit, 10)


# Predict Function for regsubsets()
predict.regsubsets = function(object, newdata, id,...) {
  form = as.formula(object$call[[2]])
  mat = model.matrix(form,newdata)
  coefi = coef(object, id=id)
  xvars = names(coefi)
  mat[,xvars]%*%coefi
}

#dimnames = list(NULL, paste(1:pred_len))
## Cross-Validation using k-fold for Best Subset Selection
k=10
set.seed(799)
folds = sample(1:k, data_len, replace = TRUE)
cv.errors = matrix(NA,k,(pred_len-1))
for (j in 1:k) {
  set.seed(799)
  bss.fit_2 = regsubsets(log_medv~.-medv, data = Boston_log[folds!=j,], 
                         nvmax = pred_len)
  for (i in 1:(pred_len-1)) {
    pred = predict(bss.fit_2, Boston_log[folds==j,], id=i)
    cv.errors[j,i] = mean((Boston_log$log_medv[folds==j]-pred)^2)
  }
}
cv.errors
mean.cv.errors = apply(cv.errors,2,mean)
par(mfrow=c(1,1))
plot(mean.cv.errors, type='b', xlab="Number of Predictors",
     ylab="10-Fold Cross-Validation Error")
min(mean.cv.errors)
which.min(mean.cv.errors)
points(12, mean.cv.errors[12], col ="red",cex =2, pch =20)
coef(bss.fit, 12)
mean.cv.errors

## Best names
bss_coef = coef(bss.fit, 10)
bss_coef_noint = bss_coef[bss_coef[] != 0][-1]
bss_names = names(bss_coef_noint)
bss_params_sum = paste(bss_names, collapse = "+")
bss_params_sum



##### Lasso Regression
x = model.matrix(log_medv~.-medv, data = Boston_log,)[,-1]
y = Boston_log$log_medv
grid = 10^seq(10,-2,length=100)
library(glmnet)

lasso_fit = glmnet(x,y, alpha = 1, lambda = grid)
par(mfrow=c(1,1))
plot(lasso_fit)
set.seed(799)
cv.out = cv.glmnet(x, y, alpha=1)
par(mfrow=c(1,2))
plot(cv.out)
plot(cv.out$lambda,cv.out$cvm, xlab = "Lambda", ylab = "Mean Cross-Validated Error")
bestlam = cv.out$lambda.min ; bestlam
lasso_coef = predict(lasso_fit, type = "coefficients", s=bestlam)[1:14,]
lasso_coef


## Lasso Names
lasso_coef_noint = lasso_coef[lasso_coef[] != 0][-1]
lasso_coef_noint
lasso_names = names(lasso_coef_noint)
lasso_params_sum <- paste(lasso_names, collapse = "+")
lasso_params_sum
lasso_names


## Cross-Validation using k-fold for lasso regression
k=10
set.seed(799)
folds = sample(1:k, data_len, replace = TRUE)
cv.errors_l = rep(NA,k)
for (j in 1:k) {
  x_l = model.matrix(log_medv~.-medv, data = Boston_log[folds!=j,],)[,-1]
  y_l = Boston_log[folds!=j,]$log_medv
  lasso.fit_l = glmnet(x_l,y_l, alpha = 1, lambda = grid)
  set.seed(799)
  cv.out_l = cv.glmnet(x_l, y_l, alpha=1)
  bestlam_l = cv.out_l$lambda.min
  x_l_test = model.matrix(log_medv~.-medv, data = Boston_log[folds==j,],)[,-1]
  y_l_test = Boston_log[folds==j,]$log_medv
  lasso.pred_l = predict(lasso.fit_l, s=bestlam_l, newx = x_l_test)
  cv.errors_l[j] = mean((lasso.pred_l-y_l_test)^2)
}

cv.errors_l
mean_cv.errors_l = mean(cv.errors_l)
mean_cv.errors_l



# calculate R2 and Adj-R2
predicted_values = predict(lasso_fit, s=bestlam, newx = x)
residuals = traindata_omitna$Q46_2 - predicted_values
RSS = sum(residuals^2)
TSS = sum((traindata_omitna$Q46_2 - mean(traindata_omitna$Q46_2))^2)
R.2 = 1-(RSS/TSS) ; R.2
adj.r.2 = 1-((RSS/(479-21-1))*(479-1)/TSS) ; adj.r.2
mean(residuals^2)



##### Final Model
lm.fit = lm(log_medv~crim+chas+nox+rm+dis+rad+tax+ptratio+black+lstat
            , data = Boston_log)
summary(lm.fit)
par(mfrow=c(2,2))
plot(lm.fit)
library(car)
vif(lm.fit)