Conditional AUC P-splines Estimate, three covariates case
================
Guiomar Pescador-Barrios

## Data

``` r
# Test results AD, MCI and CN respectively
yd <- ADNI$tau[ADNI$DX == 3]
ym <- ADNI$tau[ADNI$DX == 2]
yh <- ADNI$tau[ADNI$DX == 1]

# Age AD, MCI and CN respectively

xd <- ADNI$age[ADNI$DX == 3]
xm <- ADNI$age[ADNI$DX == 2]
xh <- ADNI$age[ADNI$DX == 1]

# Age AD, MCI and CN respectively
xd2 <- ADNI$MMSE[ADNI$DX == 3]
xm2 <- ADNI$MMSE[ADNI$DX == 2]
xh2 <- ADNI$MMSE[ADNI$DX == 1]

# APOE4 AD, MCI and CN respectively
xd3 <- ADNI$APOE4[ADNI$DX == 3]
xm3 <- ADNI$APOE4[ADNI$DX == 2]
xh3 <- ADNI$APOE4[ADNI$DX == 1]
```

## P-splines estimator implementation

``` r
ps_est_fun <- function(y, x1, x2, x3, x1_pred, x2_pred, x3_pred) {
  
  fit <- gam(y ~ s(x1, bs = "ps") +  s(x2, bs = "ps") + factor(x3))
  
  df_pred <- data.frame(x1 = x1_pred, x2 = x2_pred, x3 = x3_pred)
  
  mean_fitted_values <- fit$fitted.values
  mean_pred <- predict(fit, df_pred)
  
  y_wr <- log((y - mean_fitted_values)^2)
  
  fit_wr <- gam(y_wr ~ s(x1, bs = "ps") +  s(x2, bs = "ps") + factor(x3))
  
  # Step 3
  cov <- exp(fit_wr$fitted.values)
  theta_hat <- sum(((y - mean_fitted_values)^2) * cov)/sum((cov)^2)
  
  sigma2_fitted_values <- theta_hat*exp(fit_wr$fitted.values)
  
  sigma2_pred <- theta_hat * exp(predict(fit_wr, df_pred))
  
  return(list("mu_fitted" = mean_fitted_values, "sigma2_fitted" = sigma2_fitted_values,
              "mu_pred" = mean_pred, "sigma2_pred" = sigma2_pred))
}
```

``` r
roc_ps <- function(yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, x1_pred, x2_pred, x3_pred) {
  
  fit_h <- ps_est_fun(y = yh, x1 = x1h, x2 = x2h, x3 = x3h,
                      x1_pred = x1_pred, x2_pred = x2_pred, x3_pred = x3_pred)
  fit_d <- ps_est_fun(y = yd, x1 = x1d, x2 = x2d, x3 = x3d, 
                      x1_pred = x1_pred, x2_pred = x2_pred, x3_pred = x3_pred)
  
  mu_h <- fit_h$mu_pred
  sigma_h <- sqrt(fit_h$sigma2_pred)
  
  mu_d <- fit_d$mu_p
  sigma_d <- sqrt(fit_d$sigma2_pred)
  
  roc_cov_est <- matrix(0, nrow = length(p), ncol = length(x1_pred))
  auc_est <- numeric(length(x1_pred))
  for(j in 1:length(x1_pred)){
    roc_cov_est[, j] <- 1 - pnorm(((mu_h[j] - mu_d[j])/sigma_d[j]) + (sigma_h[j]/sigma_d[j])*qnorm(1 - p))
    auc_est[j] <- sum(roc_cov_est[, j])/length(p)
  }
  
  return(list("fit_h" = fit_h, "fit_d" = fit_d,
              "rocs" = roc_cov_est, "auc" = auc_est, 
              "mu_h" = mu_h, "mu_d" = mu_d, 
              "sigma_h" = sigma_h, "sigma_d" = sigma_d))
}
```

## Bootstrap intervals implementation

``` r
boot_fun <- function(b, yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, x1_pred, x2_pred, x3_pred) {
  B <- b
  auc_est_boot <- matrix(0, nrow = length(x1_pred), ncol = B)
  
  for(l in 1:B){
    ind_h <- sample(1:length(yh), size = length(yh), replace = TRUE)
    yh_boot <- yh[ind_h]  
    x1h_boot <- x1h[ind_h]
    x2h_boot <- x2h[ind_h]
    x3h_boot <- x3h[ind_h]
    
    ind_d <- sample(1:length(yd), size = length(yd), replace = TRUE)  
    yd_boot <- yd[ind_d]  
    x1d_boot <- x1d[ind_d]
    x2d_boot <- x2d[ind_d]
    x3d_boot <- x3d[ind_d]
    
    aux <- roc_ps(yd = yd_boot, x1d = x1d_boot, x2d = x2d_boot, x3d = x3d_boot,
                  yh = yh_boot, x1h = x1h_boot, x2h = x2h_boot, x3h = x3h_boot,
                  p = p, x1_pred = x1_pred, x2_pred = x2_pred, x3_pred = x3_pred)
    
    auc_est_boot[, l] <- aux$auc
  }
  
  auc_boot_l <- apply(auc_est_boot, 1, quantile, prob = 0.025)
  auc_boot_u <- apply(auc_est_boot, 1, quantile, prob = 0.975)
  
  return(list("auc_boot_l" = auc_boot_l, "auc_boot_u" = auc_boot_u))
  
}
```

``` r
boot_res_fun <- function(b, yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, 
                         x1_pred, x2_pred, x3_pred, roc_original_sample) {
  B <- b
  auc_est_boot_res <- matrix(0, nrow = length(x1_pred), ncol = B)
  
  std_res_d_original_sample <- (yd - roc_original_sample$fit_d$mu_fitted)/sqrt(roc_original_sample$fit_d$sigma2_fitted)
  std_res_h_original_sample <- (yh - roc_original_sample$fit_h$mu_fitted)/sqrt(roc_original_sample$fit_h$sigma2_fitted)
  
  fit_d <- ps_est_fun(y = yd, x1 = x1d, x2 = x2d, x3 = x3d, 
                      x1_pred = x1d, x2_pred = x2d, x3_pred = x3d)
  
  mu_d <- fit_d$mu_pred
  sigma_d <- sqrt(fit_d$sigma2_pred)
  
  fit_h <- ps_est_fun(y = yh, x1 = x1h, x2 = x2h, x3 = x3h,
                      x1_pred = x1h, x2_pred = x2h, x3_pred = x3h)
  
  mu_h <- fit_h$mu_pred
  sigma_h <- sqrt(fit_h$sigma2_pred)
  
  for(l in 1:B){
    std_res_d_boot <- sample(std_res_d_original_sample, length(yd), replace = TRUE)
    yd_boot <- mu_d + sigma_d*std_res_d_boot
    
    std_res_h_boot <- sample(std_res_h_original_sample, length(yh), replace = TRUE)
    yh_boot <- mu_h + sigma_h*std_res_h_boot
    
    aux <- roc_ps(yd = yd_boot, x1d = x1d, x2d = x2d, x3d = x3d,
                  yh = yh_boot, x1h = x1h, x2h = x2h, x3h = x3h,
                  p = p, x1_pred = x1_pred, x2_pred = x2_pred, x3_pred = x3_pred)
    
    auc_est_boot_res[, l] <- aux$auc
  }
  
  auc_boot_res_l <- apply(auc_est_boot_res, 1, quantile, prob = 0.025)
  auc_boot_res_u <- apply(auc_est_boot_res, 1, quantile, prob = 0.975)
  
  return(list("auc_boot_res_l" = auc_boot_res_l, "auc_boot_res_u" = auc_boot_res_u))
  
}
```

### Plotting AUC surfaces

``` r
plot_cov_fun_1 <- function(yd, x1d, x2d, x3d, yh, x1h, x2h, x3h,
                          x1_pred, x2_pred, ncat, cat) {
  # Define sequence p
  p <- seq(0, 1, len = 101)
  b <- 100
  
  for (i in n_cat) {
    x3_pred <- rep(i, length(x1_pred))
    
    roc_original_sample <- roc_ps(yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, x1_pred, x2_pred, x3_pred)
    
    matrix_auc <- matrix(roc_original_sample$auc, nrow = length(x1_pred), ncol = length(x1_pred))
    persp(x1_pred, x2_pred, matrix_auc,theta = 30, phi = 20, col = "pink",
          xlab="Age", ylab="MMSE", zlab="AUC", sub = cat[i+1])
  }
  
}
```

### Example

``` r
# Prediction values and levels of categorical covariate
x_p <- seq(65,85, by = 0.5)
x2_p <- seq(20,30, by = 0.25)
n_cat <- 0:2

# Plotting set-up
CAT <- c("APOE4 0", "APOE4 1", "APOE4 2")
par(cex.axis=1.5, cex.lab=1.2, cex.main=2, cex.sub=1.5, mfrow=c(1,3),  mar = c(5.5, 4.5, 6.5, 2.5))
plot_cov_fun_1(yd, xd, xd2, xd3, yh, xh, xh2, xh3, x_p, x2_p, n_cat, CAT)
title("Age-/MMSE/APOE4-specific AUC surface AD vs. CN", outer = T, line=-3)
```

![](README_figs/README-unnamed-chunk-7-1.png)<!-- -->

## Plotting for specific values of MMSE

``` r
plot_cov_fun_2 <- function(yd, x1d, x2d, x3d, yh, x1h, x2h, x3h,
                           x1_pred, x2_pred, ncat, cat) {
  # Define sequence p
  p <- seq(0, 1, len = 101)
  b <- 100

  for (i in n_cat) {
    x3_pred <- rep(i, length(x1_pred))
    
    roc_original_sample <- roc_ps(yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, x1_pred, x2_pred, x3_pred)
    
    boot_auc <- boot_fun(b, yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, x1_pred, x2_pred, x3_pred)
    boot_res_auc <- boot_res_fun(b, yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p,
     x1_pred, x2_pred, x3_pred, roc_original_sample)
    
    plot(x1_pred, roc_original_sample$auc, ylim = c(0, 1),
         type = "l", ylab = "AUC", xlab = "Age")
     lines(x1_pred, boot_auc$auc_boot_l, lty = 2, lwd = 1.5, col = "blue")
     lines(x1_pred, boot_auc$auc_boot_u, lty = 2, lwd = 1.5, col = "blue")

    lines(x1_pred, boot_res_auc$auc_boot_res_l, lty = 2, lwd = 1.5, col = "red")
    lines(x1_pred, boot_res_auc$auc_boot_res_u, lty = 2, lwd = 1.5, col = "red")
    text(82,0.1, cat[i+1], cex=1)
  }
}
```

### Example

``` r
par(mfrow = c(1,1))
n_cat <- 0
x_p <- seq(65,85, by = 0.5)
x2_p <- rep(26, length(x_p))
plot_cov_fun_2(yd, xd, xd2, xd3, yh, xh, xh2, xh3, x_p, x2_p, n_cat, CAT)
title("Age-specific AUC surface AD vs. CN, MMSE = 28")
```

![](README_figs/README-unnamed-chunk-9-1.png)<!-- -->

### Plotting for specific values of age

``` r
plot_cov_fun_3 <- function(yd, x1d, x2d, x3d, yh, x1h, x2h, x3h,
                           x1_pred, x2_pred, ncat, cat) {
  # Define sequence p
  p <- seq(0, 1, len = 101)
  b <- 100
  
  
  for (i in n_cat) {
    x3_pred <- rep(i, length(x1_pred))
    
    roc_original_sample <- roc_ps(yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, x1_pred, x2_pred, x3_pred)
    
    boot_auc <- boot_fun(b, yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p, x1_pred, x2_pred, x3_pred)
    boot_res_auc <- boot_res_fun(b, yd, x1d, x2d, x3d, yh, x1h, x2h, x3h, p,
                                 x1_pred, x2_pred, x3_pred, roc_original_sample)
    
    plot(x2_pred, roc_original_sample$auc, ylim = c(0, 1),
         type = "l", ylab = "AUC", xlab = "MMSE")
    lines(x2_pred, boot_auc$auc_boot_l, lty = 2, lwd = 1.5, col = "blue")
    lines(x2_pred, boot_auc$auc_boot_u, lty = 2, lwd = 1.5, col = "blue")
    
    lines(x2_pred, boot_res_auc$auc_boot_res_l, lty = 2, lwd = 1.5, col = "red")
    lines(x2_pred, boot_res_auc$auc_boot_res_u, lty = 2, lwd = 1.5, col = "red")
    text(28,0.1, cat[i+1], cex=1)
  }
}
```

### Example

``` r
par(mfrow = c(1,1))
n_cat <- 0
x2_p <- seq(22,30, by = 0.25)
x_p <- rep(65, length(x2_p))
plot_cov_fun_3(yd, xd, xd2, xd3, ym, xm, xm2, xm3, x_p, x2_p, n_cat, CAT)
title("Age-APOE4-MMSE-specific AUC, AD vs. MCI, Age = 65")
```

![](README_figs/README-unnamed-chunk-11-1.png)<!-- -->