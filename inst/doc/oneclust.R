## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  fig.asp = 10 / 16,
  eval = TRUE
)

## ---- message=FALSE-----------------------------------------------------------
library("oneclust")

## -----------------------------------------------------------------------------
df_levels <- sim_postcode_levels(nlevels = 500, seed = 42)
train <- sim_postcode_samples(df_levels, n = 100000, threshold = 3000, prob = c(0.2, 0.1), seed = 43)
test <- sim_postcode_samples(df_levels, n = 100000, threshold = 3000, prob = c(0.2, 0.1), seed = 44)

## -----------------------------------------------------------------------------
head(df_levels)
head(train)
head(test)

## -----------------------------------------------------------------------------
k <- 32
level_hist <- table(train$postcode)
level_new <- oneclust(level_hist, k)$cluster
feature_tr <- level_new[match(train$postcode, names(level_hist))] %>%
  as.character() %>%
  ordered(levels = as.character(1:k))

## -----------------------------------------------------------------------------
op <- par(las = 1)
plot(feature_tr, train$label, lty = 0, xlab = "Cluster", ylab = "Label")
abline(h = 0.2, col = cud(1))
abline(h = 0.1, col = cud(2))
par(op)

## -----------------------------------------------------------------------------
sum(train$is_rare)
sum(table(feature_tr)[1:5])

## -----------------------------------------------------------------------------
feature_te <- level_new[match(test$postcode, names(level_hist))] %>%
  as.character() %>%
  ordered(levels = as.character(1:k))

## -----------------------------------------------------------------------------
op <- par(las = 1)
plot(feature_te, test$label, lty = 0, xlab = "Cluster", ylab = "Label")
abline(h = 0.2, col = cud(1))
abline(h = 0.1, col = cud(2))
par(op)

## -----------------------------------------------------------------------------
sum(test$is_rare)
sum(table(feature_te)[1:5])

## -----------------------------------------------------------------------------
set.seed(42)
n <- 100
i <- 1:n
y <- (i > 20 & i < 30) + 5 * (i > 50 & i < 70) + rnorm(n, sd = 0.1)

## -----------------------------------------------------------------------------
# # If genlasso is available:
# out <- genlasso::fusedlasso1d(y)
out <- readRDS("out.rds")

## -----------------------------------------------------------------------------
# beta1 <- coef(out, lambda = 1.5)$beta
beta1 <- readRDS("beta1.rds")
plot(beta1)
abline(h = 0)

## -----------------------------------------------------------------------------
# beta2 <- genlasso::softthresh(out, lambda = 1.5, gamma = 1)
beta2 <- readRDS("beta2.rds")
grp <- as.integer(beta2 != 0) + 1L
plot(beta2, col = cud(grp))
abline(h = 0)
legend("topleft", legend = c("Zero", "Non-zero"), col = cud(unique(grp)), pch = 1)

## -----------------------------------------------------------------------------
cl1 <- oneclust(beta1, k = 2)$cluster
plot(beta1, col = cud(cl1))
abline(h = 0)
legend("topleft", legend = paste("Cluster", unique(cl1)), col = cud(unique(cl1)), pch = 1)

## -----------------------------------------------------------------------------
cl2 <- oneclust(beta1, k = 3)$cluster
plot(beta1, col = cud(cl2))
abline(h = 0)
legend("topleft", legend = paste("Cluster", unique(cl2)), col = cud(unique(cl2)), pch = 1)

## -----------------------------------------------------------------------------
cl3 <- oneclust(beta1, k = 5, sort = FALSE)$cluster
plot(beta1, col = cud(cl3))
abline(h = 0)
legend("topleft", legend = paste("Cluster", unique(cl3)), col = cud(unique(cl3)), pch = 1)

## -----------------------------------------------------------------------------
x <- seq(0, 1, len = 1024)
pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wdt <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)

psignal <- numeric(length(x))
for (i in seq(along = pos)) {
  psignal <- psignal + hgt[i] / (1 + abs((x - pos[i]) / wdt[i]))^4
}

plot(psignal, type = "l")

## -----------------------------------------------------------------------------
cl <- oneclust(psignal, k = 2)
plot(psignal, type = "h", col = cud(cl$cluster))
legend("topright", legend = paste("Cluster", unique(cl$cluster)), col = cud(unique(cl$cluster)), lty = 1)

## -----------------------------------------------------------------------------
cl <- oneclust(psignal, k = 4)
plot(psignal, type = "h", col = cud(cl$cluster + 2))
legend("topright", legend = paste("Cluster", unique(cl$cluster)), col = cud(unique(cl$cluster + 2)), lty = 1)

## -----------------------------------------------------------------------------
cl <- oneclust(psignal, k = 6, sort = FALSE)
plot(psignal, type = "h", col = cud(cl$cluster))
legend("topright", legend = paste("Cluster", unique(cl$cluster)), col = cud(unique(cl$cluster)), lty = 1, cex = 0.8)

