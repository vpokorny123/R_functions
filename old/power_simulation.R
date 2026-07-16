library(nnet)
library(MASS)   # for mvrnorm

# Generate N x P continuous predictors ~ multivariate normal (mean 0, SD 1)
gen_continuous_X <- function(N, P=11, rho = 0, sd = 1) {
  if (P == 1) {
    X <- matrix(rnorm(N, 0, sd), ncol=1)
  } else {
    Sigma <- matrix(rho, nrow=P, ncol=P)
    diag(Sigma) <- 1
    X <- mvrnorm(N, mu = rep(0, P), Sigma = Sigma)
  }
  colnames(X) <- paste0("X", 1:P)
  as.data.frame(scale(X))  # center & scale -> mean 0, sd 1
}

# Turn betas (3 x P) and X into probabilities and sample outcomes
sim_multinom_cont <- function(N, base_probs, betas, X) {
  K <- length(base_probs)
  stopifnot(K == 4)
  stopifnot(nrow(betas) == K-1, ncol(betas) == ncol(X))
  Xmat <- as.matrix(X)
  eta <- Xmat %*% t(betas)                    # N x 3
  base_ratio <- base_probs[-1] / base_probs[1]
  denom <- 1 + rowSums(sweep(exp(eta), 2, base_ratio, "*"))
  p1 <- 1/denom
  pk_rest <- sweep(exp(eta), 1, denom, "/") * rep(base_ratio, each=nrow(eta))
  Pmat <- cbind(p1, pk_rest)
  Y <- apply(Pmat, 1, function(p) sample.int(4,1,prob=p))
  data.frame(y = factor(Y), X)
}

# Power function: test_predictors is numeric vector of column indices to test
power_multinom_cont <- function(N, base_probs, betas, P=11, test_predictors = 1,
                                nsim=500, alpha=0.05, rho=0, trace=FALSE) {
  pvals <- numeric(nsim)
  for (i in 1:nsim) {
    X <- gen_continuous_X(N, P=P, rho=rho)
    dat <- sim_multinom_cont(N, base_probs, betas, X)
    fit_full <- suppressWarnings(multinom(y ~ ., data=dat, trace=FALSE))
    # build reduced formula by removing tested predictors
    remove_names <- paste0("X", test_predictors)
    rhs_terms <- setdiff(names(dat)[names(dat) != "y"], remove_names)
    if (length(rhs_terms) == 0) {
      fit_red <- suppressWarnings(multinom(y ~ 1, data=dat, trace=FALSE))
    } else {
      fit_red <- suppressWarnings(multinom(as.formula(paste("y ~", paste(rhs_terms, collapse = " + "))),
                                           data=dat, trace=FALSE))
    }
    LRT <- 2*(logLik(fit_full) - logLik(fit_red))
    df <- length(test_predictors) * (nrow(betas))
    pvals[i] <- 1 - pchisq(as.numeric(LRT), df=df)
    if (trace && i %% 50 == 0) message("sim ", i)
  }
  mean(pvals < alpha)
}

# --------------------------
# Example usage (replace base_probs/betas as desired)
set.seed(2025)
P <- 11
base_probs <- c(0.25, 0.25, 0.25, 0.25)
# create betas: rows = 3 logits, cols = 11 predictors (log-odds per 1 SD)
                
betas <- matrix(c(log(.76),0,0, #KB
                  log(1.41),0,0, # PRL
                  log(1.33),0,0, # SWS
                  0,0,0, # HR
                  log(1.29),0,0, #CFT
                  log(.54),log(.49),log(.59), #MFT
                  0,0,log(.75), #ebbinghaus
                  0,0,0, #CH
                  0,0,0, #EEFRT
                  0,0,0, #DD
                  log(1.25),log(1.25),log(1.25) #GLIAT
                  ), nrow=3, ncol=P)
# make predictor 1 the main exposure (moderate effects)
#betas[,1] <- c(log(.4), log(.4), log(.4))
# small background effects on a couple predictors
#betas[1,2] <- log(.4)
#betas[2,3] <- log(.4)

# quick check: power for predictor 1 at N = 600
#power_single <- power_multinom_cont(N=600, base_probs=base_probs, betas=betas, P=P,
#                                    test_predictors = 11, nsim=400, alpha=0.05, rho=0)
#cat("power (predictor 1) at N=600:", round(power_single,3), "\n")

# grid over candidate sample sizes
Ns <- c(800, 1440, 1800)
res <- sapply(Ns, function(n) power_multinom_cont(n, base_probs, betas, P=P,
                                                  test_predictors=11, nsim=300, alpha=0.05, rho=0))
data.frame(N=Ns, power=round(res,3))

