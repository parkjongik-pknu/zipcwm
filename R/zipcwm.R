# ============================================================
# zip_cwm.R
# ============================================================


library(mvtnorm) # mutivariate gaussian dist. for covariates
library(cluster) # k-medoids

zipcwm <- function(X, Z, Y, 
                   init_method = c("kmeans", "kmedoids", "hierarchical"),
                   max_iter=100, tol=1e-6, irls_max_iter=20) {
  
  init_method <- match.arg(init_method)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  Y <- as.numeric(Y)
  N <- length(Y)
  K <- 2
  
  # include intercept
  X_in <- cbind(1, X)
  Z_in <- cbind(1, Z)
  Pin_x <- ncol(X_in)
  Pin_z <- ncol(Z_in)
  
  # # combination of covariates X and Z
  All_Covariates <- cbind(X, Z[, !(colnames(Z) %in% colnames(X)), drop=FALSE])
  Pall <- ncol(All_Covariates)
  init_data <- All_Covariates
  
  # ------------------------------------------------------------
  # initialization(K=2)
  # ------------------------------------------------------------
  set.seed(123)
  initial_clusters <- NULL
  
  if (init_method == "kmeans") {
    km <- kmeans(init_data, centers = K, nstart = 20)
    initial_clusters <- km$cluster} 
  
  else if (init_method == "kmedoids") {
    pam_res <- pam(init_data, k = K)
    initial_clusters <- pam_res$clustering} 
  
  else if (init_method == "hierarchical") {
    dist_mat <- dist(init_data)
    hc <- hclust(dist_mat, method = "ward.D2")
    initial_clusters <- cutree(hc, k = K)
  }
  
  pi_k    <- rep(1/K, K)
  mu_k    <- list()
  sigma_k <- list()
  beta_k  <- matrix(0, Pin_x, K)
  gamma_k <- matrix(0, Pin_z, K)
  
  for(k in 1:K) {
    idx <- which(initial_clusters == k)
    if(length(idx) < 5) idx <- sample(1:N, 10) 
    
    mu_k[[k]] <- colMeans(All_Covariates[idx, , drop=FALSE])
    sigma_k[[k]] <- cov(All_Covariates[idx, , drop=FALSE]) + diag(1e-5, Pall)
    
    beta_k[, k]  <- coef(glm(Y[idx] ~ X[idx, ], family = poisson))
    gamma_k[, k] <- coef(glm(as.numeric(Y[idx] == 0) ~ Z[idx, ], family = binomial))
  }
  
  ll_history <- numeric()
  
  # ------------------------------------------------------------
  # em algorithm
  # ------------------------------------------------------------
  for(iter in 1:max_iter) {
    
    # --- E-step ---
    tau <- matrix(0, N, K) # posterior probability of cwm
    z_hat_2_list <- matrix(0, N, K) 
    
    for(k in 1:K) {
      f_covs <- dmvnorm(All_Covariates, mean = mu_k[[k]], sigma = sigma_k[[k]]) # phi_p(x;mu_k, sigma_k)
      
      mu_y <- exp(X_in %*% beta_k[, k])
      pi_inv <- 1 / (1 + exp(-(Z_in %*% gamma_k[, k])))
      
      f_y_zip <- ifelse(Y == 0, # phi_1(y; beta_k, gamma_k)
                        pi_inv + (1 - pi_inv) * exp(-mu_y),
                        (1 - pi_inv) * (exp(-mu_y) * mu_y^Y / (factorial(Y) + 1e-15)))
      
      tau[, k] <- pi_k[k] * f_covs * f_y_zip # 분자 : pi_k * phi_p * phi_1
      
      z_h2 <- rep(1, N) # poisson probability to update beta
      z_h2[Y == 0] <- ((1 - pi_inv[Y == 0]) * exp(-mu_y[Y == 0])) / (f_y_zip[Y == 0] + 1e-15)
      z_hat_2_list[, k] <- z_h2
    }
    
    # 분모 : sum_{j=1}^{2} pi_j * phi_p * phi_1
    row_sums <- rowSums(tau)
    curr_ll <- sum(log(row_sums + 1e-15))
    tau <- tau / (row_sums + 1e-15)
    
    ll_history <- c(ll_history, curr_ll)
    if(iter > 1 && abs(curr_ll - ll_history[iter-1]) < tol) break
    
    # --- M-step ---
    for(k in 1:K) {
      W_k <- tau[, k]
      N_k <- sum(W_k)
      
      pi_k[k] <- N_k / N
      mu_k[[k]] <- colSums(W_k * All_Covariates) / N_k
      diff <- sweep(All_Covariates, 2, mu_k[[k]])
      sigma_k[[k]] <- (t(diff) %*% (W_k * diff)) / N_k + diag(1e-6, Pall)
      
      # Beta Update (Count Part) - X_in 사용
      z_h2 <- z_hat_2_list[, k]
      for(s in 1:irls_max_iter) {
        mu_curr <- exp(X_in %*% beta_k[, k])
        W_beta <- as.vector(W_k * z_h2 * mu_curr)
        W_beta <- pmax(W_beta, 1e-10)
        score_b <- t(X_in) %*% (W_k * z_h2 * (Y - mu_curr))
        info_b  <- t(X_in) %*% (W_beta * X_in)
        step_b <- solve(info_b + diag(1e-8, Pin_x)) %*% score_b
        beta_k[, k] <- beta_k[, k] + as.vector(step_b)
        if(max(abs(step_b)) < 1e-6) break
      }
      
      # Gamma Update (Inflation Part) - Z_in 사용
      z_h1 <- 1 - z_h2
      for(s in 1:irls_max_iter) {
        pi_curr <- 1 / (1 + exp(-(Z_in %*% gamma_k[, k])))
        W_gamma <- as.vector(W_k * pi_curr * (1 - pi_curr))
        W_gamma <- pmax(W_gamma, 1e-10)
        score_g <- t(Z_in) %*% (W_k * (z_h1 - pi_curr))
        info_g  <- t(Z_in) %*% (W_gamma * Z_in)
        step_g <- solve(info_g + diag(1e-8, Pin_z)) %*% score_g
        gamma_k[, k] <- gamma_k[, k] + as.vector(step_g)
        if(max(abs(step_g)) < 1e-6) break
      }
    }
  }
  
  return(list(pi = pi_k, mu = mu_k, sigma = sigma_k, 
              beta = beta_k, gamma = gamma_k, loglik = ll_history,
              init_method = init_method))
}

