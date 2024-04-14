# Fast SparCC-derived approximation of absolute abundance variances
Absolute_abn_var <- function(OTU.abn.log, var_min = 1e-5){
  p <- ncol(OTU.abn.log)
  x_var <- apply(OTU.abn.log, 2, var)
  log_sum <- rowSums(OTU.abn.log)
  t <- p*x_var + sum(x_var) - 2*apply(OTU.abn.log, 2, function(x) cov(x, log_sum))
  M <- matrix(1,p,p)
  diag(M) <- p-1
  alpha <- solve(M, t)
  alpha[alpha < var_min] <- var_min
  return(alpha)
}

# SparCEV_base_internal <- function(OTU.abn, phenotype, pseudo_count = 1, variances = NULL, idx_include = rep(T, ncol(OTU.abn)), var_min = 1e-5){
#   p <- ncol(OTU.abn)
#   Rc <- sum(idx_include)
#   if(Rc < 2){
#     warning("Too few correlations were below the selected threshold. Consider increasing the value of cutoff.")
#     return(rep(NA, p))
#   }
#   # Estimating relative abundances
#   OTU_TSS <- (OTU.abn+pseudo_count)/rowSums(OTU.abn+pseudo_count)
#   log_OTU <- log(OTU_TSS)
#   log_sum <- rowSums( log_OTU[,idx_include] )
#   
#   if(is.null(variances)){
#     variances <- Absolute_abn_var(log_OTU, var_min = var_min)
#   }
#   
#   # Approximation of correlations between abundances and phenotype
#   log_ratio_cov <- apply(log_OTU, 2, function(x) cov( Rc*x - log_sum, phenotype ) )
#   phenotype_sd <- sd(phenotype)
#   cors <- log_ratio_cov/(sqrt(variances)*phenotype_sd*(Rc - as.numeric(idx_include)) )
#   cors[abs(cors)>1] <- sign(cors[abs(cors)>1])
#   return(cors)
# }

SparCEV_base_internal <- function(OTU.abn, phenotype, pseudo_count = 1,
                                  variances = NULL, idx_include, var_min){
  p <- ncol(OTU.abn)
  Rc <- sum(idx_include)
  if(Rc < 2){
    warning("Too few correlations were below the selected threshold. Consider increasing the value of t.")
    return(rep(NA, p))
  }
  # Estimating relative abundances
  OTU_TSS <- (OTU.abn+pseudo_count)/rowSums(OTU.abn+pseudo_count)
  log_OTU <- log(OTU_TSS)
  log_sum <- rowSums( log_OTU[,idx_include] )

  if(is.null(variances)){
    variances <- Absolute_abn_var(log_OTU, var_min = var_min)
  }

  # Approximation of correlations between abundances and phenotype
  log_ratio_cov <- apply(log_OTU, 2, function(x) cov( Rc*x - log_sum, phenotype ) )
  phenotype_sd <- sd(phenotype)
  cors <- log_ratio_cov/(sqrt(variances)*phenotype_sd*(Rc - as.numeric(idx_include)) )
  cors[abs(cors)>1] <- sign(cors[abs(cors)>1])
  return(cors)
}

#' Estimation of Cross-correlations between OTU abundances (or components of other
#' compositinal datasets) and a non-compositional varible. This method is appropriate
#' when the average covariance is close to zero.
#' @param OTU.abn An OTU table with observarions as rows and OTUs columns
#' (or a similar compositional dataset).
#' @param phenotype Continuous non-compositional phenotypic variable.
#' @param pseudo_count The pseudo-count added to all read counts to avoid division by zero.
#' @param var_min In the event that a variance estimate of an OTU is negative,
#' it is replaced by this number.
#' @param Find_m Toggles whether or not the permutation threshold, m, should be 
#' calculated and returned. All correlations above m are likely to be real signals.
#' @param B_m The number of bootstraps to use when calculating m. Only used when Find_m = TRUE.
#' @param cores The number of cores used for finding m.
#' @return A vector containing the correlations between the OTU abundances and the phenotypic variable.
#' @export
SparCEV_base <- function(OTU.abn, phenotype, pseudo_count = 1, var_min = 1e-5,
                         cores = 1, Find_m = TRUE, B_m = 100){
  if(nrow(OTU.abn) != length(phenotype)){
    stop("The number of rows in OTU.abn and the length of phenotype mustb be identical.")
  }

  p <- ncol(OTU.abn)
  # Calling internal version of SparCEV_base
  cors <- SparCEV_base_internal(OTU.abn = OTU.abn,
                                phenotype = phenotype,
                                pseudo_count = pseudo_count,
                                var_min = var_min,
                                idx_include = rep(T, p))
  
  if(Find_m){
    if(cores == 1){
      bt_list <- lapply(1:B_m,
                        function(x){
                          idx <- sample(1:nrow(OTU.abn))
                          SparCEV_base_internal(OTU.abn = OTU.abn,
                                                phenotype = phenotype[idx],
                                                pseudo_count = pseudo_count,
                                                var_min = var_min,
                                                idx_include = rep(T, p))})
    }
    if(cores > 1){
      bt_list <- parallel::mclapply(1:B_m,
                                    function(x){
                                      idx <- sample(1:nrow(OTU.abn))
                                      SparCEV_base_internal(OTU.abn = OTU.abn,
                                                            phenotype = phenotype[idx],
                                                            pseudo_count = pseudo_count,
                                                            var_min = var_min,
                                                            idx_include = rep(T, p))
                                    }, mc.cores = cores )
    }
    m <- mean(unlist(lapply(bt_list, function(x) max(abs(x)))))
    Correlated <- abs(cors) > m
  }
  
  if(!Find_m){
    m <- NULL
    Correlated <- NULL
  }

  out <- list(cor = cors,
              m = m,
              Correlated = Correlated)
  
  return(out)
}

SparCEV_iter_inner <- function(OTU.abn, phenotype, pseudo_count = 1, var_min = 1e-5, iter, t){
  p <- ncol(OTU.abn)
  
  OTU_TSS <- (OTU.abn+pseudo_count)/rowSums(OTU.abn+pseudo_count)
  log_OTU <- log(OTU_TSS)
  alpha <- Absolute_abn_var(log_OTU, var_min = var_min)
  
  idx_include_current <- rep(T, p)
  idx_include_old <- rep(F, p)
  current_iter <- 1
  while(!all(idx_include_current == idx_include_old) & current_iter <= iter){
    Spar_cor <- SparCEV_base_internal(OTU.abn = OTU.abn,
                                      phenotype = phenotype,
                                      pseudo_count = pseudo_count,
                                      variances = alpha,
                                      var_min = var_min,
                                      idx_include = idx_include_current)
    if(all(is.na(Spar_cor))){
      out <- list(cor = rep(NA, p), Final_idx = idx_include_old)
      return(out)
    }
    idx_include_old <- idx_include_current
    idx_include_current <- abs(Spar_cor) < t
    current_iter <- current_iter + 1
  }

  out <- list(cor = Spar_cor,
              Final_idx = idx_include_old)
  
  return(out)
}

#' Estimation of Cross-correlations between OTU abundances (or components of other
#' compositional variables) and a non-compositional variable.
#' @param OTU.abn An OTU table with observarions as rows and OTUs columns
#' (or a similar compositional dataset).
#' @param phenotype Continuous non-compositional phenotypic variable.
#' @param pseudo_count The pseudo-count added to all read counts to avoid division by zero.
#' @param var_min In the event that a variance estimate of an OTU is negative,
#' it is replaced by this number.
#' @param t A cutoff below which a correlation is considered weak. Only
#' OTUs weakly correlated with the phenotypic variable are included in the next
#' iteration. When t = NULL, it is chosen via bootstrap.
#' @param B_t The number of bootstraps to use when selecting t. Only used when t = NULL.
#' @param t_quant The quantile used under bootstrap for selecting t. Only used whern t = NULL.
#' @param Find_m Toggles whether or not the permutation threshold, m, should be 
#' calculated and returned. All correlations above m are likely to be real signals.
#' @param B_m The number of bootstraps to use when calculating m. Only used when Find_m = TRUE.
#' @param cores The number of cores used for selecting t and m.
#' @param iter The maximum number of iteration carried out.
#' @return A list containing the following: cor: a vector containing the
#' correlations between the OTU abundances and the phenotypic variable, Final_idx:
#' A logical vector indicating whether or not an OTU was included in the final iteration,
#' t: The value of t (either user specified or returned by bootstrap), m: the value of
#' m (NULL if Find_m = FALSE), Correlated: A logical vector indicating if an OTU's estimated
#' absolute correlation coefficient is above m.
#' @export
SparCEV <- function(OTU.abn, phenotype, pseudo_count = 1, var_min = 1e-5, iter = 20,
                    t = NULL, B_t = 100, t_quant = 0.8, cores = 1, Find_m = TRUE, B_m = 100){
  p <- ncol(OTU.abn)

  if(is.null(colnames(OTU.abn))){
    colnames(OTU.abn) <- paste0("OTU", 1:p)
  }

  if(is.null(t)){
    if(cores == 1){
      bt_list <- lapply(1:B_t,
                        function(x){
                          idx <- sample(1:nrow(OTU.abn))
                          SparCEV_base(OTU.abn = OTU.abn,
                                       phenotype = phenotype[idx],
                                       pseudo_count = pseudo_count,
                                       var_min = var_min,
                                       Find_m = F)$cor})
    }
    if(cores > 1){
      bt_list <- parallel::mclapply(1:B_t,
                                    function(x){
                                      idx <- sample(1:nrow(OTU.abn))
                                      SparCEV_base(OTU.abn = OTU.abn,
                                                   phenotype = phenotype[idx],
                                                   pseudo_count = pseudo_count,
                                                   var_min = var_min,
                                                   Find_m = F)$cor
                                    }, mc.cores = cores )
    }
    t <- mean(unlist(lapply(bt_list, function(x) quantile(abs(x), t_quant))))
  }

  Spar_run <- SparCEV_iter_inner(OTU.abn = OTU.abn,
                                 phenotype = phenotype,
                                 pseudo_count = pseudo_count,
                                 var_min = var_min,
                                 iter = iter,
                                 t = t)

  if(Find_m){
    if(cores == 1){
      bt_list <- lapply(1:B_m,
                        function(x){
                          idx <- sample(1:nrow(OTU.abn))
                          SparCEV_iter_inner(OTU.abn = OTU.abn,
                                             phenotype = phenotype[idx],
                                             pseudo_count = pseudo_count,
                                             var_min = var_min,
                                             iter = iter,
                                             t = t)$cor })
    }
    if(cores > 1){
      bt_list <- parallel::mclapply(1:B_m,
                                    function(x){
                                      idx <- sample(1:nrow(OTU.abn))
                                      SparCEV_iter_inner(OTU.abn = OTU.abn,
                                                         phenotype = phenotype[idx],
                                                         pseudo_count = pseudo_count,
                                                         var_min = var_min,
                                                         iter = iter,
                                                         t = t)$cor
                                    }, mc.cores = cores )
    }
    
    m <- mean(unlist(lapply(bt_list, function(x) max(abs(x)))))
    Correlated <- abs(Spar_run$cor) > m
  }

  if(!Find_m){
    m <- NULL
    Correlated <- NULL
  }

  out <- list(cor = Spar_run$cor,
              Final_idx = Spar_run$Final_idx,
              t = t,
              m = m,
              Correlated = Correlated)
  return(out)
}
