# Fast SparCC-derived approximation of absolute abundance variances
Absolute_abn_var <- function(OTU.abn.log, var_min = 1e-5){
  p <- ncol(OTU.abn.log)
  x_var <- apply(OTU.abn.log, 2, var)
  log_sum <- rowSums(OTU.abn.log)
  t <- p*x_var + sum(x_var) - 2*apply(OTU.abn.log, 2, function(x) cov(x, log_sum))
  M <- matrix(1,p,p) + diag(p-2, p, p)
  alpha <- solve(M, t)
  alpha[alpha < var_min] <- var_min
  return(alpha)
}

SparXCC_base_internal <- function(OTU.abn, gene.expr, pseudo_count = 1, var_min = 1e-5,
                                  OTU_include = rep(T, ncol(OTU.abn)),
                                  Gene_include = rep(T, ncol(gene.expr))){
  D1 <- t(apply(OTU.abn+pseudo_count, 1, function(x) x/sum(x)))
  D2 <- t(apply(gene.expr+pseudo_count, 1, function(x) x/sum(x)))

  log.OTU <- log(D1); p <- ncol(log.OTU)
  log.genes <- log(D2); q <- ncol(log.genes)

  OTU_SLR <- t( apply(log.OTU, 1, function(x) x - mean(x[OTU_include]) ) )
  Gene_SLR <- t( apply(log.genes, 1, function(x) x - mean(x[Gene_include]) ) )

  R <- sum(OTU_include)
  S <- sum(Gene_include)

  t_ik <- R*S*cov(OTU_SLR, Gene_SLR)

  # Basis variances
  alpha <- Absolute_abn_var(log.OTU)
  beta <- Absolute_abn_var(log.genes)

  alpha_mat <- matrix(alpha, p, q)
  beta_mat <- matrix(beta, p, q, byrow = T)

  # constructing matrices that contain |R\{i}| and |S\{k}|
  R_mat <- matrix(R - OTU_include, p, q)
  S_mat <- matrix(S - Gene_include, p, q, byrow = T)

  # Computing correlation approximations
  rho <- t_ik/( R_mat*S_mat*sqrt(alpha_mat)*sqrt(beta_mat) )
  rho[abs(rho)>1] <- sign(rho[abs(rho)>1])
  rownames(rho) <- colnames(OTU.abn)
  colnames(rho) <- colnames(gene.expr)

  return(rho)
}

#' Estimation of Cross-correlations between OTU abundances (or components of other
#' compositional variables) and gene expressions (or components of other
#' compositional variables). This method is appropriate
#' when the average covariance is close to zero.
#' @param OTU.abn An OTU table with observarions as rows and OTUs columns
#' (or a similar compositional dataset).
#' @param gene.expr A table of gene expressions with observarions as rows and OTUs columns
#' (or a similar compositional dataset).
#' @param pseudo_count The pseudo-count added to all read counts to avoid division by zero.
#' @param var_min In the event that a variance estimate of an OTU is negative,
#' it is replaced by this number.
#' @param Find_m Toggles whether or not the permutation threshold, m, should be 
#' calculated and returned. All correlations above m are likely to be real signals.
#' @param B_m The number of bootstraps to use when calculating m. Only used when Find_m = TRUE.
#' @param cores The number of cores used for finding m.
#' @return A list that contains cor, cross-correlation between the OTU abundances and gene expressions, m, the dynamic permutation threshold, and cor_above_m, a logical matrix indicating whether the absoulte cross-correlation between the corresponding OTU and gene is above m.
#' @export
SparXCC_base <- function(OTU.abn, gene.expr, pseudo_count = 1, var_min = 1e-5,
                         Find_m = TRUE, B_m = 100, cores = 1){
  if(nrow(OTU.abn) != nrow(gene.expr)){
    stop("OTU.abn and gene.expr must contain the same number of replicates.")
  }
  
  p <- ncol(OTU.abn)
  q <- ncol(gene.expr)
  n <- nrow(OTU.abn)
  
  rho <-   SparXCC_base_internal(OTU.abn = OTU.abn,
                                 gene.expr = gene.expr,
                                 pseudo_count = pseudo_count,
                                 var_min = var_min,
                                 OTU_include = rep(T, p),
                                 Gene_include = rep(T, q))
  
  if(Find_m){
    if(cores == 1){
      bt_list <- lapply(1:B_m,
                        function(x){
                          SparXCC_base_internal(OTU.abn = OTU.abn[sample(1:n),],
                                                gene.expr = gene.expr[sample(1:n),],
                                                pseudo_count = pseudo_count,
                                                var_min = var_min,
                                                OTU_include = rep(T, p),
                                                Gene_include = rep(T, q))})
    }
    if(cores > 1){
      bt_list <- parallel::mclapply(1:B_m,
                                    function(x){
                                      SparXCC_base_internal(OTU.abn = OTU.abn[sample(1:n),],
                                                            gene.expr = gene.expr[sample(1:n),],
                                                            pseudo_count = pseudo_count,
                                                            var_min = var_min,
                                                            OTU_include = rep(T, p),
                                                            Gene_include = rep(T, q))
                                    }, mc.cores = cores )
    }
    m <- mean(unlist(lapply(bt_list, function(x) max(abs(x)))))
    Correlated <- abs(rho) > m
  }
  
  if(!Find_m){
    m <- NULL
    Correlated <- NULL
  }
  
  out <- list(cor = rho, 
              m = m, 
              Cor_above_m = Correlated)
  
  return(out)
}

SparXCC_iter_internal <- function(OTU.abn, gene.expr, pseudo_count = 1, var_min = 1e-5,
                                  iter = 10, t1, t2){
  if(nrow(OTU.abn) != nrow(gene.expr)){
    stop("OTU.abn and gene.expr must contain the same number of replicates.")
  }
  
  p <- ncol(OTU.abn)
  q <- ncol(gene.expr)
  n <- nrow(OTU.abn)
  
  if(is.null(colnames(OTU.abn))){
    colnames(OTU.abn) <- paste0("OTU", 1:ncol(OTU.abn))
  }
  if(is.null(colnames(gene.expr))){
    colnames(gene.expr) <- paste0("Gene", 1:ncol(gene.expr))
  }
  
  OTU_include_current <- rep(T, ncol(OTU.abn))
  OTU_include_old <- rep(F, ncol(OTU.abn))
  Gene_include_current <- rep(T, ncol(gene.expr))
  Gene_include_old <- rep(F, ncol(gene.expr))
  current_iter <- 1
  while(any(OTU_include_current != OTU_include_old) &
        any(Gene_include_current != Gene_include_old) &
        current_iter <= iter &
        sum(OTU_include_current) > 2 &
        sum(Gene_include_current) > 2){
    Spar_cor <- SparXCC_base_internal(OTU.abn = OTU.abn,
                                      gene.expr = gene.expr,
                                      pseudo_count = pseudo_count,
                                      var_min = var_min,
                                      OTU_include = OTU_include_current,
                                      Gene_include = Gene_include_current)
    if(all(is.na(Spar_cor))){
      return( matrix(NA, ncol(OTU.abn), ncol(gene.expr)) )
    }
    OTU_include_old <- OTU_include_current
    Gene_include_old <- Gene_include_current
    OTU_include_current <- rowMeans(abs(Spar_cor)) < t1
    Gene_include_current <- colMeans(abs(Spar_cor)) < t2
    current_iter <- current_iter + 1
  }
  
  if(current_iter == 2){
    warning("No iteration took place, since all Variables of both datasets had low average absolute correlation. \nConsider choosing a smaller t1 or t2.")
  }
  
  if(sum(OTU_include_current) < 2 | sum(Gene_include_current) < 2){
    warning("Iterative procedure was terminated prematurely due to insufficient low-correlation pairs. \nConsider changing t1 or t2 to a larger value.")
  }
  
  out <- list(cor = Spar_cor,
              Final_idx_OTU = OTU_include_old,
              Final_idx_Gene = Gene_include_old)
  return(out)
}

#' Estimation of Cross-correlations between OTU abundances (or components of other
#' compositional variables) and gene expressions (or components of other
#' compositional variables).
#' @param OTU.abn An OTU table with observarions as rows and OTUs columns
#' (or a similar compositional dataset).
#' @param gene.expr A table of gene expressions with observarions as rows and OTUs columns
#' (or a similar compositional dataset).
#' @param pseudo_count The pseudo-count added to all read counts to avoid division by zero.
#' @param var_min In the event that a variance estimate of an OTU is negative,
#' it is replaced by this number.
#' @param t1 A cutoff to be used in the iterative procedure. Only
#' OTUs whose average absolute correlation with the gene expressions are included in the next
#' iteration. When t1 = NULL, it is chosen via bootstrap.
#' @param t2 A cutoff to be used in the iterative procedure. Only
#' genes whose average absolute correlation with the OTU abundances are included in the next
#' iteration. When t2 = NULL, it is chosen via bootstrap.
#' @param B_t The number of bootstraps to use when selecting t1 and t2. Only used when t1 = NULL and/or t2 = NULL.
#' @param t1_quant The quantile used under bootstrap for selecting t1. Only used whern t1 = NULL.
#' @param t2_quant The quantile used under bootstrap for selecting t2. Only used whern t2 = NULL.
#' @param Find_m Toggles whether or not the permutation threshold, m, should be 
#' calculated and returned. All correlations above m are likely to be real signals.
#' @param B_m The number of bootstraps to use when calculating m. Only used when Find_m = TRUE.
#' @param cores The number of cores used for selecting t and m.
#' @param iter The maximum number of iteration carried out.
#' @return A list containing the following: cor: a matrix containing the
#' correlations between the OTU abundances and gene expressions, Final_idx_OTU:
#' A logical vector indicating whether or not an OTU was included in the final iteration,
#' Final_idx_gene: A logical vector indicating whether or not a gene was included in the final iteration,
#' t1: The value of t1 (either user specified or returned by bootstrap),
#' t2: The value of t2 (either user specified or returned by bootstrap),
#' m: the value of m (NULL if Find_m = FALSE), Correlated: A logical matrix 
#' indicating if the absolute correlation coefficient estimate between and OTU
#'  and a gene is above m.
#' @export
SparXCC <- function(OTU.abn, gene.expr, pseudo_count = 1, var_min = 1e-5,
                    iter = 10, t1 = NULL, t2 = NULL, B_t = 100, t1_quant = 0.8,
                    t2_quant = 0.8, Find_m = TRUE, B_m = 100, cores = 1){
  if(nrow(OTU.abn) != nrow(gene.expr)){
    stop("OTU.abn and gene.expr must contain the same number of replicates.")
  }

  p <- ncol(OTU.abn)
  q <- ncol(gene.expr)
  n <- nrow(OTU.abn)

  if(is.null(colnames(OTU.abn))){
    colnames(OTU.abn) <- paste0("OTU", 1:ncol(OTU.abn))
  }
  if(is.null(colnames(gene.expr))){
    colnames(gene.expr) <- paste0("Gene", 1:ncol(gene.expr))
  }

  if( is.null(t1) | is.null(t2) ){
    if(cores == 1){
      S <- lapply(1:B_t,
                  function(x){
                    S <- SparXCC_base(OTU.abn = OTU.abn[sample(1:n),],
                                      gene.expr = gene.expr[sample(1:n),],
                                      pseudo_count = pseudo_count,
                                      var_min = var_min,
                                      Find_m = F)$cor
                    list(OTU = rowMeans(abs(S)), Gene = colMeans(abs(S)))
                    })
    }
    if(cores > 1){
      S <- parallel::mclapply(1:B_t,
                              function(x){
                                S <- SparXCC_base(OTU.abn = OTU.abn[sample(1:n),],
                                                  gene.expr = gene.expr[sample(1:n),],
                                                  pseudo_count = pseudo_count,
                                                  var_min = var_min,
                                                  Find_m = F)$cor
                                list(OTU = rowMeans(abs(S)), Gene = colMeans(abs(S)))
                              }, mc.cores = cores)
    }
    if(is.null(t1)){
      t1 <- mean(unlist(lapply(S, function(x) quantile(x$OTU, t1_quant) )))
    }
    if(is.null(t2)){
      t2 <- mean(unlist(lapply(S, function(x) quantile(x$Gene, t2_quant) )))
    }
  }
  
  Spar_run <- SparXCC_iter_internal(OTU.abn = OTU.abn,
                                    gene.expr = gene.expr,
                                    pseudo_count = pseudo_count,
                                    var_min = var_min,
                                    iter = iter,
                                    t1 = t1,
                                    t2 = t2)

  if(Find_m){
    if(cores == 1){
      bt_list <- lapply(1:B_m,
                        function(x){
                          SparXCC_iter_internal(OTU.abn = OTU.abn[sample(1:n),],
                                                gene.expr = gene.expr[sample(1:n),],
                                                pseudo_count = pseudo_count,
                                                var_min = var_min,
                                                iter = iter,
                                                t1 = t1,
                                                t2 = t2)$cor
                                    })
    }
    if(cores > 1){
      bt_list <- parallel::mclapply(1:B_m,
                                    function(x){
                                      SparXCC_iter_internal(OTU.abn = OTU.abn[sample(1:n),],
                                                            gene.expr = gene.expr[sample(1:n),],
                                                            pseudo_count = pseudo_count,
                                                            var_min = var_min,
                                                            iter = iter,
                                                            t1 = t1,
                                                            t2 = t2)$cor
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
              Final_idx_OTU = Spar_run$Final_idx_OTU,
              Final_idx_Gene = Spar_run$Final_idx_Gene,
              t1 = t1,
              t2 = t2,
              m = m,
              Correlated = Correlated)
  return(out)
}
