# Simulation study (under mixing)

require("cubature") # adaptive multidimensional integration (e.g., Cuhre)
require("doFuture") # registers future-based backends for foreach (%dopar%)
require("fs") # cross-platform filesystem ops (paths, files, dirs)
require("future.batchtools") # runs futures via batchtools (Slurm templates/resources)
require("ksm") # SPD kernel density tools (kdens_symmat, integrate_spd, bandwidth_optim)
require("parallel") # base R parallelism (PSOCK clusters, detectCores)

# Define the list of libraries to load on each cluster node

libraries_to_load <- c(
  "cubature",
  "doFuture",
  "fs",
  "future.batchtools",
  "parallel",
  "ksm"
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "combo",
  "cores_per_node",
  "criteria",
  "iw_logsumexp",
  "iw_lmgamma",
  "iw_spd_array",
  "iw_prepare",
  "make_kdens_inverseWishart",
  "kdens_inverseWishart",
  "bandwidth_optim_inverseWishart",
  "ISE",
  "IAE",
  "kernels",
  "libraries_to_load",
  "ISE",
  "Mmod",
  "K",
  "models",
  "nobs",
  "path",
  "replication",
  "resources_list",
  "RR",
  "setup_parallel_cluster",
  "simu_rWAR",
  "simu_fdens_WAR",
  "Smod",
  "timing",
  "vars_to_export"
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function(export_env = parent.frame()) {
  num_cores <- as.integer(cores_per_node)
  cl <- parallel::makeCluster(num_cores)
  
  # Export the list of libraries to the worker nodes
  clusterExport(cl, varlist = "libraries_to_load", envir = export_env)
  
  # Load necessary libraries on each cluster node
  invisible(clusterEvalQ(cl, {
    lapply(libraries_to_load, library, character.only = TRUE)
  }))
  
  # Export all necessary objects, functions, and parameters to the worker nodes
  clusterExport(cl, varlist = vars_to_export, envir = export_env)
  
  return(cl) # Return the cluster object
}

# Initialize missing placeholders without replacing existing libraries or functions

invisible(
  lapply(
    vars_to_export[
      !(vars_to_export %in% c("vars_to_export", "setup_parallel_cluster"))
    ],
    function(x) {
      if (!exists(x, envir = .GlobalEnv, inherits = FALSE)) {
        assign(x, NULL, envir = .GlobalEnv)
      }
    }
  )
)

########################################
## Inverse Wishart KDE and selectors  ##
########################################

# These helpers implement the manuscript's variable-interchanged kernel:
# f_IW(S) = mean_i K_{1/b + d + 1, b X_i}(S).
# ksm currently supports only Wishart, smlnorm and smnorm kernels.
# To reuse just these functions in another simulation script, evaluate this
# section without running the experiment below.

iw_logsumexp <- function(z) {
  zmax <- max(z)
  if (!is.finite(zmax)) return(zmax)
  zmax + log(sum(exp(z - zmax)))
}

iw_lmgamma <- function(a, d) {
  d * (d - 1) * log(pi) / 4 + sum(lgamma(a - (seq_len(d) - 1) / 2))
}

iw_spd_array <- function(x) {
  if (is.matrix(x)) x <- array(x, dim = c(nrow(x), ncol(x), 1L))
  if (length(dim(x)) != 3L || dim(x)[1] != dim(x)[2] || dim(x)[3] < 1L ||
      any(!is.finite(x))) stop("Expected a finite d x d x n array of SPD matrices.")
  for (i in seq_len(dim(x)[3])) {
    Xi <- x[, , i, drop = FALSE]
    dim(Xi) <- dim(x)[1:2]
    if (!isSymmetric(Xi, tol = 1e-10)) stop("Matrices must be symmetric.")
    chol(Xi) # also validates positive definiteness
  }
  x
}

iw_prepare <- function(x, overlap = FALSE) {
  x <- iw_spd_array(x)
  d <- dim(x)[1]
  n <- dim(x)[3]
  xvec <- matrix(x, nrow = d * d, ncol = n)
  inverses <- matrix(NA_real_, nrow = d * d, ncol = n)
  logdet <- numeric(n)
  for (i in seq_len(n)) {
    Xi <- matrix(x[, , i], nrow = d, ncol = d)
    Ri <- chol(Xi)
    logdet[i] <- 2 * sum(log(diag(Ri)))
    inverses[, i] <- as.vector(chol2inv(Ri))
  }
  out <- list(x = x, d = d, n = n, xvec = xvec,
              inverses = inverses, logdet = logdet)
  if (overlap) {
    logdet_sum <- matrix(NA_real_, nrow = n, ncol = n)
    for (i in seq_len(n)) {
      for (j in seq_len(i)) {
        Zij <- matrix(x[, , i] + x[, , j], nrow = d, ncol = d)
        logdet_sum[i, j] <- logdet_sum[j, i] <- 2 * sum(log(diag(chol(Zij))))
      }
    }
    out$logdet_sum <- logdet_sum
  }
  out
}

# Construct once per sample and bandwidth before numerical integration.
# The returned function accepts evaluation matrices x and an optional log flag.
make_kdens_inverseWishart <- function(xs, b) {
  if (length(b) != 1L || !is.finite(b) || b <= 0) stop("b must be positive.")
  sample <- iw_prepare(xs)
  d <- sample$d
  alpha <- 1 / (2 * b)
  a <- alpha + (d + 1) / 2
  log_weights <- -a * sample$logdet -
    d * a * log(2 * b) - iw_lmgamma(a, d)
  log_n <- log(sample$n)

  function(x, log = FALSE) {
    x <- iw_spd_array(x)
    if (dim(x)[1] != d) stop("x and xs must have the same matrix dimension.")
    neval <- dim(x)[3]
    xvec <- matrix(x, nrow = d * d, ncol = neval)
    logdet <- vapply(seq_len(neval), function(i) {
      2 * sum(base::log(diag(chol(matrix(x[, , i], nrow = d, ncol = d)))))
    }, numeric(1))
    # Rows are evaluation points; columns are observations in the KDE.
    log_kernel <- outer(alpha * logdet, log_weights, "+") -
      crossprod(xvec, sample$inverses) / (2 * b)
    log_density <- apply(log_kernel, 1L, iw_logsumexp) - log_n
    if (log) log_density else exp(log_density)
  }
}

# x: evaluation matrices; xs: observations. All densities use Lebesgue measure.
kdens_inverseWishart <- function(x, xs, b, log = FALSE) {
  make_kdens_inverseWishart(xs, b)(x, log = log)
}

bandwidth_optim_inverseWishart <- function(
    x, criterion = c("lscv", "lcv"), tol = 1e-4,
    h = 1L, bounds = c(1e-4, 10)
) {
  criterion <- match.arg(criterion)
  if (length(bounds) != 2L || any(!is.finite(bounds))) stop("Invalid bandwidth bounds.")
  bounds <- sort(bounds)
  if (bounds[1] <= 0 || bounds[2] <= bounds[1]) stop("Invalid bandwidth bounds.")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be positive.")
  if (length(h) != 1L || !is.finite(h) || h < 1 || h != as.integer(h)) {
    stop("h must be a positive integer.")
  }
  cache <- iw_prepare(x, overlap = criterion == "lscv")
  d <- cache$d
  n <- cache$n
  keep <- abs(outer(seq_len(n), seq_len(n), "-")) >= h
  nkeep <- rowSums(keep)
  if (any(nkeep == 0L)) stop("h leaves at least one validation point without observations.")
  trace_mat <- crossprod(cache$xvec, cache$inverses)
  logdet_pair <- outer(cache$logdet, cache$logdet, "+")
  cdim <- (d + 1) / 2
  rdim <- d * cdim

  objective <- function(b) {
    alpha <- 1 / (2 * b)
    a <- alpha + cdim
    log_kernel <- outer(alpha * cache$logdet, -a * cache$logdet, "+") -
      trace_mat / (2 * b) - d * a * log(2 * b) - iw_lmgamma(a, d)
    log_kernel[!keep] <- -Inf
    log_cv <- apply(log_kernel, 1L, iw_logsumexp) - log(nkeep)
    if (criterion == "lcv") return(-mean(log_cv))

    # This squared-density integral equals the overlap for the Wishart KDE.
    log_overlap <- iw_lmgamma(2 * alpha + cdim, d) -
      2 * iw_lmgamma(a, d) - rdim * log(2 * b) +
      alpha * logdet_pair - (2 * alpha + cdim) * cache$logdet_sum
    log_l2 <- iw_logsumexp(as.vector(log_overlap)) - 2 * log(n)
    # Match ksm's LSCV convention: pool all retained validation pairs.
    log_cross <- log(2) + iw_logsumexp(log_kernel[keep]) - log(sum(nkeep))

    # A strictly increasing signed-log transform preserves the LSCV minimizer
    # and avoids overflow in the difference of two positive sums.
    if (log_l2 == log_cross) return(0)
    sign_q <- if (log_l2 > log_cross) 1 else -1
    log_abs_q <- max(log_l2, log_cross) +
      log(-expm1(-abs(log_l2 - log_cross)))
    sign_q * (max(log_abs_q, 0) + log1p(exp(-abs(log_abs_q))))
  }

  # Match ksm::bandwidth_optim, including its endpoint-bound expansions.
  opt <- stats::optimize(objective, interval = bounds, tol = tol)$minimum
  if ((opt - bounds[1]) < 1e-3) {
    opt <- stats::optimize(objective, interval = c(1e-8, bounds[1] * 1.1),
                           tol = tol)$minimum
  }
  if ((bounds[2] - opt) < 1e-3) {
    opt <- stats::optimize(objective, interval = c(bounds[2] * 0.9, 1e4),
                           tol = tol)$minimum
  }
  opt
}

#' Simulation study serial dependence simulation from WAR(1)
#' @param n sample size
#' @param Mmod model for autoregressive coefficients matrix
#' @param Smod covariance matrix of VAR components
#' @param K degrees of freedom
#' @return a cube of dimension 2 by 2 by \code{n}.

simu_rWAR <- function(
    n,
    Mmod = c("M1", "M2", "M3"),
    Smod = c("S1", "S2", "S3"),
    K = 5L
) {
  Mmod <- match.arg(Mmod)
  Smod <- match.arg(Smod)
  n <- as.integer(n)
  K <- as.integer(K)
  stopifnot(n >= 1, K >= 2) # otherwise matrix not invertible / degenerate
  if (Mmod == "M1") {
    M <- matrix(
      c(0.9, 0, 1, 0),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else if (Mmod == "M2") {
    M <- matrix(
      c(0.3, -0.3, -0.3, 0.3),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else {
    M <- diag(rep(0.5, 2))
  }
  if (Smod == "S1") {
    S <- matrix(c(1, 0.9, 0.9, 1), nrow=2)
  } else if (Smod == "S2") {
    S <- matrix(c(1, 0.95, 0.95, 1), nrow=2)
  } else {
    S <- matrix(c(1, 0.99, 0.99, 1), nrow=2)
  }
  ksm::rWAR(
    n = n,
    M = M,
    S = S,
    K = K
  )
}

#' Simulation study asymptotic density for dependent
#' @param x array of positive definite matrices of size 2 by 2 by \code{n}.
#' @param Mmod model for autoregressive coefficients matrix
#' @param Smod covariance matrix of VAR components
#' @param K degrees of freedom
#' @return a vector of density functions

simu_fdens_WAR <- function(
    x,
    Mmod = c("M1", "M2", "M3"),
    Smod = c("S1", "S2", "S3"),
    K = 5L
) {
  Mmod <- match.arg(Mmod)
  Smod <- match.arg(Smod)
  K <- as.integer(K)
  stopifnot(K >= 2) # otherwise matrix not invertible / degenerate
  if (Mmod == "M1") {
    M <- matrix(
      c(0.9, 0, 1, 0),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else if (Mmod == "M2") {
    M <- matrix(
      c(0.3, -0.3, -0.3, 0.3),
      byrow = TRUE,
      nrow = 2,
      ncol = 2
    )
  } else {
    M <- diag(rep(0.5, 2))
  }
  if (Smod == "S1") {
    S <- matrix(c(1, 0.9, 0.9, 1), nrow=2)
  } else if (Smod == "S2") {
    S <- matrix(c(1, 0.95, 0.95, 1), nrow=2)
  } else {
    S <- matrix(c(1, 0.99, 0.99, 1), nrow=2)
  }
  
  Sigma_inf <- ksm::Riccati(M = M, S = S)$solution
  c(ksm::dWishart(x, df = K, S = Sigma_inf, log = FALSE))
}

#' Integrated squared error of kernel density estimator for symmetric matrices
#'
#' Given a sample \code{x} from a target density, and the optimal bandwidth, compute the integrated squared error via numerical integration.
#' @param x sample of symmetric matrix observations from which to build the kernel density kernel
#' @param bandwidth double for the bandwidth of the kernel
#' @param fdens density of the target function from which points are drawn.
#' @param kernel string, one of \code{Wishart}, \code{inverseWishart}, \code{smlnorm} (log-Gaussian) or \code{smnorm} (Gaussian).
#' @param tol double for tolerance of numerical integral
#' @param lb lower bound for integration range of eigenvalues
#' @param ub upper bound for integration range of eigenvalues
#' @param return double with value of the integrated squared error
#' @param ... additional parameters, currently ignored

ISE <- function(
    x,
    bandwidth,
    fdens,
    kernel = c("Wishart", "inverseWishart", "smlnorm", "smnorm"),
    tol = 1e-3,
    lb = 0,
    ub = Inf,
    neval = 1e6L,
    method = c("suave", "hcubature"),
    ...
) {
  method <- match.arg(
    method[1],
    choices = c("hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas")
  )
  kernel <- match.arg(kernel)
  dim <- ncol(x)
  stopifnot(
    dim %in% c(2L, 3L),
    length(lb) == 1L,
    length(ub) == 1L,
    lb >= 0,
    ub > lb,
    tol > 0,
    is.numeric(bandwidth) & length(bandwidth) == 1L,
    is.finite(bandwidth),
    bandwidth > 0
  )
  # Reuse the sample determinants and inverses throughout numerical integration.
  iw_density <- if (kernel == "inverseWishart") {
    make_kdens_inverseWishart(xs = x, b = bandwidth)
  } else {
    NULL
  }
  integrand <- function(vars, kernel) {
    # Construct SPD matrix S
    if (length(vars) == 3L) {
      S <- ksm::rotation_scaling(vars[1], vars[2:3])
      jacobian_value <- abs(vars[2] - vars[3]) / 4
    } else {
      S <- ksm::rotation_scaling(vars[1:3], vars[4:6])
      jacobian_value <- abs(vars[4] - vars[5]) *
        abs(vars[4] - vars[6]) *
        abs(vars[5] - vars[6]) *
        sin(vars[2]) /
        24
    }
    S <- array(S, dim = c(nrow(S), ncol(S), 1))
    
    # Compute the squared difference between the kernel and the target density.
    density_estimate <- if (kernel == "inverseWishart") {
      iw_density(S)
    } else {
      ksm::kdens_symmat(
        x = S,
        xs = x,
        b = bandwidth,
        kernel = kernel,
        log = FALSE
      )
    }
    diff_squared <- (density_estimate - fdens(S))^2
    
    # Return the product of diff_squared and the Jacobian factor
    return(diff_squared * jacobian_value)
  }
  if (dim == 2L) {
    # Define the limits of integration
    lower_limit <- c(0, rep(lb, length.out = dim))
    upper_limit <- c(2 * pi, rep(ub, length.out = dim))
  } else if (dim == 3L) {
    lower_limit <- c(rep(0, 3), rep(lb, length.out = dim))
    upper_limit <- c(2 * pi, pi, 2 * pi, rep(ub, length.out = dim))
  }
  result <- cubature::cubintegrate(
    f = integrand,
    lower = lower_limit,
    upper = upper_limit,
    method = method,
    relTol = tol,
    maxEval = neval,
    kernel = kernel
  )
  
  # print(result)
  # Return the result of the integral * 10^6
  return(1e6 * result$integral)
}

##################
## Set the path ##
##################

path <- getwd()
# setwd(path)

##############################
## Parallelization on nodes ##
##############################

cores_per_node <- 64 # number of cores for each node in the super-computer

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "240G",
  walltime = "12:00:00",
  nodes = 1
  # Omit 'partition' to let SLURM choose
)

##############################
## Hyper-parameters         ##
##############################

RR <- 1:1024 # replications
nobs <- c(100L, 200L, 300L)
Mmod <- c("M1", "M2", "M3")
Smod <- c("S1", "S2", "S3")
kernels <- c("Wishart", "inverseWishart", "smlnorm", "smnorm")
criteria <- c("lscv", "lcv")
K <- 5

###############
## Main code ##
###############

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# Disable the check for random number generation misuse in doFuture
options(doFuture.rng.onMisuse = "ignore")

# Register the doFuture parallel backend
doFuture::registerDoFuture()

# Tweak the batchtools_slurm with the custom template and resources
myslurm <- future::tweak(
  batchtools_slurm,
  template = "batchtools.slurm.dependent.tmpl",
  resources = resources_list
)

# Set the plan for future
future::plan(list(myslurm, multisession))

# Create empty data frames to store the results
raw_results <- data.frame(
  nobs = integer(),
  Mmod = integer(),
  Smod = integer(),
  kernel = character(),
  criterion = character(),
  ISE = numeric(),
  bandwidth = numeric(),
  timing = numeric()
)

# Capture the start time
start_time <- Sys.time()

# Split RR into node-sized blocks
RR_blocks <- split(RR, ceiling(seq_along(RR) / cores_per_node))

# For each over blocks (one Slurm job per block); inside each, parallelize reps across node cores
res <- foreach::foreach(
  block = RR_blocks,
  .combine = "rbind",
  .export = vars_to_export,
  .packages = libraries_to_load
) %dopar%
  {
    # spin up a per-node PSOCK cluster
    cl <- setup_parallel_cluster()
    
    block_results <- data.frame(
      nobs = integer(),
      Mmod = integer(),
      Smod = integer(),
      kernel = character(),
      criterion = character(),
      ISE = numeric(),
      bandwidth = numeric(),
      timing = numeric()
    )
    
    # each r in this block runs in parallel across the node's cores
    r_list <- parallel::parLapply(cl, X = block, fun = function(r_local) {
      set.seed(r_local)
      
      out_rows_all <- vector(
        "list",
        length(nobs) * length(Mmod) * length(Smod)
      )
      idx_out <- 0L
      
      for (i in seq_along(nobs)) {
        for (j in seq_along(Mmod)) {
          for (k in seq_along(Smod)) {
            local_rows <- list()
            
            xs <- simu_rWAR(
              n = nobs[i],
              Mmod = Mmod[j],
              Smod = Smod[k],
              K = K # degrees of freedom must be at least d
            )
            
            idx_local <- 0L
            for (l in seq_along(kernels)) {
              for (crit in criteria) {
                start_time <- Sys.time()
                band <- try(
                  if (kernels[l] == "inverseWishart") {
                    bandwidth_optim_inverseWishart(
                      x = xs,
                      criterion = crit,
                      h = ceiling(nobs[i]^0.25) # adjustment for serial dependence
                    )
                  } else {
                    ksm:::bandwidth_optim(
                      x = xs,
                      criterion = crit,
                      h = ceiling(nobs[i]^0.25),
                      kernel = kernels[l]
                    )
                  },
                  silent = FALSE
                )
                
                if (inherits(band, "try-error") || !is.numeric(band) ||
                    length(band) != 1L || !is.finite(band) || band <= 0) {
                  band <- NA_real_
                  try_ISE <- NA_real_
                } else {
                  try_ISE <- try(
                    ISE(
                      x = xs,
                      bandwidth = band,
                      fdens = function(x) {
                        simu_fdens_WAR(x = x, Mmod = Mmod[j], Smod = Smod[k], K = K)
                      },
                      kernel = kernels[l],
                      lb = 0,
                      ub = Inf,
                      tol = 1e-2,
                      method = "hcubature"
                    ),
                    silent = FALSE
                  )
                  if (inherits(try_ISE, "try-error")) {
                    try_ISE <- NA_real_
                  }
                }
                timing <- as.numeric(difftime(
                  Sys.time(),
                  start_time,
                  units = "secs"
                ))
                
                idx_local <- idx_local + 1L
                local_rows[[idx_local]] <- data.frame(
                  nobs = nobs[i],
                  Mmod = j,
                  Smod = k,
                  kernel = kernels[l],
                  criterion = crit,
                  ISE = try_ISE,
                  bandwidth = band,
                  timing = timing
                )
              }
            }
            
            idx_out <- idx_out + 1L
            out_rows_all[[idx_out]] <- do.call(rbind, local_rows)
          }
        }
      }
      
      do.call(rbind, out_rows_all)
    })
    
    block_results <- rbind(block_results, do.call(rbind, r_list))
    
    # tear down the per-node cluster
    try(parallel::stopCluster(cl), silent = TRUE)
    
    block_results
  }

# Combine results from all nodes
raw_results <- res

# Stop parallel execution
future::plan(sequential)

# Calculate the duration in minutes
elapsed_time_minutes <- as.numeric(difftime(
  Sys.time(),
  start_time,
  units = "mins"
))
print(paste("Elapsed time:", round(elapsed_time_minutes, 2), "minutes"))

# Save the raw results to a CSV file
raw_output_file <- file.path(path, "raw_ISE_dep.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_dep.csv")

#########################
## Process the results ##
#########################

# Build a summary table by (nobs, Mmod, Smod, kernel, criterion)

# Guard against accidental factor coercion
raw_results$kernel <- as.character(raw_results$kernel)
raw_results$criterion <- as.character(raw_results$criterion)

# Split by grouping keys
.grp <- interaction(
  raw_results$nobs,
  raw_results$Mmod,
  raw_results$Smod,
  raw_results$kernel,
  raw_results$criterion,
  drop = TRUE
)

summ_list <- lapply(split(raw_results, .grp), function(df) {
  data.frame(
    nobs = df$nobs[1],
    Mmod = df$Mmod[1],
    Smod = df$Smod[1],
    kernel = df$kernel[1],
    criterion = df$criterion[1],
    # ISE stats
    mean_ISE = mean(df$ISE, na.rm = TRUE),
    sd_ISE = sd(df$ISE, na.rm = TRUE),
    median_ISE = median(df$ISE, na.rm = TRUE),
    IQR_ISE = IQR(df$ISE, na.rm = TRUE),
    mean_time_min = mean(df$timing, na.rm = TRUE) / 60,
    stringsAsFactors = FALSE
  )
})

summary_results <- do.call(rbind, summ_list)
row.names(summary_results) <- NULL

# Save the summary to CSV next to the raw file
summary_output_file <- file.path(path, "summary_dep.csv")
write.csv(summary_results, summary_output_file, row.names = FALSE)

print("Summary results saved to summary_dep.csv")
