library(interval)

Y <- readRDS("ADNI_Y_Final.rds")
X_full <- readRDS("ADNI_X_Final.rds")
align_to_half_month <- function(x) { round(x * 2) / 2 }

Y$L_i <- align_to_half_month(Y$L_i)
Y$R_i <- ifelse(is.infinite(Y$R_i), Inf, align_to_half_month(Y$R_i))

conflict_idx <- which(Y$R_i <= Y$L_i)
if (length(conflict_idx) > 0) {
  Y$R_i[conflict_idx] <- Y$L_i[conflict_idx] + 0.5
}

snp_freq <- colMeans(X_full, na.rm = TRUE)
snp_var  <- snp_freq * (1 - snp_freq)

maf_pass <- snp_freq >= 0.01 & snp_freq <= 0.99
X_maf_temp  <- X_full[, maf_pass, drop = FALSE]
X_phase2 <- X_maf_temp
storage.mode(X_phase2) <- "integer"
X_phase2[is.na(X_phase2)] <- 0L
snp_names_p2 <- colnames(X_phase2)
l         <- Y$L_i
r         <- Y$R_i

npmle_on_grid <- function(l, r, supp) {
  fit <- initcomputeMLE(l, r)
  im  <- fit$intmap

  rpts <- if (nrow(im) == 2L) as.numeric(im[2, ]) else as.numeric(im[, 2])
  pf   <- as.numeric(fit$pf)
  p    <- numeric(length(supp))
  keep <- which(pf > 0 & is.finite(rpts))
  for (k in keep) {
    j <- which.min(abs(supp - rpts[k]))
    p[j] <- p[j] + pf[k]
  }
  tot <- sum(p)
  if (tot <= 0) {
    attr(p, "mass_beyond_support") <- 1
    return(p)
  }
  p <- p / tot
  attr(p, "mass_beyond_support") <- 1 - tot
  p
}

precompute_logrank_globals <- function(l, r) {
  n      <- length(l)
  supp   <- sort(unique(c(l[is.finite(l)], r[is.finite(r)])))
  p_hat  <- npmle_on_grid(l, r, supp)
  m_tail <- attr(p_hat, "mass_beyond_support")
  p_hat  <- as.numeric(p_hat)
  if (m_tail > 1e-8) {
    supp  <- c(supp, max(supp) + 1)
    p_hat <- c(p_hat * (1 - m_tail), m_tail)
  }
  active  <- p_hat > 1e-14
  supp_a  <- supp[active]
  p_hat_a <- p_hat[active] / sum(p_hat[active])
  Ja      <- length(supp_a)
  alpha   <- matrix(0L, n, Ja)
  for (i in seq_len(n)) {
    li <- l[i]; ri <- r[i]
    if      ( is.finite(li) &&  is.finite(ri)) alpha[i, supp_a >  li & supp_a <= ri] <- 1L
    else if ( is.finite(li) && !is.finite(ri)) alpha[i, supp_a >  li               ] <- 1L
    else if (!is.finite(li) &&  is.finite(ri)) alpha[i,                supp_a <= ri] <- 1L
  }
  denom <- as.vector(alpha %*% p_hat_a)
  denom[denom < 1e-14] <- 1e-14
  w     <- sweep(alpha * matrix(p_hat_a, n, Ja, byrow = TRUE), 1, denom, "/")
  d_j   <- colSums(w)
  n_j   <- rev(cumsum(rev(d_j)))
  ok    <- n_j > 1e-9
  c_i   <- as.vector(w %*% (ok - cumsum(ifelse(ok, d_j / n_j, 0))))
  list(supp = supp_a, p_hat = p_hat_a, alpha = alpha,
       w = w, d_j = d_j, n_j = n_j, c_i = c_i, J = Ja, n = n, tail_mass = m_tail)
}

screen_logrank_raw <- function(no, x_mat, g) {
  x_k <- x_mat[, no]
  if (sum(x_k) < 2 || sum(1 - x_k) < 2) return(0)
  U_n     <- sum(x_k * g$c_i) / g$n
  p_hat_k <- mean(x_k)

  abs(U_n) / sqrt(p_hat_k * (1 - p_hat_k) + 1e-8)
}

logrank_globals_adni <- precompute_logrank_globals(l, r)
p_total_p2 <- ncol(X_phase2)
score_lr_p2 <- sapply(seq_len(p_total_p2), screen_logrank_raw,
                      x_mat = X_phase2, g = logrank_globals_adni)
rank_lr_p2   <- rank(-score_lr_p2)

ord   <- order(score_lr_p2, decreasing = TRUE)[1:60]
top60 <- data.frame(
  Rank        = rank_lr_p2[ord],
  SNP         = snp_names_p2[ord],
  rsID        = sub("_[^_]+$", "", snp_names_p2[ord]),
  Score       = round(score_lr_p2[ord], 6),
  CarrierFreq = round(colMeans(X_phase2[, ord]), 4)
)
cat("Number of SNPs:", ncol(X_phase2), "\n")
print(top60, row.names = FALSE)
write.csv(top60[1:50, ], "ADNI_top50_genomewide_R.csv", row.names = FALSE)
writeLines(top60$rsID[1:50], "ADNI_top50_rsID.txt")
