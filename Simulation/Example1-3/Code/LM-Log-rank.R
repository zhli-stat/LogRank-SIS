library(Icens)
library(interval)
library(MASS)
library(dplyr)
library(energy)
library(parallel)
library(doParallel)

n_sim <- 200

output_dir <- file.path("Simulation", "Example1-3", "Intermediate", "LM_Log-rank")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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
  alpha <- matrix(0L, n, Ja)
  for (i in seq_len(n)) {
    li <- l[i]; ri <- r[i]
    if      ( is.finite(li) &&  is.finite(ri)) alpha[i, supp_a >  li & supp_a <= ri] <- 1L
    else if ( is.finite(li) && !is.finite(ri)) alpha[i, supp_a >  li               ] <- 1L
    else if (!is.finite(li) &&  is.finite(ri)) alpha[i,                supp_a <= ri] <- 1L
  }
  denom <- as.vector(alpha %*% p_hat_a)
  denom[denom < 1e-14] <- 1e-14
  w   <- sweep(alpha * matrix(p_hat_a, n, Ja, byrow = TRUE), 1, denom, "/")
  d_j <- colSums(w)
  n_j <- rev(cumsum(rev(d_j)))
  ok  <- n_j > 1e-9
  c_i <- as.vector(w %*% (ok - cumsum(ifelse(ok, d_j / n_j, 0))))
  list(supp = supp_a, p_hat = p_hat_a, alpha = alpha,
       w = w, d_j = d_j, n_j = n_j, c_i = c_i, J = Ja, n = n, tail_mass = m_tail)
}

screen_logrank_raw <- function(no, s, g) {
  x_k <- s[, no]
  if (sum(x_k) < 2 || sum(1 - x_k) < 2) return(0)
  U_n     <- sum(x_k * g$c_i) / g$n
  p_hat_k <- mean(x_k)

  abs(U_n) / sqrt(p_hat_k * (1 - p_hat_k) + 1e-8)
}

screen_npmle1 <- function(i,l,r,s,basis){
  every = data.frame(l=l,r=r,V3=s[,i])
  part1 = every[every$V3==0,c('l','r')]
  part2 = every[every$V3==1,c('l','r')]
  a <- initcomputeMLE(part1$l,part1$r)
  b <- initcomputeMLE(part2$l,part2$r)
  df1 <- data.frame(cbind(t(a$intmap),a$pf))
  df2 <- data.frame(cbind(t(b$intmap),b$pf))

  if (length(a$pf) == 1 && a$pf == 1){
    data2 <- transs(df2)
    tmp2 <- merge(data2, basis, all.y = TRUE)
    sta <- cumsum(tmp2$X3)
  } else if (length(b$pf) == 1 && b$pf == 1){
    data1 <- transs(df1)
    tmp1 <- merge(data1, basis, all.y = TRUE)
    sta <- cumsum(tmp1$X3)
  } else {
    data1 <- transs(df1)
    data2 <- transs(df2)
    tmp1 <- merge(data1, basis, all.y = TRUE)
    tmp2 <- merge(data2, basis, all.y = TRUE)
    tmp <- merge(tmp1, tmp2, by = 'X1')
    tmp[is.na(tmp$X3.x), c('X3.x')] <- 0
    tmp[is.na(tmp$X3.y), c('X3.y')] <- 0
    sta <- cumsum(tmp$X3.x) - cumsum(tmp$X3.y)
  }
  return(sum(abs(sta)))
}

screen_npmle2 <- function(i,l,r,s,basis){
  every = data.frame(l=l,r=r,V3=s[,i])
  part1 = every[every$V3==0,c('l','r')]
  part2 = every[every$V3==1,c('l','r')]
  a <- initcomputeMLE(part1$l,part1$r)
  b <- initcomputeMLE(part2$l,part2$r)
  df1 <- data.frame(cbind(t(a$intmap),a$pf))
  df2 <- data.frame(cbind(t(b$intmap),b$pf))

  if (length(a$pf) == 1 && a$pf == 1){
    data2 <- transs(df2)
    tmp2 <- merge(data2, basis, all.y = TRUE)
    sta <- cumsum(tmp2$X3)
  } else if (length(b$pf) == 1 && b$pf == 1){
    data1 <- transs(df1)
    tmp1 <- merge(data1, basis, all.y = TRUE)
    sta <- cumsum(tmp1$X3)
  } else {
    data1 <- transs(df1)
    data2 <- transs(df2)
    tmp1 <- merge(data1, basis, all.y = TRUE)
    tmp2 <- merge(data2, basis, all.y = TRUE)
    tmp <- merge(tmp1, tmp2, by = 'X1')
    tmp[is.na(tmp$X3.x), c('X3.x')] <- 0
    tmp[is.na(tmp$X3.y), c('X3.y')] <- 0
    sta <- cumsum(tmp$X3.x) - cumsum(tmp$X3.y)
  }
  return(sum(abs(sta)**2))
}

screen_npmle3 <- function(i,l,r,s,basis){
  every = data.frame(l=l,r=r,V3=s[,i])
  part1 = every[every$V3==0,c('l','r')]
  part2 = every[every$V3==1,c('l','r')]
  a <- initcomputeMLE(part1$l,part1$r)
  b <- initcomputeMLE(part2$l,part2$r)
  df1 <- data.frame(cbind(t(a$intmap),a$pf))
  df2 <- data.frame(cbind(t(b$intmap),b$pf))

  if (length(a$pf) == 1 && a$pf == 1){
    data2 <- transs(df2)
    tmp2 <- merge(data2, basis, all.y = TRUE)
    sta <- cumsum(tmp2$X3)
  } else if (length(b$pf) == 1 && b$pf == 1){
    data1 <- transs(df1)
    tmp1 <- merge(data1, basis, all.y = TRUE)
    sta <- cumsum(tmp1$X3)
  } else {
    data1 <- transs(df1)
    data2 <- transs(df2)
    tmp1 <- merge(data1, basis, all.y = TRUE)
    tmp2 <- merge(data2, basis, all.y = TRUE)
    tmp <- merge(tmp1, tmp2, by = 'X1')
    tmp[is.na(tmp$X3.x), c('X3.x')] <- 0
    tmp[is.na(tmp$X3.y), c('X3.y')] <- 0
    sta <- cumsum(tmp$X3.x) - cumsum(tmp$X3.y)
  }
  return(max(abs(sta)))
}

int_split <- function(i,df){
  tmp1 <- as.numeric(df[i,][1])
  tmp2 <- as.numeric(df[i,][2])
  if (is.infinite(tmp1)) {start<- tmp2-1}
  else {start <- as.numeric(df[i,][1])}
  if (is.infinite(tmp2)) {end <- tmp1+1}
  else {end <- as.numeric(df[i,][2])}
  difs <- as.numeric(end-start-1)
  return (data.frame(cbind(start+c(0:difs),start+1+c(0:difs),rep(df[i,][3]/(difs+1),difs+1))))
}

transs <- function(df){
  dff <- df$X2-df$X1
  subs <- df[dff==1,c('X1','X2','X3')]  
  no <- which(dff!=1)
  if (length(no)==0) {res <- subs}
  else if (length(no)==1) {res <- data.frame(rbind(subs,int_split(no[1],df=df)))
  res$X1 <- as.numeric(res$X1)
  res$X2 <- as.numeric(res$X2)
  res$X3 <- as.numeric(res$X3)}
  else if (length(no)>=2) {ress <- data.frame(rbind(subs,int_split(no[1],df=df)))
  res <- data.frame(rbind(ress,int_split(no[2],df=df)))
  res$X1 <- as.numeric(res$X1)
  res$X2 <- as.numeric(res$X2)
  res$X3 <- as.numeric(res$X3)
  }
  return (res)
}

screen_ks<- function(i,y_mid,s){
  every = data.frame(cbind(y_mid,s[,i]))
  part1 = every[every$V2==0,c('y_mid')]
  part2 = every[every$V2==1,c('y_mid')]
  return (as.numeric(ks.test(part1,part2)$statistic))
}

screen_dc <- function(i,y_mid,s){
  return (dcor(y_mid,s[,i]))
}

screen_t <- function(i,y_mid,s){
  return (abs(cor(y_mid,s[,i])))
}

Fk<-function(X0,x) {
  Fk=c()
  for (i in 1:length(x)){
    Fk[i]=sum(X0<=x[i])/length(X0)
  }
  return(Fk)
}

Fkr<-function(X0,Y,yr,x) {
  Fkr=c()
  ind_yr=(Y==yr)
  for (i in 1:length(x)){
    Fkr[i]=sum((X0<=x[i])*ind_yr)/sum(ind_yr)
  }
  return(Fkr)
}

MV<-function(Xk,Y) {
  Fk0 <- Fk(Xk,Xk)
  Yr <- unique(Y)
  MVr <- c()
  for (r in 1:length(Yr)) {
    MVr[r] <- (sum(Y==Yr[r])/length(Y))*mean((Fkr(Xk,Y,Yr[r],Xk)-Fk0)^2)
  }
  MV <- sum(MVr)
  return(MV)
}

screen_mv <- function(i,s,y_mid){
  return (abs(MV(Y=s[,i],Xk=y_mid)))
}

Simu_Multi_Norm<-function(x_len, sd = sd, pho =pho){
  V <- matrix(data = NA, nrow = x_len, ncol = x_len)
  for(i in 1:x_len){
    for(j in 1:x_len){
      V[i,j] <- pho^abs(i-j)
    }
  }
  V<-(sd^2) * V
  return(V)
}

sigma1 <- Simu_Multi_Norm(x_len = 1000, sd = 1, pho = 0.2)
sigma2 <- Simu_Multi_Norm(x_len = 1000, sd = 1, pho = 0.5)

num_cores <- detectCores() - 2
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, varlist = c(
  "screen_npmle1", "screen_npmle2", "screen_npmle3",
  "npmle_on_grid", "precompute_logrank_globals", "screen_logrank_raw",
  "screen_ks", "screen_dc", "screen_mv", "screen_t",
  "Fk", "Fkr", "MV", "transs", "int_split",
  "initcomputeMLE", "ictest",
  "mvrnorm", "rnorm", "rpois", "rchisq",
  "rank", "quantile", "floor", "is.infinite", "is.finite",
  "sigma1", "sigma2", "n_sim"
))

clusterEvalQ(cl, {
  library(Icens)
  library(interval)
  library(MASS)
  library(energy)
  library(dplyr)
})

run_simulation_parallel_lm <- function(sim_indices, pp, nn, rat, balanced = TRUE) {
  results <- foreach(kk = sim_indices, .combine = cbind, 
                     .packages = c("Icens", "interval", "MASS", "energy", "dplyr"),
                     .errorhandling = "remove") %dopar% {

                       set.seed(kk)

                       tryCatch({
                         x <- mvrnorm(nn, mu = rep(0,1000), pp)

                         if (balanced) {
                           s <- x > 0  
                         } else {
                           s <- matrix(FALSE, nrow = nn, ncol = 1000)
                           for(j in 1:1000) {
                             threshold <- quantile(x[,j], 0.3)
                             s[,j] <- x[,j] > threshold
                           }
                         }

                         prob <- 0.2
                         kk_func <- function(n){
                           return(2 * c((log(200)/sqrt(200)) + (abs(rnorm(sd=1, n=n)))))
                         }

                         beta = c(1 * kk_func(3), rep(0, 15), -1 * kk_func(2), rep(0, 980))
                         y <- s %*% beta + rnorm(sd=1, n=200)

                         l = floor(y - rpois(n=200, lambda=2))
                         r = floor(y + 2 * rchisq(n=200, df=3) + 1)

                         if (sum(r <= l) > 0) {
                           r[r <= l] <- r[r <= l] + 1
                         }

                         tot <- prob * nn
                         nos <- sample(1:nn, size = tot)
                         tt <- floor(rat * tot)

                         r[nos[1:tt]] <- Inf
                         if (tt < tot) {
                           l[nos[(tt+1):tot]] <- -Inf
                         }

                         transform <- function(i, l, r){
                           if (is.infinite(l[i])) { return(r[i]) }
                           else if (is.infinite(r[i])) { return(l[i]) }
                           else { return(0.5 * (l[i] + r[i])) }
                         }
                         transform2 <- function(i, l, r){
                           if (is.infinite(l[i])) { return(r[i]) }
                           else if (is.infinite(r[i])) { return(l[i]) }
                           else { return(2/3 * l[i] + 1/3 * r[i]) }
                         }
                         transform3 <- function(i, l, r){
                           if (is.infinite(l[i])) { return(r[i]) }
                           else if (is.infinite(r[i])) { return(l[i]) }
                           else { return(1/3 * l[i] + 2/3 * r[i]) }
                         }

                         y_mid <- sapply(1:nn, transform, l=l, r=r)
                         y_3 <- sapply(1:nn, transform2, l=l, r=r)
                         y_6 <- sapply(1:nn, transform3, l=l, r=r)

                         intervals <- initcomputeMLE(l, r)
                         max_tmp <- max(intervals$intmap[2,])
                         min_tmp <- min(intervals$intmap[1,])
                         if (is.finite(max_tmp)) { end <- max_tmp } 
                         else { end <- max(intervals$intmap[1,]) + 1 }
                         if (is.finite(min_tmp)) { start <- min_tmp } 
                         else { start <- min(intervals$intmap[2,]) - 1 }

                         X1 <- start:(end-1)
                         X2 <- (start+1):end
                         basis <- cbind(X1, X2)

                         time_results <- numeric(13)

                         a1 <- Sys.time()
                         s_npmle1 <- sapply(1:1000, screen_npmle1, l=l, r=r, s=s, basis=basis)
                         time_results[1] <- as.numeric(Sys.time() - a1, units = "secs")

                         a2 <- Sys.time()
                         s_npmle2 <- sapply(1:1000, screen_npmle2, l=l, r=r, s=s, basis=basis)
                         time_results[2] <- as.numeric(Sys.time() - a2, units = "secs")

                         a3 <- Sys.time()
                         s_npmle3 <- sapply(1:1000, screen_npmle3, l=l, r=r, s=s, basis=basis)
                         time_results[3] <- as.numeric(Sys.time() - a3, units = "secs")

                         a_lr <- Sys.time()
                         lr_globals <- precompute_logrank_globals(l, r) 
                         s_logrank <- sapply(1:1000, screen_logrank_raw, s=s, g=lr_globals) 
                         time_results[4] <- as.numeric(Sys.time() - a_lr, units = "secs")

                         a4 <- Sys.time()
                         s_ks <- sapply(1:1000, screen_ks, y_mid=y_mid, s=s)
                         time_results[5] <- as.numeric(Sys.time() - a4, units = "secs")

                         a5 <- Sys.time()
                         s_dc <- sapply(1:1000, screen_dc, y_mid=y_mid, s=s)
                         time_results[6] <- as.numeric(Sys.time() - a5, units = "secs")

                         a6 <- Sys.time()
                         s_mv <- sapply(1:1000, screen_mv, y_mid=y_mid, s=s)
                         time_results[7] <- as.numeric(Sys.time() - a6, units = "secs")

                         a7 <- Sys.time()
                         s_ks3 <- sapply(1:1000, screen_ks, y_mid=y_3, s=s)
                         time_results[8] <- as.numeric(Sys.time() - a7, units = "secs")

                         a8 <- Sys.time()
                         s_dc3 <- sapply(1:1000, screen_dc, y_mid=y_3, s=s)
                         time_results[9] <- as.numeric(Sys.time() - a8, units = "secs")

                         a9 <- Sys.time()
                         s_mv3 <- sapply(1:1000, screen_mv, y_mid=y_3, s=s)
                         time_results[10] <- as.numeric(Sys.time() - a9, units = "secs")

                         a10 <- Sys.time()
                         s_ks6 <- sapply(1:1000, screen_ks, y_mid=y_6, s=s)
                         time_results[11] <- as.numeric(Sys.time() - a10, units = "secs")

                         a11 <- Sys.time()
                         s_dc6 <- sapply(1:1000, screen_dc, y_mid=y_6, s=s)
                         time_results[12] <- as.numeric(Sys.time() - a11, units = "secs")

                         a12 <- Sys.time()
                         s_mv6 <- sapply(1:1000, screen_mv, y_mid=y_6, s=s)
                         time_results[13] <- as.numeric(Sys.time() - a12, units = "secs")

                         a <- (1001 - rank(s_npmle1))[c(1,2,3,19,20)]
                         b <- (1001 - rank(s_npmle2))[c(1,2,3,19,20)]
                         c <- (1001 - rank(s_npmle3))[c(1,2,3,19,20)]

                         lr_rank <- (1001 - rank(s_logrank))[c(1,2,3,19,20)]

                         d <- (1001 - rank(s_ks))[c(1,2,3,19,20)]
                         e <- (1001 - rank(s_dc))[c(1,2,3,19,20)]
                         f <- (1001 - rank(s_mv))[c(1,2,3,19,20)]

                         d3 <- (1001 - rank(s_ks3))[c(1,2,3,19,20)]
                         e3 <- (1001 - rank(s_dc3))[c(1,2,3,19,20)]
                         f3 <- (1001 - rank(s_mv3))[c(1,2,3,19,20)]

                         d6 <- (1001 - rank(s_ks6))[c(1,2,3,19,20)]
                         e6 <- (1001 - rank(s_dc6))[c(1,2,3,19,20)]
                         f6 <- (1001 - rank(s_mv6))[c(1,2,3,19,20)]

                         return(c(a, b, c, lr_rank, d, e, f, d3, e3, f3, d6, e6, f6, time_results))

                       }, error = function(e) {
                         warning(paste("Simulation", kk, "failed:", e$message))
                         return(rep(NA, 78))
                       })
                     }

  return(results)
}

cat("Running LM parallel simulations...\n")
cat("Replications:", n_sim, "\n")
cat("Cores used:", num_cores, "\n")

res1  <- run_simulation_parallel_lm(1:n_sim, pp=sigma1, nn=200, rat=0.3, balanced=TRUE)
res2  <- run_simulation_parallel_lm(1:n_sim, pp=sigma2, nn=200, rat=0.3, balanced=TRUE)
res4  <- run_simulation_parallel_lm(1:n_sim, pp=sigma1, nn=200, rat=0.7, balanced=TRUE)
res5  <- run_simulation_parallel_lm(1:n_sim, pp=sigma2, nn=200, rat=0.7, balanced=TRUE)
res7  <- run_simulation_parallel_lm(1:n_sim, pp=sigma1, nn=200, rat=0.3, balanced=FALSE)
res8  <- run_simulation_parallel_lm(1:n_sim, pp=sigma2, nn=200, rat=0.3, balanced=FALSE)
res10 <- run_simulation_parallel_lm(1:n_sim, pp=sigma1, nn=200, rat=0.7, balanced=FALSE)
res11 <- run_simulation_parallel_lm(1:n_sim, pp=sigma2, nn=200, rat=0.7, balanced=FALSE)

stopCluster(cl)

write.csv(res1,  file.path(output_dir, "Full_res1_Time_LR.csv"))
write.csv(res2,  file.path(output_dir, "Full_res2_Time_LR.csv"))
write.csv(res4,  file.path(output_dir, "Full_res4_Time_LR.csv"))
write.csv(res5,  file.path(output_dir, "Full_res5_Time_LR.csv"))
write.csv(res7,  file.path(output_dir, "Full_res7_Time_LR.csv"))
write.csv(res8,  file.path(output_dir, "Full_res8_Time_LR.csv"))
write.csv(res10, file.path(output_dir, "Full_res10_Time_LR.csv"))
write.csv(res11, file.path(output_dir, "Full_res11_Time_LR.csv"))

cat("All LM simulations complete. Results saved to:", output_dir, "\n")

calc_mms <- function(res_matrix, rows) {
  if (is.null(res_matrix) || nrow(res_matrix) < 78) return("NA")
  mms_rep <- apply(res_matrix[rows, ], 2, max, na.rm = TRUE)
  return(sprintf("%.1f(%.1f)", median(mms_rep, na.rm = TRUE), IQR(mms_rep, na.rm = TRUE) / 1.34))
}

calc_time <- function(res_list, time_row) {
  times <- c(as.numeric(res_list$c1_b[time_row, ]),
             as.numeric(res_list$c1_u[time_row, ]),
             as.numeric(res_list$c2_b[time_row, ]),
             as.numeric(res_list$c2_u[time_row, ]))
  return(sprintf("%.2f(%.2f)", mean(times, na.rm = TRUE), sd(times, na.rm = TRUE)))
}

generate_lm_table <- function(list_02, list_05) {
  methods_info <- list(
    list(name="Log-Rank",    rank_rows=16:20, time_row=69),
    list(name="ADD-SIS",     rank_rows=1:5,   time_row=66),
    list(name="DC-SIS (M1)", rank_rows=41:45, time_row=74),
    list(name="KF (M1)",     rank_rows=36:40, time_row=73),
    list(name="MV-SIS (M1)", rank_rows=46:50, time_row=75),
    list(name="DC-SIS (M2)", rank_rows=26:30, time_row=71),
    list(name="KF (M2)",     rank_rows=21:25, time_row=70),
    list(name="MV-SIS (M2)", rank_rows=31:35, time_row=72),
    list(name="DC-SIS (M3)", rank_rows=56:60, time_row=77),
    list(name="KF (M3)",     rank_rows=51:55, time_row=76),
    list(name="MV-SIS (M3)", rank_rows=61:65, time_row=78)
  )
  build_rho_df <- function(res_list, rho_label) {
    df <- data.frame()
    for (minfo in methods_info) {
      df <- rbind(df, data.frame(
        Rho              = rho_label,
        Method           = minfo$name,
        Case1_Balanced   = calc_mms(res_list$c1_b, minfo$rank_rows),
        Case1_Unbalanced = calc_mms(res_list$c1_u, minfo$rank_rows),
        Case2_Balanced   = calc_mms(res_list$c2_b, minfo$rank_rows),
        Case2_Unbalanced = calc_mms(res_list$c2_u, minfo$rank_rows),
        Time_Sec         = calc_time(res_list, minfo$time_row)
      ))
    }
    return(df)
  }
  final_table <- rbind(build_rho_df(list_02, "0.2"), build_rho_df(list_05, "0.5"))
  final_table$Rho[duplicated(final_table$Rho)] <- ""
  return(final_table)
}

list_02 <- list(c1_b = res1, c1_u = res7,  c2_b = res4, c2_u = res10)
list_05 <- list(c1_b = res2, c1_u = res8,  c2_b = res5, c2_u = res11)

lm_table <- generate_lm_table(list_02, list_05)
print(lm_table, row.names = FALSE)
