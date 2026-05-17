n  <- 200L
d1 <- floor(n / log(n))
d2 <- 2L * d1
d3 <- 3L * d1

method_rows <- list(
  "Log-Rank"    = 16:20,
  "ADD-SIS"     = 1:5,
  "KF (M1)"     = 36:40,
  "DC-SIS (M1)" = 41:45,
  "MV-SIS (M1)" = 46:50,
  "KF (M2)"     = 21:25,
  "DC-SIS (M2)" = 26:30,
  "MV-SIS (M2)" = 31:35,
  "KF (M3)"     = 51:55,
  "DC-SIS (M3)" = 56:60,
  "MV-SIS (M3)" = 61:65
)

read_res <- function(path) {
  if (!file.exists(path)) { warning("File not found: ", path); return(NULL) }
  as.matrix(read.csv(path, row.names = 1))
}

calc_prop <- function(res_mat, rows) {
  if (is.null(res_mat) || nrow(res_mat) < max(rows)) return(rep(NA_real_, 8))
  mat   <- res_mat[rows, , drop = FALSE]
  valid <- colSums(!is.na(mat)) == nrow(mat)
  if (!any(valid)) return(rep(NA_real_, 8))
  mat <- mat[, valid, drop = FALSE]
  Pk  <- rowMeans(mat <= d1)
  Pa1 <- mean(apply(mat <= d1, 2, all))
  Pa2 <- mean(apply(mat <= d2, 2, all))
  Pa3 <- mean(apply(mat <= d3, 2, all))
  round(c(Pk, Pa1, Pa2, Pa3), 3)
}

make_row <- function(bal_label, case_label, rho_label, method_name, prop_vals) {
  data.frame(
    Balanced = bal_label, Case = case_label, Rho = rho_label, Method = method_name,
    P1 = prop_vals[1], P2 = prop_vals[2], P3 = prop_vals[3],
    P19 = prop_vals[4], P20 = prop_vals[5],
    Pa1 = prop_vals[6], Pa2 = prop_vals[7], Pa3 = prop_vals[8]
  )
}

build_prop_table <- function(res_grid, rho_vals) {
  all_rows <- list()
  for (bal in c("Balanced", "Unbalanced")) {
    bal_key <- if (bal == "Balanced") "b" else "u"
    for (case_num in 1:2) {
      case_label <- paste("Case", case_num)
      for (rho in rho_vals) {
        key     <- paste(bal_key, paste0("c", case_num),
                         paste0("r", gsub("\\.", "", rho)), sep = "_")
        res_mat <- res_grid[[key]]
        for (m_idx in seq_along(method_rows)) {
          all_rows[[length(all_rows) + 1]] <- make_row(
            bal_label   = if (case_num == 1 && rho == rho_vals[1] && m_idx == 1) bal else "",
            case_label  = if (rho == rho_vals[1] && m_idx == 1) case_label else "",
            rho_label   = if (m_idx == 1) rho else "",
            method_name = names(method_rows)[m_idx],
            prop_vals   = calc_prop(res_mat, method_rows[[m_idx]])
          )
        }
      }
    }
  }
  do.call(rbind, all_rows)
}


lm_dir <- file.path("Simulation", "Example1-3", "Intermediate", "LM_Log-rank")
lm_res_grid <- list(
  b_c1_r02 = read_res(file.path(lm_dir, "Full_res1_Time_LR.csv")),
  b_c1_r05 = read_res(file.path(lm_dir, "Full_res2_Time_LR.csv")),
  b_c2_r02 = read_res(file.path(lm_dir, "Full_res4_Time_LR.csv")),
  b_c2_r05 = read_res(file.path(lm_dir, "Full_res5_Time_LR.csv")),
  u_c1_r02 = read_res(file.path(lm_dir, "Full_res7_Time_LR.csv")),
  u_c1_r05 = read_res(file.path(lm_dir, "Full_res8_Time_LR.csv")),
  u_c2_r02 = read_res(file.path(lm_dir, "Full_res10_Time_LR.csv")),
  u_c2_r05 = read_res(file.path(lm_dir, "Full_res11_Time_LR.csv"))
)
lm_prop_table <- build_prop_table(lm_res_grid, rho_vals = c("0.2", "0.5"))
cat("\n=== LM Example (d1=", d1, ", d2=", d2, ", d3=", d3, ") ===\n\n", sep = "")
print(lm_prop_table, row.names = FALSE)


glm_dir <- file.path("Simulation", "Example1-3", "Intermediate", "GLM_Log-rank")
glm_res_grid <- list(
  b_c1_r02 = read_res(file.path(glm_dir, "Full_res1_Time_LR.csv")),
  b_c1_r05 = read_res(file.path(glm_dir, "Full_res2_Time_LR.csv")),
  b_c1_r08 = read_res(file.path(glm_dir, "Full_res3_Time_LR.csv")),
  b_c2_r02 = read_res(file.path(glm_dir, "Full_res4_Time_LR.csv")),
  b_c2_r05 = read_res(file.path(glm_dir, "Full_res5_Time_LR.csv")),
  b_c2_r08 = read_res(file.path(glm_dir, "Full_res6_Time_LR.csv")),
  u_c1_r02 = read_res(file.path(glm_dir, "Full_res7_Time_LR.csv")),
  u_c1_r05 = read_res(file.path(glm_dir, "Full_res8_Time_LR.csv")),
  u_c1_r08 = read_res(file.path(glm_dir, "Full_res9_Time_LR.csv")),
  u_c2_r02 = read_res(file.path(glm_dir, "Full_res10_Time_LR.csv")),
  u_c2_r05 = read_res(file.path(glm_dir, "Full_res11_Time_LR.csv")),
  u_c2_r08 = read_res(file.path(glm_dir, "Full_res12_Time_LR.csv"))
)
glm_prop_table <- build_prop_table(glm_res_grid, rho_vals = c("0.2", "0.5", "0.8"))
cat("\n=== GLM Example ===\n\n")
print(glm_prop_table, row.names = FALSE)


cox_dir <- file.path("Simulation", "Example1-3", "Intermediate", "COX_Log-rank")
cox_res_grid <- list(
  b_c1_r02 = read_res(file.path(cox_dir, "Full_res1_Time_LR.csv")),
  b_c1_r05 = read_res(file.path(cox_dir, "Full_res2_Time_LR.csv")),
  b_c1_r08 = read_res(file.path(cox_dir, "Full_res3_Time_LR.csv")),
  b_c2_r02 = read_res(file.path(cox_dir, "Full_res4_Time_LR.csv")),
  b_c2_r05 = read_res(file.path(cox_dir, "Full_res5_Time_LR.csv")),
  b_c2_r08 = read_res(file.path(cox_dir, "Full_res6_Time_LR.csv")),
  u_c1_r02 = read_res(file.path(cox_dir, "Full_res7_Time_LR.csv")),
  u_c1_r05 = read_res(file.path(cox_dir, "Full_res8_Time_LR.csv")),
  u_c1_r08 = read_res(file.path(cox_dir, "Full_res9_Time_LR.csv")),
  u_c2_r02 = read_res(file.path(cox_dir, "Full_res10_Time_LR.csv")),
  u_c2_r05 = read_res(file.path(cox_dir, "Full_res11_Time_LR.csv")),
  u_c2_r08 = read_res(file.path(cox_dir, "Full_res12_Time_LR.csv"))
)
cox_prop_table <- build_prop_table(cox_res_grid, rho_vals = c("0.2", "0.5", "0.8"))
cat("\n=== COX Example ===\n\n")
print(cox_prop_table, row.names = FALSE)


out_dir <- file.path("Simulation", "Example1-3", "Tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(lm_prop_table,  file.path(out_dir, "Prop_Table_LM.csv"),  row.names = FALSE)
write.csv(glm_prop_table, file.path(out_dir, "Prop_Table_GLM.csv"), row.names = FALSE)
write.csv(cox_prop_table, file.path(out_dir, "Prop_Table_COX.csv"), row.names = FALSE)
cat("\nSaved to:", out_dir, "\n")
