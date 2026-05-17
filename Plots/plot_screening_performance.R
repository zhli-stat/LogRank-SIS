library(ggplot2)
library(dplyr)

# rho = 0.2, Case 1 (rat = 0.3); single-point methods use M1 imputation.
# Row indices in the 78-row result matrix for each method's active-predictor ranks.
method_rows <- list(
  "Log-Rank (Proposed)" = 16:20,
  "ADD-SIS"             = 1:5,
  "DC-SIS"              = 41:45,
  "KF"                  = 36:40,
  "MV-SIS"              = 46:50
)

read_res <- function(path) {
  if (!file.exists(path)) { warning("File not found: ", path); return(NULL) }
  as.matrix(read.csv(path, row.names = 1))
}

calc_mms_stats <- function(res_mat, rows) {
  if (is.null(res_mat) || nrow(res_mat) < max(rows)) return(c(NA_real_, NA_real_))
  mms <- apply(res_mat[rows, , drop = FALSE], 2, max, na.rm = TRUE)
  mms <- mms[is.finite(mms)]
  if (length(mms) == 0) return(c(NA_real_, NA_real_))
  c(median(mms), IQR(mms) / 1.34)
}

lm_dir  <- file.path("Simulation", "Example1-3", "Intermediate", "LM_Log-rank")
glm_dir <- file.path("Simulation", "Example1-3", "Intermediate", "GLM_Log-rank")
cox_dir <- file.path("Simulation", "Example1-3", "Intermediate", "COX_Log-rank")

scenarios <- list(
  list(model = "Linear",  setting = "Balanced",
       res = read_res(file.path(lm_dir,  "Full_res1_Time_LR.csv"))),
  list(model = "Linear",  setting = "Unbalanced",
       res = read_res(file.path(lm_dir,  "Full_res7_Time_LR.csv"))),
  list(model = "Poisson", setting = "Balanced",
       res = read_res(file.path(glm_dir, "Full_res1_Time_LR.csv"))),
  list(model = "Poisson", setting = "Unbalanced",
       res = read_res(file.path(glm_dir, "Full_res7_Time_LR.csv"))),
  list(model = "Cox",     setting = "Balanced",
       res = read_res(file.path(cox_dir, "Full_res1_Time_LR.csv"))),
  list(model = "Cox",     setting = "Unbalanced",
       res = read_res(file.path(cox_dir, "Full_res7_Time_LR.csv")))
)

rows_list <- list()
for (sc in scenarios) {
  for (m_name in names(method_rows)) {
    stats <- calc_mms_stats(sc$res, method_rows[[m_name]])
    rows_list[[length(rows_list) + 1]] <- data.frame(
      Model      = sc$model,
      Setting    = sc$setting,
      Method     = m_name,
      Median_MMS = stats[1],
      RSD        = stats[2]
    )
  }
}
plot_data <- do.call(rbind, rows_list)

y_axis_order <- rev(c("Log-Rank (Proposed)", "ADD-SIS", "DC-SIS", "MV-SIS", "KF"))
plot_data$Method  <- factor(plot_data$Method,  levels = y_axis_order)
plot_data$Model   <- factor(plot_data$Model,   levels = c("Linear", "Poisson", "Cox"))
plot_data$Setting <- factor(plot_data$Setting, levels = c("Balanced", "Unbalanced"))

legend_display_order <- c("Log-Rank (Proposed)", "ADD-SIS", "DC-SIS", "MV-SIS", "KF")

p <- ggplot(plot_data, aes(x = Median_MMS, y = Method)) +

  geom_segment(aes(xend = Median_MMS + RSD, yend = Method, color = Method),
               linewidth = 1.2, alpha = 0.75) +

  geom_point(aes(color = Method, shape = Method), size = 4.5, fill = "white", stroke = 1.5) +

  facet_grid(Setting ~ Model, scales = "free_x") +

  scale_color_manual(
    values = c(
      "Log-Rank (Proposed)" = "#D50000",
      "ADD-SIS"             = "#1565C0",
      "DC-SIS"              = "#2E7D32",
      "MV-SIS"              = "#8E24AA",
      "KF"                  = "#F57C00"
    ),
    breaks = legend_display_order
  ) +

  scale_shape_manual(
    values = c(19, 15, 17, 18, 1),
    breaks = legend_display_order
  ) +

  scale_y_discrete(
    expand = expansion(mult = 0.1),
    labels = c(
      "Log-Rank (Proposed)" = "Log-Rank\n(Proposed)",
      "ADD-SIS" = "ADD-SIS",
      "DC-SIS"  = "DC-SIS",
      "MV-SIS"  = "MV-SIS",
      "KF"      = "KF"
    )
  ) +

  theme_bw(base_size = 18, base_family = "serif") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
    axis.text.y  = element_text(size = 16, color = "black", face = "bold", lineheight = 0.9),
    axis.text.x  = element_text(size = 15, color = "black", face = "bold"),

    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.8),
    strip.text.x = element_text(face = "bold", size = 18, margin = margin(t = 10, b = 10)),
    strip.text.y = element_text(face = "bold", size = 18, angle = 270,
                                margin = margin(t = 12, r = 12, b = 12, l = 12)),

    panel.spacing    = unit(0.2, "lines"),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5, linetype = "dashed"),

    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 16, face = "bold"),
    legend.margin     = margin(t = 5, b = 0),
    legend.box.margin = margin(t = -2, b = 0),
    legend.key.size   = unit(1.2, "cm")
  ) +

  labs(
    x = "Median Minimum Model Size (MMS) with Robust Standard Deviation (RSD)",
    y = NULL
  )

print(p)

ggsave(filename = "Screening_Performance_2x3_LargeFont.pdf",
       plot = p, width = 13, height = 8.5, device = "pdf")
