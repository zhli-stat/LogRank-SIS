library(ggplot2)
library(patchwork)
library(Icens)

fig_input_dir <- file.path("CaseStudy", "Output", "Salary")
data_dir      <- file.path("CaseStudy", "Data")
if (!file.exists(file.path(fig_input_dir, "salary_scores.csv")))
  stop("salary_scores.csv not found. Knit the salary Rmd first and run this ",
       "script from the repository root.")

scores_fig <- read.csv(file.path(fig_input_dir, "salary_scores.csv"),
                       stringsAsFactors = FALSE, fileEncoding = "UTF-8")
cfg_fig    <- read.csv(file.path(fig_input_dir, "salary_fig_config.csv"),
                       stringsAsFactors = FALSE, fileEncoding = "UTF-8")

dt_fig <- read.csv(file.path(data_dir, "WordMatrix_Filter.csv"))
s_fig  <- subset(dt_fig, select = -c(X1, l, r))
l_fig  <- as.numeric(dt_fig$l)
r_fig  <- as.numeric(dt_fig$r)
r_fig[is.infinite(r_fig)] <- Inf
word_names <- scores_fig$word
stopifnot(length(word_names) == ncol(s_fig))

GLOBAL_MAX_X <- cfg_fig$GLOBAL_MAX_X[1]
d            <- cfg_fig$d[1]
words_cfg <- lapply(seq_len(nrow(cfg_fig)), function(i)
  list(word = cfg_fig$word[i], label = cfg_fig$label[i], tag = cfg_fig$tag[i]))

emicm_to_df <- function(em, group_label) {
  pf   <- as.numeric(em$pf)
  rpts <- as.numeric(em$intmap[2L, ])

  if (length(pf) != length(rpts)) {
    min_len <- min(length(pf), length(rpts))
    pf      <- pf[1:min_len]
    rpts    <- rpts[1:min_len]
  }

  ok   <- is.finite(rpts) & !is.na(rpts) & !is.na(pf)
  pf   <- pf[ok]
  rpts <- rpts[ok]

  if (length(pf) == 0) {
    return(data.frame(x = c(0, GLOBAL_MAX_X), cdf = c(0, 0),
                      group = group_label, stringsAsFactors = FALSE))
  }

  tbl       <- tapply(pf, rpts, sum)
  uniq_pts  <- as.numeric(names(tbl))
  ord       <- order(uniq_pts)
  uniq_pts  <- uniq_pts[ord]
  cum_mass  <- cumsum(as.numeric(tbl)[ord])

  x_vals <- c(0, uniq_pts, max(GLOBAL_MAX_X, max(uniq_pts, na.rm=TRUE)))
  y_vals <- c(0, cum_mass, max(cum_mass))

  data.frame(x = x_vals, cdf = y_vals, group = group_label, stringsAsFactors = FALSE)
}

pull_scores_and_ranks <- function(idx) {
  list(add_sc = scores_fig$add_sc[idx], lr_sc = scores_fig$lr_sc[idx],
       add_rk = scores_fig$add_rk[idx], lr_rk = scores_fig$lr_rk[idx])
}

fmt_rank <- function(rk, d_thresh) {
  if (is.na(rk)) return("   N/A    ")
  if (rk <= d_thresh) return(sprintf(" Rank:%3d ", rk))
  return("Unselected")
}

make_panel <- function(cfg) {
  word_str    <- cfg$word
  word_label  <- cfg$label
  tag         <- cfg$tag

  idx <- which(word_names == word_str)[1]

  if(is.na(idx)) {
    return(ggplot() + theme_void() + labs(title = paste(tag, word_label, "[Word Not Found]")))
  }

  x_k   <- as.integer(s_fig[[idx]] == TRUE | s_fig[[idx]] == "TRUE" | s_fig[[idx]] == 1)
  n_pos <- sum(x_k == 1)
  n_neg <- sum(x_k == 0)

  if (n_pos < 2 || n_neg < 2) {
    return(ggplot() + theme_void() + labs(title = paste(tag, word_label, "[Insufficient Data]")))
  }

  em_pos <- suppressMessages(suppressWarnings(EMICM(cbind(l_fig[x_k == 1], r_fig[x_k == 1]))))
  em_neg <- suppressMessages(suppressWarnings(EMICM(cbind(l_fig[x_k == 0], r_fig[x_k == 0]))))

  lbl_pos <- sprintf("With keyword")
  lbl_neg <- sprintf("Without keyword (Baseline)")

  df_all <- rbind(emicm_to_df(em_pos, lbl_pos), emicm_to_df(em_neg, lbl_neg))

  d_val <- if (exists("d")) d else 60
  info  <- pull_scores_and_ranks(idx)

  ann_txt <- if (!is.na(info$add_sc) && !is.na(info$lr_sc)) {
    sprintf("ADD-SIS : %9.3f (%s)\nLog-rank: %9.3f (%s)", 
            info$add_sc, fmt_rank(info$add_rk, d_val), 
            info$lr_sc,  fmt_rank(info$lr_rk, d_val))
  } else {
    "Score/Rank unavailable"
  }

  pal  <- setNames(c("#B2182B", "#2166AC"), c(lbl_pos, lbl_neg)) 
  ltys <- setNames(c("solid",   "dashed"),  c(lbl_pos, lbl_neg))
  szs  <- setNames(c(1.4,       1.1),       c(lbl_pos, lbl_neg)) 

  p <- ggplot(df_all, aes(x = x, y = cdf, colour = group, linetype = group, linewidth = group)) +
    geom_step(direction = "hv", alpha = 0.9) +
    annotate("text", x = GLOBAL_MAX_X * 0.96, y = 0.05, label = ann_txt,
             hjust = 1, vjust = 0, size = 6.5, family = "mono", fontface = "bold", 
             colour = "#222222", lineheight = 1.1) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_linetype_manual(values = ltys, name = NULL) +
    scale_linewidth_manual(values = szs, name = NULL) +
    scale_x_continuous(
      name   = "Posted Salary (yuan / month)",
      breaks = seq(0, 40000, by = 10000),
      labels = c("0", "10k", "20k", "30k", "40k"), expand = c(0, 0)
    ) +
    scale_y_continuous(
      name   = "Cumulative Probability",
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"), expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(0, GLOBAL_MAX_X), ylim = c(0, 1.05)) +
    labs(title = paste(tag, word_label)) +
    theme_bw(base_size = 20) +
    theme(
      plot.title        = element_text(size = 22, face = "bold", margin = margin(b = 5)),
      axis.title        = element_text(size = 20, face = "bold"),
      axis.text         = element_text(size = 18, face = "bold", colour = "black"),
      legend.position   = "bottom", 
      legend.key.width  = unit(3, "cm"),
      legend.text       = element_text(size = 18, face = "bold"),
      panel.grid.minor  = element_blank(), 
      panel.grid.major.x= element_blank(),
      panel.grid.major.y= element_line(colour = "#E8E8E8", linewidth = 0.4),
      panel.border      = element_rect(colour = "black", linewidth = 1.2), 
      plot.margin       = margin(5, 5, 5, 5)
    )

  return(p)
}

panels <- lapply(words_cfg, make_panel)
names(panels) <- sapply(words_cfg, `[[`, "tag")

fig_out <- (panels[["(a)"]] | panels[["(b)"]]) /
  (panels[["(c)"]] | panels[["(d)"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box.margin = margin(t = 5, b = 5))

print(fig_out)

ggsave(
  filename = "salary_analysis_plot.pdf", 
  plot = fig_out, 
  width = 14,
  height = 11,
  device = cairo_pdf,
  dpi = 300
)