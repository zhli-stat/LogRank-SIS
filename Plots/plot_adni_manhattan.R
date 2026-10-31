library(ggplot2)

# PATH: run from the repository root. Adjust bim_file if your PLINK .bim/.map lives elsewhere.
bim_file <- file.path("CaseStudy", "Code", "ADNI_1_GWAS_Plink", "ADNI_cluster_01_forward_757LONI.bim")
out_png  <- file.path("Plots", "ADNI_Manhattan_LogrankSIS.png")

if (!exists("score_lr_p2") || !exists("snp_names_p2")) {
  old_wd <- getwd()
  setwd(file.path("CaseStudy", "Code"))
  source("get_top50_genomewide.R")
  setwd(old_wd)
}

map <- read.table(bim_file, header = FALSE, colClasses = "character")
if (ncol(map) >= 6) {
  map <- data.frame(chr = map[[1]], snp = map[[2]], bp = as.numeric(map[[4]]))
} else if (ncol(map) == 4) {
  map <- data.frame(chr = map[[1]], snp = map[[2]], bp = as.numeric(map[[4]]))
} else stop("Unrecognised map file: expected 6 (.bim) or 4 (.map) columns.")

rs  <- sub("_[^_]+$", "", snp_names_p2)
idx <- match(rs, map$snp)
cat("SNPs matched to the map file:", sum(!is.na(idx)), "of", length(rs), "\n")
if (mean(!is.na(idx)) < 0.99) stop("Fewer than 99% of SNPs matched; check bim_file.")

chr <- toupper(map$chr[idx])
chr[chr %in% c("23", "25", "XY")] <- "X"
keep <- !is.na(idx) & chr %in% c(as.character(1:22), "X")
cat("Plotted:", sum(keep), "| not plotted (unmatched, Y, MT):", sum(!keep), "\n")

d <- data.frame(snp = rs, chr = factor(chr, levels = c(as.character(1:22), "X")),
                bp = map$bp[idx], score = score_lr_p2,
                rank = rank(-score_lr_p2))[keep, ]
d <- d[order(d$chr, d$bp), ]

len    <- tapply(d$bp, d$chr, max)
gap    <- 0.008 * sum(len)
offset <- c(0, cumsum(as.numeric(len) + gap))[seq_along(len)]
names(offset) <- names(len)
d$x    <- d$bp + offset[as.character(d$chr)]
d$band <- factor(as.integer(d$chr) %% 2)
axis_x <- offset + len / 2

thr     <- sort(score_lr_p2, decreasing = TRUE)[50]
top     <- d[d$rank <= 50, ]
bg      <- d[d$rank > 50, ]
genes   <- merge(data.frame(snp = c("rs10821702", "rs11951879"),
                            gene = c("ANK3", "PDE4D")), d)
x_right <- max(d$x)

col_top  <- "#2a78d6"
greys    <- c("#d3d2cd", "#a3a29c")
ink      <- "#52514e"
x_lab    <- x_right + 0.012 * x_right

p <- ggplot() +
  geom_point(data = bg, aes(x, score, colour = band), size = 0.22, shape = 16) +
  scale_colour_manual(values = greys, guide = "none") +
  geom_segment(aes(x = 0, xend = x_right, y = thr, yend = thr),
               colour = ink, linewidth = 0.3) +
  geom_point(data = top, aes(x, score), colour = col_top, size = 0.9, shape = 16) +
  annotate("text", x = x_lab, y = thr, hjust = 0, vjust = 0.5, size = 2.3,
           colour = ink, family = "serif", lineheight = 0.9,
           label = "Top 50\n(d = 50)") +
  scale_x_continuous(breaks = axis_x, labels = names(len),
                     expand = expansion(mult = c(0.004, 0.08))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(x = "Chromosome", y = expression("Utility " * hat(omega)[k])) +
  theme_classic(base_size = 8, base_family = "serif") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x  = element_text(size = 5.5, colour = ink),
        axis.text.y  = element_text(size = 6.5, colour = ink),
        axis.line    = element_line(linewidth = 0.3, colour = ink),
        axis.ticks.y = element_line(linewidth = 0.3, colour = ink))

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(data = genes, aes(x, score, label = gene),
                                    family = "serif", fontface = "italic", size = 2.4, nudge_y = 0.012,
                                    segment.size = 0.25, segment.colour = ink, min.segment.length = 0,
                                    seed = 1)
} else {
  p <- p + annotate("text", x = genes$x, y = genes$score, vjust = -0.8, size = 2.4,
                    family = "serif", fontface = "italic", label = genes$gene)
}

ggsave(out_png, p, width = 7.2, height = 2.6, units = "in", dpi = 600)
cat("Saved:", normalizePath(out_png), "\n")
cat("Threshold (d = 50):", round(thr, 6), "\n")
