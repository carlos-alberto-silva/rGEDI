# Package-only README example. Uses the original bundled Level 2A fixture.
# Run from the repository root: Rscript readme/build-package-filtering.R
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else library(rGEDI)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

main <- function() {
  fixture <- system.file("extdata", "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip", package = "rGEDI")
  stopifnot(nzchar(fixture))
  scratch <- tempfile("rgedi-package-example-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  file <- unzip(fixture, exdir = scratch)
  h5 <- readLevel2A(file)
  on.exit(close(h5), add = TRUE)
  m <- getLevel2AM(h5)
  power <- c("BEAM0101", "BEAM0110", "BEAM1000", "BEAM1011")
  labels <- c("All", "Quality = 1", "+ No degradation", "+ Sensitivity >= 0.95", "+ Night", "+ Power beam")
  tables <- list(m)
  tables[[2]] <- tables[[1]][quality_flag == 1]
  tables[[3]] <- tables[[2]][degrade_flag == 0]
  tables[[4]] <- tables[[3]][sensitivity >= 0.95 & sensitivity <= 1]
  tables[[5]] <- tables[[4]][solar_elevation < 0]
  tables[[6]] <- tables[[5]][beam %in% power]
  counts <- data.table(step = factor(labels, levels = labels), n = vapply(tables, nrow, integer(1)))
  distributions <- rbindlist(list(
    m[, .(rh98, group = "All footprints")],
    tables[[6]][, .(rh98, group = "Retained footprints")]
  ))
  distributions[, group := factor(group, levels = c("All footprints", "Retained footprints"))]
  theme_package <- theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"), legend.position = "none")
  a <- ggplot(counts, aes(step, n)) +
    geom_col(fill = "#267F82", width = .68) +
    geom_text(aes(label = n), vjust = -.45, size = 4) +
    scale_y_continuous(expand = expansion(mult = c(0, .13))) +
    labs(title = "A  Apply successive footprint filters", x = NULL, y = "Footprints (n)") +
    theme_package + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  b <- ggplot(distributions[is.finite(rh98)], aes(group, rh98, fill = group)) +
    geom_boxplot(width = .48, outlier.shape = NA, linewidth = .5) +
    geom_point(position = position_jitter(width = .13, height = 0, seed = 24), size = 1.3, alpha = .4) +
    scale_fill_manual(values = c("#B9C9D2", "#267F82")) +
    labs(title = "B  Inspect the retained RH98 distribution", x = NULL, y = "RH98 (m)") +
    theme_package
  p <- a + b + plot_layout(widths = c(1.4, 1)) +
    plot_annotation(title = "Filter GEDI footprint products with rGEDI",
      subtitle = "Bundled Level 2A subset | 18 April 2019 | no external data required",
      caption = "The bundled subset is already quality-screened and acquired at night; some filters leave its count unchanged.",
      theme = theme(plot.title = element_text(face = "bold", size = 18),
                    plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10, hjust = 0)))
  dir.create("readme", showWarnings = FALSE)
  ggsave("readme/fig-package-filtering.png", p, width = 11.6, height = 5.5, dpi = 300, bg = "white")
  ggsave("readme/fig-package-filtering.pdf", p, width = 11.6, height = 5.5, bg = "white")
  print(counts)
}
set.seed(24)
main()
