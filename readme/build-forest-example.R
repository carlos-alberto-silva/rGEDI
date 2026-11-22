# Redraw the saved forest-model validation with shared axes and density colours.
library(data.table)
library(ggplot2)
library(patchwork)
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1L] else Sys.getenv("RGEDI_FOREST_OUTPUT",
                                                  unset = file.path(tempdir(), "rGEDI-forest"))
if (!file.exists(file.path(outdir, "predictions.csv"))) {
  source("inst/scripts/forest_biomass_workflow.R")
}
pred <- fread(file.path(outdir, "predictions.csv"))
stats <- pred[, .(n = .N, r2 = 1 - sum((predicted - observed)^2) / sum((observed - mean(observed))^2),
                  rmse = sqrt(mean((predicted - observed)^2)), mae = mean(abs(predicted - observed)),
                  bias = mean(predicted - observed)), by = stage]
ink <- "#18262D"
theme_set(theme_minimal(base_size = 10) +
          theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                strip.text = element_text(face = "bold")))
save_plate <- function(plot, name, width, height) {
  ggplot2::ggsave(file.path("readme", paste0(name, ".png")), plot,
                 width = width, height = height, dpi = 400, bg = "white")
}
pred[, stage := factor(stage, c("Training", "Spatial and 2022 holdout"))]
stats[, stage := factor(stage, levels(pred$stage))]
max_value <- ceiling(max(pred$observed, pred$predicted) / 50) * 50
stats[, label := sprintf("n = %s\nR\u00b2 = %.2f\nRMSE = %.1f Mg/ha\nBias = %.1f Mg/ha",
                         format(n, big.mark = ","), r2, rmse, bias)]
# Density at the observed points gives each plotted point its actual prediction pair.
pred[, density := MASS::kde2d(observed, predicted, n = 100, lims = c(0, max_value, 0, max_value))$z[
  cbind(pmax(1, pmin(100, 1 + floor(99 * observed / max_value))),
        pmax(1, pmin(100, 1 + floor(99 * predicted / max_value))))], by = stage]
p <- ggplot(pred, aes(observed, predicted)) +
  geom_abline(aes(slope = 1, intercept = 0, linetype = "1:1"), colour = "#6B7375", linewidth = .55) +
  geom_point(aes(colour = density), size = .95, alpha = .8) +
  geom_smooth(aes(linetype = "Linear fit"), method = "lm", formula = y ~ x, se = FALSE,
              colour = "#9D326F", linewidth = .65) +
  geom_text(data = stats, aes(x = 15, y = max_value - 15, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.0, colour = ink) +
  facet_wrap(vars(stage), nrow = 1) +
  scale_colour_viridis_c(option = "magma", name = "Point density", trans = "sqrt") +
  scale_linetype_manual(values = c("1:1" = "dashed", "Linear fit" = "solid"), name = NULL) +
  coord_equal(xlim = c(0, max_value), ylim = c(0, max_value), expand = FALSE) +
  labs(x = "Observed GEDI AGBD (Mg/ha)", y = "Predicted AGBD (Mg/ha)") +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        legend.key.width = grid::unit(12, "mm"), panel.spacing = grid::unit(5, "mm"))
save_plate(p, "fig-forest-validation-full-range", 8.8, 4.45)
outside <- pred[observed > 750 | predicted > 750, .N]
focused <- p + coord_equal(xlim = c(0, 750), ylim = c(0, 750), expand = FALSE) +
  geom_text(data = stats, aes(x = 18, y = 735, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3, colour = ink) +
  labs(caption = paste(outside, "training observation exceeds the displayed limits.",
                       "Statistics use all observations. The companion figure shows the full range.")) +
  theme(plot.caption = element_text(size = 8, hjust = 0))
save_plate(focused, "fig-gedi-model-validation", 8.8, 4.65)

