################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2025

################################################################################

library(tidyverse)
library(rio)
library(mgcv)
library(ggpubr)
library(dplyr)

rm(list = ls())
setwd("E:/Research_Projects/SENECA_network_energy")

# FUNCTIONS -----
partialRsq <- function(fullmodel, redmodel) {
  sse.full <- sum((fullmodel$y - fullmodel$fitted.values)^2)
  sse.red <- sum((redmodel$y - redmodel$fitted.values)^2)
  partialRsq <- (sse.red - sse.full) / sse.red
  return(partialRsq)
}
fit_gam <- function(i, dataset) {
  # Full model
  mod_energy <- mgcv::gam(
    scale(dataset[, 6 + i]) ~ s(Age_Cog, k = 3) +
      scale(Handness) + Gender + scale(TIV) + scale(FD),
    data = dataset,
    method = "REML"
  )

  # Reduced model
  mod_energy_reduced <- mgcv::gam(
    scale(dataset[, 6 + i]) ~
      scale(Handness) + Gender + scale(TIV) + scale(FD),
    data = dataset,
    method = "REML"
  )

  # Extract F and p-values for the smooth term
  model_summary <<- summary(mod_energy)
  f_value <- model_summary$s.table[1, "F"] # F-statistic for smooth term
  edf <- model_summary$s.table[1, "edf"] # Effective degrees of freedom for smooth term
  p_value <- model_summary$s.table[1, "p-value"] # p-value for smooth term
  partial_r2 <- partialRsq(mod_energy, mod_energy_reduced)

  results <<- rbind(
    results,
    data.frame(
      element = colnames(dataset)[6 + i],
      F_value = f_value,
      edf = edf,
      P_value = p_value,
      partial_R2 = partial_r2,
      stringsAsFactors = FALSE
    )
  )
}

plot_energy <- function(i, model, dataset, add.intercept = TRUE) {
  df <- model$model
  mod.intercept <- model$coefficients["(Intercept)"]
  pterms <- predict(model, type = "terms", se.fit = TRUE)

  if (add.intercept == TRUE) {
    pterms.fit <- pterms$fit + mod.intercept
  } else {
    pterms.fit <- pterms$fit
  }
  pterms.sefit <- pterms$se.fit

  colnames(pterms.fit) <- gsub(x = colnames(pterms.fit), pattern = "s\\(", replacement = "") %>%
    gsub(pattern = "\\)", replacement = "")
  colnames(pterms.sefit) <- gsub(x = colnames(pterms.sefit), pattern = "s\\(", replacement = "") %>%
    gsub(pattern = "\\)", replacement = "")

  pterms.df <- data.frame(pterms.fit) %>%
    dplyr::select(Age_Cog) %>%
    plyr::rename(c("Age_Cog" = "fit")) %>%
    cbind(data.frame(pterms.sefit) %>%
      dplyr::select(Age_Cog) %>%
      plyr::rename(c("Age_Cog" = "se.fit"))) %>%
    mutate(
      upr = fit + 1.96 * se.fit,
      lwr = fit - 1.96 * se.fit
    )

  partial.residuals <- data.frame(pterms.fit) %>%
    mutate(across(.cols = everything(), .fns = function(x) {
      x + resid(model)
    })) %>%
    cbind(rawdata = df[, "Age_Cog"])

  plot.df <- cbind(partial.residuals, pterms.df)

  # Check the sign of the trajectory
  save_sign <- plot.df %>% arrange(rawdata)
  if (min(plot.df$fit) < 0){save_sign$fit_positive = plot.df$fit + min(plot.df$fit)} else {save_sign$fit_positive = plot.df$fit}
  if((save_sign$fit_positive[dim(save_sign)[1]] - save_sign$fit_positive[1]) < 0){
    sign[[i]] <<- -1
  } else {sign[[i]] <<- 1}
  
  
  # default plot
  ggplot(plot.df, aes(x = rawdata, y = fit)) +
    geom_hline(yintercept = 0, lty = "dashed") +
    geom_point(size = 2, colour = "gray56") +
    geom_ribbon(aes(x = rawdata, y = fit, ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = .5) +
    geom_line(aes(x = rawdata, y = fit), inherit.aes = FALSE) +
    xlab("Age (in years)") +
    ylab(paste0(colnames(dataset)[6 + i], "\n (a.u.)")) +
    scale_x_continuous(breaks = seq(20, 90, 10)) +
    scale_y_continuous(breaks = round(
      seq(
        min(plot.df$fit), max(plot.df$fit),
        (max(plot.df$fit) - min(plot.df$fit)) / 3
      ),
      digits = 3
    )) +
    theme_pubr(base_size = 14)
}
plot_energy_derv <- function(i, model, this_font_size = 12) {
  # get model derivatives
  derv <- gratia::derivatives(model, interval = "confidence", unconditional = F, order = 1L, eps = 1) # from gratia. "confidence" for point-wise intervals

  # add significance variable (true or false)
  derv <- derv %>%
    mutate(sig = !(0 > .lower_ci & 0 < .upper_ci)) # derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  # new variable with only significant derivatives (non-sig. ones are set to 0)
  derv$sig_deriv <- derv$.derivative * derv$sig

  # print changes range if significant
  if (all(derv$sig_deriv == 0)) {
    cat(sprintf("\n No significant change in G%s \n", i))
  } else {
    cat(sprintf("\nSig change: %1.2f - %1.2f\n", min(derv$Age_Cog[derv$sig == T]), max(derv$Age_Cog[derv$sig == T])))
  }
  
  # plot change
  derv[derv == 0] <- NA
  derv_save <<- derv 
  age_of_onset[[i]] <<- derv[,c("Age_Cog", "sig_deriv")] %>% na.omit() %>% .[1,1]
  
  d <- ggplot(data = derv) +
    geom_tile(aes(x = Age_Cog, y = .5, fill = sig_deriv)) +
    scale_fill_gradient(
      low = "darkblue", high = "darkorange", na.value = "white",
      limits = c(min(derv$sig_deriv), max(derv$sig_deriv))
    )

  d +
    coord_cartesian(xlim = c(18, 90)) +
    labs(x = "", fill = "") +
    scale_x_continuous(breaks = seq(20, 90, 10)) +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = this_font_size),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      text = element_text(size = this_font_size),
      legend.text = element_text(size = this_font_size),
      legend.title = element_text(size = this_font_size),
      axis.title = element_text(size = this_font_size),
      legend.key.width = unit(1, "cm"),
      legend.position = "right",
      plot.margin = unit(c(0, 0, 0.5, 0), "cm"), # Top, left,Bottom, right
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black", size = 1.5),
      axis.line.y = element_line(colour = "black", size = 1.5),
      axis.ticks.length = unit(.25, "cm"),
      axis.text = element_text(size = 50)
    ) +
    guides(fill = guide_colorbar(
      ticks = T,
      ticks.linewidth = 1,
      ticks.colour = "black",
      draw.ulim = T,
      frame.colour = "black",
      frame.linetype = 1,
      frame.linewidth = 1,
      reverse = T,
      direction = "horizontal",
      title.position = "top"
    )) +
    geom_rect(aes(ymin = 0, ymax = 1, xmin = min(Age_Cog), xmax = max(Age_Cog)), color = "black", fill = "white", alpha = 0)
  
}
main_plotting_function <- function(i, dataset, save_label) {
  network_name <- colnames(dataset)[6 + i]
  print(network_name)

  # Full model
  mod_energy <<- mgcv::gam(
    dataset[, 6 + i] ~ s(Age_Cog, k = 3) +
      scale(Handness) + Gender + scale(TIV) + scale(FD),
    data = dataset,
    method = "REML"
  )
  print(summary(mod_energy))

  # Plot residuals
  resid_plot <- plot_energy(i, model = mod_energy, dataset = dataset, add.intercept = TRUE)
  print(resid_plot)

  # Plot derivatives
  dev_plot <- plot_energy_derv(i, model = mod_energy)
  print(dev_plot)

  final_plot <- cowplot::plot_grid(
    rel_heights = c(12, 2),
    plotlist = list(resid_plot, dev_plot),
    align = "v", axis = "lr", ncol = 1
  )

  ggsave(
    plot = final_plot,
    filename = paste0(
      "./output/Glasser360/network_energy/",save_label,"/",
      colnames(dataset[6 + i]), "_",
      save_label, "_energy.pdf"
    ),
    device = "pdf",
    width = 320, height = 290,
    units = "mm"
  )
}

# Whole-brain energy ----------------------------------------------------------
df_energy <- rio::import("./output/Glasser360/Glasser360_df_energy.csv")[, c(1, 3:7, 2)]

# Equivalence testing
bayes_mod <- brms::brm(
  whole_brain_energy ~ Age_Cog +
    scale(Handness) + Gender + scale(TIV) + scale(FD),
  data = df_energy
)

bayestestR::rope_range(bayes_mod)
mod_eq <- bayestestR::equivalence_test(bayes_mod, ci = .89, range = c(-0.001699776, 0.001699776))
print(mod_eq, digits = 5)
graphics.off()
plot(mod_eq) + theme_pubr()


# Synergistic-Redundant Energy ----
df_energy_within_networks <- rio::import("./output/Glasser360/Glasser360_df_energy_networks.csv")
df_energy_between_networks <- rio::import("./output/Glasser360/Glasser360_df_energy_between_networks.csv")

between_set <- df_energy_between_networks[,7:ncol(df_energy_between_networks)]
within_set <- df_energy_within_networks[,7:ncol(df_energy_within_networks)]

# Combine sets and adjust for negative values
union_set <- cbind(within_set, between_set)
adjustment <- abs(min(union_set, na.rm = TRUE)) + 10e-6  # Ensure non-negative values

# Compute synergistic and redundant energy
if (min(union_set) < 0) {
  synergistic_energy <- if (ncol(as.matrix(between_set)) > 1) { # If we consider multiple pairs
    rowMeans(between_set + adjustment)
  } else {
    between_set + adjustment
  }
  redundant_energy <- rowMeans(within_set + adjustment)
}

synergy_index <-cbind(df_energy_between_networks[,1:6], 
                      synergy_index = 
                        1/ncol(within_set) * # Normalization
                        (synergistic_energy/redundant_energy)
)
colnames(synergy_index)[7] <- "Synergistic-Redundant Energy"

results <- data.frame(
  element = character(),
  F_value = numeric(),
  edf = numeric(),
  P_value = numeric(),
  partial_R2 = numeric(),
  stringsAsFactors = FALSE
)

lapply(X = seq(dim(synergy_index)[2] - 6), FUN = fit_gam, dataset = synergy_index)
results["P_value (FDR)"] <- p.adjust(results[, 4], method = "fdr")
results <- results %>% dplyr::filter(results$`P_value (FDR)` < .05)

new_df <- cbind(
  synergy_index[, 1:6],
  synergy_index %>% dplyr::select(any_of(results$element))
)

sign <- list()
age_of_onset <- list()
main_plotting_function(i = 1, dataset = new_df, save_label = "")


# Acceleration between age 45 and age 65
derv_save$sig_deriv[67]/derv_save$sig_deriv[39]
# Acceleration between age 45 and age 75
derv_save$sig_deriv[81]/derv_save$sig_deriv[39]




# FIGURE 2A
df_plot <- cbind(df_energy, syn = synergy_index$`Synergistic-Redundant Energy`)
df_plot %>%
  pivot_longer(c(whole_brain_energy, syn), names_to = "Energies", values_to = "value") %>%
  group_by(Energies) %>%
  mutate(value = as.numeric(scale(value))) %>%
  ggplot(aes(Age_Cog, value, color = Energies)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.1, size = 2) +
  geom_smooth(linewidth = 2.5, method = "gam", formula = y ~ s(x, k = 3), alpha = .1) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(-0.5, 0.7, 0.2)) +
  coord_cartesian(ylim = c(-0.5, 0.5), xlim = c(20, 90)) +
  scale_color_manual(
    values = c("darkorange", 
               "#08519C"),
    labels = c("Synergistic-Redundant Energy", "Whole-brain energy (homeostatis)")
  ) +
  theme_pubclean(base_size = 16) +
  labs(
    x = "Age (years)", 
    y = "Normalized energies (U)"
  ) + 
  annotate(
    "text", x = 20, y = .05, 
    label = "not significant", 
    hjust = 0, vjust = 0, 
    size = 5, color = "#08519C"
  )



# BETWEEN Canonical RSN ---------------------------------------------------

df_energy_between_networks <- rio::import("./output/Glasser360/Glasser360_df_energy_between_networks.csv")

results <- data.frame(
  element = character(),
  F_value = numeric(),
  edf = numeric(),
  P_value = numeric(),
  partial_R2 = numeric(),
  stringsAsFactors = FALSE
)

lapply(X = seq(dim(df_energy_between_networks)[2] - 6), FUN = fit_gam, dataset = df_energy_between_networks)
results["P_value (FDR)"] <- p.adjust(results[, 4], method = "fdr")
results <- results %>% dplyr::filter(results$`P_value (FDR)` < .05, 
                                     results$partial_R2>.01)

# Plot only the elements which are significant after FDR correction
new_df <- cbind(
  df_energy_between_networks[, 1:6],
  df_energy_between_networks %>% dplyr::select(any_of(results$element))
)

sign <- list()
age_of_onset <- list()
lapply(
  X = seq(dim(results)[1]),
  FUN = main_plotting_function,
  dataset = new_df,
  save_label = "between_network"
)


sign_unlisted <- do.call(rbind, sign) %>% as.data.frame()
age_of_onset_unlisted <- do.call(rbind, age_of_onset) %>% as.data.frame()
sign_results <- results %>% cbind(., sign_unlisted, age_of_onset_unlisted) %>% 
  mutate(signed_partial_R2 = partial_R2*V1) %>% arrange(desc(partial_R2))

library(knitr)
library(kableExtra)

# Create and display the kable table
table_html <- sign_results %>% head(4) %>%  
  dplyr::select(element, Age_Cog, F_value, edf, signed_partial_R2, `P_value (FDR)`) %>%
  kable(
    format = "html",
    digits = 2,
    col.names = c("Region", "Age of onset", "F (age)", "edf", "signed partial R2", "pFDR (q < 0.05)")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive")
  ) %>%
  add_header_above(c(" " = 1, "Age-related change in between-network energy" = 5)) %>%
  row_spec(0, bold = TRUE) %>%
  print()



# WITHIN Canonical RSN -----------------------------------------------------------
df_energy_within_networks <- rio::import("./output/Glasser360/Glasser360_df_energy_networks.csv")

results <- data.frame(
  element = character(),
  F_value = numeric(),
  edf = numeric(),
  P_value = numeric(),
  partial_R2 = numeric(),
  stringsAsFactors = FALSE
)

lapply(X = seq(dim(df_energy_within_networks)[2] - 6), FUN = fit_gam, dataset = df_energy_within_networks)
results["P_value (FDR)"] <- p.adjust(results[, 4], method = "fdr")
results <- results %>% dplyr::filter(results$`P_value (FDR)` < .05, 
                                     results$partial_R2>.01)

# Plot only the elements which are significant after FDR correction
new_df <- cbind(
  df_energy_within_networks[, 1:6],
  df_energy_within_networks %>% dplyr::select(any_of(results$element))
)
sign <- list()
age_of_onset <- list()
lapply(
  X = seq(dim(results)[1]),
  FUN = main_plotting_function,
  dataset = new_df,
  save_label = "within_network"
)

sign_unlisted <- do.call(rbind, sign) %>% as.data.frame()
age_of_onset_unlisted <- do.call(rbind, age_of_onset) %>% as.data.frame()
sign_results <- results %>% cbind(., sign_unlisted, age_of_onset_unlisted) %>% 
  mutate(signed_partial_R2 = partial_R2*V1) %>% arrange(desc(partial_R2))


library(knitr)
library(kableExtra)

# Create and display the kable table
table_html <- sign_results %>% head(3) %>%  
  dplyr::select(element, Age_Cog, F_value, edf, signed_partial_R2, `P_value (FDR)`) %>%
  kable(
    format = "html",
    digits = 2,
    col.names = c("Region", "Age of onset", "F (age)", "edf", "signed partial R2", "pFDR (q < 0.05)")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive")
  ) %>%
  add_header_above(c(" " = 1, "Age-related change in within-network energy" = 5)) %>%
  row_spec(0, bold = TRUE) %>%
  print()

# Region Energy at the whole-brain level ----
library(ggsegGlasser)

region_energy <- rio::import("./output/Glasser360/Glasser360_region_energy.csv")[, 2:361] %>%
  colMeans() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(new_label = gsub("\\.", "-", rowname))
colnames(region_energy) <- c("old_label", "Region Energy", "new_label")

atlas <- ggsegGlasser::glasser %>%
  as.data.frame() %>%
  mutate(new_label = substring(label, 4)) %>%
  left_join(., region_energy, by = "new_label") %>%
  as_brain_atlas() %>%
  as_tibble()

# Average instability index
ggplot() +
  geom_brain(
    atlas = atlas,
    mapping = aes(fill = `Region Energy`),
    position = position_brain(hemi ~ side),
    color = "black",
    size = 0.25,
    show.legend = T
  ) +
  theme_void() +
  scale_fill_viridis_c(option = "D")

gradients <- RcppCNPy::npyLoad("../SENECA_gradients/output/Glasser360/Glasser360_mean_aligned_grad.npy") %>%
  as.data.frame()
mod_grad <- lm(`Region Energy` ~ V1 + V2 + V3,
  data = cbind(region_energy, gradients)
)
summary(mod_grad)
effectsize::eta_squared(model = mod_grad, alternative = "two.sided")

# Region energy BETWEEN DMN-FPN ---------------------------------------------------

df_energy_DMN_FPN <- rio::import("./output/Glasser360/region_contribution_to_DMN_FPN_between_energy.csv")

results <- data.frame(
  element = character(),
  F_value = numeric(),
  edf = numeric(),
  P_value = numeric(),
  partial_R2 = numeric(),
  stringsAsFactors = FALSE
)

lapply(X = seq(dim(df_energy_DMN_FPN)[2] - 6), FUN = fit_gam, dataset = df_energy_DMN_FPN)
results["P_value (FDR)"] <- p.adjust(results[, 4], method = "fdr")
results <- results %>% dplyr::filter(results$`P_value (FDR)` < .05,
                                     results$partial_R2>.01)

# Plot only the elements which are significant after FDR correction
new_df <- cbind(
  df_energy_DMN_FPN[, 1:6],
  df_energy_DMN_FPN %>% dplyr::select(any_of(results$element))
)

sign <- list()
age_of_onset <- list()
lapply(
  X = seq(dim(results)[1]),
  FUN = main_plotting_function,
  dataset = new_df,
  save_label = "DMN_FPN_synergy"
)

sign_unlisted <- do.call(rbind, sign) %>% as.data.frame()
age_of_onset_unlisted <- do.call(rbind, age_of_onset) %>% as.data.frame()
# Only R AVI has a U-shaped trajectory,  R PreS decline but stabilizes
# All other regions decline linearly or accelerate
sign_results <- results %>% cbind(., sign_unlisted, age_of_onset_unlisted) %>% 
  mutate(signed_partial_R2 = partial_R2*V1) %>% arrange(desc(partial_R2))


# Perform k-means clustering
set.seed(123)
sign_results$kmeans_result <- kmeans(sign_results$edf, centers = 2)$cluster

# Linear trajectories
atlas <- ggsegGlasser::glasser %>%
  as.data.frame() %>%
  mutate(new_label = substring(label, 4)) %>%
  left_join(.,
    sign_results %>% 
      dplyr::filter(kmeans_result == 1) %>%
      mutate(new_label = gsub("\\.", "-", element)),
    by = "new_label"
  ) %>%
  as_brain_atlas() %>%
  as_tibble()

ggplot() +
  geom_brain(
    atlas = atlas,
    mapping = aes(fill = signed_partial_R2),
    position = position_brain(hemi ~ side),
    color = "black",
    size = 0.25,
    show.legend = T
  ) +
  theme_void() +
  scale_fill_gradient2()

atlas <- ggsegGlasser::glasser %>%
  as.data.frame() %>%
  mutate(new_label = substring(label, 4)) %>%
  left_join(.,
            sign_results %>% 
              dplyr::filter(kmeans_result == 2) %>%
              mutate(new_label = gsub("\\.", "-", element)),
            by = "new_label"
  ) %>%
  as_brain_atlas() %>%
  as_tibble()

ggplot() +
  geom_brain(
    atlas = atlas,
    mapping = aes(fill = signed_partial_R2),
    position = position_brain(hemi ~ side),
    color = "black",
    size = 0.25,
    show.legend = T
  ) +
  theme_void() +
  scale_fill_gradient2()


# Tables and Correlations with MNI coordinates ----

HCP_coords <- rio::import("./code/HCPex_CAB_NP/HCP-MMP1_MNI_coord.csv")[, c(2, 4:5, 7, 10:13)] %>% as.data.frame()
HCP_coords$region[HCP_coords$region == "H"] <- "Hipp"

results_prep <- sign_results %>% mutate(
  new_label = gsub("\\.", "-", element),
  LR = substring(element, 1, 1),
  region = substring(new_label, 3)
)

HCP_composition <- rio::import("./output/Glasser360/Glasser360_RSN_composition.csv")[,2:3] %>% as.data.frame() %>% dplyr::rename(element=Region)
results_merge <- left_join(results_prep, HCP_coords, by = c("region", "LR")) %>% 
  left_join(., HCP_composition, by = "element")

results_save_early <- results_merge %>%
  dplyr::arrange(desc(abs(signed_partial_R2))) %>%
  dplyr::filter(kmeans_result == 1)
results_save_midlife <- results_merge %>%
  dplyr::arrange(desc(abs(signed_partial_R2))) %>%
  dplyr::filter(kmeans_result == 2)

library(knitr)
library(kableExtra)

# Create and display the kable table
table_html <- results_save_early %>% head(10) %>%  
  dplyr::select(element, Age_Cog, F_value, edf, signed_partial_R2, `P_value (FDR)`, cortex, network) %>%
  kable(
    format = "html",
    digits = 2,
    col.names = c("Region", "Age of onset", "F (age)", "edf", "signed partial R2", "pFDR (q < 0.05)", "Location", "Network")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive")
  ) %>%
  add_header_above(c(" " = 1, "Top 10 regions disrupting DMN-FPN dynamics in early adulthood" = 7)) %>%
  row_spec(0, bold = TRUE) %>%
  print() # Print the table

table_html <- results_save_midlife %>% head(10) %>%  
  dplyr::select(element, Age_Cog, F_value, edf, signed_partial_R2, `P_value (FDR)`, cortex, network) %>%
  kable(
    format = "html",
    digits = 2,
    col.names = c("Region", "Age of onset", "F (age)", "edf", "signed partial R2", "pFDR (q < 0.05)", "Location", "Network")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive")
  ) %>%
  add_header_above(c(" " = 1, "Top 10 regions disrupting DMN-FPN dynamics in midlife" = 7)) %>%
  row_spec(0, bold = TRUE) %>%
  print() # Print the table

# Posterior-to-Anterior Y - MNI
cor.test(results_merge$kmeans_result * results_merge$partial_R2, results_merge$`y-cog`)
# Inferior-to-Superior Z - MNI
cor.test(results_merge$kmeans_result * results_merge$partial_R2, results_merge$`z-cog`)
# medial-lateral - MNI
medial <- (results_merge$`x-cog` - median(results_merge$`x-cog`))^2
cor.test(results_merge$kmeans_result * results_merge$partial_R2, medial)

# Region energy within DMN -------------------------------------------------

DMN_region_energy <- rio::import("./output/Glasser360/Glasser360_DMN_region_energy.csv") %>% 
  dplyr::select(
    colnames(.)[c(1:6)],
    "L_31pd", "R_31pd", "L_31pv", "R_31pv", "L_7m", "R_7m", 
    "L_v23ab", "R_v23ab", "L_d23ab",  "R_d23ab", 
    "L_31a", # Homolog region is 60% FPN and 25% DMN
    "L_23d", "R_23d", "L_POS1", "R_POS1")

results <- data.frame(
  element = character(),
  F_value = numeric(),
  edf = numeric(),
  P_value = numeric(),
  partial_R2 = numeric(),
  stringsAsFactors = FALSE
)

lapply(X = seq(dim(DMN_region_energy)[2] - 6), FUN = fit_gam, dataset = DMN_region_energy)
results["P_value (FDR)"] <- p.adjust(results[, 4], method = "fdr")
results <- results %>% dplyr::filter(results$`P_value (FDR)` < .05,
                                     results$partial_R2>.01)

# Plot only the elements which are significant after FDR correction
new_df <- cbind(
  DMN_region_energy[, 1:6],
  DMN_region_energy %>% dplyr::select(any_of(results$element))
)

sign <- list()
age_of_onset <- list()
lapply(
  X = seq(dim(results)[1]),
  FUN = main_plotting_function,
  dataset = new_df,
  save_label = "within_DMN"
)

sign_unlisted <- do.call(rbind, sign) %>% as.data.frame()
age_of_onset_unlisted <- do.call(rbind, age_of_onset) %>% as.data.frame()
sign_results <- results %>% cbind(., sign_unlisted, age_of_onset_unlisted) %>% 
  mutate(signed_partial_R2 = partial_R2*V1) %>% arrange(desc(partial_R2))

# Perform k-means clustering
# set.seed(123)
# sign_results$kmeans_result <- kmeans(sign_results$edf, centers = 2)$cluster

library(knitr)
library(kableExtra)

# All trajectories
atlas <- ggsegGlasser::glasser %>%
  as.data.frame() %>%
  mutate(new_label = substring(label, 4)) %>%
  left_join(.,
            sign_results %>% 
              mutate(new_label = gsub("\\.", "-", element)),
            by = "new_label"
  ) %>%
  as_brain_atlas() %>%
  as_tibble()

ggplot() +
  geom_brain(
    atlas = atlas,
    mapping = aes(fill = signed_partial_R2),
    position = position_brain(hemi ~ side),
    color = "black",
    size = 0.25,
    show.legend = T
  ) +
  theme_void() +
  scale_fill_gradient2()



# Tables and Correlations with MNI coordinates ----

library(knitr)
library(kableExtra)

# Create and display the kable table
table_html <- results_merge %>%  
  dplyr::select(element, Age_Cog, F_value, edf, signed_partial_R2, `P_value (FDR)`, cortex, network) %>%
  kable(
    format = "html",
    digits = 2,
    col.names = c("Region", "Age of onset", "F (age)", "edf", "signed partial R2", "pFDR (q < 0.05)", "Location", "Network")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive")
  ) %>%
  add_header_above(c(" " = 1, "PCC regions contributing to age-related energy changes within the DMN" = 7)) %>%
  row_spec(0, bold = TRUE) %>%
  print() 

# Left-to-Right X - MNI
cor.test(results_merge$edf * results_merge$signed_partial_R2, results_merge$`x-cog`)


# Region energy within FPN -------------------------------------------------

df_instability_index <- rio::import("./output/Glasser360/Glasser360_FPN_region_energy.csv")

results <- data.frame(
  element = character(),
  F_value = numeric(),
  edf = numeric(),
  P_value = numeric(),
  partial_R2 = numeric(),
  stringsAsFactors = FALSE
)

lapply(X = seq(dim(df_instability_index)[2] - 6), FUN = fit_gam, dataset = df_instability_index)
results["P_value (FDR)"] <- p.adjust(results[, 4], method = "fdr")
results <- results %>% dplyr::filter(results$`P_value (FDR)` < .05,
                                     results$partial_R2>.01)

# Plot only the elements which are significant after FDR correction
new_df <- cbind(
  df_instability_index[, 1:6],
  df_instability_index %>% dplyr::select(any_of(results$element))
)

# Check which significant regions increase or decrease their instability index with age
sign <- list()
lapply(
  X = seq(dim(results)[1]),
  FUN = main_plotting_function,
  dataset = new_df,
  save_label = "within_FPN_region"
)

sign_unlisted <- do.call(rbind, sign) %>% as.data.frame()
results <- results %>% cbind(., sign_unlisted) %>% 
  mutate(signed_partial_R2 = partial_R2*V1)

# All trajectories
library(ggseg)
library(ggsegGlasser)
atlas <- ggsegGlasser::glasser %>%
  as.data.frame() %>%
  mutate(new_label = substring(label, 4)) %>%
  left_join(.,
            results %>%
              dplyr::filter(edf < 1.75) %>% 
              mutate(new_label = gsub("\\.", "-", element)),
            by = "new_label"
  ) %>%
  as_brain_atlas() %>%
  as_tibble()

ggplot() +
  geom_brain(
    atlas = atlas,
    mapping = aes(fill = signed_partial_R2),
    position = position_brain(hemi ~ side),
    color = "black",
    size = 0.25,
    show.legend = T
  ) +
  theme_void() +
  scale_fill_viridis_c(option = "A")
