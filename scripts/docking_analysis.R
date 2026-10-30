library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(purrr)
library(stringr)
library(forcats)
library(patchwork)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(tibble)
library(pheatmap)
library(RColorBrewer)

# --- Output directories (repository layout) ---
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/docking", recursive = TRUE, showWarnings = FALSE)

dir.create("outputs_tables", showWarnings = FALSE)
dir.create("outputs_plots", showWarnings = FALSE)

show_and_save <- function(plot, filename, width, height, dpi = 300) {
  print(plot)
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

############################################################
# Read and unify data
############################################################

vina_11 <- read_excel("data/raw/docking/VINA_11.xlsx")
vina_12 <- read_excel("data/raw/docking/VINA_12.xlsx")
vina_21 <- read_excel("data/raw/docking/VINA_21.xlsx")

vina_all <- bind_rows(
  vina_11,
  vina_12,
  vina_21
)


############################################################
# Read descriptors calculated with RDKit
############################################################

pk_raw <- read_excel("data/raw/docking/PK.xlsx")

pk_unique <- pk_raw %>%
  select(
    Antibiotic_clean,
    SMILES,
    MW,
    LogP,
    TPSA,
    HBD,
    HBA,
    RotBonds,
    HeavyAtoms,
    FormalCharge,
    AromRings,
    FractionSP3
  ) %>%
  distinct(Antibiotic_clean, .keep_all = TRUE)

vina_all_pk <- vina_all %>%
  left_join(pk_unique, by = "Antibiotic_clean")

missing_pk <- vina_all_pk %>%
  filter(is.na(SMILES)) %>%
  distinct(Antibiotic_clean) %>%
  arrange(Antibiotic_clean)

print(missing_pk)


############################################################
# 3. Ligand Efficiency calculation
############################################################

vina_all_pk <- vina_all_pk %>%
  mutate(
    Ligand_Efficiency_row = energy / HeavyAtoms
  )


############################################################
# Rename systems to 1-1, 1-2 and 2-1
############################################################

vina_all_pk <- vina_all_pk %>%
  mutate(
    system = recode(
      as.character(system),
      "11" = "1-1",
      "12" = "1-2",
      "21" = "2-1"
    )
  )


############################################################
# Top scoring by antibiotic and system
############################################################

top_scoring_df <- vina_all_pk %>%
  group_by(system, Antibiotic_clean) %>%
  slice_min(
    order_by = energy,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    top_scoring = energy,
    Ligand_Efficiency_top = top_scoring / HeavyAtoms
  ) %>%
  select(
    system,
    frame,
    grid,
    pose,
    Antibiotic_clean,
    Antibiotic_class,
    top_scoring,
    Ligand_Efficiency_top,
    SMILES,
    MW,
    LogP,
    TPSA,
    HBD,
    HBA,
    RotBonds,
    HeavyAtoms,
    FormalCharge,
    AromRings,
    FractionSP3
  ) %>%
  arrange(system, top_scoring)


write_xlsx(top_scoring_df, "results/tables/docking/top_scoring_df.xlsx")

############################################################
# Ranking by docking 
############################################################

dock_summary <- vina_all_pk %>%
  group_by(system, Antibiotic_clean, Antibiotic_class) %>%
  summarise(
    n_poses = n(),
    min_energy = min(energy, na.rm = TRUE),
    mean_energy = mean(energy, na.rm = TRUE),
    median_energy = median(energy, na.rm = TRUE),
    sd_energy = sd(energy, na.rm = TRUE),
    q25_energy = quantile(energy, 0.25, na.rm = TRUE),
    q75_energy = quantile(energy, 0.75, na.rm = TRUE),
    mean_LE = mean(energy / HeavyAtoms, na.rm = TRUE),
    min_LE = min(energy / HeavyAtoms, na.rm = TRUE),
    MW = first(MW),
    LogP = first(LogP),
    TPSA = first(TPSA),
    HBD = first(HBD),
    HBA = first(HBA),
    RotBonds = first(RotBonds),
    HeavyAtoms = first(HeavyAtoms),
    FormalCharge = first(FormalCharge),
    AromRings = first(AromRings),
    FractionSP3 = first(FractionSP3),
    .groups = "drop"
  )


############################################################
# Frequency top 5%, 10% y 20%
############################################################

top_frequency <- vina_all_pk %>%
  group_by(system, frame, grid) %>%
  mutate(
    rank_energy = rank(energy, ties.method = "first"),
    percentile_rank = percent_rank(energy),
    top_5 = percentile_rank <= 0.05,
    top_10 = percentile_rank <= 0.10,
    top_20 = percentile_rank <= 0.20
  ) %>%
  ungroup() %>%
  group_by(system, Antibiotic_clean) %>%
  summarise(
    freq_top5 = mean(top_5, na.rm = TRUE),
    freq_top10 = mean(top_10, na.rm = TRUE),
    freq_top20 = mean(top_20, na.rm = TRUE),
    .groups = "drop"
  )

dock_summary <- dock_summary %>%
  left_join(top_frequency, by = c("system", "Antibiotic_clean"))


write_xlsx(dock_summary, "results/tables/docking/dock_summary.xlsx")

############################################################
# Affinity by family
############################################################

family_summary <- dock_summary %>%
  group_by(system, Antibiotic_class) %>%
  summarise(
    n_antibiotics = n(),
    median_min_energy = median(min_energy, na.rm = TRUE),
    mean_min_energy = mean(min_energy, na.rm = TRUE),
    median_LE = median(min_LE, na.rm = TRUE),
    median_freq_top10 = median(freq_top10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(system, median_min_energy)

write_xlsx(family_summary, "results/tables/docking/family_summary.xlsx")

############################################################
# Fig. 1
############################################################

family_order <- dock_summary %>%
  group_by(Antibiotic_class) %>%
  summarise(
    med_energy = median(mean_energy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(med_energy) %>%
  pull(Antibiotic_class)

dock_summary2 <- dock_summary %>%
  mutate(
    Antibiotic_class = factor(Antibiotic_class, levels = family_order)
  )

right_marginal <- dock_summary2 %>%
  group_by(Antibiotic_class) %>%
  summarise(
    total_antibiotics = n_distinct(Antibiotic_clean),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = total_antibiotics, y = Antibiotic_class)) +
  geom_col(fill = "grey75", color = "grey30", linewidth = 0.3) +
  labs(x = "# antibiotics", y = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(t = 42, r = 5, b = 5, l = 0)
  )

p_family_box <- ggplot(
  dock_summary2,
  aes(x = Antibiotic_class, y = mean_energy, fill = system)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    color = "grey20",
    alpha = 0.85,
    linewidth = 0.45
  ) +
  facet_wrap(~system, scales = "free_x") +
  coord_flip() +
  labs(
    x = "Antibiotic family",
    y = "Mean Vina Energy (kcal/mol)"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20, face = "bold", color = "black"),
    axis.text = element_text(size = 17, color = "black"),
    strip.text = element_text(size = 18, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

p_family_final <- p_family_box + right_marginal +
  plot_layout(widths = c(5, 1))

show_and_save(
  plot = p_family_final,
  filename = "results/figures/Fig1.png",
  width = 13,
  height = 8,
  dpi = 300
)

write_xlsx(dock_summary2, "results/tables/docking/Table_Fig1.xlsx")

############################################################
# Fig. S1 + table
############################################################
library(gt)
library(scales)

pk_vars <- c(
  "MW",
  "LogP",
  "TPSA",
  "HBD",
  "HBA",
  "RotBonds",
  "HeavyAtoms",
  "FormalCharge",
  "AromRings",
  "FractionSP3"
)


calculate_spearman <- function(data, descriptor, response = "mean_energy") {
  
  x <- data[[descriptor]]
  y <- data[[response]]
  
  complete_rows <- complete.cases(x, y)
  
  x <- x[complete_rows]
  y <- y[complete_rows]
  
  # Evitar errores cuando una variable no presenta variación
  if (
    length(x) < 3 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      tibble(
        rho = NA_real_,
        p_value = NA_real_,
        n = length(x)
      )
    )
  }
  
  test_result <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    rho = unname(test_result$estimate),
    p_value = test_result$p.value,
    n = length(x)
  )
}



spearman_results_long <- map_dfr(
  pk_vars,
  function(variable_name) {
    
    dock_summary %>%
      group_by(system) %>%
      group_modify(
        ~ calculate_spearman(
          data = .x,
          descriptor = variable_name,
          response = "mean_energy"
        )
      ) %>%
      ungroup() %>%
      mutate(
        descriptor = variable_name,
        .before = 1
      )
  }
) %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(
    descriptor,
    system
  )

spearman_rho_wide <- spearman_results_long %>%
  select(
    descriptor,
    system,
    rho
  ) %>%
  pivot_wider(
    names_from = system,
    values_from = rho
  )

spearman_p_wide <- spearman_results_long %>%
  select(
    descriptor,
    system,
    p_value
  ) %>%
  pivot_wider(
    names_from = system,
    values_from = p_value
  )

spearman_n_wide <- spearman_results_long %>%
  select(
    descriptor,
    system,
    n
  ) %>%
  pivot_wider(
    names_from = system,
    values_from = n
  )

write_xlsx(
  list(
    complete_results = spearman_results_long,
    spearman_rho = spearman_rho_wide,
    p_values = spearman_p_wide,
    sample_size = spearman_n_wide
  ),
  "results/tables/docking/Table_S1_spearman_descriptors_vs_mean_energy.xlsx"
)

# fig S1
heavy_atoms_spearman_stats <- dock_summary %>%
  group_by(system) %>%
  group_modify(
    ~ calculate_spearman(
      data = .x,
      descriptor = "HeavyAtoms",
      response = "mean_energy"
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_label = case_when(
      is.na(p_value) ~ "p = NA",
      p_value < 0.001 ~ "p < 0.001",
      TRUE ~ paste0(
        "p = ",
        formatC(
          p_value,
          format = "f",
          digits = 3
        )
      )
    ),
    
    correlation_label = paste0(
      "Spearman ",
      "\u03C1 = ",
      sprintf("%.2f", rho),
      "\n",
      p_label
    )
  )


p_heavy_atoms_energy <- ggplot(
  dock_summary,
  aes(
    x = HeavyAtoms,
    y = mean_energy
  )
) +
  geom_point(
    alpha = 0.70,
    size = 2.3
  ) +
  
  # Linear fit shown only to visualize the general trend.
  # Statistical inference is based on Spearman correlation.
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.18
  ) +
  
  geom_text(
    data = heavy_atoms_spearman_stats,
    aes(
      x = Inf,
      y = Inf,
      label = correlation_label
    ),
    inherit.aes = FALSE,
    hjust = 1.12,
    vjust = 1.25,
    size = 4.2,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ system,
    nrow = 1
  ) +
  
  labs(
    x = "Number of heavy atoms",
    y = "Mean Vina energy (kcal/mol)"
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    axis.title = element_text(
      size = 17,
      face = "bold",
      color = "black"
    ),
    
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    strip.text = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    panel.grid.minor = element_blank()
  )


show_and_save(
  plot = p_heavy_atoms_energy,
  filename = "results/figures/FigS1.png",
  width = 12,
  height = 5.5,
  dpi = 600
)


# Table 
descriptor_labels <- c(
  "MW" = "Molecular weight",
  "LogP" = "LogP",
  "TPSA" = "Topological polar surface area",
  "HBD" = "Hydrogen-bond donors",
  "HBA" = "Hydrogen-bond acceptors",
  "RotBonds" = "Rotatable bonds",
  "HeavyAtoms" = "Heavy atoms",
  "FormalCharge" = "Formal charge",
  "AromRings" = "Aromatic rings",
  "FractionSP3" = "Fraction sp3"
)

spearman_heatmap_table <- spearman_rho_wide %>%
  mutate(
    descriptor = recode(
      descriptor,
      !!!descriptor_labels
    )
  )



gt_spearman_table <- spearman_heatmap_table %>%
  
  gt(
    rowname_col = "descriptor"
  ) %>%
  
  tab_header(
    title = md(
      "**Spearman correlations between molecular descriptors and mean Vina energy**"
    ),
    subtitle = md(
      "Values correspond to Spearman's rho within each Alg–Chi nanoparticle composition"
    )
  ) %>%
  
  cols_label(
    `1-1` = "1-1",
    `1-2` = "1-2",
    `2-1` = "2-1"
  ) %>%
  
  tab_stubhead(
    label = "Molecular descriptor"
  ) %>%
  
  fmt_number(
    columns = any_of(
      c(
        "1-1",
        "1-2",
        "2-1"
      )
    ),
    decimals = 2
  ) %>%
  
  data_color(
    columns = any_of(
      c(
        "1-1",
        "1-2",
        "2-1"
      )
    ),
    
    method = "numeric",
    
    palette = c(
      "#2166AC",
      "#67A9CF",
      "#D1E5F0",
      "#FFFFFF",
      "#FDDBC7",
      "#EF8A62",
      "#B2182B"
    ),
    
    domain = c(
      -1,
      1
    ),
    
    na_color = "grey90"
  ) %>%
  
  tab_style(
    style = cell_text(
      weight = "bold"
    ),
    locations = cells_column_labels(
      everything()
    )
  ) %>%
  
  tab_style(
    style = cell_text(
      weight = "bold"
    ),
    locations = cells_stub()
  ) %>%
  
  tab_style(
    style = cell_text(
      align = "center",
      weight = "bold"
    ),
    locations = cells_body(
      columns = any_of(
        c(
          "1-1",
          "1-2",
          "2-1"
        )
      )
    )
  ) %>%
  
  tab_options(
    table.font.size = px(14),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(13),
    
    table.border.top.width = px(1.5),
    table.border.bottom.width = px(1.5),
    
    column_labels.border.top.width = px(1),
    column_labels.border.bottom.width = px(1),
    
    data_row.padding = px(7),
    
    table.width = pct(100)
  ) %>%
  
  tab_source_note(
    source_note = md(
      paste0(
        "**Interpretation:** because more negative Vina energies represent stronger predicted affinity, ",
        "a negative rho indicates that increasing the descriptor is associated with more favorable docking energy, ",
        "whereas a positive rho indicates an association with less favorable energy."
      )
    )
  )

gt_spearman_table

gtsave(
  gt_spearman_table,
  filename =
    "results/figures/Table_S1.png",
  expand = 5,
  vwidth = 1400
)


############################################################
# 9.3 Plot familia + marginal (mean LE) - Fig. S2
############################################################

family_order <- dock_summary %>%
  group_by(Antibiotic_class) %>%
  summarise(
    med_energy = median(mean_LE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(med_energy) %>%
  pull(Antibiotic_class)

dock_summary2 <- dock_summary %>%
  mutate(
    Antibiotic_class = factor(Antibiotic_class, levels = family_order)
  )

right_marginal <- dock_summary2 %>%
  group_by(Antibiotic_class) %>%
  summarise(
    total_antibiotics = n_distinct(Antibiotic_clean),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = total_antibiotics, y = Antibiotic_class)) +
  geom_col(fill = "grey75", color = "grey30", linewidth = 0.3) +
  labs(x = "# antibiotics", y = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(t = 42, r = 5, b = 5, l = 0)
  )

p_family_box <- ggplot(
  dock_summary2,
  aes(x = Antibiotic_class, y = mean_LE, fill = system)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    color = "grey20",
    alpha = 0.85,
    linewidth = 0.45
  ) +
  facet_wrap(~system, scales = "free_x") +
  coord_flip() +
  labs(
    x = "Antibiotic family",
    y = "Ligand Efficiency (kcal/mol per heavy atom)"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20, face = "bold", color = "black"),
    axis.text = element_text(size = 17, color = "black"),
    strip.text = element_text(size = 18, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

p_family_final <- p_family_box + right_marginal +
  plot_layout(widths = c(5, 1))

show_and_save(
  plot = p_family_final,
  filename = "results/figures/FigS2.png",
  width = 13,
  height = 8,
  dpi = 300
)


write_xlsx(dock_summary2, "results/tables/docking/Table_FigS2.xlsx")


############################################################
# 10. Priorización robusta de candidatos para MD
#
# Estrategia:
# 1) conservar la mejor pose por frame, grid y antibiótico;
# 2) resumir el comportamiento entre frames y grids;
# 3) calcular frecuencia de aparición en los percentiles superiores;
# 4) generar rankings normalizados;
# 5) construir un score de consenso.
############################################################


############################################################
# 10.1 Mejor pose por frame × grid × antibiótico
############################################################

best_pose_frame_grid <- vina_all_pk %>%
  filter(
    !is.na(energy),
    !is.na(HeavyAtoms),
    HeavyAtoms > 0
  ) %>%
  group_by(
    system,
    frame,
    grid,
    Antibiotic_clean,
    Antibiotic_class
  ) %>%
  slice_min(
    order_by = energy,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    best_LE = energy / HeavyAtoms
  )

write_xlsx(
  best_pose_frame_grid,
  "results/tables/docking/best_pose_per_frame_grid.xlsx"
)


############################################################
# 10.2 Ranking dentro de cada frame y grid
#
# Una energía más negativa representa mejor afinidad.
# Un LE más negativo representa mayor eficiencia.
############################################################

best_pose_frame_grid_ranked <- best_pose_frame_grid %>%
  group_by(system, frame, grid) %>%
  mutate(
    percentile_energy_frame = percent_rank(energy),
    percentile_LE_frame = percent_rank(best_LE),
    
    top5_energy_frame = percentile_energy_frame <= 0.05,
    top10_energy_frame = percentile_energy_frame <= 0.10,
    top20_energy_frame = percentile_energy_frame <= 0.20,
    
    top5_LE_frame = percentile_LE_frame <= 0.05,
    top10_LE_frame = percentile_LE_frame <= 0.10,
    top20_LE_frame = percentile_LE_frame <= 0.20
  ) %>%
  ungroup()


############################################################
# 10.3 Resumen robusto por antibiótico y sistema
############################################################

candidate_metrics <- best_pose_frame_grid_ranked %>%
  group_by(
    system,
    Antibiotic_clean,
    Antibiotic_class
  ) %>%
  summarise(
    n_frame_grid = n(),
    
    # Afinidad absoluta
    min_best_energy = min(energy, na.rm = TRUE),
    mean_best_energy = mean(energy, na.rm = TRUE),
    median_best_energy = median(energy, na.rm = TRUE),
    sd_best_energy = sd(energy, na.rm = TRUE),
    se_best_energy = sd_best_energy / sqrt(n_frame_grid),
    q25_best_energy = quantile(energy, 0.25, na.rm = TRUE),
    q75_best_energy = quantile(energy, 0.75, na.rm = TRUE),
    iqr_best_energy = IQR(energy, na.rm = TRUE),
    
    # Ligand efficiency
    min_best_LE = min(best_LE, na.rm = TRUE),
    mean_best_LE = mean(best_LE, na.rm = TRUE),
    median_best_LE = median(best_LE, na.rm = TRUE),
    sd_best_LE = sd(best_LE, na.rm = TRUE),
    se_best_LE = sd_best_LE / sqrt(n_frame_grid),
    
    # Frecuencia de aparición favorable
    freq_top5_energy = mean(top5_energy_frame, na.rm = TRUE),
    freq_top10_energy = mean(top10_energy_frame, na.rm = TRUE),
    freq_top20_energy = mean(top20_energy_frame, na.rm = TRUE),
    
    freq_top5_LE = mean(top5_LE_frame, na.rm = TRUE),
    freq_top10_LE = mean(top10_LE_frame, na.rm = TRUE),
    freq_top20_LE = mean(top20_LE_frame, na.rm = TRUE),
    
    # Propiedades moleculares
    SMILES = first(SMILES),
    MW = first(MW),
    LogP = first(LogP),
    TPSA = first(TPSA),
    HBD = first(HBD),
    HBA = first(HBA),
    RotBonds = first(RotBonds),
    HeavyAtoms = first(HeavyAtoms),
    FormalCharge = first(FormalCharge),
    AromRings = first(AromRings),
    FractionSP3 = first(FractionSP3),
    
    .groups = "drop"
  )


############################################################
# 10.4 Rankings normalizados dentro de cada sistema
#
# Cada score varía entre 0 y 1:
# 1 = mejor comportamiento
# 0 = peor comportamiento
############################################################

candidate_rankings <- candidate_metrics %>%
  group_by(system) %>%
  mutate(
    # Energías negativas: percent_rank entrega 0 al valor más negativo.
    score_energy = 1 - percent_rank(mean_best_energy),
    score_LE = 1 - percent_rank(mean_best_LE),
    
    # Mayor frecuencia = mejor
    score_frequency = percent_rank(freq_top10_energy),
    
    # Menor dispersión = mayor consistencia
    score_consistency = 1 - percent_rank(sd_best_energy),
    
    # Ranking ordinal
    rank_energy = min_rank(mean_best_energy),
    rank_LE = min_rank(mean_best_LE),
    rank_frequency = min_rank(desc(freq_top10_energy)),
    rank_consistency = min_rank(sd_best_energy),
    
    percentile_energy = percent_rank(mean_best_energy),
    percentile_LE = percent_rank(mean_best_LE),
    percentile_frequency = percent_rank(desc(freq_top10_energy)),
    percentile_consistency = percent_rank(sd_best_energy),
    
    # Score de consenso propuesto
    consensus_score =
      0.40 * score_energy +
      0.25 * score_LE +
      0.25 * score_frequency +
      0.10 * score_consistency,
    
    rank_consensus = min_rank(desc(consensus_score))
  ) %>%
  ungroup()


############################################################
# 10.5 Clasificación de perfiles de candidatos
############################################################

candidate_rankings <- candidate_rankings %>%
  mutate(
    top10_energy = percentile_energy <= 0.10,
    top20_energy = percentile_energy <= 0.20,
    
    top10_LE = percentile_LE <= 0.10,
    top20_LE = percentile_LE <= 0.20,
    
    top10_frequency = rank_frequency <= ceiling(
      0.10 * n_distinct(Antibiotic_clean)
    ),
    
    top20_frequency = rank_frequency <= ceiling(
      0.20 * n_distinct(Antibiotic_clean)
    )
  ) %>%
  mutate(
    candidate_profile = case_when(
      top20_energy & top20_LE & top20_frequency ~
        "High-confidence: energy + LE + frequency",
      
      top10_energy & top20_LE ~
        "High absolute affinity",
      
      top10_LE & top20_energy ~
        "High ligand efficiency",
      
      top10_frequency & top20_energy ~
        "Highly recurrent",
      
      top10_energy ~
        "Energy-driven candidate",
      
      top10_LE ~
        "LE-driven candidate",
      
      TRUE ~
        "Background"
    )
  )


############################################################
# 10.6 Corregir definición de top frequency
#
# La sección anterior genera los rankings correctamente,
# pero aquí dejamos los cortes de frecuencia basados
# directamente en percentiles dentro de cada sistema.
############################################################

candidate_rankings <- candidate_rankings %>%
  group_by(system) %>%
  mutate(
    percentile_freq_top10 =
      percent_rank(desc(freq_top10_energy)),
    
    top10_frequency =
      percentile_freq_top10 <= 0.10,
    
    top20_frequency =
      percentile_freq_top10 <= 0.20
  ) %>%
  ungroup() %>%
  mutate(
    candidate_profile = case_when(
      top20_energy & top20_LE & top20_frequency ~
        "High-confidence: energy + LE + frequency",
      
      top10_energy & top20_LE ~
        "High absolute affinity",
      
      top10_LE & top20_energy ~
        "High ligand efficiency",
      
      top10_frequency & top20_energy ~
        "Highly recurrent",
      
      top10_energy ~
        "Energy-driven candidate",
      
      top10_LE ~
        "LE-driven candidate",
      
      TRUE ~
        "Background"
    )
  )


write_xlsx(
  candidate_rankings %>%
    arrange(system, desc(consensus_score)),
  "results/tables/docking/candidate_rankings_all.xlsx"
)


############################################################
# 11. Mejores familias según energía y LE
############################################################

family_candidate_summary <- candidate_rankings %>%
  group_by(system, Antibiotic_class) %>%
  summarise(
    n_antibiotics = n_distinct(Antibiotic_clean),
    
    mean_family_energy =
      mean(mean_best_energy, na.rm = TRUE),
    
    median_family_energy =
      median(mean_best_energy, na.rm = TRUE),
    
    mean_family_LE =
      mean(mean_best_LE, na.rm = TRUE),
    
    median_family_LE =
      median(mean_best_LE, na.rm = TRUE),
    
    mean_family_frequency =
      mean(freq_top10_energy, na.rm = TRUE),
    
    median_family_consensus =
      median(consensus_score, na.rm = TRUE),
    
    n_high_confidence =
      sum(
        candidate_profile ==
          "High-confidence: energy + LE + frequency",
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  group_by(system) %>%
  mutate(
    rank_family_energy =
      min_rank(mean_family_energy),
    
    rank_family_LE =
      min_rank(mean_family_LE),
    
    rank_family_consensus =
      min_rank(desc(median_family_consensus)),
    
    top3_family_energy =
      rank_family_energy <= 3,
    
    top3_family_LE =
      rank_family_LE <= 3
  ) %>%
  ungroup()


write_xlsx(
  family_candidate_summary %>%
    arrange(system, rank_family_consensus),
  "results/tables/docking/family_candidate_summary_robust.xlsx"
)


############################################################
# 11.1 Añadir información de las mejores familias
#      a cada candidato
############################################################

candidate_rankings <- candidate_rankings %>%
  left_join(
    family_candidate_summary %>%
      select(
        system,
        Antibiotic_class,
        rank_family_energy,
        rank_family_LE,
        rank_family_consensus,
        top3_family_energy,
        top3_family_LE
      ),
    by = c("system", "Antibiotic_class")
  )


############################################################
# 12. Definir pool de candidatos para MD
#
# Se incorporan:
# - top por consenso;
# - top por energía;
# - top por LE;
# - alta recurrencia;
# - pertenencia a familias prioritarias.
############################################################

candidate_pool <- candidate_rankings %>%
  mutate(
    selected_by_consensus = rank_consensus <= 10,
    selected_by_energy = rank_energy <= 10,
    selected_by_LE = rank_LE <= 10,
    selected_by_frequency = rank_frequency <= 10,
    
    selected_by_family =
      top3_family_energy |
      top3_family_LE,
    
    candidate_pool_flag =
      selected_by_consensus |
      selected_by_energy |
      selected_by_LE |
      selected_by_frequency |
      (
        selected_by_family &
          (
            top20_energy |
              top20_LE |
              top20_frequency
          )
      )
  ) %>%
  filter(candidate_pool_flag) %>%
  arrange(system, desc(consensus_score))


write_xlsx(
  candidate_pool,
  "results/tables/docking/candidate_pool_for_MD.xlsx"
)


############################################################
# 13. Propuesta automática de candidatos por sistema
#
# Selecciona:
# - 2 mejores por energía;
# - 2 mejores por LE;
# - 2 mejores por consenso;
# - 1 más recurrente;
#
# Después elimina duplicados.
############################################################

selected_energy <- candidate_rankings %>%
  group_by(system) %>%
  slice_min(
    order_by = mean_best_energy,
    n = 2,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(selection_reason = "Top mean Vina energy")


selected_LE <- candidate_rankings %>%
  group_by(system) %>%
  slice_min(
    order_by = mean_best_LE,
    n = 2,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(selection_reason = "Top mean ligand efficiency")


selected_consensus <- candidate_rankings %>%
  group_by(system) %>%
  slice_max(
    order_by = consensus_score,
    n = 2,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(selection_reason = "Top consensus score")


selected_frequency <- candidate_rankings %>%
  group_by(system) %>%
  slice_max(
    order_by = freq_top10_energy,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(selection_reason = "Highest top-10% frequency")


automatic_MD_selection_long <- bind_rows(
  selected_energy,
  selected_LE,
  selected_consensus,
  selected_frequency
)


############################################################
# 13.1 Consolidar razones cuando un antibiótico
#      es seleccionado por varios criterios
############################################################

automatic_MD_selection <- automatic_MD_selection_long %>%
  group_by(system, Antibiotic_clean) %>%
  summarise(
    Antibiotic_class = first(Antibiotic_class),
    
    selection_reason = paste(
      sort(unique(selection_reason)),
      collapse = "; "
    ),
    
    mean_best_energy = first(mean_best_energy),
    sd_best_energy = first(sd_best_energy),
    mean_best_LE = first(mean_best_LE),
    freq_top10_energy = first(freq_top10_energy),
    consensus_score = first(consensus_score),
    rank_energy = first(rank_energy),
    rank_LE = first(rank_LE),
    rank_frequency = first(rank_frequency),
    rank_consensus = first(rank_consensus),
    candidate_profile = first(candidate_profile),
    
    MW = first(MW),
    LogP = first(LogP),
    TPSA = first(TPSA),
    HeavyAtoms = first(HeavyAtoms),
    FormalCharge = first(FormalCharge),
    
    .groups = "drop"
  ) %>%
  arrange(system, desc(consensus_score))


write_xlsx(
  automatic_MD_selection,
  "results/tables/docking/automatic_MD_selection_by_system.xlsx"
)


############################################################
# 14. Análisis de candidatos universales
############################################################

universal_candidates <- candidate_rankings %>%
  group_by(Antibiotic_clean, Antibiotic_class) %>%
  summarise(
    n_systems = n_distinct(system),
    
    mean_consensus_all =
      mean(consensus_score, na.rm = TRUE),
    
    min_consensus_all =
      min(consensus_score, na.rm = TRUE),
    
    mean_energy_all =
      mean(mean_best_energy, na.rm = TRUE),
    
    mean_LE_all =
      mean(mean_best_LE, na.rm = TRUE),
    
    mean_frequency_all =
      mean(freq_top10_energy, na.rm = TRUE),
    
    max_rank_consensus =
      max(rank_consensus, na.rm = TRUE),
    
    systems_top20_energy =
      sum(top20_energy, na.rm = TRUE),
    
    systems_top20_LE =
      sum(top20_LE, na.rm = TRUE),
    
    systems_top20_frequency =
      sum(top20_frequency, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  filter(n_systems == 3) %>%
  mutate(
    universal_high_confidence =
      systems_top20_energy == 3 &
      systems_top20_LE >= 2 &
      systems_top20_frequency >= 2,
    
    universal_score =
      0.60 * mean_consensus_all +
      0.40 * min_consensus_all,
    
    rank_universal =
      min_rank(desc(universal_score))
  ) %>%
  arrange(rank_universal)


write_xlsx(
  universal_candidates,
  "results/tables/docking/universal_candidates.xlsx"
)


############################################################
# 15. Análisis de selectividad entre proporciones
#
# Se considera como proporción preferida aquella con
# mean_best_energy más negativa.
#
# selectivity_gap:
# diferencia entre el segundo mejor y el mejor sistema.
# Un valor positivo y grande representa mayor selectividad.
############################################################

selectivity_long <- candidate_rankings %>%
  select(
    system,
    Antibiotic_clean,
    Antibiotic_class,
    mean_best_energy,
    se_best_energy,
    mean_best_LE,
    consensus_score
  )


selectivity_ranked <- selectivity_long %>%
  group_by(Antibiotic_clean, Antibiotic_class) %>%
  arrange(mean_best_energy, .by_group = TRUE) %>%
  mutate(
    within_antibiotic_rank =
      row_number(),
    
    best_system =
      first(system),
    
    best_energy =
      first(mean_best_energy),
    
    best_se =
      first(se_best_energy),
    
    second_best_energy =
      nth(mean_best_energy, 2),
    
    second_best_se =
      nth(se_best_energy, 2),
    
    selectivity_gap =
      second_best_energy - best_energy,
    
    selectivity_se =
      sqrt(best_se^2 + second_best_se^2),
    
    selectivity_z =
      selectivity_gap / selectivity_se
  ) %>%
  ungroup()


selective_candidates <- selectivity_ranked %>%
  filter(within_antibiotic_rank == 1) %>%
  mutate(
    statistically_supported_selectivity =
      !is.na(selectivity_z) &
      selectivity_z >= 1.96,
    
    practically_relevant_selectivity =
      selectivity_gap >= 0.50,
    
    selective_candidate =
      statistically_supported_selectivity &
      practically_relevant_selectivity,
    
    selectivity_class = case_when(
      selective_candidate ~
        "Supported selective candidate",
      
      practically_relevant_selectivity ~
        "Energetically selective, uncertain variability",
      
      statistically_supported_selectivity ~
        "Statistically distinguishable, small effect",
      
      TRUE ~
        "Non-selective"
    )
  ) %>%
  arrange(desc(selectivity_gap))


write_xlsx(
  selective_candidates,
  "results/tables/docking/selective_candidates.xlsx"
)


############################################################
# 16. Tabla maestra para decisión final de MD
############################################################

master_MD_decision_table <- candidate_rankings %>%
  left_join(
    universal_candidates %>%
      select(
        Antibiotic_clean,
        universal_score,
        rank_universal,
        universal_high_confidence
      ),
    by = "Antibiotic_clean"
  ) %>%
  left_join(
    selective_candidates %>%
      select(
        Antibiotic_clean,
        best_system,
        selectivity_gap,
        selectivity_z,
        selective_candidate,
        selectivity_class
      ),
    by = "Antibiotic_clean"
  ) %>%
  mutate(
    preferred_system =
      best_system,
    
    MD_priority = case_when(
      universal_high_confidence &
        rank_consensus <= 10 ~
        "Priority 1: universal high-confidence",
      
      selective_candidate &
        system == preferred_system &
        rank_consensus <= 20 ~
        "Priority 1: composition-selective",
      
      candidate_profile ==
        "High-confidence: energy + LE + frequency" ~
        "Priority 2: multiparametric",
      
      rank_consensus <= 10 ~
        "Priority 2: consensus",
      
      top10_energy | top10_LE ~
        "Priority 3: single-metric",
      
      TRUE ~
        "Not prioritized"
    )
  ) %>%
  arrange(
    system,
    factor(
      MD_priority,
      levels = c(
        "Priority 1: universal high-confidence",
        "Priority 1: composition-selective",
        "Priority 2: multiparametric",
        "Priority 2: consensus",
        "Priority 3: single-metric",
        "Not prioritized"
      )
    ),
    desc(consensus_score)
  )


write_xlsx(
  list(
    all_candidates = master_MD_decision_table,
    automatic_selection = automatic_MD_selection,
    universal_candidates = universal_candidates,
    selective_candidates = selective_candidates,
    family_summary = family_candidate_summary
  ),
  "results/tables/docking/master_MD_candidate_selection.xlsx"
)


############################################################
# 17. Fig2
############################################################

scatter_cutoffs <- candidate_rankings %>%
  group_by(system) %>%
  summarise(
    energy_cutoff =
      quantile(mean_best_energy, 0.20, na.rm = TRUE),
    
    LE_cutoff =
      quantile(mean_best_LE, 0.20, na.rm = TRUE),
    
    .groups = "drop"
  )


plot_scatter_data <- candidate_rankings %>%
  left_join(
    scatter_cutoffs,
    by = "system"
  ) %>%
  mutate(
    selection_zone = case_when(
      mean_best_energy <= energy_cutoff &
        mean_best_LE <= LE_cutoff ~
        "Top energy and LE",
      
      mean_best_energy <= energy_cutoff ~
        "Top energy only",
      
      mean_best_LE <= LE_cutoff ~
        "Top LE only",
      
      TRUE ~
        "Other"
    ),
    
    label_candidate =
      rank_consensus <= 8 |
      candidate_profile ==
      "High-confidence: energy + LE + frequency"
  )


p_energy_LE <- ggplot(
  plot_scatter_data,
  aes(
    x = mean_best_energy,
    y = mean_best_LE
  )
) +
  geom_vline(
    aes(xintercept = energy_cutoff),
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey45"
  ) +
  geom_hline(
    aes(yintercept = LE_cutoff),
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey45"
  ) +
  geom_point(
    aes(
      fill = selection_zone,
      size = freq_top10_energy
    ),
    shape = 21,
    color = "grey20",
    alpha = 0.85,
    stroke = 0.35
  ) +
  geom_text_repel(
    data = plot_scatter_data %>%
      filter(label_candidate),
    aes(label = Antibiotic_clean),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  facet_wrap(~system, scales = "free") +
  scale_size_continuous(
    name = "Frequency in\ntop 10%",
    range = c(2, 8),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Mean best Vina energy (kcal/mol)",
    y = "Mean ligand efficiency\n(kcal/mol per heavy atom)",
    fill = "Candidate profile"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(
      size = 17,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      size = 13,
      color = "black"
    ),
    strip.text = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(size = 11)
  )


show_and_save(
  plot = p_energy_LE,
  filename =
    "results/figures/Fig2.png",
  width = 15,
  height = 6.5,
  dpi = 600
)

write_xlsx(plot_scatter_data, "results/tables/docking/Table_Fig2.xlsx")


############################################################
# Fig 3
############################################################

MD_selected_keys <- automatic_MD_selection %>%
  distinct(
    system,
    Antibiotic_clean
  ) %>%
  mutate(
    selected_for_MD = TRUE
  )

top_consensus_plot_data <- candidate_rankings %>%
  group_by(system) %>%
  slice_max(
    order_by = consensus_score,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  left_join(
    MD_selected_keys,
    by = c(
      "system",
      "Antibiotic_clean"
    )
  ) %>%
  mutate(
    selected_for_MD =
      replace_na(selected_for_MD, FALSE),
    
    # Crear una etiqueta única para cada antibiótico y sistema.
    # Esto permite ordenar independientemente dentro de los paneles.
    antibiotic_system = paste(
      Antibiotic_clean,
      system,
      sep = "___"
    )
  )



antibiotic_levels <- top_consensus_plot_data %>%
  arrange(
    system,
    consensus_score
  ) %>%
  pull(antibiotic_system) %>%
  unique()


top_consensus_plot_data <- top_consensus_plot_data %>%
  mutate(
    antibiotic_system = factor(
      antibiotic_system,
      levels = antibiotic_levels
    )
  )



x_min_consensus <- min(
  top_consensus_plot_data$consensus_score,
  na.rm = TRUE
)

x_max_consensus <- max(
  top_consensus_plot_data$consensus_score,
  na.rm = TRUE
)

x_range_consensus <-
  x_max_consensus - x_min_consensus

# Evitar margen igual a cero en caso de valores muy similares
x_margin <- max(
  0.015,
  0.08 * x_range_consensus
)

x_left <- max(
  0,
  x_min_consensus - x_margin
)

x_right <- min(
  1,
  x_max_consensus + 2.2 * x_margin
)

segment_start <- x_left

top_consensus_plot_data <- top_consensus_plot_data %>%
  mutate(
    Family_group = case_when(
      Antibiotic_class == "Tetracyclines" ~ "Tetracyclines",
      Antibiotic_class == "Beta-lactams" ~ "Beta-lactams",
      Antibiotic_class == "Streptogramins" ~ "Streptogramins",
      Antibiotic_class == "Quinolones" ~ "Quinolones",
      TRUE ~ "Others"
    ),
    Family_group = factor(
      Family_group,
      levels = c(
        "Tetracyclines",
        "Beta-lactams",
        "Streptogramins",
        "Quinolones",
        "Others"
      )
    )
  )

family_palette <- c(
  "Tetracyclines" = "#1F78B4",
  "Beta-lactams" = "#33A02C",
  "Streptogramins" = "#E31A1C",
  "Quinolones" = "#FF7F00",
  "Others" = "#6B6B6B"
)


top_consensus_plot_data <- top_consensus_plot_data %>%
  mutate(
    star_x = consensus_score + 0.50 * x_margin
  )


p_consensus_improved <- ggplot(
  top_consensus_plot_data,
  aes(
    x = consensus_score,
    y = antibiotic_system
  )
) +
  
  geom_segment(
    aes(
      xend = consensus_score,
      yend = antibiotic_system
    ),
    x = segment_start,
    linewidth = 0.65,
    color = "grey70"
  ) +
  
  geom_point(
    aes(
      size = freq_top10_energy,
      fill = Family_group
    ),
    shape = 21,
    color = "grey15",
    stroke = 0.45,
    alpha = 0.95
  ) +
  
  geom_point(
    data = top_consensus_plot_data %>%
      filter(selected_for_MD),
    aes(
      x = star_x,
      y = antibiotic_system
    ),
    inherit.aes = FALSE,
    shape = 8,
    size = 3.2,
    stroke = 0.8,
    color = "black"
  ) +
  
  facet_wrap(
    ~system,
    scales = "free_y",
    nrow = 1
  ) +
  
  scale_x_continuous(
    name = "Consensus score"
  ) +
  
  scale_y_discrete(
    name = NULL,
    labels = function(x) {
      sub("___.*$", "", x)
    }
  ) +
  
  scale_size_continuous(
    name = "Frequency in\ntop 10%",
    range = c(2.8, 8),
    limits = c(0, 1),
    breaks = c(0.20, 0.40, 0.60, 0.80),
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  scale_fill_manual(
    name = "Antibiotic family",
    values = family_palette
  ) +
  
  coord_cartesian(
    xlim = c(x_left, x_right),
    clip = "off"
  ) +
  
  labs(
    caption = paste0(
      "\u2605 Preselected candidate for molecular dynamics. ",
      "Point size represents recurrence in the top 10% ",
      "across frame-grid combinations."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(
      size = 17,
      face = "bold",
      color = "black"
    ),
    
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    strip.text = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "grey88",
      linewidth = 0.45
    ),
    
    panel.spacing = unit(
      1.2,
      "lines"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    legend.text = element_text(
      size = 10,
      color = "black"
    ),
    
    plot.caption = element_text(
      size = 10.5,
      color = "grey25",
      hjust = 0,
      margin = margin(t = 10)
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )



show_and_save(
  plot = p_consensus_improved,
  filename =
    "results/figures/Fig3.png",
  width = 17,
  height = 9,
  dpi = 600
)

write_xlsx(top_consensus_plot_data, "results/tables/docking/Table_Fig3.xlsx")
  
############################################################
# Fig 4
############################################################

top_consensus_figure_table <- top_consensus_plot_data %>%
  select(
    system,
    Antibiotic_clean,
    Antibiotic_class,
    mean_best_energy,
    mean_best_LE,
    freq_top10_energy,
    sd_best_energy,
    score_energy,
    score_LE,
    score_frequency,
    score_consistency,
    consensus_score,
    rank_consensus,
    selected_for_MD
  ) %>%
  arrange(
    system,
    rank_consensus
  )

write_xlsx(
  top_consensus_figure_table,
  "results/tables/docking/top15_consensus_figure_data.xlsx"
)

############################################################
# 19. Heatmap: consensus scores across systems
############################################################

heatmap_antibiotics <- candidate_rankings %>%
  group_by(Antibiotic_clean) %>%
  summarise(
    maximum_consensus = max(consensus_score, na.rm = TRUE),
    mean_consensus = mean(consensus_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(maximum_consensus)) %>%
  slice_head(n = 30) %>%
  pull(Antibiotic_clean)


heatmap_df <- candidate_rankings %>%
  filter(
    Antibiotic_clean %in% heatmap_antibiotics
  ) %>%
  select(
    Antibiotic_clean,
    system,
    consensus_score
  ) %>%
  pivot_wider(
    names_from = system,
    values_from = consensus_score
  ) %>%
  column_to_rownames("Antibiotic_clean") %>%
  as.matrix()


heatmap_df <- heatmap_df[
  ,
  c("1-1", "1-2", "2-1"),
  drop = FALSE
]


heatmap_palette <- colorRampPalette(
  rev(
    brewer.pal(
      9,
      "RdYlBu"
    )
  )
)(100)


heatmap_breaks <- seq(
  min(
    heatmap_df,
    na.rm = TRUE
  ),
  max(
    heatmap_df,
    na.rm = TRUE
  ),
  length.out = 101
)


p_heatmap <- pheatmap(
  heatmap_df,
  
  scale = "none",
  
  cluster_rows = TRUE,
  
  cluster_cols = FALSE,
  
  color = heatmap_palette,
  
  breaks = heatmap_breaks,
  
  border_color = "white",
  
  fontsize = 12,
  
  fontsize_row = 10,
  
  fontsize_col = 15,
  
  angle_col = 0,
  
  main = "Consensus ranking across Alg–Chi compositions",
  
  treeheight_row = 60,
  
  treeheight_col = 0,
  
  cellwidth = 75,
  
  cellheight = 26,
  
  legend = TRUE,
  
  legend_breaks = pretty(
    range(
      heatmap_df,
      na.rm = TRUE
    ),
    n = 5
  ),
  
  silent = TRUE
)



main_index <- which(
  p_heatmap$gtable$layout$name == "main"
)

p_heatmap$gtable$grobs[[main_index]]$gp <- grid::gpar(
  fontsize = 17,
  font = 2,
  col = "black"
)



png(
  filename =
    "results/figures/FigS3.png",
  width = 3300,
  height = 4500,
  res = 300,
  bg = "white"
  )

grid::grid.newpage()

grid::grid.draw(
  p_heatmap$gtable
)

dev.off()

write_xlsx(candidate_rankings, "results/tables/docking/Table_FigS3.xlsx")

############################################################
# 20. Composition-selective candidates
#
# X       = diferencia entre el mejor y segundo mejor sistema
# Color   = composición Alg-Chi preferida
# Forma   = selectividad estadística y prácticamente apoyada
#
# selective_candidate = TRUE cuando:
#   selectivity_z >= 1.96
#   AND
#   selectivity_gap >= 0.50 kcal/mol
############################################################


############################################################
# 20.1 Parámetros del análisis y de la figura
############################################################

selectivity_threshold <- 0.50
z_threshold <- 1.96
n_candidates_plot <- 25


############################################################
# 20.2 Preparar los datos
############################################################

top_selective_plot <- selective_candidates %>%
  filter(
    !is.na(selectivity_gap),
    selectivity_gap > 0,
    !is.na(best_system)
  ) %>%
  
  slice_max(
    order_by = selectivity_gap,
    n = n_candidates_plot,
    with_ties = FALSE
  ) %>%
  
  mutate(
    # Ordenar según el gap de selectividad:
    # los valores mayores aparecerán arriba.
    Antibiotic_clean = fct_reorder(
      Antibiotic_clean,
      selectivity_gap
    ),
    
    # Mostrar únicamente las composiciones realmente presentes.
    best_system = factor(
      best_system,
      levels = c(
        "1-2",
        "2-1"
      )
    ),
    
    # Etiqueta más legible para la leyenda.
    selectivity_status = factor(
      if_else(
        selective_candidate,
        "Supported",
        "Not supported"
      ),
      levels = c(
        "Not supported",
        "Supported"
      )
    ),
    
    # Valor que se mostrará junto a los candidatos selectivos.
    gap_label = if_else(
      selective_candidate,
      sprintf("%.2f", selectivity_gap),
      NA_character_
    )
  )


############################################################
# 20.3 Comprobaciones de los datos
############################################################

print(
  top_selective_plot %>%
    count(
      best_system,
      selectivity_status
    )
)

print(
  top_selective_plot %>%
    filter(selective_candidate) %>%
    select(
      Antibiotic_clean,
      Antibiotic_class,
      best_system,
      selectivity_gap,
      selectivity_z,
      consensus_score
    )
)


############################################################
# 20.4 Paleta de composiciones
#
# Azul = 1-2
# Verde = 2-1
############################################################

composition_palette <- c(
  "1-2" = "green",
  "2-1" = "red"
)


############################################################
# 20.5 Datos que recibirán etiquetas
############################################################

selective_label_data <- top_selective_plot %>%
  filter(
    selective_candidate,
    !is.na(gap_label)
  )


############################################################
# 20.6 Límites dinámicos del eje X
############################################################

maximum_gap <- max(
  top_selective_plot$selectivity_gap,
  na.rm = TRUE
)

# Espacio adicional a la derecha para las etiquetas.
x_upper_limit <- max(
  0.70,
  maximum_gap + 0.07
)

# Posición vertical de la anotación del umbral.
threshold_annotation_y <- length(
  levels(top_selective_plot$Antibiotic_clean)
) + 0.25


############################################################
# 20.7 Construcción de la figura
############################################################

p_selectivity_improved <- ggplot(
  top_selective_plot,
  aes(
    x = selectivity_gap,
    y = Antibiotic_clean
  )
) +
  
  # Líneas horizontales del lollipop.
  geom_segment(
    aes(
      x = 0,
      xend = selectivity_gap,
      yend = Antibiotic_clean
    ),
    color = "grey70",
    linewidth = 0.70
  ) +
  
  # Umbral práctico de selectividad.
  geom_vline(
    xintercept = selectivity_threshold,
    linetype = "dashed",
    linewidth = 0.85,
    color = "grey25"
  ) +
  
  # Puntos:
  # relleno = formulación preferida;
  # forma = selectividad apoyada o no apoyada.
  geom_point(
    aes(
      fill = best_system,
      shape = selectivity_status
    ),
    size = 4.5,
    color = "black",
    stroke = 0.80,
    alpha = 1
  ) +
  
  # Etiquetar exclusivamente los candidatos selectivos.
  geom_text_repel(
    data = selective_label_data,
    aes(
      label = gap_label
    ),
    direction = "x",
    nudge_x = 0.025,
    hjust = 0,
    size = 4,
    color = "black",
    segment.color = "grey45",
    segment.linewidth = 0.4,
    box.padding = 0.30,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  
  # Anotación del umbral.
  annotate(
    geom = "text",
    x = selectivity_threshold + 0.01,
    y = threshold_annotation_y,
    label = expression(
      Delta * E == 0.50
    ),
    hjust = 0,
    vjust = 1,
    size = 4.2,
    color = "grey15"
  ) +
  
  # Escala de colores de las formulaciones.
  scale_fill_manual(
    name = "Preferred\ncomposition",
    values = composition_palette,
    drop = TRUE
  ) +
  
  # Círculo = no selectivo;
  # triángulo = selectivo.
  scale_shape_manual(
    name = "Composition\nselectivity",
    values = c(
      "Not supported" = 21,
      "Supported" = 24
    ),
    labels = c(
      "Not supported",
      "Supported"
    )
  ) +
  
  # Forzar que la leyenda de color muestre puntos rellenos
  # con los colores reales.
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        size = 4.5,
        color = "black",
        stroke = 0.8,
        alpha = 1
      )
    ),
    
    shape = guide_legend(
      order = 2,
      override.aes = list(
        fill = "white",
        size = 4.5,
        color = "black",
        stroke = 0.8,
        alpha = 1
      )
    )
  ) +
  
  coord_cartesian(
    xlim = c(
      0,
      x_upper_limit
    ),
    clip = "off"
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      x_upper_limit,
      by = 0.10
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  labs(
    x = expression(
      "Selectivity gap, " *
        Delta * "E (kcal/mol)"
    ),
    
    y = NULL,
    
    caption = paste0(
      "Triangles denote antibiotics with statistically supported ",
      "and practically relevant composition selectivity ",
      "(z \u2265 ",
      z_threshold,
      "; \u0394E \u2265 ",
      selectivity_threshold,
      " kcal/mol). Labels indicate \u0394E."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_text(
      size = 17,
      face = "bold",
      color = "black"
    ),
    
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    strip.text = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    panel.grid.minor = element_blank(),
    
    # Ajustes específicos del lollipop
    panel.grid.major.y = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "grey88",
      linewidth = 0.45
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    legend.text = element_text(
      size = 10,
      color = "black"
    ),
    
    legend.key.height = unit(
      0.65,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.30,
      "cm"
    ),
    
    legend.box = "vertical",
    
    plot.caption = element_text(
      size = 10.5,
      color = "grey25",
      hjust = 0,
      margin = margin(t = 10)
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )


############################################################
# 20.8 Mostrar y guardar como PNG
############################################################

show_and_save(
  plot = p_selectivity_improved,
  filename =
    "results/figures/Fig4.png",
  width = 11.5,
  height = 9.5,
  dpi = 600
)

write_xlsx(top_selective_plot, "results/tables/docking/Table_Fig4.xlsx")

############################################################
# 20.10 Tabla asociada exactamente a la figura
############################################################

top_selective_figure_table <- top_selective_plot %>%
  mutate(
    Antibiotic_clean =
      as.character(Antibiotic_clean),
    
    best_system =
      as.character(best_system),
    
    rank_selectivity_gap =
      min_rank(
        desc(selectivity_gap)
      )
  ) %>%
  arrange(
    rank_selectivity_gap
  ) %>%
  select(
    rank_selectivity_gap,
    
    Antibiotic_clean,
    Antibiotic_class,
    
    best_system,
    best_energy,
    best_se,
    
    second_best_energy,
    second_best_se,
    
    selectivity_gap,
    selectivity_se,
    selectivity_z,
    
    statistically_supported_selectivity,
    practically_relevant_selectivity,
    selective_candidate,
    selectivity_class,
    
    mean_best_energy,
    mean_best_LE,
    consensus_score
  )


write_xlsx(
  top_selective_figure_table,
  "results/tables/docking/top_composition_selective_candidates_figure_data.xlsx"
)


selective_candidates %>%
  count(best_system, selective_candidate)

