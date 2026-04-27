rm(list=ls())

require('tidyverse')
require('gtsummary')
require('sjmisc')
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(lmtp)




report_results <- function(main_list, Outcome) {
  
  control_DR <- lapply(main_list, function(x) summary(x$tmle_control$density_ratios))
  trt_DR <- lapply(main_list, function(x) summary(x$tmle_trt$density_ratios))
  Effect_Mdf <- lapply(main_list, function(x) summary(x$effect_mdf))
  Effect_Contrast <- lapply(main_list, function(x) x$tmle_contrast$vals)
  
  if (Outcome=="CJP") {
    # 1. Extract the data from elements 2 through 7
    plot_df <- map_df(2:7, function(i) {
      data.frame(
        month_label = paste0("Month ", i - 1),
        theta = Effect_Contrast[[i]]$theta,
        low   = Effect_Contrast[[i]]$conf.low,
        high  = Effect_Contrast[[i]]$conf.high,
        order = i - 1 # Used to ensure chronological order on the x-axis
      )
    })
    # Determine significance for color mapping
    plot_df$significant <- (plot_df$low > 0 & plot_df$high > 0) | (plot_df$low < 0 & plot_df$high < 0)
    
    # 2. Create the Forest Plot / Timeline Plot
    p<- ggplot(plot_df, aes(x = reorder(month_label, order), y = theta, color = significant)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + # Reference line (no effect)
      #geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, size = 0.8, color = "#2c3e50") +
      geom_errorbar(aes(ymin = low, ymax = high), width = 0.05, size = 1) +
      
      #geom_point(size = 4, color = "#e67e22") + # Bold color for the point estimate
      geom_point() + # Bold color for the point estimate
      
      scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#2c3e50")) +
      
      labs(
        #title = "Timeline of Estimated Effects (Consequences)",
        #subtitle = "Monthly contrasts relative to baseline (Months 1-6)",
        x = "Observation Period",
        y = "Parameter estimate"#,caption = "Error bars represent 95% Confidence Intervals"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
  } else if (Outcome=="CJS") {
    
    p <- list()
    
    # 2. Create an empty list to store the 4 plots
    plot_list <- list()
    
    for (i in target_indices) {
      # Extract results for the specific index i
      res <- Effect_Contrast[[i]]
      
      # Create a small data frame for this plot
      # i no longer defines "month", so we use a generic Assessment label
      plot_df <- data.frame(
        label = paste0("Assessment ", i),
        theta = res$theta,
        low   = res$conf.low,
        high  = res$conf.high
      )
      
      # Determine significance (CI excludes 0)
      # TRUE if both low and high are on the same side of zero
      is_significant <- (plot_df$low > 0 & plot_df$high > 0) | (plot_df$low < 0 & plot_df$high < 0)
      
      # Set color: red if significant, black if not
      point_color <- if (is_significant) "#c0392b" else "#2c3e50"
      
      # 3. Create the plot for this index
      p[[i]] <- ggplot(plot_df, aes(x = label, y = theta)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
        geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, size = 1, color = "#2c3e50") +
        geom_point(size = 4, color = point_color) +
        labs(
          title = paste("Effect Estimate - Assessment", i),
          subtitle = paste("Point estimate with 95% Confidence Intervals"),
          x = NULL,
          y = expression(theta),
          caption = "Error bars indicate 95% CI; Red indicates significance (CI excludes 0)"
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
      
    }
  }
  
  
  return (list(control_DR=control_DR,
               trt_DR=trt_DR,
               Effect_Mdf=Effect_Mdf,
               Effect_Contrast=Effect_Contrast,
               p=p))
}
# CJP ----
load("CJP.RData")

## CJP Main graph ----
# 1. Helper function to safely extract the first element's contrast data
extract_first_theta <- function(target_list, label) {
  # Accessing the first element of the list as requested
  res <- target_list[[1]]$tmle_contrast$vals
  
  print(target_list[[1]]$effect_mdf %>% summary)
  
  # Logic: Significant if the interval [low, high] does NOT cross 0
  is_sig <- !(res$conf.low <= 0 & res$conf.high >= 0)
  
  data.frame(
    Analysis = label,
    theta = res$theta,
    low = res$conf.low,
    high = res$conf.high,
    pval = res$p.value,
    sig_status = ifelse(is_sig, "Sig", "NonSig")  )
}

# 2. Combine the first elements of all four lists
plot_df <- bind_rows(
  extract_first_theta(target_list=ITT_CC_CJP, label="ITT: CC"),
  extract_first_theta(ITT_IPCW_CJP, "ITT: IPCW"),
  extract_first_theta(PP_IPCW_CJP, "PP: IPCW"),
  extract_first_theta(PP_STD_CJP, "PP: CC")
)

# Ensure the facets appear in the specified order
plot_df$Analysis <- factor(plot_df$Analysis, levels = c(
   "ITT: IPCW", "ITT: CC","PP: IPCW", "PP: CC"
))

# 3. Create the faceted comparison plot
ggplot(plot_df, aes(x = Analysis, y = theta, color=sig_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  #geom_errorbar(aes(ymin = low, ymax = high), width = 0.1, size = 1.1) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.05, size = 1) +
  
  #geom_point(size = 4.5, color = "#2980b9") +
  geom_point() +
  
  # Map colors: Red for significant, Black for non-significant
  scale_color_manual(values = c("Sig" = "red", "NonSig" = "black")) +
  # Use free_x to ensure each facet only shows its own label
  facet_wrap(~Analysis, scales = "free_x", nrow = 2) +
  labs(
    #title = "Sensitivity Analysis: Cost-Consequence Estimate Consistency",
    #subtitle = "Comparison of theta across ITT and Per-Protocol estimation methods",
    x = NULL, 
    y = "Parameter estimate"
    #y = expression(paste("Primary Outcome Effect (", theta, ")"))
    #caption = "ITT: Intention-to-Treat; PP: Per-Protocol; IPCW: Inverse Probability of Censoring Weighting"
  ) +
  theme_bw() +
  theme(
    legend.position = 'none',
    strip.background = element_blank(),
    #strip.background = element_rect(fill = "#f7f9f9"),
    strip.text = element_text(face = "bold", size = 14, hjust = 0),
    axis.text.x = element_blank(), # Hidden because the facet title carries the label
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )


## CJP Secondary anal ----
ITT_CC_CJP_R <- report_results(ITT_CC_CJP, Outcome="CJP"); ITT_CC_CJP_R$Effect_Mdf[[1]]; ITT_CC_CJP_R$p
ITT_IPCW_CJP_R <- report_results(ITT_IPCW_CJP, Outcome="CJP"); ITT_IPCW_CJP_R$Effect_Mdf[[1]]; ITT_IPCW_CJP_R$p
PP_IPCW_CJP_R <- report_results(PP_IPCW_CJP, Outcome="CJP"); PP_IPCW_CJP_R$Effect_Mdf[[1]]; PP_IPCW_CJP_R$p
PP_STD_CJP_R <- report_results(PP_STD_CJP, Outcome="CJP"); PP_STD_CJP_R$Effect_Mdf[[1]]; PP_STD_CJP_R$p


# CJS ----
ClaCos_anal <- FALSE
load("CJS_t1.RData")
load("CJS_t2.RData")
all_outcomes_CJS <- c("SERS_SC_TOT", "ISMI_SC_TOT", "STORI_CROISSANCE", 
                      "STAIA_TOT", "STAIB_TOT")
if (ClaCos_anal) {
  all_outcomes_CJS <- c( "ACSO", "ACSO_H", 
                         "MASC_NB", 
                         "AIHQ_HB_NB","AIHQ_ATTRIB_RESP_NB", "AIHQ_AB_NB",  
                         "PERSO_FLU_PONC", "PERSO_INTER_LIB", "PERSO_INTER_INDIC", "PERSO_INTER_TOTAL", "PERSO_CONV_SOC", 
                         "TREF_BR_TOTAL_NB","TREF_SD_TOTAL_NB")
  load("CJS_t1_ClaCos.RData")
  load("CJS_t2_ClaCos.RData")
}
# CJS Define the target indices 
target_indices <- seq_along(all_outcomes_CJS)
target_outcomes <- all_outcomes_CJS[target_indices]

report_final_CJS <- function(main_list_t1, main_list_t2) {
  
  main_list_t1_R <- report_results(main_list = main_list_t1, Outcome="CJS"); 
  main_list_t2_R <- report_results(main_list = main_list_t2, Outcome="CJS");
  
  for (i in seq_along(target_indices)) {
    # print(target_outcomes[i])
    # print(main_list_t1_R$Effect_Mdf[[i]])
    # print(main_list_t2_R$Effect_Mdf[[i]])
  }
  
  combined_plot_list <- list()
  
  for (i in seq_along(target_indices)) {
    idx <- target_indices[i]
    
    # 2. Extract data from t1 and t2
    # Assuming Effect_Contrast_t1 and _t2 are the source lists for your objects
    df_t1 <- data.frame(
      which_outcome = all_outcomes_CJS[idx],
      theta = main_list_t1_R$Effect_Contrast[[idx]]$theta,
      low   = main_list_t1_R$Effect_Contrast[[idx]]$conf.low,
      high  = main_list_t1_R$Effect_Contrast[[idx]]$conf.high,
      pval  = main_list_t1_R$Effect_Contrast[[idx]]$p.value,
      
      model = "Immediate post-intervention"
    )
    
    df_t2 <- data.frame(
      which_outcome = all_outcomes_CJS[idx],
      theta = main_list_t2_R$Effect_Contrast[[idx]]$theta,
      low   = main_list_t2_R$Effect_Contrast[[idx]]$conf.low,
      high  = main_list_t2_R$Effect_Contrast[[idx]]$conf.high,
      pval  = main_list_t2_R$Effect_Contrast[[idx]]$p.value,
      
      model = "6 months post-intervention"
    )
    print(df_t1);print(df_t2);
    
    # Bind them together
    plot_df <- rbind(df_t1, df_t2)
    plot_df$model <- factor(plot_df$model , levels = c(
      "Immediate post-intervention", "6 months post-intervention"
    ))
    
    # Determine significance for color mapping
    plot_df$significant <- (plot_df$low > 0 & plot_df$high > 0) | (plot_df$low < 0 & plot_df$high < 0)
    
    # 3. Build the Plot
    p <- ggplot(plot_df, aes(x = model, y = theta, color = significant)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
      geom_errorbar(aes(ymin = low, ymax = high), width = 0.05, size = 1) +
      geom_point() +
      # Manual color scale: Red for significant, Black/DarkBlue for not
      scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#2c3e50")) +
      labs(
        title = target_outcomes[idx],
        #subtitle = "T1 vs T2 Estimates on Shared Y-Axis",
        x = NULL,
        y = "Parameter estimate"
        #y = expression(theta)
        #caption = "Red indicates 95% CI excludes 0"
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    
    # 4. Store and Save
    combined_plot_list[[i]] <- p
  }
  
  return(combined_plot_list)
  
}

plots_ITT_CC_CJS <- report_final_CJS(main_list_t1 = ITT_CC_CJS_t1,
                                     main_list_t2 = ITT_CC_CJS_t2)
plots_ITT_IPCW_CJS <- report_final_CJS(ITT_IPCW_CJS_t1,ITT_IPCW_CJS_t2)
plots_PP_IPCW_CJS <- report_final_CJS(PP_IPCW_CJS_t1,PP_IPCW_CJS_t2)
plots_PP_STD_CJS <- report_final_CJS(PP_STD_CJS_t1,PP_STD_CJS_t2)

# 1. Store your 4 result objects in a named list for easier iteration
all_results <- list(
  ITT_CC   = plots_ITT_CC_CJS,
  ITT_IPCW = plots_ITT_IPCW_CJS,
  PP_IPCW  = plots_PP_IPCW_CJS,
  PP_STD   = plots_PP_STD_CJS
)

# 2. Create a list to store the 4 "Master" combined plots
master_comparison_plots <- list()

# 3. Loop through the 4 plot positions (1, 2, 3, 4) in each list
for (i in seq_along(target_indices)) {
  
  # Extract the i-th plot from each of the four analysis sets
  p_cc   <- all_results$ITT_CC[[i]]   + ggtitle("ITT: CC")
  p_itt_i <- all_results$ITT_IPCW[[i]] + ggtitle("ITT: IPCW")
  p_pp_i  <- all_results$PP_IPCW[[i]]  + ggtitle("PP: IPCW")
  p_pp_s  <- all_results$PP_STD[[i]]   + ggtitle("PP: CC")
  
  # 4. Combine into a 2x2 grid using patchwork syntax
  # (p1 | p2) / (p3 | p4) creates two rows of two
  library(grid)
  y_label <- wrap_elements(
    panel = textGrob(
      "Parameter estimate", 
      rot = 90#, gp = gpar(fontsize = 12, fontface = "bold")
    )
  )
  
  master_plot <- (p_itt_i | p_cc) / (p_pp_i | p_pp_s)
  
  # 5. Add universal styling and labels
  master_plot_clean <- (y_label + (master_plot & labs(y = NULL)))
  master_comparison_plots[[i]] <- master_plot_clean + 
    plot_layout(guides = "collect",widths = c(1, 20))  &
    # plot_annotation(
    #   left = grid::textGrob(
    #     expression(paste("Effect Size (", theta, ")")), 
    #     rot = 90, 
    #     gp = grid::gpar(fontsize = 12, fontface = "bold")
    #   )
    #   title =target_outcomes[i],
    #   #subtitle = "Comparing ITT and PP estimates across different confounding adjustment methods",
    #   caption = "Red points indicate statistical significance (95% CI excludes 0)",
    #   tag_levels = 'A' # Adds A, B, C, D labels to the subplots
       
    theme(plot.title = element_text(size = 14, face = "bold")) # Apply to all subplots
  }
master_comparison_plots


#ITT_MI_CJP,ITT_MI_CJS_t1,ITT_MI_CJS_t2,


# As Treated=crumble ----
library(crumble)
#ASSIDUITE_SOIN, ETAPE, MOTIF
#use crumble and mediation analysis
#df.full$ASSIDUITE_SOIN <- as.factor(df.full$ASSIDUITE_SOIN)
df.full %>% 
  mutate(ASSIDUITE_SOIN=as.factor(ASSIDUITE_SOIN)) %>%
  # select only the variables you want to summarize
  select(BRAS, ASSIDUITE_SOIN , all_of(starts_with(all_outcomes))) %>% 
  tbl_summary(
    by = BRAS,
    # Use 'continuous' type for ratios to get Mean (SD) or Median (IQR)
    type = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ "{mean} ({sd})",
    digits = all_continuous() ~ 2,
    missing = "always" # Optional: hides the 'Unknown' row if you have NAs
  ) %>%
  add_p() %>% # Optional: adds a p-value to compare ratios between groups
  bold_labels()

df.full_AsTreated <- df.full_imputed %>% complete(1) %>% 
  mutate(cens=df.full$cens,
         !!sym(outcome) := df.full[[outcome]],
         ASSIDUITE_SOIN_crit = df.full[["ASSIDUITE_SOIN_crit"]],
         ASSIDUITE_SOIN = df.full[["ASSIDUITE_SOIN"]],
         ASSID_gp = case_when(ASSIDUITE_SOIN >= 7 & ASSIDUITE_SOIN <=10 ~ "7_10",
                              TRUE ~ as.factor(ASSIDUITE_SOIN))
  )
df.full_AsTreated$ASSID_gp <- as.factor(df.full_AsTreated$ASSID_gp)
table(df.full_AsTreated %>% filter(cens==1) %>% select(ASSID_gp) %>% pull)
set.seed(123)
tmle_control <- lmtp::lmtp_tmle(data = df.full_AsTreated,
                                trt = "BRAS_ASSID",
                                outcome = outcome, 
                                baseline = baseline,
                                outcome_type = "continuous",
                                cens="cens",
                                shift = static_binary_off, 
                                folds = nfolds,
                                learners_trt=learners_trt,
                                learners_outcome=learners_outcome
)
set.seed(123)
tmle_trt <- lmtp::lmtp_tmle(data = df.full_cens,
                            trt = "BRAS",
                            outcome = outcome, 
                            baseline = baseline,
                            outcome_type = "continuous",
                            cens="cens", 
                            shift = static_binary_on,
                            folds=nfolds,
                            learners_trt=learners_trt,
                            learners_outcome=learners_outcome)
lmtp_contrast(tmle_trt, ref=tmle_control)


# # boxplot(01036KR avec 0.5 dans RCS, M1 arrêt maladie donc absence ++)
# library(ggpubr)
# ggboxplot(
#   data = df.full.pp,
#   x = "BRAS",
#   y = "RATIO_H_T",
#   color = "BRAS",
#   fill = "BRAS",
#   palette = c("#00AFBB", "#FC4E07"),
#   legend = "none",
#   add = "jitter",
#   alpha = 0.3
# ) +
#   stat_compare_means(
#     method = "wilcox.test",
#     label.y = max(df.full.pp$RATIO_H_T, na.rm = T)*1.05    
#   ) +
#   labs(
#     title = "Ratio du nombre d'heures effectuées sur le nombre d'heures prévues",
#     x = NULL,
#     y = "RATIO_H_T"
#   ) -> res.plot
# 
# ggsave(
#   filename = "comp_ratio_PerProtocol.png",
#   plot = res.plot,
#   width = 8, height = 8, dpi = 300
# )

# 27112025 : Analayses complémentaire Nicolas
# df.full %>% 
#   tbl_summary(include = RATIO_H_T, by = BRAS) 
# 
# df.full %>% 
#   tbl_strata(
#     strata = ACSO_H_t1_cat,
#     .tbl_fun =
#       ~ .x |>
#       tbl_summary(by = "BRAS",include = "RATIO_H_T",  missing = "no", type = RATIO_H_T ~ "continuous") %>% 
#       add_n(),
#     .header = "**{strata}**, N = {n}"
#   ) 