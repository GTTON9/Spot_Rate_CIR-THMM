plot_single_viterbi <- function(V_simulated, nstates, Gamma, kappa, theta, sigma, Reg_chain, subtitle_text) {
  
  
  prob <- forward_filter(V_simulated, nstates, Gamma, kappa, theta, sigma)
  states_estimate <- prob$states_estimate
  a <- prob$posterior
  
  time_index <- 1:ncol(a)
  
  prob_data <- data.frame(
    Time = time_index,
    Prob_State1 = a[1,],
    Prob_State2 = 1 - a[1,]
  )
  
  
  state_data <- data.frame(
    Time = time_index,
    Estimated_State = states_estimate - 1, 
    True_State = as.factor(Reg_chain[1:ncol(a)])
  )
  
  prob_data_long <- pivot_longer(
    prob_data,
    cols = starts_with("Prob_"),
    names_to = "Probability_Type",
    values_to = "Probability"
  )
  
  max_time <- max(state_data$Time)
  regime_start_times <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>%
    rename(Time_Start = Time)
  
  
  accuracy <- mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State))) * 100
  
  plot <- ggplot() +
    geom_rect(
      data = regime_start_times,
      aes(xmin = Time_Start, xmax = Time_End, ymin = -0.05, ymax = 1.05, fill = True_State),
      alpha = 0.2, 
      inherit.aes = FALSE
    ) +
    geom_line(
      data = prob_data_long, 
      aes(x = Time, y = Probability, color = Probability_Type),
      linewidth = 0.9
    ) +
    geom_point(
      data = state_data,
      aes(x = Time, y = Estimated_State, shape = "Estimated State"), 
      size = 1,
      color = "black",
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"), 
      labels = c("0" = "Regime 1\n(Calm)", "1" = "Regime 2\n(Turbulent)"),
      name = "True Regime"
    ) +
    scale_color_manual(
      values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
      labels = c("Prob_State1" = expression(paste(italic(P), "(", bold(V)[paste("1:", italic(t))], ", ", italic(R)[italic(t)], " = 1 | ", bold(theta), ")")),
                 "Prob_State2" = expression(paste(italic(P), "(", bold(V)[paste("1:", italic(t))], ", ", italic(R)[italic(t)], " = 2 | ", bold(theta), ")"))),
      name = "Filtered\nProbability"
    ) +
    scale_shape_manual(
      values = c("Estimated State" = 1), 
      labels = c("Estimated State" = ""),
      name = "Regime Estimate"
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      name = "Probability / State"
    ) +
    labs(
      title = sprintf("Accuracy: %.2f%%", accuracy),
      subtitle = subtitle_text,
      x = "Time Step"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # 加大标题
      plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
      
      legend.position = "bottom",         
      legend.box = "vertical",            
      legend.margin = margin(t = 10),     
      legend.spacing.y = unit(0.2, "cm"), 
      
      legend.text = element_text(size = 12),   
      legend.title = element_text(size = 11, face = "bold"),
      legend.key.size = unit(1.2, "lines")      
    ) +
    guides(
      fill = guide_legend(order = 1, ncol = 2),  
      color = guide_legend(order = 2, ncol = 1),  
      shape = guide_legend(order = 3, ncol = 1)
    )
  
  return(plot)
}





baum_welch_single_analysis <- function(
    X_series,
    Reg_chain,
    date_sequence,
    nstates = 2,
    runs = 10,
    dt = 1/252,
    interval_step = 10,
    subtitle_text = "X series"
) {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  cat("=================================================================\n")
  cat("Baum-Welch Single Analysis\n")
  cat("=================================================================\n\n")
  
  # ===== 1) Fit BW/HMM =====
  my_data_df <- data.frame(Date = date_sequence, Var = X_series)
  
  series_control <- Heston_set_controls(
    states = nstates, sdds = "Heston", date_column = "Date",
    file = my_data_df, data_column = "Var", logreturns = FALSE,
    from = date_sequence[1], to = date_sequence[length(date_sequence)],
    runs = runs
  )
  
  data_hmm   <- prepare_data(series_control)
  model_hmm  <- Heston_fit_model(data_hmm)
  final_model <- decode_states_heston(model_hmm)
  
  param_hat <- parUncon2par_heston(
    final_model$estimate, series_control,
    FALSE, numerical_safeguard = TRUE
  )
  
  states_estimate <- final_model$decoding
  
  cat("  BW fit done.\n\n")
  
  # ===== 2) Single Viterbi/Filtered probability plot =====
  plot_single_viterbi <- function(V_simulated, param, Reg_chain, subtitle_text) {
    
    prob <- forward_filter(V_simulated, nstates, param$Gamma,
                           param$kappa, param$theta, param$sigma)
    states_est <- prob$states_estimate
    a <- prob$posterior
    
    time_index <- 1:ncol(a)
    
    prob_data <- data.frame(
      Time = time_index,
      Prob_State1 = a[1,],
      Prob_State2 = 1 - a[1,]
    )
    
    state_data <- data.frame(
      Time = time_index,
      Estimated_State = states_est - 1,
      True_State = as.factor(Reg_chain[1:ncol(a)])
    )
    
    prob_data_long <- tidyr::pivot_longer(
      prob_data,
      cols = starts_with("Prob_"),
      names_to = "Probability_Type",
      values_to = "Probability"
    )
    
    max_time <- max(state_data$Time)
    regime_start_times <- state_data %>%
      mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
      filter(change | row_number() == 1) %>%
      dplyr::select(Time, True_State) %>%
      mutate(Time_End = lead(Time, default = max_time + 1)) %>%
      rename(Time_Start = Time)
    
    accuracy <- mean(state_data$Estimated_State ==
                       as.numeric(as.character(state_data$True_State))) * 100
    
    ggplot() +
      geom_rect(
        data = regime_start_times,
        aes(xmin = Time_Start, xmax = Time_End, ymin = -0.05, ymax = 1.05, fill = True_State),
        alpha = 0.2, inherit.aes = FALSE
      ) +
      geom_line(
        data = prob_data_long,
        aes(x = Time, y = Probability, color = Probability_Type),
        linewidth = 0.9
      ) +
      geom_point(
        data = state_data,
        aes(x = Time, y = Estimated_State, shape = "Estimated State"),
        size = 1, color = "black", alpha = 0.8
      ) +
      scale_fill_manual(
        values = c("0" = "skyblue", "1" = "salmon"),
        labels = c("0" = "Regime 1\n(Calm)", "1" = "Regime 2\n(Turbulent)"),
        name = "True Regime"
      ) +
      scale_color_manual(
        values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
        labels = c("Prob_State1" = "P(State=1|data)",
                   "Prob_State2" = "P(State=2|data)"),
        name = "Filtered\nProbability"
      ) +
      scale_shape_manual(values = c("Estimated State" = 1), name = "Regime Estimate") +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      labs(
        title = sprintf("%s | Accuracy: %.2f%%", subtitle_text, accuracy),
        x = "Time Step",
        y = "Probability / State"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        legend.box = "vertical"
      )
  }
  
  viterbi_plot <- plot_single_viterbi(
    V_simulated = X_series,
    param = param_hat,
    Reg_chain = Reg_chain,
    subtitle_text = subtitle_text
  )
  
  # ===== 3) Confidence interval plot =====
  ci_result <- plot_cir_confidence_improved(
    true_vol = X_series,
    param = param_hat,
    states_estimate = states_estimate,
    Reg_chain = Reg_chain,
    dt = dt,
    interval_step = interval_step,
    subtitle_text = subtitle_text
  )
  
  cat("=================================================================\n")
  cat("Analysis Complete!\n")
  cat("=================================================================\n")
  
  return(list(
    params = param_hat,
    states_estimate = states_estimate,
    viterbi_plot = viterbi_plot,
    ci_plot = ci_result$plot,
    ci_data = ci_result$ci_data
  ))
}