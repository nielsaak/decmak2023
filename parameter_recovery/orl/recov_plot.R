# More informative plotting - code-courtesy of Lasse - fanx!
recov_plot <- function(true, infer, plot_lab, plot_col, lim = FALSE, lim_val = 1) {
  
  # library(ggplot2)
  
  if (lim == TRUE) {  
    df <- data.frame(true, infer)
    
    pl <- ggplot(df, aes(x = true,
                         y = infer,
                         color = plot_col)) + #Setting aesthetics for plot
      geom_point() + #Giving points a color each
      #scale_linetype(name="perfect recovery") +
      geom_abline(intercept=0, slope=1, linetype=2) +
      geom_smooth(method = "lm", se = T, formula = "y ~ x") +
      theme_minimal() + #Setting theme
      xlab(plot_lab[1]) + #Setting x label
      ylab(plot_lab[2]) + #Setting y label
      labs(color = "") + #Setting legend title
      ggtitle(paste0("'", plot_lab[2], "' compared to '", plot_lab[1], "'")) +
      theme(legend.position="none") +
      coord_cartesian(xlim = c(0, lim_val), ylim = c(0,lim_val)) 
  }
  else {
    df <- data.frame(true, infer)
    
    pl <- ggplot(df, aes(x = true,
                         y = infer,
                         color = plot_col)) + #Setting aesthetics for plot
      geom_point() + #Giving points a color each
      #scale_linetype(name="perfect recovery") +
      geom_abline(intercept=0, slope=1, linetype=2) +
      geom_smooth(method = "lm", se = T, formula = "y ~ x") +
      theme_minimal() + #Setting theme
      xlab(plot_lab[1]) + #Setting x label
      ylab(plot_lab[2]) + #Setting y label
      labs(color = "") + #Setting legend title
      ggtitle(paste0("'", plot_lab[2], "' compared to '", plot_lab[1], "'")) +
      theme(legend.position="none")
  }
  
  return(pl)
  
}