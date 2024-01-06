# More informative plotting - code-courtesy of Lasse - fanx!
recov_plot <- function(true, infer, plot_lab_1, plot_lab_2, plot_col, title) {
  
  # library(ggplot2)
  
  df <- data.frame(true, infer)
  
  pl <- ggplot(df, aes(x = true,
                       y = infer,
                       color = plot_col)) + #Setting aesthetics for plot
    geom_point() + #Giving points a color each
    #scale_linetype(name="perfect recovery") +
    geom_abline(intercept=0, slope=1, linetype=2) +
    geom_smooth(method = "lm", se = T, formula = "y ~ x") +
    theme_minimal() + #Setting theme
    xlab(plot_lab_1) + #Setting x label
    ylab(plot_lab_2) + #Setting y label
    labs(color = "") + #Setting legend title
    ggtitle(title) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))
  
  return(pl)
  
}