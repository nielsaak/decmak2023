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
    stat_cor(method = "pearson", aes(label = ..r.label..), label) +
    theme_minimal() + #Setting theme
    xlab(plot_lab_1) + #Setting x label
    ylab(plot_lab_2) + #Setting y label
    labs(color = "") + #Setting legend title
    ggtitle(title) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))
  
  return(pl)
  
}

# # Load ggplot2
# library(ggplot2)
# 
# # Your data
# data <- data.frame(x = ..., y = ...)
# 
# # Linear regression
# model <- lm(y ~ x, data = data)
# 
# # Get R-squared value
# r_squared <- summary(model)$r.squared
# 
# # Create the plot
# ggplot(data, aes(x = x, y = y)) +
#   geom_point() +
#   geom_abline(slope = coef(model)[2], intercept = coef(model)[1]) +
#   annotate("text", x = your_x, y = your_y, label = paste("R^2 =", round(r_squared, 2)), parse = TRUE)
