
kneedle_method = function(x, y, pval, threshold = 1, alpha = 0.05, n_try = 3) {
  spline_fit = smooth.spline(x, y)
  y_smooth = predict(spline_fit)$y
  
  x_norm = (x - min(x)) / (max(x) - min(x))
  y_norm = (y_smooth - min(y_smooth)) / (max(y_smooth) - min(y_smooth))
  
  rotated_diff = (y_norm - x_norm) / sqrt(2)  # Rotate differences by 45 degrees

  elbow_index_top = order(rotated_diff, decreasing = TRUE)
  elbow_index_top = elbow_index_top[1:min(n_try, length(elbow_index_top))]

  elbow_index = NA
  for (i in elbow_index_top){
    if ((i >= threshold) & (pval[i] <= alpha)) {
      elbow_index = i
      break
    }
  }
  
  return(elbow_index)
}

