inv_trans_theta=function(data){
  theta_inv=data$theta
  theta_inv[which(theta_inv<0)]=theta_inv[which(theta_inv<0)]+2*pi
  data$theta <- theta_inv
  return(data)
}
trans_theta=function(data){
  theta_trans <- data$theta
  theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  data$theta <- theta_trans
  return(data)
}

# max_min_feat=function(x,tt_thres_maxmin=2.5,l=5){
#   library(dplyr)
#   maxs <- as.numeric(splus2R::peaks(x$theta, span = l))
#   mins <- as.numeric(splus2R::peaks(-x$theta, span = l))
# 
#   # Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
#   maxs[which(x$theta>(tt_thres_maxmin))] <- 0
#   mins[which(x$theta<(-tt_thres_maxmin))] <- 0
# 
#   #
#   x$maxs=maxs
#   x$mins=mins
# 
#   find_closest <- function(t_val, t_events, theta_events) {
#     if (length(t_events) == 0) return(NA)  # No events found
#     diffs <- abs(t_events - t_val)  # Compute absolute time differences
#     closest_index <- which.min(diffs)  # Find index of closest event
#     return(theta_events[closest_index])  # Return corresponding theta value
#   }
# 
#   # Extract time and values of local maxima and minima
#   t_maxs <- x$t[x$maxs == 1]
#   theta_maxs <- x$theta[x$maxs == 1]
# 
#   t_mins <- x$t[x$mins == 1]
#   theta_mins <- x$theta[x$mins == 1]
# 
#   # Compute the closest local max and min for each row
#   x <- x %>%
#     mutate(
#       value_max = sapply(t, find_closest, t_maxs, theta_maxs),
#       value_min = sapply(t, find_closest, t_mins, theta_mins)
#     )
# 
#   # # TP indicator
#   # x$I_TP=as.numeric(x$value_max*x$value_min>0 & x$value_max>x$value_min)
#   #
#   # # HS indicator
#   # x$I_HS=as.numeric(x$value_max*x$value_min<0&x$value_max<x$value_min)
#   #
#   # # QS indicator
#   # x$I_QS=as.numeric(x$value_max*x$value_min<0 & x$value_max>x$value_min)
#   #
#   # # CP indicator
#   # x$I_CP=as.numeric(x$value_max*x$value_min>0&x$value_max<=x$value_min)
#   #
#   #
#   # x$mean_osc=x$I_HS*(x$value_min+x$value_max+2*pi)/2+
#   #   x$I_QS*(x$value_min+x$value_max)/2+
#   #   x$I_TP*(x$value_min+x$value_max)/2
# 
#   x=subset(x,select=-c(maxs,mins))
# 
#   return(x)
# }

max_min_feat=function(x,var_name="theta",tt_thres_maxmin=2.5,l=5){
  library(dplyr)
  y=x[[var_name]]
  maxs <- as.numeric(splus2R::peaks(y, span = l))
  mins <- as.numeric(splus2R::peaks(-y, span = l))
  
  # Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
  maxs[which(y>(tt_thres_maxmin))] <- 0
  mins[which(y<(-tt_thres_maxmin))] <- 0
  
  #
  x$maxs=maxs
  x$mins=mins
  
  find_closest <- function(t_val, t_events, theta_events) {
    if (length(t_events) == 0) return(NA)  # No events found
    diffs <- abs(t_events - t_val)  # Compute absolute time differences
    closest_index <- which.min(diffs)  # Find index of closest event
    return(theta_events[closest_index])  # Return corresponding theta value
  }
  
  # Extract time and values of local maxima and minima
  t_maxs <- x$t[x$maxs == 1]
  theta_maxs <- y[x$maxs == 1]
  
  t_mins <- x$t[x$mins == 1]
  theta_mins <- y[x$mins == 1]
  
  # Compute the closest local max and min for each row
  x <- x %>%
    mutate(
      value_max = sapply(t, find_closest, t_maxs, theta_maxs),
      value_min = sapply(t, find_closest, t_mins, theta_mins)
    )
  
  if(var_name=="theta"){
    # TP indicator
      x$I_TP=as.numeric(x$value_max*x$value_min>0 & x$value_max>x$value_min)

      # # HS indicator
      x$I_HS=as.numeric(x$value_max*x$value_min<0&x$value_max<x$value_min)

      # # QS indicator
      x$I_QS=as.numeric(x$value_max*x$value_min<0 & x$value_max>x$value_min)

      # # CP indicator
      x$I_CP=as.numeric(x$value_max*x$value_min>0&x$value_max<=x$value_min)
  }
  
  old <- c("value_max", "value_min")
  new <- paste0(old, "_", var_name)
  names(x)[names(x) %in% old] <- new
  
  matched_cols <- grepl(var_name, names(x))
  matched_cols[1]=T
  matched_cols[which(names(x)==var_name)]=F
  
  if(var_name=="theta"){
    matched_cols[(length(matched_cols)-3):length(matched_cols)] <- TRUE
  }
  
  return(x[ , matched_cols])
}



# compute_feat=function(x,l,ma=T,tt_thres_maxmin=2.5){
#   
#   # First diff
#   x$dtheta=c(NA,diff(x$theta))
#   
#   # Moving sd of theta
#   x$sd_theta=zoo::rollapply(x$theta, l, sd, fill=NA)
#   x$sd_dtheta=zoo::rollapply(x$dtheta, l, sd, fill=NA)
#   
#   # Moving av and sd of a
#   
#   if(ma){
#     x$ma_theta=zoo::rollapply(x$theta, l, mean, fill=NA)
#   }
#   
#   return(x)
#   
# }

compute_feat <- function(x,
                         var_name,
                         l,
                         ma_flag   = TRUE,
                         sd_flag=TRUE) {
  library(zoo)
  # pull out the vector
  v <- x[[var_name]]
  
  # 1) first difference
  d  <- c(NA, diff(v))
  x[[paste0("d", var_name)]] <- d
  
  # 2) rolling SD of the original and the diff
  if(sd_flag){
    x[[paste0("sd_",   var_name)]] <- rollapply(v, l, sd, fill = NA)
    x[[paste0("sd_d",  var_name)]] <- rollapply(d, l, sd, fill = NA)
  }
  
  # 3) optional moving average of the original
  if (ma_flag) {
    x[[paste0("ma_",  var_name)]] <- rollapply(v, l, mean, fill = NA)
  }
  
  matched_cols <- grepl(var_name, names(x))
  matched_cols[1]=T
  matched_cols[which(names(x)==var_name)]=F
  
  return(x[ , matched_cols])
}
