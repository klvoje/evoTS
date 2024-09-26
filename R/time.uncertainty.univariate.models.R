#' @title Assessing the effects of uncertainty in sample times on relative model fit and estimated parameters.
#'
#' @description This function investigates whether and how uncertainty in sample times affects model parameters. It also allows for the comparison of relative model fit between two or more models given uncertainty in sample times.
#'
#' @param y an univariate evoTS object.
#' 
#' @param model one or more models to be investigated, "Stasis", URW", "GRW", "OU.fixed", "OUBM", "decel", "accel".
#'
#' @param mode to model uncertainty in time, the default method is "uniform", which draws a random number between time.min and time.max. Alternatively, if "normal" is selected, uncertainty at each time point is drawn from a normal distribution with the original time point as the mean and sd.trait as the standard deviation of the distribution.
#'
#' @param time.min the vector defining the lower limits of the uniform distribution when mode = "uniform" (default). The default value is NULL, in which case the lower limit is set to 10% of the time vector in the paleoTS object.
#' 
#' @param time.max the vector defining the upper limits of the uniform distribution when mode = "uniform" (default). The default value is NULL, in which case the upper limit is set to 10% of the time vector in the paleoTS object.
#' 
#' @param sd.time value defining the standard deviation of the normal distribution of uncertainty in time when mode = "normal". The default is NULL, in which case the value is set to 0.05.
#' 
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param iter the number of times the optimization method is run with different time vectors. The default is 100.
#'
#' @details If only one model is defined using the model argument, the output includes only the mean parameter estimates and their standard deviations based on the iterated model fitting. If more than one model is defined (e.g., model = c("URW", "GRW")), the output also shows the number of times each model had the lowest AICc. Note that the Stasis model is time-independent, meaning that different sample times do not affect the estimated model parameters or the AICc. However, the Stasis model can still be included in the model argument, allowing for comparisons of its relative fit to the other models.
#'
#' @return The output is organized as a list. The component mean.parameter.values contains the average parameter values from all iterations, while sd.parameters includes the standard deviations for each model parameter. The component n.best.AICc indicates the number of times the model achieved the lowest (best) fit according to AICc.  
#'
#' @author Kjetil Lysne Voje
#'
#'@export
#'
#'@examples
#'
#'## Generate a paleoTS object by simulating a time series
#'y<-paleoTS:::sim.GRW(15,0.04,0.08)
#'
#'##Evaluate the effect of uncertainty in the timing of samples.
#'time.uncertainty.univariate.models(y, model=c("URW", "GRW"))
#'

time.uncertainty.univariate.models<-function(y, mode = "uniform", time.min = NULL, time.max = NULL, sd.time = NULL, model = c("URW","GRW"), pool = TRUE, iter=100){
  
  if (is.null(time.min) == TRUE) time.min<-y$tt*0.9
  if (is.null(time.max) == TRUE) time.max<-y$tt*1.1
  
  if (is.null(sd.time) == TRUE) sd.time<-0.05
  
  # Create lists and assign names
  result_list <- lapply(model, function(name) {
    vector(mode = 'list', length = iter)
  })
  
  # Assign model names to each list
  names(result_list) <- paste0("result.list_", model)
  
  
  # Creating a new time vector based on assumtions about the error
    for (i in 1:iter){
      new.time<-rep(NA, length(y$tt))
      
      # Assuming uniform error
      if (mode == "uniform") {
      for (j in 1:length(y$tt))
        {
        new.time[j]<-runif(1, min = time.min[j], max = time.max[j])
        }  
      }
      
      # Assuming normal error
      if (mode == "normal") {
        for (j in 1:length(y$tt))
        {
          new.time[j]<-rnorm(1, mean = y$mm[j], sd = sd.time)
        }  
      }

      # Check if samples are ordered correctly given the new time vector
      if (all(diff(new.time) < 0) == TRUE) {
        tmp_data_matrix<-cbind(y$mm,y$vv,y$nn,new.time)
        tmp_sorted_matrix <- tmp_data_matrix[order(tmp_data_matrix[, 4]),]
        y_new_tt<-paleoTS:::as.paleoTS(mm = tmp_sorted_matrix[,1], vv = tmp_sorted_matrix[,2], nn = tmp_sorted_matrix[,3], tt = tmp_sorted_matrix[,4])
      } 
      
      if (all(diff(new.time) < 0) == FALSE)  {y_new_tt<-paleoTS:::as.paleoTS(mm = y$mm, vv = y$vv, nn = y$nn, tt = new.time)}
      
      # Fitting models based on new time vector
      if ("Stasis" %in% model) {result_list$result.list_Stasis[[i]]<-paleoTS:::opt.joint.Stasis(y_new_tt, pool = pool)
      }
      if ("URW" %in% model) {result_list$result.list_URW[[i]]<-paleoTS:::opt.joint.URW(y_new_tt, pool = pool)
      }
      if ("GRW" %in% model) {result.list$result.list_GRW[[i]]<-paleoTS:::opt.joint.GRW(y_new_tt, pool = pool)
      }
      if ("OU.fixed" %in% model) {result.list$result.list_OU.fixed[[i]]<-paleoTS:::opt.joint.OU(y_new_tt, pool = pool)
      }
      if ("OUBM" %in% model) {result.list$result.list_OUBM[[i]]<-opt.joint.OUBM(y_new_tt, pool = pool)
      }
      if ("accel" %in% model) {result.list$result.list_accel[[i]]<-opt.joint.accel(y_new_tt, pool = pool)
      }
      if ("decel" %in% model) {result.list$result.list_decel[[i]]<-opt.joint.decel(y_new_tt, pool = pool)
      }
      
    }

  #Structuring output in lists
  if ("Stasis" %in% model) {Stasis_aicc_values <- sapply(result_list$result.list_Stasis, function(x) {x$AICc})} else {Stasis_aicc_values<-NULL; } 
  if ("Stasis" %in% model) {Stasis_logl_values <- sapply(result_list$result.list_Stasis, function(x) {x$logL});
  Stasis_par_values <- sapply(result_list$result.list_Stasis, function(x) {x$parameters});
  Stasis.mean.par<-apply(Stasis_par_values,1,mean); Stasis.SD.par<-apply(Stasis_par_values,1,sd) 
  }
  
  if ("URW" %in% model) {URW_aicc_values <- sapply(result_list$result.list_URW, function(x) {x$AICc})} else {URW_aicc_values<-NULL; } 
  if ("URW" %in% model) {URW_logl_values <- sapply(result_list$result.list_URW, function(x) {x$logL});
  URW_par_values <- sapply(result_list$result.list_URW, function(x) {x$parameters});
  URW.mean.par<-apply(URW_par_values,1,mean); URW.SD.par<-apply(URW_par_values,1,sd) 
  }
  
  if ("GRW" %in% model) {GRW_aicc_values <- sapply(result.list$result.list_GRW, function(x) {x$AICc})} else {GRW_aicc_values<-NULL; } 
  if ("GRW" %in% model) {GRW_logl_values <- sapply(result.list$result.list_GRW, function(x) {x$logL});
  GRW_par_values <- sapply(result.list$result.list_GRW, function(x) {x$parameters});
  GRW.mean.par<-apply(GRW_par_values,1,mean); GRW.SD.par<-apply(GRW_par_values,1,sd)
  }
  
  if ("OU.fixed" %in% model) {OU.fixed_aicc_values <- sapply(result.list$result.list_OU.fixed, function(x) {x$AICc})} else {OU.fixed_aicc_values<-NULL; }
  if ("OU.fixed" %in% model) {OU.fixed_logl_values <- sapply(result.list$result.list_OU.fixed, function(x) {x$logL});
  OU.fixed_par_values <- sapply(result.list$result.list_OU.fixed, function(x) {x$parameters});
  OU.fixed.mean.par<-apply(OU.fixed_par_values,1,mean); OU.fixed.SD.par<-apply(OU.fixed_par_values,1,sd)
  }
  
  if ("OUBM" %in% model) {OUBM_aicc_values <- sapply(result.list$result.list_OUBM, function(x) {x$AICc})} else {OUBM_aicc_values<-NULL; }
  if ("OUBM" %in% model) {OUBM_logl_values <- sapply(result.list$result.list_OUBM, function(x) {x$logL});
  OUBM_par_values <- sapply(result.list$result.list_OUBM, function(x) {x$parameters});
  OUBM.mean.par<-apply(OUBM_par_values,1,mean); OUBM.SD.par<-apply(OUBM_par_values,1,sd)
  }
  
  if ("accel" %in% model) {accel_aicc_values <- sapply(result.list$result.list_accel, function(x) {x$AICc})} else {accel_aicc_values<-NULL; }
  if ("accel" %in% model) {accel_logl_values <- sapply(result.list$result.list_accel, function(x) {x$logL});
  accel_par_values <- sapply(result.list$result.list_accel, function(x) {x$parameters});
  accel.mean.par<-apply(accel_par_values,1,mean); accel.SD.par<-apply(accel_par_values,1,sd)
  }
  
  if ("decel" %in% model) {decel_aicc_values <- sapply(result.list$result.list_decel, function(x) {x$AICc})} else {decel_aicc_values<-NULL; }
  if ("decel" %in% model) {decel_logl_values <- sapply(result.list$result.list_decel, function(x) {x$logL});
  decel_par_values <- sapply(result.list$result.list_decel, function(x) {x$parameters});
  decel.mean.par<-apply(decel_par_values,1,mean); decel.SD.par<-apply(decel_par_values,1,sd)
  }
  
  # Combine vectors into a list
  vectors <- list(Stasis_aicc_values, URW_aicc_values, GRW_aicc_values, OU.fixed_aicc_values, OUBM_aicc_values, decel_aicc_values, accel_aicc_values)
  
  # Filter out NULL values
  non_null_vectors <- Filter(Negate(is.null), vectors)
  
  if (length(model) >1) {
  # Combine non-NULL vectors into a matrix
  comparison <- do.call(rbind, non_null_vectors)

  # Find the minimum AICc value for each model run 
  min_values <- apply(comparison, 2, min)
  
  # Create a matrix to count how many times the smallest AICc appears at each position
  counts_matrix <- matrix(0, nrow = length(model), ncol = length(min_values))
  
  # Count occurrences of the minimum AICc values
  for (i in 1:length(min_values)) {
    counts_matrix[, i] <- as.integer(comparison[, i] == min_values[i])
  }
  }
  
  # Preparing output
  w.list <- vector(mode='list', length(model))
  names(w.list)<-model
 i<-0
  if ("Stasis" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(Stasis.mean.par)
    w.list[[i]]$sd.parameters<-(Stasis.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
   
  if ("URW" %in% model) {
    i<-i+1
  w.list[[i]]$mean.parameter.values<-(URW.mean.par)
  w.list[[i]]$sd.parameters<-(URW.SD.par)
  if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("GRW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(GRW.mean.par)
    w.list[[i]]$sd.parameters<-(GRW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("OU.fixed" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(OU.fixed.mean.par)
    w.list[[i]]$sd.parameters<-(OU.fixed.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("OUBM" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(OUBM.mean.par)
    w.list[[i]]$sd.parameters<-(OUBM.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("decel" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(decel.mean.par)
    w.list[[i]]$sd.parameters<-(decel.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("accel" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(accel.mean.par)
    w.list[[i]]$sd.parameters<-(accel.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  return(w.list)

}
