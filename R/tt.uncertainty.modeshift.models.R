#' @title Assessing the effects of uncertainty in sample times on relative model fit and estimated parameters in mode shift models.
#'
#' @description This function investigates whether and how uncertainty in sample times affects model parameters in mode shift models. It also allows for the comparison of relative model fit between two or more models given uncertainty in sample times.
#'
#' @param y an univariate evoTS object.
#' 
#' @param model one or more mode shift models to be investigated.
#' 
#' @param shift.point The sample that split the time-series into two segments. 
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
#'y<-paleoTS:::sim.GRW(20,0.4,0.08)
#'
#'##Evaluate the effect of uncertainty in the timing of samples.
#'tt.uncertainty.modeshift.models(y, model=c("URW.Stasis", "URW.URW"), shift.point = 10, iter = 10)
#'

tt.uncertainty.modeshift.models<-function(y, mode = "uniform", time.min = NULL, time.max = NULL, sd.time = NULL, model = c("Stasis.Stasis", "URW.URW"), shift.point = NULL, pool = TRUE, iter=20){
  
  if (is.null(shift.point) == TRUE) print("You need to define a shift-point")
  
  if (is.null(time.min) == TRUE) time.min<-y$tt*0.9
  if (is.null(time.max) == TRUE) time.max<-y$tt*1.1
  
  if (is.null(sd.time) == TRUE) sd.time<-0.05
  
  # Create lists and assign names
  result.list <- lapply(model, function(name) {
    vector(mode = 'list', length = iter)
  })
  
  # Assign model names to each list
  names(result.list) <- paste0("result.list_", model)
  
  
  # Creating a new time vector based on assumptions about the error
  for (i in 1:iter){
    new.time<-rep(NA, length(y$tt))
    
    # Assuming uniform error
    if (mode == "uniform") {
      for (j in 1:length(y$tt))
      {
        new.time[j]<-stats::runif(1, min = time.min[j], max = time.max[j])
      }  
    }
    
    # Assuming normal error
    if (mode == "normal") {
      for (j in 1:length(y$tt))
      {
        new.time[j]<-stats::rnorm(1, mean = y$mm[j], sd = sd.time)
      }  
    }
    
    # Check if samples are ordered correctly given the new time vector
    if (all(diff(new.time) < 0) == TRUE) {
      tmp_data_matrix<-cbind(y$mm,y$vv,y$nn,new.time)
      tmp_sorted_matrix <- tmp_data_matrix[order(tmp_data_matrix[, 4]),]
      y_new_tt<-paleoTS::as.paleoTS(mm = tmp_sorted_matrix[,1], vv = tmp_sorted_matrix[,2], nn = tmp_sorted_matrix[,3], tt = tmp_sorted_matrix[,4])
    } 
    
    if (all(diff(new.time) < 0) == FALSE)  {y_new_tt<-paleoTS::as.paleoTS(mm = y$mm, vv = y$vv, nn = y$nn, tt = new.time)}
    
    GG <-shift.point
    GG<-as.matrix(GG)
    nc <- ncol(GG)
    gg <- shift2gg(GG[, nc], ns <- length(y_new_tt$mm))
    
    # Fitting models based on new time vector
    if ("Stasis.Stasis" %in% model) {result.list$result.list_Stasis.Stasis[[i]]<-opt.joint.Stasis.Stasis(y_new_tt, gg, pool = pool) 
    }
    if ("Stasis.URW" %in% model) {result.list$result.list_Stasis.URW[[i]]<-opt.joint.Stasis.URW(y_new_tt, gg, pool = pool) 
    }
    if ("Stasis.GRW" %in% model) {result.list$result.list_Stasis.GRW[[i]]<-opt.joint.Stasis.GRW(y_new_tt, gg, pool = pool) 
    }
    if ("Stasis.OU" %in% model) {result.list$result.list_Stasis.OU[[i]]<-opt.joint.Stasis.OU(y_new_tt, gg, pool = pool) 
    }
    if ("URW.Stasis" %in% model) {result.list$result.list_URW.Stasis[[i]]<-opt.joint.URW.Stasis(y_new_tt, gg, pool = pool) 
    }
    if ("GRW.Stasis" %in% model) {result.list$result.list_GRW.Stasis[[i]]<-opt.joint.GRW.Stasis(y_new_tt, gg, pool = pool) 
    }
    if ("OU.Stasis" %in% model) {result.list$result.list_OU.Stasis[[i]]<-opt.joint.OU.Stasis(y_new_tt, gg, pool = pool) 
    }
    if ("URW.URW" %in% model) {result.list$result.list_URW.URW[[i]]<-opt.joint.URW.URW(y_new_tt, gg, pool = pool) 
    }
    if ("URW.GRW" %in% model) {result.list$result.list_URW.GRW[[i]]<-opt.joint.URW.GRW(y_new_tt, gg, pool = pool) 
    }
    if ("URW.OU" %in% model) {result.list$result.list_URW.OU[[i]]<-opt.joint.URW.OU(y_new_tt, gg, pool = pool) 
    }
    if ("OU.URW" %in% model) {result.list$result.list_OU.URW[[i]]<-opt.joint.OU.URW(y_new_tt, gg, pool = pool) 
    }
    if ("GRW.URW" %in% model) {result.list$result.list_GRW.URW[[i]]<-opt.joint.GRW.URW(y_new_tt, gg, pool = pool) 
    }
    if ("GRW.OU" %in% model) {result.list$result.list_GRW.OU[[i]]<-opt.joint.GRW.OU(y_new_tt, gg, pool = pool) 
    }
    if ("GRW.GRW" %in% model) {result.list$result.list_GRW.GRW[[i]]<-opt.joint.GRW.GRW(y_new_tt, gg, pool = pool) 
    }
    if ("OU.GRW" %in% model) {result.list$result.list_OU.GRW[[i]]<-opt.joint.OU.GRW(y_new_tt, gg, pool = pool) 
    }
    if ("OU.OU" %in% model) {result.list$result.list_OU.OU[[i]]<-opt.joint.OU.OU(y_new_tt, gg, pool = pool) 
    }
  }
  
  #Structuring output in lists
  if ("Stasis.Stasis" %in% model) {Stasis.Stasis_aicc_values <- sapply(result.list$result.list_Stasis.Stasis, function(x) {x$AICc})} else {Stasis.Stasis_aicc_values<-NULL; } 
  if ("Stasis.Stasis" %in% model) {Stasis.Stasis_logl_values <- sapply(result.list$result.list_Stasis.Stasis, function(x) {x$logL});
  Stasis.Stasis_par_values <- sapply(result.list$result.list_Stasis.Stasis, function(x) {x$parameters});
  if (is.matrix(Stasis.Stasis_par_values) == TRUE) { 
    Stasis.Stasis.mean.par<-apply(Stasis.Stasis_par_values, 1, function(x) mean(x, na.rm = TRUE))
    Stasis.Stasis.SD.par<-apply(Stasis.Stasis_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    Stasis.Stasis_par_values<-do.call(rbind, Stasis.Stasis_par_values)
    Stasis.Stasis.mean.par<-apply(Stasis.Stasis_par_values, 2, function(x) mean(x, na.rm = TRUE))
    Stasis.Stasis.SD.par<-apply(Stasis.Stasis_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("Stasis.URW" %in% model) {Stasis.URW_aicc_values <- sapply(result.list$result.list_Stasis.URW, function(x) {x$AICc})} else {Stasis.URW_aicc_values<-NULL; } 
  if ("Stasis.URW" %in% model) {Stasis.URW_logl_values <- sapply(result.list$result.list_Stasis.URW, function(x) {x$logL});
  Stasis.URW_par_values <- sapply(result.list$result.list_Stasis.URW, function(x) {x$parameters}, simplify =FALSE);
  if (is.matrix(Stasis.URW_par_values) == TRUE) { 
    Stasis.URW.mean.par<-apply(Stasis.URW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    Stasis.URW.SD.par<-apply(Stasis.URW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    Stasis.URW_par_values<-do.call(rbind, Stasis.URW_par_values)
    Stasis.URW.mean.par<-apply(Stasis.URW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    Stasis.URW.SD.par<-apply(Stasis.URW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("Stasis.GRW" %in% model) {Stasis.GRW_aicc_values <- sapply(result.list$result.list_Stasis.GRW, function(x) {x$AICc})} else {Stasis.GRW_aicc_values<-NULL; } 
  if ("Stasis.GRW" %in% model) {Stasis.GRW_logl_values <- sapply(result.list$result.list_Stasis.GRW, function(x) {x$logL});
  Stasis.GRW_par_values <- sapply(result.list$result.list_Stasis.GRW, function(x) {x$parameters});
  if (is.matrix(Stasis.GRW_par_values) == TRUE) { 
    Stasis.GRW.mean.par<-apply(Stasis.GRW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    Stasis.GRW.SD.par<-apply(Stasis.GRW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    Stasis.GRW_par_values<-do.call(rbind, Stasis.GRW_par_values)
    Stasis.GRW.mean.par<-apply(Stasis.GRW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    Stasis.GRW.SD.par<-apply(Stasis.GRW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("Stasis.OU" %in% model) {Stasis.OU_aicc_values <- sapply(result.list$result.list_Stasis.OU, function(x) {x$AICc})} else {Stasis.OU_aicc_values<-NULL; } 
  if ("Stasis.OU" %in% model) {Stasis.OU_logl_values <- sapply(result.list$result.list_Stasis.OU, function(x) {x$logL});
  Stasis.OU_par_values <- sapply(result.list$result.list_Stasis.OU, function(x) {x$parameters});
  if (is.matrix(Stasis.OU_par_values) == TRUE) { 
    Stasis.OU.mean.par<-apply(Stasis.OU_par_values, 1, function(x) mean(x, na.rm = TRUE))
    Stasis.OU.SD.par<-apply(Stasis.OU_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    Stasis.OU_par_values<-do.call(rbind, Stasis.OU_par_values)
    Stasis.OU.mean.par<-apply(Stasis.OU_par_values, 2, function(x) mean(x, na.rm = TRUE))
    Stasis.OU.SD.par<-apply(Stasis.OU_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("URW.Stasis" %in% model) {URW.Stasis_aicc_values <- sapply(result.list$result.list_URW.Stasis, function(x) {x$AICc})} else {URW.Stasis_aicc_values<-NULL; } 
  if ("URW.Stasis" %in% model) {URW.Stasis_logl_values <- sapply(result.list$result.list_URW.Stasis, function(x) {x$logL});
  URW.Stasis_par_values <- sapply(result.list$result.list_URW.Stasis, function(x) {x$parameters});
  if (is.matrix(URW.Stasis_par_values) == TRUE) { 
  URW.Stasis.mean.par<-apply(URW.Stasis_par_values, 1, function(x) mean(x, na.rm = TRUE))
  URW.Stasis.SD.par<-apply(URW.Stasis_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
  URW.Stasis_par_values<-do.call(rbind, URW.Stasis_par_values)
  URW.Stasis.mean.par<-apply(URW.Stasis_par_values, 2, function(x) mean(x, na.rm = TRUE))
  URW.Stasis.SD.par<-apply(URW.Stasis_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("GRW.Stasis" %in% model) {GRW.Stasis_aicc_values <- sapply(result.list$result.list_GRW.Stasis, function(x) {x$AICc})} else {GRW.Stasis_aicc_values<-NULL; } 
  if ("GRW.Stasis" %in% model) {GRW.Stasis_logl_values <- sapply(result.list$result.list_GRW.Stasis, function(x) {x$logL});
  GRW.Stasis_par_values <- sapply(result.list$result.list_GRW.Stasis, function(x) {x$parameters});
  if (is.matrix(GRW.Stasis_par_values) == TRUE) { 
    GRW.Stasis.mean.par<-apply(GRW.Stasis_par_values, 1, function(x) mean(x, na.rm = TRUE))
    GRW.Stasis.SD.par<-apply(GRW.Stasis_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    GRW.Stasis_par_values<-do.call(rbind, GRW.Stasis_par_values)
    GRW.Stasis.mean.par<-apply(GRW.Stasis_par_values, 2, function(x) mean(x, na.rm = TRUE))
    GRW.Stasis.SD.par<-apply(GRW.Stasis_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("OU.Stasis" %in% model) {OU.Stasis_aicc_values <- sapply(result.list$result.list_OU.Stasis, function(x) {x$AICc})} else {OU.Stasis_aicc_values<-NULL; } 
  if ("OU.Stasis" %in% model) {OU.Stasis_logl_values <- sapply(result.list$result.list_OU.Stasis, function(x) {x$logL});
  OU.Stasis_par_values <- sapply(result.list$result.list_OU.Stasis, function(x) {x$parameters});
  if (is.matrix(OU.Stasis_par_values) == TRUE) { 
    OU.Stasis.mean.par<-apply(OU.Stasis_par_values, 1, function(x) mean(x, na.rm = TRUE))
    OU.Stasis.SD.par<-apply(OU.Stasis_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    OU.Stasis_par_values<-do.call(rbind, OU.Stasis_par_values)
    OU.Stasis.mean.par<-apply(OU.Stasis_par_values, 2, function(x) mean(x, na.rm = TRUE))
    OU.Stasis.SD.par<-apply(OU.Stasis_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("URW.URW" %in% model) {URW.URW_aicc_values <- sapply(result.list$result.list_URW.URW, function(x) {x$AICc})} else {URW.URW_aicc_values<-NULL; } 
  if ("URW.URW" %in% model) {URW.URW_logl_values <- sapply(result.list$result.list_URW.URW, function(x) {x$logL});
  URW.URW_par_values <- sapply(result.list$result.list_URW.URW, function(x) {x$parameters});
  if (is.matrix(URW.URW_par_values) == TRUE) { 
    URW.URW.mean.par<-apply(URW.URW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    URW.URW.SD.par<-apply(URW.URW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    URW.URW_par_values<-do.call(rbind, URW.URW_par_values)
    URW.URW.mean.par<-apply(URW.URW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    URW.URW.SD.par<-apply(URW.URW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("URW.GRW" %in% model) {URW.GRW_aicc_values <- sapply(result.list$result.list_URW.GRW, function(x) {x$AICc})} else {URW.GRW_aicc_values<-NULL; } 
  if ("URW.GRW" %in% model) {URW.GRW_logl_values <- sapply(result.list$result.list_URW.GRW, function(x) {x$logL});
  URW.GRW_par_values <- sapply(result.list$result.list_URW.GRW, function(x) {x$parameters});
  if (is.matrix(URW.GRW_par_values) == TRUE) { 
    URW.GRW.mean.par<-apply(URW.GRW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    URW.GRW.SD.par<-apply(URW.GRW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    URW.GRW_par_values<-do.call(rbind, URW.GRW_par_values)
    URW.GRW.mean.par<-apply(URW.GRW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    URW.GRW.SD.par<-apply(URW.GRW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("URW.OU" %in% model) {URW.OU_aicc_values <- sapply(result.list$result.list_URW.OU, function(x) {x$AICc})} else {URW.OU_aicc_values<-NULL; } 
  if ("URW.OU" %in% model) {URW.OU_logl_values <- sapply(result.list$result.list_URW.OU, function(x) {x$logL});
  URW.OU_par_values <- sapply(result.list$result.list_URW.OU, function(x) {x$parameters});
  if (is.matrix(URW.OU_par_values) == TRUE) { 
    URW.OU.mean.par<-apply(URW.OU_par_values, 1, function(x) mean(x, na.rm = TRUE))
    URW.OU.SD.par<-apply(URW.OU_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    URW.OU_par_values<-do.call(rbind, URW.OU_par_values)
    URW.OU.mean.par<-apply(URW.OU_par_values, 2, function(x) mean(x, na.rm = TRUE))
    URW.OU.SD.par<-apply(URW.OU_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("OU.URW" %in% model) {OU.URW_aicc_values <- sapply(result.list$result.list_OU.URW, function(x) {x$AICc})} else {OU.URW_aicc_values<-NULL; } 
  if ("OU.URW" %in% model) {OU.URW_logl_values <- sapply(result.list$result.list_OU.URW, function(x) {x$logL});
  OU.URW_par_values <- sapply(result.list$result.list_OU.URW, function(x) {x$parameters});
  if (is.matrix(OU.URW_par_values) == TRUE) { 
    OU.URW.mean.par<-apply(OU.URW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    OU.URW.SD.par<-apply(OU.URW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    OU.URW_par_values<-do.call(rbind, OU.URW_par_values)
    OU.URW.mean.par<-apply(OU.URW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    OU.URW.SD.par<-apply(OU.URW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("GRW.URW" %in% model) {GRW.URW_aicc_values <- sapply(result.list$result.list_GRW.URW, function(x) {x$AICc})} else {GRW.URW_aicc_values<-NULL; } 
  if ("GRW.URW" %in% model) {GRW.URW_logl_values <- sapply(result.list$result.list_GRW.URW, function(x) {x$logL});
  GRW.URW_par_values <- sapply(result.list$result.list_GRW.URW, function(x) {x$parameters});
  if (is.matrix(GRW.URW_par_values) == TRUE) { 
    GRW.URW.mean.par<-apply(GRW.URW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    GRW.URW.SD.par<-apply(GRW.URW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    GRW.URW_par_values<-do.call(rbind, GRW.URW_par_values)
    GRW.URW.mean.par<-apply(GRW.URW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    GRW.URW.SD.par<-apply(GRW.URW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("GRW.OU" %in% model) {GRW.OU_aicc_values <- sapply(result.list$result.list_GRW.OU, function(x) {x$AICc})} else {GRW.OU_aicc_values<-NULL; } 
  if ("GRW.OU" %in% model) {GRW.OU_logl_values <- sapply(result.list$result.list_GRW.OU, function(x) {x$logL});
  GRW.OU_par_values <- sapply(result.list$result.list_GRW.OU, function(x) {x$parameters});
  if (is.matrix(GRW.OU_par_values) == TRUE) { 
    GRW.OU.mean.par<-apply(GRW.OU_par_values, 1, function(x) mean(x, na.rm = TRUE))
    GRW.OU.SD.par<-apply(GRW.OU_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    GRW.OU_par_values<-do.call(rbind, GRW.OU_par_values)
    GRW.OU.mean.par<-apply(GRW.OU_par_values, 2, function(x) mean(x, na.rm = TRUE))
    GRW.OU.SD.par<-apply(GRW.OU_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("GRW.GRW" %in% model) {GRW.GRW_aicc_values <- sapply(result.list$result.list_GRW.GRW, function(x) {x$AICc})} else {GRW.GRW_aicc_values<-NULL; } 
  if ("GRW.GRW" %in% model) {GRW.GRW_logl_values <- sapply(result.list$result.list_GRW.GRW, function(x) {x$logL});
  GRW.GRW_par_values <- sapply(result.list$result.list_GRW.GRW, function(x) {x$parameters});
  if (is.matrix(GRW.GRW_par_values) == TRUE) { 
    GRW.GRW.mean.par<-apply(GRW.GRW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    GRW.GRW.SD.par<-apply(GRW.GRW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    GRW.GRW_par_values<-do.call(rbind, GRW.GRW_par_values)
    GRW.GRW.mean.par<-apply(GRW.GRW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    GRW.GRW.SD.par<-apply(GRW.GRW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("OU.GRW" %in% model) {OU.GRW_aicc_values <- sapply(result.list$result.list_OU.GRW, function(x) {x$AICc})} else {OU.GRW_aicc_values<-NULL; } 
  if ("OU.GRW" %in% model) {OU.GRW_logl_values <- sapply(result.list$result.list_OU.GRW, function(x) {x$logL});
  OU.GRW_par_values <- sapply(result.list$result.list_OU.GRW, function(x) {x$parameters});
  if (is.matrix(OU.GRW_par_values) == TRUE) { 
    OU.GRW.mean.par<-apply(OU.GRW_par_values, 1, function(x) mean(x, na.rm = TRUE))
    OU.GRW.SD.par<-apply(OU.GRW_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    OU.GRW_par_values<-do.call(rbind, OU.GRW_par_values)
    OU.GRW.mean.par<-apply(OU.GRW_par_values, 2, function(x) mean(x, na.rm = TRUE))
    OU.GRW.SD.par<-apply(OU.GRW_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  #Structuring output in lists
  if ("OU.OU" %in% model) {OU.OU_aicc_values <- sapply(result.list$result.list_OU.OU, function(x) {x$AICc})} else {OU.OU_aicc_values<-NULL; } 
  if ("OU.OU" %in% model) {OU.OU_logl_values <- sapply(result.list$result.list_OU.OU, function(x) {x$logL});
  OU.OU_par_values <- sapply(result.list$result.list_OU.OU, function(x) {x$parameters});
  if (is.matrix(OU.OU_par_values) == TRUE) { 
    OU.OU.mean.par<-apply(OU.OU_par_values, 1, function(x) mean(x, na.rm = TRUE))
    OU.OU.SD.par<-apply(OU.OU_par_values,1, function(x) stats::sd(x, na.rm = TRUE))
  } else { 
    OU.OU_par_values<-do.call(rbind, OU.OU_par_values)
    OU.OU.mean.par<-apply(OU.OU_par_values, 2, function(x) mean(x, na.rm = TRUE))
    OU.OU.SD.par<-apply(OU.OU_par_values,2, function(x) stats::sd(x, na.rm = TRUE))
  }
  }
  
  
  # Combine vectors into a list
  vectors <- list(Stasis.Stasis_aicc_values, Stasis.URW_aicc_values, Stasis.GRW_aicc_values, Stasis.OU_aicc_values, URW.Stasis_aicc_values,
                  GRW.Stasis_aicc_values, OU.Stasis_aicc_values,URW.URW_aicc_values, URW.GRW_aicc_values, URW.OU_aicc_values,OU.URW_aicc_values,
                  GRW.URW_aicc_values, GRW.OU_aicc_values, GRW.GRW_aicc_values, OU.GRW_aicc_values, OU.OU_aicc_values)
                  
                  
  # Filter out NULL values
  non_null_vectors <- Filter(Negate(is.null), vectors)
  
  if (length(model) >1) {
    # Combine non-NULL vectors into a matrix
    comparison <- do.call(rbind, non_null_vectors)
    
    complete_cols <- apply(comparison, 2, function(x) all(!is.na(x)))
    
    # Subset the matrix to keep only complete columns
    comparison <- comparison[, complete_cols]
    
    
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
  if ("Stasis.Stasis" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(Stasis.Stasis.mean.par)
    w.list[[i]]$sd.parameters<-(Stasis.Stasis.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("Stasis.URW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(Stasis.URW.mean.par)
    w.list[[i]]$sd.parameters<-(Stasis.URW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("Stasis.GRW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(Stasis.GRW.mean.par)
    w.list[[i]]$sd.parameters<-(Stasis.GRW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("Stasis.OU" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(Stasis.OU.mean.par)
    w.list[[i]]$sd.parameters<-(Stasis.OU.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("URW.Stasis" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(URW.Stasis.mean.par)
    w.list[[i]]$sd.parameters<-(URW.Stasis.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("GRW.Stasis" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(GRW.Stasis.mean.par)
    w.list[[i]]$sd.parameters<-(GRW.Stasis.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("OU.Stasis" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(OU.Stasis.mean.par)
    w.list[[i]]$sd.parameters<-(OU.Stasis.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("URW.URW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(URW.URW.mean.par)
    w.list[[i]]$sd.parameters<-(URW.URW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("URW.GRW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(URW.GRW.mean.par)
    w.list[[i]]$sd.parameters<-(URW.GRW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  
  if ("URW.OU" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(URW.OU.mean.par)
    w.list[[i]]$sd.parameters<-(URW.OU.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("OU.URW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(OU.URW.mean.par)
    w.list[[i]]$sd.parameters<-(OU.URW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("GRW.URW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(GRW.URW.mean.par)
    w.list[[i]]$sd.parameters<-(GRW.URW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("GRW.OU" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(GRW.OU.mean.par)
    w.list[[i]]$sd.parameters<-(GRW.OU.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("GRW.GRW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(GRW.GRW.mean.par)
    w.list[[i]]$sd.parameters<-(GRW.GRW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("OU.GRW" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(OU.GRW.mean.par)
    w.list[[i]]$sd.parameters<-(OU.GRW.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  if ("OU.OU" %in% model) {
    i<-i+1
    w.list[[i]]$mean.parameter.values<-(OU.OU.mean.par)
    w.list[[i]]$sd.parameters<-(OU.OU.SD.par)
    if (length(model) >1) {w.list[[i]]$n.best.AICc<-sum(counts_matrix[i,])}
  }
  
  return(w.list)
  
}
