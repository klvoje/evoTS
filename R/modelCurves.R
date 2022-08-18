# Function from paleoTS (License:	GPL-3)
modelCurves <-
function(x, w, np=500)
# returns list of model means, upper and lower 95% probability envelopes
{
  ee<- ii<- array(dim=np)  # set up arrays
  mn<- w$modelName
  mp<- w$par
  x0<- ifelse(w$method=="AD", x$mm[1], w$par["anc"])
  ttp<- seq(x$tt[1], x$tt[length(x$tt)], length.out=np)
  #ttp<- x$tt
  #ttp<- sort(c(seq(x$tt[1], x$tt[length(x$tt)], length.out=np), x$tt))
  #print(ttp)
  
  # list of models implemented
  okModels<- c("URW", "GRW", "Stasis", "OU", paste("Punc-", 1:100, sep=""))
  if (mn %in% okModels){  
	  if(mn=="URW"){ ee<- rep(x0, np); vv<- mp["vstep"]*ttp}  
	  if(mn=="GRW"){ ee<- x0+mp["mstep"]*ttp; vv<- mp["vstep"]*ttp}
	  if(mn=="Stasis"){ ee<- rep(mp["theta"], np); vv<- rep(mp["omega"],np) }
	  if(mn=="OU"){ if(!is.null(x$start.age)) tto<- x$start.age-ttp
	  				else	tto<- ttp
	  				ee<- mp["theta"] * (1-exp(-mp["alpha"]*tto)) + mp["anc"]*exp(-mp["alpha"]*tto)
	  				vv<- (mp["vstep"]/(2*mp["alpha"])) * (1-exp(-2*mp["alpha"]*tto)) }
	  if(grepl("Punc", mn)){
	  		# handle time shifts
	  		sp<- w$par[grep("shift", names(w$par))]  # shift samples
	  		st<- x$tt[sp]  # shift times
	  		ng<- length(sp)+1
	  		ggt<- cut(ttp, breaks=c(min(x$tt)-1, st, max(x$tt)+1), right=TRUE, labels=1:ng)
	
	
	  		# extract needed paramaters
	  		th<- w$par[grep("theta", names(w$par))]  # theta estimates
	  		om<- w$par[grep("omega", names(w$par))]  # omega estimates
	  		if(length(om)==1)	om<- rep(om, ng)	 # make into vector if needed
			
			# get ee and vv
			ee<- th[ggt]
			vv<- om[ggt]  		
		}  
	} else {ee<-NA; vv<- NA; warning(paste("modelFit argument not implemented for model", mn, "\n"))}
	
	
  if (!is.null(x$start.age))	tto<- x$start.age - ttp  else tto<- ttp 
  res<- list(tt=ttp, ee=ee, ll=ee-1.96*sqrt(vv), uu=ee+1.96*sqrt(vv))	
  return(res)	
}
