do.multivariate <- function(reads, N, densNames) {
  #reads:  matrix [T*Nmod] ---> read from memory
  #N = number of combinatorial states that the user wants to use
  #T = number of bins
  T = nrow(reads)
  Nmod = ncol(reads)
  #TODO this here is not true!!! The reads will have other columns at the beginning, with the chromosome, the bins, etc...!!!
  modifications = colnames(reads)#modifications = colnames(reads)[,4:ncol(reads)]

  print(paste("Nmod = ",Nmod,", seqlength = ",T,", modifications = "))
  print(modifications)

  uniresults = list()
  #read the univariate results. They are the output from the univariate part, called uni_modification.RData. 
  #If we don't run the univariate: we may want to consider giving the user the option to upload its own results here!!!!!!
  for (modification in modifications) {
	fname = paste("uni_", modification, ".RData", sep="")
	if (file.exists(fname)) {
	  name = load(fname) 
	  print(paste("loading results from file for modification ",modification))
	  uniresults[[modification]] = get(name)
	} else {
	  #TODO: this is the default densityNames. Can we change it???
	  print(paste("doing univariate for modification ",modification))
	  posterior = univariate.hmm(reads[,modification], densityNames=densNames, max.it=500, eps=1e-2) 
	  #TODO: what happens if we run different chromosomes???????????????????????????????
	  # the result is a posterior matrix with attr "params"
	  uniresults[[modification]] = posterior
	  save(posterior, file=fname)
	}
  }

  params = lapply(uniresults, attr, "params")
  
  #this is hard coded for only NB or Z. If that's not the case, the density parameters would be very different!!!!!!!
  densityType = array(NA, dim=c(Nmod,2))
  for(imod in 1:Nmod){
	for(um in 1:2){
		if(attr(uniresults[[imod]],"params")$distributions[um,"w"]==0) densityType[imod,um]="NB"
		if(attr(uniresults[[imod]],"params")$distributions[um,"w"]>0) densityType[imod,um]="Z"
		#if(params[[imod]]$distributions[um,3]==0) densityType[imod,um]="NB"
		#if(params[[imod]]$distributions[um,3]>0) densityType[imod,um]="Z"
	}
  }
  
  enrichment = NULL
  for(imod in 1:Nmod){
	enrichment = cbind(enrichment, apply(uniresults[[imod]], 1, which.max ) )
  }
  enrichment = enrichment - 1

  CombState = rep(0,T)
  for(imod in 1:Nmod) {
	CombState = CombState + 2^(Nmod-imod)*enrichment[,imod]
  }
  #RANK
  s = sort(table(CombState),decreasing=TRUE)
  states = as.numeric(names(s))

  print("densityType for each mark:")
  print(densityType)
  print("combinatorial states ranking:")
  print(s)

  lib = require(gsl)
  if( lib == FALSE ){
	install.packages("gsl")
	require(gsl)
  }
  
  print("computing correlation matrices ... \n")
  MAX = max(reads)
  temp = matrix(rep(NA,c(MAX+1)*Nmod*2),ncol=Nmod*2) #rows:T, columns: 2*Nmod --> for each modification, first unmodified then modified
  for(imod in 1:Nmod)
  {
	for(um in 1:2)
	{
		if(enc2utf8(densityType[imod,um]) == enc2utf8("Z"))
		{
			r = attr(uniresults[[imod]],"params")$distributions[um,"r"]
			p = attr(uniresults[[imod]],"params")$distributions[um,"p"]
			w = attr(uniresults[[imod]],"params")$distributions[um,"w"]
			for(t in 1:c(MAX+1))
			{
				#u = pzinba (reads[imod,t], size = r[imod,um], mu = mu[imod,um], beta = w[imod,um]) #CDF, (x, size, mu, beta, lower.tail = TRUE)
				Ot = t-1 #reads[t,imod]
				lGamma1plusRplusX=lgamma(1+r+Ot)
				lGammaR=lgamma(r)
				lGamma2plusX=lgamma(2+Ot)
				lHyper=log( hyperg_2F1( 1, 1+r+Ot, 2+Ot, 1-p ) )#log(gsl_sf_hyperg_2F1(1, 1+r[imod,um]+Ot, 2+Ot, 1-p[imod,um]));
				lppowert=(1+Ot)*log(1-p)
				lppowerr=r*log(p)
				u = 1 - exp( log(1-w) + lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX )
				temp[t,(imod - 1) * 2 + um] = qnorm(u) #z[t,(imod - 1) * 2 + um] = qnorm(u) #inverse normal
				if(temp[t,(imod - 1) * 2 + um] == Inf){
				   if(max(reads[,modifications[imod]]) >=t){
					if(um==0) cat("for reads=", t, ", imod=", imod, ", unmodified component ---> copula transformation would be= ", temp[t,(imod-1) *2+um], "\n")
					if(um==1) cat("for reads=", t, ", imod=", imod, ", modified component ---> copula transformation would be = ", temp[t,(imod-1)*2+um], "\n")					
				   }
				   temp[t,(imod-1)*2+um] = temp[t-1, (imod-1)*2+um]
				}
			}
		}
		if(enc2utf8(densityType[imod,um]) == enc2utf8("NB"))
		{
			r = attr(uniresults[[imod]],"params")$distributions[um,"r"]
			p = attr(uniresults[[imod]],"params")$distributions[um,"p"]
			for(t in 1:c(MAX+1))
			{
				Ot = t-1 #reads[t,imod]
				lGamma1plusRplusX=lgamma(1+r+Ot);
				lGammaR=lgamma(r);
				lGamma2plusX=lgamma(2+Ot);
				#f21hyper {BMS}
				lHyper=log( hyperg_2F1( 1, 1+r+Ot, 2+Ot, 1-p ) )#f21hyper( 1, 1+r+Ot, 2+Ot, 1-p ) ) #gsl_sf_hyperg_2F1(1, 1+r+Ot, 2+Ot, 1-p));
				lppowert=(1+Ot)*log(1-p);
				lppowerr=r*log(p);
				u = 1 - exp( lppowerr + lppowert + lGamma1plusRplusX + lHyper - lGammaR - lGamma2plusX );
				temp[t,(imod - 1) * 2 + um] = qnorm(u) #z[t,(imod - 1) * 2 + um] = qnorm(u) #inverse normal
				if(temp[t,(imod - 1) * 2 + um] == Inf){
				   if(max(reads[,modifications[imod]]) >=t){
					if(um==0) cat("for reads=", t, ", imod=", imod, ", unmodified component ---> copula transformation would be= ", temp[t,(imod-1) *2+um], "\n")
					if(um==1)1+1# cat("for reads=", t, ", imod=", imod, ", modified component ---> copula transformation would be = ", temp[t,(imod-1)*2+um], "\n") 
				   }
				   temp[t,(imod-1)*2+um] = temp[t-1,(imod-1)*2+um]
				}
			}
		}
	}
  }

  z = matrix(rep(NA,T*Nmod*2),ncol=Nmod*2) #rows:T, columns: 2*Nmod --> for each modification, first unmodified then modified
  for(imod in 1:Nmod)
  {
	for(um in 1:2)
	{
		for(t in 1:T)
		{
			z[t,(imod-1)*2+um] = temp[c(reads[t,imod]+1),(imod-1)*2+um]
		}
	}
  }
	

  #correlationMatrix = array(NA, dim=c(N,Nmod,Nmod))
  correlationMatrixInverse = array(NA, dim=c(Nmod,Nmod,N))
  deter = array(NA, dim=N)
  for(i in 1:N)
  {
	iN = states[i]
	temp = iN# iN-1 ???????????????????
	#cat(iN)
	#cat("\t")
	#cat(temp)
	enrich = NULL
	for(imod in 1:Nmod){
		if(temp>=1) enrich = c( temp%%2 ,enrich)
		else enrich = c(0, enrich)
		if(temp > 0) temp = (temp-enrich[1]) / 2
		if(temp<0) break
	}
	#cat("\n")
	#cat(enrich)
	#cat("\n\n")
	enrichTF = NULL
	for(imod in 1:Nmod)
	{
		if (enrich[imod] == 0) enrichTF = c( enrichTF, c(TRUE, FALSE))
		if (enrich[imod] == 1) enrichTF = c( enrichTF, c(FALSE, TRUE))
	}
	some = NULL
	for(imod in 1:c(Nmod*2))
	{
		if(enrichTF[imod] == TRUE) some = c(some, imod)
	}
	#correlationMatrix[i,,] = cor(z[which(c(CombState+1)==iN),some]) #correlationMatrix[i] corresponds to state iN = states[i]!!!!
	corMat = cor(z[which(c(CombState)==iN),some])
cat("CORRELATION MATRIX state ",i,"\n")
print(corMat)
	correlationMatrixInverse[,,i] = solve(corMat) #correlationMatrix[i] corresponds to state iN = states[i]!!!!
cat("CORRELATION MATRIX INV state ",i,"\n")
print(correlationMatrixInverse[,,i])
	deter[i] = det( corMat )
  }
  
  print("doing multivariate")
  multiresult = multivariate.hmm(reads, correlationMatrixInverse, deter, params, states)
  return(multiresult)
  
}

