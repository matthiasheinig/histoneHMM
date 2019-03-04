table2EMResult <- function(tab, prop=NULL) {
  emfit = new("EMResult")
  for (i in 1:nrow(tab)) {
    emfit@models[[i]] = new("Zinba")
    emfit@models[[i]]@mu = tab[i,"mu"]
    emfit@models[[i]]@size = tab[i,"size"]
    emfit@models[[i]]@beta = tab[i,"beta"]
  }
  if (is.null(prop)) {
    prop = tab[,"prop"]
  }
  emfit@proportions = prop
  return(emfit)
}

threecomponent.to.zinba.mixture <- function(trunc.mix) {
  unmodified.idx = which.min(sapply(trunc.mix@models[2:3], function(x) {
    if (inherits(x, "TruncatedDistribution")) {
      return(x@distribution@mu)
    } else {
      return(x@mu)
    }
  }))
  modified.idx = setdiff(1:2, unmodified.idx)
  ## transform to a regular zinba mixture 
  zero.prop = trunc.mix@proportions[1]
  unmodified.prop = trunc.mix@proportions[1 + unmodified.idx] + zero.prop
  modified.prop = trunc.mix@proportions[1 + modified.idx]
  unmodified = trunc.mix@models[[1 + unmodified.idx]]
  if (inherits(unmodified, "TruncatedDistribution")) {
    unmodified = unmodified@distribution
  }
  unmodified = new("Zinba", mu=unmodified@mu, size=unmodified@size, beta=zero.prop / unmodified.prop)
  modified = trunc.mix@models[[1 + modified.idx]]
  if (inherits(modified, "TruncatedDistribution")) {
    modified = modified@distribution
  }
  models = list(unmodified, modified)
  zinba.mix = new("EMResult", models=models, proportions=c(unmodified.prop, modified.prop))
  return(zinba.mix)
}

run.univariate.hmm <- function(fname, data=NULL, n.expr.bins=5, em=FALSE, chrom=NULL, maxq=1-1e-4, redo=FALSE, model.constructor=NULL, baum.welch=FALSE, plot="pdf") {

  if (is.null(data)) {
    data = read.csv(fname, sep="\t", stringsAsFactors=F)
  }

  max = round(quantile(data$signal, p=maxq))
  data$signal[data$signal > max] = max
  
  if (is.null(chrom)) {
    chrom = unique(data$chrom)
  }

  if (!em) {
    if (!baum.welch) {
      outname = gsub(".txt$", "-posterior.txt", fname)
      pfile = gsub(".txt$", "-zinba-params.txt", fname)
      if (!file.exists(pfile) || redo) {
        ## fit zero inflated negative binomials for highest / lowest expressed
        high = data$signal[data$bin.type == paste("expressed", n.expr.bins, sep=".") & data$chrom %in% chrom]
        low = data$signal[data$bin.type == paste("expressed", 1, sep=".") & data$chrom %in% chrom]
        
        zinba.fit.low = fitzinba(low)
        zinba.fit.high = fitzinba(high)
        
        zinba.params = rbind(highly.expressed=zinba.fit.high,
          lowly.expressed=zinba.fit.low)
        zinba.params = cbind(zinba.params, prop=c(length(high), length(low)) / (length(high) + length(low)))
        write.table(zinba.params, file=pfile, sep="\t", quote=F)
  
        zinba.mix = table2EMResult(zinba.params)
        save(zinba.mix, file=gsub(".txt$", "-zinba-params.RData", fname))

        if (plot != "none") {
          if (plot == "pdf") {
            pdf(file=gsub(".txt$", "-zinba-fit.pdf", fname))
          } else if (plot == "svg") {
            svg(file=gsub(".txt$", "-zinba-fit.svg", fname))
          }
          plotDensity(model=zinba.mix, x=data$signal, xlim=c(0, max), add=FALSE, col="black", alpha=1)
          dev.off()
        }
      } else {
        zinba.params = read.csv(pfile, sep="\t")
      }
    } else {
      ## baum welch here
      outname = gsub(".txt$", "-bw-posterior.txt", fname)
      pfile = gsub(".txt$", "-zinba-params-bw.txt", fname)
      if (!file.exists(pfile) || redo) {
        signal = data$signal[data$chrom %in% chrom]
        signal = signal[signal < max] ## remove the point of truncation
        ## initialize like aaron does
        m = mean(signal)
        v = var(signal)
        v1 = v
        v2 = v
        if (m >= v) {
          v1 = m + 1
          v2 = m + 2
        }
        ## zero size and mu means this is the zero inflation component
        mu = c(0, m, m + 1)
        size = c(0, m^2 / (v1 - m), (m + 1)^2 / (v2 - m - 1))
        post = hmm.baumwelch.negbinom(signal, mu=mu, size=size, eps=0.1)
        ## determine which is the high and low occupancy state
        mu = attr(post, "param")$mu
        size = attr(post, "param")$size
        high.idx = which.max(mu)
        low.idx = setdiff(2:3, high.idx)
        idx = c(low.idx, high.idx)
        ## get the proportions from the posteriors
        prop = colSums(post) / nrow(post)
        ## put the zero inflation into the low occupancy component
        prop2 = prop[idx]
        prop2[1] = prop2[1] + prop[1]
        zinba.params = cbind(mu=mu[idx], size=size[idx], beta=c(prop[1] * 1 / prop2[1], 0), prop=prop2)
        rownames(zinba.params) = c("highly.expressed", "lowly.expressed")
        write.table(zinba.params, file=pfile, sep="\t", quote=F)
  
        zinba.mix = table2EMResult(zinba.params)
        save(zinba.mix, file=gsub(".txt$", "-zinba-params-bw.RData", fname))
        
        pdf(file=gsub(".txt$", "-zinba-bwfit.pdf", fname))
        plotDensity(model=zinba.mix, x=data$signal, xlim=c(0, max), add=FALSE, col="black", alpha=1)
        dev.off()
      } else {
        zinba.params = read.csv(pfile, sep="\t")
      }
    }
  } else {
    ## EM here
    outname = gsub(".txt$", "-em-posterior.txt", fname)
    pfile = gsub(".txt$", "-zinba-params-em.txt", fname)
    if (!file.exists(pfile) || redo) {
      signal = data$signal[data$chrom %in% chrom]
      signal = signal[signal < max] ## remove the point of truncation
      if (is.null(model.constructor)) {
        model.constructor = c("Zero", "Negbinom", "Negbinom")
      }
      if (length(model.constructor) == 1) {
	ncomp = 2
      } else {
	ncomp = length(model.constructor)
      }
      zinba.mix = em(signal, ncomp=ncomp, model.constructor=model.constructor, eps=0.1)
      if (ncomp == 3) {
	## the three component mixture works much faster because we have a
	## c function for the M step in the EM algorithm
	zinba.mix = threecomponent.to.zinba.mixture(zinba.mix)
      }

      high.idx = which.min(sapply(zinba.mix@models, function(x) x@mu))
      low.idx = setdiff(1:2, high.idx)

      ## check if our fit is good, the modified component should have the 
      ## higher posterior until the maximum data value
      post = EMPosterior(zinba.mix, 0:max)
      lgt = post[,low.idx] >= post[,high.idx]
      equiprob = max(which(lgt)) - 1
      nruns = length(rle(lgt)$lengths)
      if (equiprob < max || nruns > 2) {
	cat("First EM result was not well behaved, using truncated distributions\n")
	if (equiprob == max) {
	  ## here we need to flip states
	  equiprob = max(which(!lgt)) - 1
	}
	## the modified state does not extend to the end
	## we rerun the estimation with truncated distributions
	m = mean(signal)
	v = var(signal)
	v1 = v
	v2 = v
	if (m >= v) {
	  v1 = m + 1
	  v2 = m + 2
	}
	models = list(zero=new("Zero"),
	  unmodified=TruncatedDistribution(new("Negbinom",
	    mu=m,
	    size=m^2 / (v1 - m)), 2 * equiprob, left=F),
	  modified  =TruncatedDistribution(new("Negbinom",
	    mu=m + 1,
	    size=(m + 1)^2 / (v2 - m - 1)), equiprob, left=T))
	trunc.mix = em(signal, 3, eps=0.1, models=models)
	zinba.mix = threecomponent.to.zinba.mixture(trunc.mix)
      }

      save(zinba.mix, file=gsub(".txt$", "-zinba-params-em.RData", fname))
      
      to.row <- function(model) {
        beta = 0
        names(beta) = "beta"
        if (class(model) == "Zinba") {
          beta = model@beta
        }
        return(c(model@size, model@mu, beta))
      }
      
      zinba.params = rbind(
        highly.expressed=to.row(zinba.mix@models[[high.idx]]),
        lowly.expressed=to.row(zinba.mix@models[[low.idx]]))
      colnames(zinba.params) = c("size", "mu", "beta")
      zinba.params = cbind(zinba.params, prop=c(zinba.mix@proportions[high.idx], zinba.mix@proportions[low.idx]))

      write.table(zinba.params, file=pfile, sep="\t", quote=F)

      if (plot != "none") {
        if (plot == "pdf") {
          pdf(file=gsub(".txt$", "-zinba-emfit.pdf", fname))
        } else if (plot == "svg") {
          svg(file=gsub(".txt$", "-zinba-emfit.svg", fname))
        }
        par(mfrow=c(2, 1))
        plotDensity(model=zinba.mix, x=signal, xlim=c(0, max), add=FALSE, col="black", alpha=1)
        title(main="Data used to fit the model")
        plotDensity(model=zinba.mix, x=data$signal, xlim=c(0, max), add=FALSE, col="black", alpha=1)
        title(main="All data")
        dev.off()
      }
    } else {
      zinba.params = read.csv(pfile, sep="\t")
    }
  }

  if (redo || !file.exists(outname)) {
    ## univariate HMM
    for (chrom in unique(data$chrom)) {
      cat(chrom, "\n")
      this.chrom = data$chrom == chrom
      posterior = zinba.hmm.posterior(data$signal[this.chrom], mu=zinba.params[,"mu"], size=zinba.params[,"size"], beta=zinba.params[,"beta"])
      colnames(posterior) = rownames(zinba.params)
      res = cbind(data[this.chrom,], posterior)
      first = (chrom == unique(data$chrom)[1])
      write.table(res, outname, sep="\t", quote=F, row.names=F, append=!first, col.names=first)
    }
  }
  posterior = read.csv(outname, sep="\t")
  return(invisible(posterior))
}


run.bivariate.hmm <- function(fname1, fname2, outdir, data1=NULL, data2=NULL, sample1=NULL, sample2=NULL, n.expr.bins=5, maxq=1-1e-4, em=FALSE, chrom=NULL, baum.welch=FALSE) {
  ## read the data
  if (is.null(data1)) {
    data1 = read.csv(fname1, sep="\t", stringsAsFactors=F)
  }
  if (is.null(data2)) {
    data2 = read.csv(fname2, sep="\t", stringsAsFactors=F)
  }
  stopifnot(nrow(data1) == nrow(data2))
  max1 = round(quantile(data1$signal, p=maxq))
  data1$signal[data1$signal > max1] = max1
  max2 = round(quantile(data2$signal, p=maxq))
  data2$signal[data2$signal > max2] = max2

  ## for the em we will not use the whole data set
  if (is.null(chrom)) {
    chrom = unique(data1$chrom)
  }


  ## change the definition of the differential states: use the posteriors
  ## of the univariate analysis instead
  legacy = FALSE
  
  if (legacy) {
    ## define states that behave the same for sample1 and sample2
    both.high = which(data1$bin.type == paste("expressed", n.expr.bins, sep=".") & data2$bin.type == paste("expressed", n.expr.bins, sep="."))
    both.low = which(data1$bin.type == "expressed.1" & data2$bin.type == "expressed.1")
    
    ## define the differential states
    not.expressed = paste("expressed", 1:2, sep=".")
    really.expressed = paste("expressed", 3:5, sep=".")
    up = which(data1$bin.type %in% not.expressed & data2$bin.type %in% really.expressed)
    down = which(data1$bin.type %in% really.expressed & data2$bin.type %in% not.expressed)  
    
  } else {
    ## use the calls from two univariate analyses
    posterior1 = run.univariate.hmm(fname1, data=data1, n.expr.bins=n.expr.bins, em=em, chrom=chrom, maxq=maxq, redo=FALSE, baum.welch=baum.welch)
    posterior2 = run.univariate.hmm(fname2, data=data2, n.expr.bins=n.expr.bins, em=em, chrom=chrom, maxq=maxq, redo=FALSE, baum.welch=baum.welch)

    ## use posterior > 0.9 to make a call (if there are no such clearly 
    ## differential regions we use less stringent definitions)
    threshold = 0.9
    while (threshold > 0.5) {
      cols = c("highly.expressed", "lowly.expressed")
      map1 = cols[apply(posterior1[,cols], 1, function(x) {
	gt = x > threshold
	if (any(gt, na.rm=T)) return(which(gt))
	else return(NA)})]
      map2 = cols[apply(posterior2[,cols], 1, function(x) {
	gt = x > threshold
	if (any(gt, na.rm=T)) return(which(gt))
	else return(NA)})]
      
      both.high = which(map1 == "highly.expressed" & map2 == "highly.expressed")
      both.low = which(map1 == "lowly.expressed" & map2 == "lowly.expressed")
      
      ## define the differential positions
      up = which(map1 == "lowly.expressed" & map2 == "highly.expressed")
      down = which(map1 == "highly.expressed" & map2 == "lowly.expressed")

      if (length(up) < 50 || length(down) < 50) {
	threshold = threshold - 0.1
      } else {
	break
      }
    }
  }

  ## fit the distributions
  fit.low = fit.zinba.copula(data1$signal[both.low], data2$signal[both.low])
  fit.high = fit.zinba.copula(data1$signal[both.high], data2$signal[both.high])
  
  ## make differential states by swapping marginals
  fit.lowx = fit.zinba.copula(data1$signal[up], data2$signal[up], marginal.x=fit.low$marginal.x, marginal.y=fit.high$marginal.y)
  fit.lowy = fit.zinba.copula(data1$signal[down], data2$signal[down], marginal.x=fit.high$marginal.x, marginal.y=fit.low$marginal.y)

  fits = list(unmod.both=fit.high, mod.x=fit.lowx, mod.y=fit.lowy, mod.both=fit.low)
  zinba.params = sapply(fits, fit.to.col)

  
  dir.create(outdir)

  if (is.null(sample1)) {
    sample1 = "sample1"
  }
  if (is.null(sample2)) {
    sample2 = "sample2"
  }
  outname = paste(sample1, "-vs-", sample2, sep="")
  write.table(zinba.params, file=file.path(outdir, paste(outname, "-zinbacopula-params.txt", sep="")), sep="\t", quote=F)

  # compute posteriors
  for (chrom in unique(data1$chrom)) {
    cat(chrom, "\n")
    this.chrom = data1$chrom == chrom
    signal1 = data1[this.chrom, "signal"]
    signal2 = data2[this.chrom, "signal"]
    signal = rbind(signal1, signal2)
    posterior.zinba = zinbacopula.hmm.posterior(signal, fits)  
    colnames(posterior.zinba) = colnames(zinba.params)
    colnames(posterior.zinba) = gsub("mod.x", sample1, colnames(posterior.zinba))
    colnames(posterior.zinba) = gsub("mod.y", sample2, colnames(posterior.zinba))
    # get the MAP state
    map.zinba = colnames(posterior.zinba)[apply(posterior.zinba, 1, which.max)]
    
    d1 = data1[this.chrom, c("signal", "bin.expr")]
    colnames(d1) = paste(colnames(d1), sample1, sep=".")
    d2 = data2[this.chrom, c("signal", "bin.expr")]
    colnames(d2) = paste(colnames(d2), sample2, sep=".")
    # and of course the coordinates of the bins
    res = cbind(data1[this.chrom, c("chrom", "start", "end", "bin.type")], d1, d2, posterior.zinba, map.zinba)
  
    first = (chrom == unique(data1$chrom)[1])
    write.table(res, file.path(outdir, paste(outname, ".txt", sep="")), sep="\t", quote=F, row.names=F, col.names=first, append=!first)
  }
  bivariate.posterior = read.csv(file.path(outdir, paste(outname, ".txt", sep="")), sep="\t", stringsAsFactors=F)
  return(invisible(bivariate.posterior))
}

plot.bivariate <- function(fname1, fname2, outdir, sample1, sample2, maxq=1-1e-3) {
  ## also check the bivariate density
  data1 = read.csv(fname1, sep="\t", stringsAsFactors=F)
  data2 = read.csv(fname2, sep="\t", stringsAsFactors=F)

  maxcount1 = quantile(data1$signal, maxq)
  data1$signal[data1$signal > maxcount1] = maxcount1
  maxcount2 = quantile(data2$signal, maxq)
  data2$signal[data2$signal > maxcount2] = maxcount2

  obs = table(factor(data1$signal, levels=0:maxcount1),
      factor(data2$signal, levels=0:maxcount2))

  zinba.copula = read.csv(paste(outdir, "/", sample1, "-vs-", sample2, "-zinbacopula-params.txt", sep=""), sep="\t")
  fits = apply(zinba.copula, 2, col.to.fit)

  ## get the mixing proportions from the posterior probabilities
  posterior = read.csv(paste(outdir, "/", sample1, "-vs-", sample2, ".txt", sep=""), sep="\t", stringsAsFactors=F)
  prop = colSums(posterior[,match(c("unmod.both", sample1, sample2, "mod.both"), colnames(posterior))]) / nrow(posterior)
  names(prop) = names(fits)
  
  ## compute the probabilities for each state
  mixture = matrix(0, nrow=maxcount1 + 1, ncol=maxcount2 + 1)
  expected = array(dim=c(maxcount1 + 1, maxcount2 + 1, length(fits)))
  for (i in 1:length(fits)) {
    fit = fits[[i]]
    expected[,,i] = sapply(0:maxcount2, function(y) sapply(0:maxcount1, function(x) dzinba.copula(x, y, fit)))
    mixture = mixture + prop[i] * expected[,,i]
  }

  ## plot everything
  cols = rev(rainbow(12, start=0, end=4/6))

  par(mfrow=c(3, 2))
  image(0:maxcount1, 0:maxcount2, log10(obs), col=cols, main="Observed", xlab="count x", ylab="count y")
  image(0:maxcount1, 0:maxcount2, log10(round(mixture * nrow(data1))), col=cols, main="Mixture", xlab="count x", ylab="count y")
  for (i in 1:length(fits)) {
    image(0:maxcount1, 0:maxcount2, log10(round(expected[,,i] * nrow(data1))), col=cols, main=names(fits)[i], xlab="count x", ylab="count y")
  }
}

get.differential.regions <- function(fname1, fname2, outdir, sample1, sample2, chrom="chr19", em=FALSE, baum.welch=FALSE, maxq=1-1e-4, cutoff=0.5) {
  dir.create(outdir)
  outname = paste(sample1, "-vs-", sample2, sep="")
  fname = paste(outdir, "/", outname, ".txt", sep="")

  do.calls = FALSE
  if (!file.exists(fname)) {
    bivariate.posterior = run.bivariate.hmm(fname1, fname2, outdir, sample1=sample1, sample2=sample2, chrom=chrom, em=em, baum.welch=baum.welch, maxq=maxq)
    do.calls = TRUE
  } else {
    if (!file.exists(paste(outdir, "/", outname, "-unmod_both.gff", sep=""))) {
      bivariate.posterior = read.csv(paste(outdir, "/", sample1, "-vs-", sample2, ".txt", sep=""), sep="\t", stringsAsFactors=F)
      do.calls = TRUE
    }
  }
  
  if (do.calls) {
    unmod.both = callRegions(bivariate.posterior, cutoff, "unmod.both", NULL)
    GR2gff(unmod.both, paste(outdir, "/", outname, "-unmod_both.gff", sep=""))
    mod.both = callRegions(bivariate.posterior, cutoff, "mod.both", NULL)
    GR2gff(mod.both, paste(outdir, "/", outname, "-mod_both.gff", sep=""))
    sample1.regions = callRegions(bivariate.posterior, cutoff, make.names(sample1), NULL)
    GR2gff(sample1.regions, paste(outdir, "/", outname, "-", sample1, ".gff", sep=""))
    sample2.regions = callRegions(bivariate.posterior, cutoff, make.names(sample2), NULL)
    GR2gff(sample2.regions, paste(outdir, "/", outname, "-", sample2, ".gff", sep=""))
  }
  
  sample1.regions = gff2GR(paste(outdir, "/", outname, "-", sample1, ".gff", sep=""))
  sample2.regions = gff2GR(paste(outdir, "/", outname, "-", sample2, ".gff", sep=""))
  unmod.both = gff2GR(paste(outdir, "/", outname, "-unmod_both.gff", sep=""))
  mod.both = gff2GR(paste(outdir, "/", outname, "-mod_both.gff", sep=""))
  hmm = c(sample1.regions, sample2.regions)
  hmm.regions = list(unmod.both=unmod.both, sample1=sample1.regions, sample2=sample2.regions, mod.both=mod.both)

  return(list(differential=hmm, regions=hmm.regions))
}
