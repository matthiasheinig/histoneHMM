#FOR 1 CHROMOSOME AS THE CURRENT POSTERIORS NOW ARE OUTPUTED FOR 1 CHROMOSOMES. Later, when we analyse more chromosomes, we just need to write a higher-level function to call this function for the different chromosomes
# Given genome coordinates of, e.g., genes, we calculate:
#1. Percentage of state overlaping the genes w.r.t state posterior
#2. Percentage of state overlaping the genes w.r.t state hard assignment (state call) and its p-value based on the hypergeometric distribution
#3. percentage of the genes overlaping a state (so, sum over all states should be 1)
#2. The fold enrichment w.r.t to state posterior

###INPUT:
#posteriors: a matrix of posteriors of states which have column names State1 State2 ...
#genome.coords: a matrix contains coordinates of the external data source with to columns for begin_coord end_coord
#bsize: bin size
#chrSize: chromosome size
genomedata.enrichment <- function(posteriors, genome.coords, bsize, chrSize){
	nBin <- nrow(posteriors)
	nState <- length(posteriors[1,])
	maxStateCol <- nState +1
	#get the state which the maximum posterior for each bin
	posteriors <- cbind(posteriors, apply(posteriors,1,which.max))

	#For a state, calculate the sum of the posterior probability over all bin and save to stateSumPosterior vector
	stateSumPosterior <- colSums(posteriors[,1:nState])
	
	#For a state, count its frequency over all bins	and save in stateFre vector
	temp <- sort(table(posteriors[,nState+1]))
	#the state in frequency-increasing order 
	state <- as.numeric(names(temp))
	stateFre <- rep(0,nState)
	for ( i in 1:nState ){
		stateFre[state[i]] = temp[i]
	}
	rm(temp)

	#For an external data source given by genome.coords, calculate the total number of bins that it overlaps by at least one base. Save to overlapAll variable
	overlapAll <- 0
	#For a state and an external data source, calculate the total sum of the posterior for the state in the bins overlapping the external data source. Save to state_dataSumPosterior vector
	state_dataSumPosterior <- rep(0,nState)
	#Also, count the number of bins overlapped by a state and this data source. Save to overlapStateFre
	overlapStateFre <- rep(0,nState)
	#Also, compute the percentage of the data source in genome.coords found in a state (hard assignment)
	perData <- rep(0,nState)
	#number of coordinates given
	nCoords = nrow(genome.coords)
	for ( iCoord in 1:nCoords ){
		 if ( genome.coords[iCoord,1] <= chrSize ){
			begin <-floor(genome.coords[iCoord,1]/bsize) + 1
			end <- 0
			if (as.numeric(genome.coords[iCoord,2]) < chrSize){
				end <- floor((genome.coords[iCoord,2]-1)/bsize) + 1
			} else {
				end <- floor(chrSize/bsize) + 1
			}
			overlapAll = overlapAll + (end-begin+1);
			#add the posterior for states in overlapping bins
                	for ( bin in begin:end )
	                {
	                	overlapStateFre[posteriors[bin,maxStateCol]] = overlapStateFre[posteriors[bin,maxStateCol]] + 1
				#curMaxState <- posteriors[bin,maxStateCol]
                        	for ( state in 1:nState ){
                        		state_dataSumPosterior[state] = state_dataSumPosterior[state] + posteriors[bin,state]
	                        }
	 		}
         	} else {
			warning <- paste("WARNING: start of the coordinate (",genome.coords[iCoord,1],",",genome.coords[iCoord,2], ") is bigger than the chromosome size! We exclude this coordinate.",sep="")
			print(warning)
		}
	}		
	#the overlap percentage w.r.t state's posteriors (i.e. percentage of a state that overlaps the external data) and fold enrichment
	overlap <- rep(0,nState)
	foldEnrichment <- rep(0,nState)
	for ( state in 1:nState ){
		overlap[state] <- state_dataSumPosterior[state]/stateSumPosterior[state]*100
		foldEnrichment[state] <- state_dataSumPosterior[state]/stateSumPosterior[state]/(overlapAll/nBin)
	}
	#Also, compute the percentage of the exteral data source in genome.coords found in a state (hard assignment)
	perData <- overlapStateFre/overlapAll*100
	#Compute the percentage of state w.r.t hard assignment overlaps the genome data and the corresponding p-value
	perState <- overlapStateFre/stateFre*100
	#compute the p-value of the overlap (hard assignment) using hypergemetric distribution
	pValue <- phyper(overlapStateFre-1,overlapAll,nBin-overlapAll,stateFre,lower.tail=F)	
	ret <- cbind(t(t(overlap)),t(t(perState)),t(t(pValue)),t(t(perData)),t(t(foldEnrichment)))
	colnames(ret,do.NULL=F)
	colnames(ret) <- c("perStatePos","perStateHard","pValue","perData","foldEnrichment")
	rownames(ret) <- c(colnames(posteriors[,1:nState]))
	
	return(ret)
}

###INPUT
#posteriors: a matrix of posteriors of states which has column names e.g. State1 State2 ....
#gene.expression: a matrix (column names is not neccessary) having 3 columns: start end expression_value
#bsize: bin size
#chrSize: chromosome size
#return a matrix of 1 column containin the average expression level for the states (rows represents states)
expression.enrichment <- function(posteriors, gene.expression, bsize, chrSize){
	#Standardized all expression values by subtracting the mean log transformed values and then dividing by the s.d. of the log transformed values
	expr <- log(gene.expression[,3])
	exprS <- (expr - mean(expr))/sd(expr)
	#attach standardized exprssion value to gene.expression
	gene.epxression <- cbind(gene.expression,t(t(exprS)))
	rm(expr)
	rm(exprS)

	#attach state calls
	nBin <- nrow(posteriors)
        nState <- length(posteriors[1,])
        maxStateCol <- nState +1
        #get the state which the maximum posterior for each bin
        posteriors <- cbind(posteriors, apply(posteriors,1,which.max))
	###ATTACH expression value to bins
	#number of coordinates given
        nCoords = nrow(gene.expression)
	binExpr <- rep(0,nBin)
	binProbeCount <- rep(0,nBin)
        for ( iCoord in 1:nCoords ){
                 if ( gene.expression[iCoord,1] <= chrSize ){
                        begin <-floor(gene.expression[iCoord,1]/bsize) + 1
                        end <- 0
                        if (as.numeric(gene.expression[iCoord,2]) < chrSize){
                                end <- floor((gene.expression[iCoord,2]-1)/bsize) + 1
                        } else {
                                end <- floor(chrSize/bsize) + 1
                        }
                        #add the posterior for states in overlapping bins
                        for ( bin in begin:end )
                        {
				binExpr[bin] = binExpr[bin]+gene.expression[iCoord,3]
				binProbeCount[bin] = binProbeCount[bin]+1
                        }
                } else {
                        warning <- paste("WARNING: start of the coordinate (",gene.expression[iCoord,1],",",gene.expression[iCoord,2], ") is bigger than the chromosome size! We exclude this coordinate.",sep="")
                        print(warning)
                }
        }
	binExpr <- binExpr/binProbeCount
	posteriors <- cbind(posteriors,t(t(binExpr)))
	rm(binExpr)
	rm(binProbeCount)
	###done ATTACH expression value to bins
	#compute average expression level of genomic interval overlaping a state	
	avgExpr <- rep(0,nState)
	for ( state in 1:nState ){
		avgExpr[state] = mean(posteriors[posteriors[,nState+1]==state,nState+2])
	}
	ret <- t(t(avgExpr))
	colnames(ret,do.NULL=F)
	colnames(ret)<-c("avgExpr")
	return(ret)
}

#Testing genomedata.enrichment
#posteriors <- read.table("example/ex2-HiddenChain_H3K27me3.txt",header=T)
#genes <- read.table("example/ex2-RefSeqGene.txt")
#chr1 <- subset(genes,genes[,1]==1)
#genome.coords <- cbind(chr1[,2:3]) 
#bsize <- 2000
#chrSize <- 267910886
#ret <- genomedata.enrichment(posteriors,genome.coords,bsize,chrSize)
#ret

if (F) {
#Testing expression.enrichment
posteriors <- read.table("example/ex1-posteriors.txt",header=T)
posteriors
tempexpr <- read.table("example/ex1-expression.txt")
chr1 <- subset(tempexpr,tempexpr[,1]=="chr1")
expr <- cbind(tempexpr[2:4])
ret <- expression.enrichment(posteriors,expr,10,25)
ret
}
