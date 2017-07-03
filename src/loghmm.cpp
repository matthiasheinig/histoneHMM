#include "hmm.h"

/* initialize the model */
LogHMM::LogHMM(double* O, int T, int N) {
    cout<<"LogHMM::LogHMM(double* O, int T, int N) T = "<<T<<" N = "<<N<<"\n";
    this->O = O;
    this->T = T;
    this->N = N;
    this->A = allocDoubleMatrix(N, N);
    this->logalpha = allocDoubleMatrix(T, N);
    this->logbeta = allocDoubleMatrix(T, N);
    this->logdensities = allocDoubleMatrix(T, N);
    this->logproba = (double*) calloc(N, sizeof(double));
    this->zn2zz = false;

cout<<"O="<<endl;
for(int t=0;t<10;t++){
cout<<"O["<<t<<"]="<<O[t]<<endl;
};
    
    /* INITIALIZATION for logproba and transition matrix*/
    /* we set initial probabilities to uniform and make self transitions more
     likely than transitions to other states */
    double self = 0.9;
    double other = (1.0 - self) / (N - 1.0);
    for (int iN=0; iN<N; iN++)
    {
        this->logproba[iN] = log((double)1/N);
        cout << "proba["<<iN<<"]= " << exp(logproba[iN]) << "; A["<<iN<<"]= ";
        for (int jN=0; jN<N; jN++)
        {
            if (iN == jN)
                A[iN][jN] = self;
            else
                A[iN][jN] = other;
            cout << A[iN][jN] << ",";
        }
        cout<<endl;
    }
    this->Nmod = 1;
}

LogHMM::LogHMM(double* O, int T, int N, int Nmod)
{
    cout<<"LogHMM::LogHMM(double* O, int T, int N) T = "<<T<<" N = "<<N<<"\n";
    this->O = O;
    this->T = T;
    this->N = N;
    this->A = allocDoubleMatrix(N, N);
    this->logalpha = allocDoubleMatrix(T, N);
    this->logbeta = allocDoubleMatrix(T, N);
    this->logdensities = allocDoubleMatrix(T, N);
    this->logproba = (double*) calloc(N, sizeof(double));
    this->zn2zz = false;
    this->Nmod = Nmod;
    
    /* INITIALIZATION for logproba and transition matrix*/
    /* we set initial probabilities to uniform and make self transitions more
     likely than transitions to other states */
    double self = 0.9;
    double other = (1.0 - self) / (N - 1.0);
    for (int iN=0; iN<N; iN++)
    {
        this->logproba[iN] = log((double)1/N);
//        cout << "proba["<<iN<<"]= " << exp(logproba[iN]) << "; A["<<iN<<"]= ";
        for (int jN=0; jN<N; jN++)
        {
            if (iN == jN)
                A[iN][jN] = self;
            else
                A[iN][jN] = other;
//            cout << A[iN][jN] << ",";
        }
//        cout<<endl;
    }
  // Somehow this does not work like this!!
  
 // cout<<"LogHMM::LogHMM(double* O, int T, int N, int Nmod) T = "<<T<<" N = "<<N<<" Nmod = "<<Nmod<<"\n";
  //LogHMM(O, T, N);
  //this->Nmod = Nmod;
}

/* clean up a model (do not touch the input data) */
LogHMM::~LogHMM() {
    freeDoubleMatrix(this->A, this->N);
    freeDoubleMatrix(this->logalpha, this->T);
    freeDoubleMatrix(this->logbeta, this->T);
    freeDoubleMatrix(this->logdensities, this->T);
    free(this->logproba);
}

/* helper functions */
double LogHMM::logSumGamma(int iN){
    double s1, s2, temp;
	double *tempvec = (double*) calloc(this->N, sizeof(double));
    int ti, j;
    s1=0.0;
    for(ti=0;ti<T-1;ti++){
        for(j=0;j<N;j++){
            tempvec[j]=this->logalpha[ti][j]+this->logbeta[ti][j];//=-inf if proba=0
        };
        temp=Max(tempvec,this->N);
        s2=0.0;
        for(j=0;j<N;j++){
            s2+=exp(this->logalpha[ti][j]+this->logbeta[ti][j]-temp);
        };
        s1+=exp(this->logalpha[ti][iN]+this->logbeta[ti][iN]-temp-log(s2));
    };
	free(tempvec);
    return log(s1);
}

double LogHMM::logSumEta(int iN, int jN) {
    int t;
    double **tempMatrix = allocDoubleMatrix(this->N, this->N);
    double temp;
    double *AnewDeno = (double*) calloc(this->T, sizeof(double));
    for(t=0;t<T-1;t++){
        for(int i=0;i<this->N;i++){
            for(int j=0;j<this->N;j++){
                tempMatrix[i][j]=this->logalpha[t][i]+log(A[i][j])+this->logdensities[t+1][j]+this->logbeta[t+1][j];
            };
        };
        temp=MaxMatrix(tempMatrix, this->N, this->N);
        AnewDeno[t]=0.0;
        for(int i=0;i<this->N;i++){
            for(int j=0;j<this->N;j++){
                AnewDeno[t]+=exp(this->logalpha[t][i]+log(A[i][j]) + this->logdensities[t+1][j]+this->logbeta[t+1][j]-temp);
            };
        };
        AnewDeno[t]=temp+log(AnewDeno[t]);
    };
    double sum1;
    sum1=0.0;
    for(int ti=0;ti<T-1;ti++){
        sum1+=exp(this->logalpha[ti][iN]+log(A[iN][jN])+this->logdensities[ti+1][jN]+this->logbeta[ti+1][jN]-AnewDeno[ti]);//-log(sum2));
    };
	free(AnewDeno);
    freeDoubleMatrix(tempMatrix,this->N);
    return log(sum1);
}


void LogHMM::forward() {
    int iN, jN, t;
    double sum, temp;
	double *tempvec = (double*) calloc(this->N, sizeof(double));
    for(iN=0;iN<this->N;iN++){
        this->logalpha[0][iN]= this->logproba[iN]+this->logdensities[0][iN];//careful: if proba[iN]=0, then logalpha=-inf
    };
    for(t=1;t<this->T;t++){
        for(iN=0;iN<this->N;iN++){
            for(jN=0;jN<this->N;jN++){//store in a vector the values for logalpha[t-1][jN]+log(A[jN][iN])
                tempvec[jN]=this->logalpha[t-1][jN]+log(A[jN][iN]);//if proba[iN]=0, then tempvec=-inf. Same problem if A[jN][iN]=0
            };
            temp=Max(tempvec,this->N);
            sum=0.0;
            for(jN=0;jN<this->N;jN++){
                sum+=exp(this->logalpha[t-1][jN]+log(A[jN][iN])-temp);
            };
            this->logalpha[t][iN]=temp+log(sum)+this->logdensities[t][iN];
if(isnan(this->logalpha[t][iN])){
cout<<"logalpha["<<t<<"]["<<iN<<"]= "<<temp<<"+log("<<sum<<")+"<<logdensities[t][iN]<<endl;
exit(1);
};
        };
    };

/*for(t=0;t<10;t++){
for(iN=0;iN<this->N;iN++){
cout<<"logalpha["<<t<<"]["<<iN<<"]="<<this->logalpha[t][iN]<<"\t";
};cout<<endl;};
for(t=T-10;t<T;t++){
for(iN=0;iN<this->N;iN++){
cout<<"logalpha["<<t<<"]["<<iN<<"]="<<this->logalpha[t][iN]<<"\t";
};cout<<endl;};
*/
	free(tempvec);
}

void LogHMM::backward() {
    int iN, jN, t;
    double sum, temp;
	double *tempvec = (double*) calloc(this->N, sizeof(double));
    for(iN=0;iN<this->N;iN++){
        this->logbeta[T-1][iN]=0.0;//=log(1)
    };
    for(t=this->T-2;t>=0;t--){
        for(iN=0;iN<this->N;iN++){
            for(jN=0;jN<this->N;jN++){
                tempvec[jN]=log(A[iN][jN])+this->logdensities[t+1][jN]+this->logbeta[t+1][jN];
            };
            temp=Max(tempvec,this->N);
            sum=0.0;
            for(jN=0;jN<this->N;jN++){
                sum+=exp(log(A[iN][jN])+this->logdensities[t+1][jN]+this->logbeta[t+1][jN]-temp);
            };
            this->logbeta[t][iN]=temp+log(sum);
        };
    };

/*for(t=0;t<10;t++){
for(iN=0;iN<this->N;iN++){
cout<<"logbeta["<<t<<"]["<<iN<<"]="<<this->logbeta[t][iN]<<"\t";
};cout<<endl;};
for(t=T-10;t<T;t++){
for(iN=0;iN<this->N;iN++){
cout<<"logbeta["<<t<<"]["<<iN<<"]="<<this->logbeta[t][iN]<<"\t";
};cout<<endl;};
*/
	free(tempvec);
}

double LogHMM::logLikelihood(){
        double sum=0.0;
        double temp;
	double *tempvec = (double*) calloc(this->N, sizeof(double));
        int iN;
        for(iN=0;iN<this->N;iN++){
                tempvec[iN]=this->logalpha[this->T-1][iN];
        };
        temp=Max(tempvec,this->N);
        for(iN=0;iN<this->N;iN++){
                sum+=exp(this->logalpha[this->T-1][iN]-temp);
        };
	free(tempvec);
        return temp+log(sum);
};


int LogHMM::baumWelch(int iterationMAX, double eps) {
    /* arbitrary initialization was done in the constructor */
    double logPold = -INFINITY;
    double logPnew, smallnumber;
    double sum, temp, temp2, probaDeno;
    double* tempvector = (double*) calloc(N, sizeof(double));
    /* useful number */
    double* logSUMofGamma = (double*) calloc(this->N, sizeof(double));
    double** Anew = allocDoubleMatrix(this->N, this->N);

    // Minh Anh, I changed the following allocation (MH)
    //double** Tauti = (double*) calloc(this->T, sizeof(double)); //weights to update Denstiy's paramaters
    double** Tauti = allocDoubleMatrix(this->N, this->T);

    int iN, jN, k, t, iteration;
    
    int countSwitch = 0;
    
    for(iteration=1; iteration <= iterationMAX; iteration++)
    {
        R_CheckUserInterrupt(); // check interupt
        cout << "iteration = " << iteration ;
        //CALCULATE AND STORE THE DENSITIES
        computeDensities();

        /* compute forward probabilities */
        forward();
        
        //CALCULATE THE LOGLIKELIHOOD
		logPnew = this->logLikelihood();
		if(isnan(logPnew)) break;
        smallnumber = logPnew - logPold; /* fabs(logPnew - logPold); */
//	cout<<"\t loglik = "<<logPnew<<endl;
        cout<<"************************************************************************************************ delta logP = " <<smallnumber<<"\n";
       
		/* check convergence */
        if(fabs(smallnumber) < eps) //it has converged
        {
            if (this->Nmod == 1) //univariate HMM, we need to check the mean of the two state 0 and 1
            {
                cout<<"Converged in univariate model, means: mu(0) = "<<this->densityFunctions[0]->getMean()<<" , mu(1) = "<<this->densityFunctions[1]->getMean()<<endl;
                if(this->densityFunctions[1]->getMean() < this->densityFunctions[0]->getMean())//wrong state order
                {
                    printf("Wrong state order\n");
                    if (this->densityFunctions[0]->getType() == NB && this->densityFunctions[1]->getType() == NB)
                    {
cout<<"means need to be exchanged, type NB/NB"<<endl;
                        printf("Swapping states\n");
                        NegativeBinomial *tempDens = new NegativeBinomial();
                        tempDens->copy(this->densityFunctions[1]); // tempDens is densifunc[1]
                        this->densityFunctions[1]->copy(this->densityFunctions[0]); 
                        this->densityFunctions[0]->copy(tempDens); 

                        //swap logproba and transition matrix
                        temp=this->logproba[0];
                        this->logproba[0]=this->logproba[1];
                        this->logproba[1]=temp;
                        temp=this->A[0][0];
                        this->A[0][0]=this->A[1][1];
                        this->A[1][1]=temp;
                        temp=this->A[0][1];
                        this->A[0][1]=this->A[1][0];
                        this->A[1][0]=temp;

                        // recompute densities and forward backward
                        computeDensities();
                        forward();
                        backward();                        
                        break;
                    }
                    else if (this->densityFunctions[0]->getType() == Z && this->densityFunctions[1]->getType() == Z)
                    {
cout<<"means need to be exchanged, type Z/Z"<<endl;
                        ZiNB *tempDens = new ZiNB();
                        tempDens->copy(this->densityFunctions[1]); // tempDens is densifunc[1]
                        this->densityFunctions[1]->copy(this->densityFunctions[0]); 
                        this->densityFunctions[0]->copy(tempDens); 

                        //swap logproba and transition matrix
                        temp=this->logproba[0];
                        this->logproba[0]=this->logproba[1];
                        this->logproba[1]=temp;
                        temp=this->A[0][0];
                        this->A[0][0]=this->A[1][1];
                        this->A[1][1]=temp;
                        temp=this->A[0][1];
                        this->A[0][1]=this->A[1][0];
                        this->A[1][0]=temp;

                        computeDensities();
                        forward();
                        backward();                
                        break;
                    }
                    else //ZN
                    {
cout << "need to restart estimation, means need to be exchanged byt type Z/NB; mu(0)="<<this->densityFunctions[0]->getMean()<<">mu(1)="<<this->densityFunctions[1]->getMean()<<endl;
                        countSwitch++;
                        if (countSwitch > 2)
                        {
                            this->zn2zz = true;
                            break;
                        }
                        logPold=-INFINITY;
                        // zinba should always be the smaller state
                        ZiNB* zinb_density = (ZiNB*)(this->densityFunctions[0]);
                        NegativeBinomial* nb_density = (NegativeBinomial*)(this->densityFunctions[1]);

                        // swap parameter R
                        temp = zinb_density->getR();
                        zinb_density->setR(nb_density->getR());
                        nb_density->setR(temp);

                        // swap parameter P
                        temp = zinb_density->getP();
                        zinb_density->setP(nb_density->getP());
                        nb_density->setP(temp);
                    }
                }
                else //the means are okay
                    break;
            }
        }
        else //not converged
        {
            logPold = logPnew;
        }
        
        /* compute backward probabilities */
        backward();
        
        /* UPDATES */
        //useful numbers
        for(iN=0; iN<this->N; iN++)
        {
            tempvector[iN]=this->logalpha[0][iN]+this->logbeta[0][iN];//is -inf if proba=0
            logSUMofGamma[iN]=logSumGamma(iN);
        }
        temp = Max(tempvector, N);
        probaDeno = 0.0;
        for(iN=0;iN<this->N;iN++)
        {
            probaDeno+=exp(this->logalpha[0][iN]+this->logbeta[0][iN]-temp);
        }
	probaDeno=temp+log(probaDeno);
	//updates
	for (iN=0; iN<this->N; iN++)
	  {   //update logproba
	    this->logproba[iN] = this->logalpha[0][iN]+this->logbeta[0][iN]-probaDeno;
	    for(jN=0; jN<this->N; jN++)
	      {   //update the A
		Anew[iN][jN]=exp(logSumEta(iN,jN)-logSUMofGamma[iN]);
	      }
	  }
	
	/* copy the new A to the model */
        for(iN=0; iN < this->N; iN++) {
	  for(jN=0; jN < this->N; jN++) {
	    this->A[iN][jN] = Anew[iN][jN];
	  }
        }
	
/*cout<<"proba=( "<<exp(this->logproba[0])<<" , "<<exp(this->logproba[1])<<" )"<<endl;
cout<<"A=( ("<<A[0][0]<<","<<A[0][1]<<") , ("<<A[1][0]<<","<<A[1][1]<<") )"<<endl;*/

	// call posterior with tauti and recompute=0
	posterior(Tauti, 0);

/*for(t=0;t<10;t++){
cout<<"tauti["<<t<<"]=("<<Tauti[0][t]<<","<<Tauti[1][t]<<")"<<endl;
};
for(t=T-2;t<T;t++){
cout<<"tauti["<<t<<"]=("<<Tauti[0][t]<<","<<Tauti[1][t]<<")"<<endl;
};*/

	for (iN=0; iN<this->N; iN++) {
//	cout<<"iN = "<<iN<<"\t";
	  //update the parameters of the distribution
	  this->densityFunctions[iN]->update(Tauti[iN], this->T);
	}
    } /* main loop end */
    
    
    //Print the last results
    cout<<endl<<endl<<"---------------------------------------------------------"<<endl;
    cout<<"| FINAL ESTIMATION RESULTS \t\t\t\t |"<<endl;
    for(iN=0;iN<N;iN++){
        if(iN==0) cout<<"| unmodified component: \t\t\t\t |"<<endl;
        if(iN==1) cout<<"| modified component: \t\t\t\t\t |"<<endl;
        cout<<"|";
        cout<<" proba["<<iN<<"]="<<exp(this->logproba[iN])<<"\t\t\t\t\t\t |"<<endl;
        cout<<"|";
        for(jN=0;jN<N;jN++){
            cout<<" A["<<iN<<"]["<<jN<<"]="<<this->A[iN][jN]<<"\t";
        };
        cout<<"\t |"<<endl;
        //cout << densityFunctions[iN] << "\t\t\t\t\t |" << endl;
        //cout << this->densityFunctions[iN];
        double curR, curP, curW;
        if (this->densityFunctions[iN]->getType() == Z)
        {
            ZiNB* temp = (ZiNB*)this->densityFunctions[iN];
            curR = temp->getR();
            curP = temp->getP();
            curW = temp->getW();
        }
        if (this->densityFunctions[iN]->getType() == NB)
        {
            NegativeBinomial* temp = (NegativeBinomial*)this->densityFunctions[iN];
            curR = temp->getR();
            curP = temp->getP();
            curW = 0;
        }
        cout<<"| r["<<iN<<"]="<<curR<<"\t\t\t\t\t\t |"<<endl;
		cout<<"| p["<<iN<<"]="<<curP<<"\t\t\t\t\t |"<<endl;
        cout<<"| w["<<iN<<"]="<<curW<<"\t\t\t\t\t\t |"<<endl;
    };
    cout<<"---------------------------------------------------------"<<endl;
    
    /* free memory */
    free(logSUMofGamma);
    freeDoubleMatrix(Anew, this->N);
    free(tempvector);
    freeDoubleMatrix(Tauti,this->N);	
    return iteration;
}

int LogHMM::estimateTransitions(int iterationMAX, double eps) {
    /* arbitrary initialization was done in the constructor */
    double logPold = -INFINITY;
    double logPnew, smallnumber;
    double sum, temp, temp2, probaDeno;
    double* tempvector = (double*) calloc(N, sizeof(double));
    /* useful number */
    double* logSUMofGamma = (double*) calloc(this->N, sizeof(double));
    double** Anew = allocDoubleMatrix(this->N, this->N);
    double** Tauti = allocDoubleMatrix(this->N, this->T); //weights to update Denstiy's paramaters
    
    int iN, jN, k, t, iteration;
    
    int countSwitch = 0;
    
    for(iteration=1; iteration <= iterationMAX; iteration++)
    {
      R_CheckUserInterrupt(); // check interupt
        cout << "iteration = " << iteration;
	//CALCULATE AND STORE THE DENSITIES
	computeDensities();
/*for(t=0;t<10;t++){
cout<<"logd=( ";
for(iN=0;iN<N;iN++){cout<<this->logdensities[t][iN]<<" ,";};
cout<<")"<<endl;
};
for(t=T-10;t<T;t++){
cout<<"logd=( ";
for(iN=0;iN<N;iN++){cout<<this->logdensities[t][iN]<<" ,";};
cout<<")"<<endl;
};*/

	/* compute forward probabilities */
        forward();

		//CALCULATE THE LOGLIKELIHOOD
		logPnew = this->logLikelihood();
		if(isnan(logPnew)) {cout<<"!!!!! loglik = "<<logPnew<<" ----> BREAK"<<endl; break;};
        smallnumber = logPnew - logPold; /* fabs(logPnew - logPold); */
	cout<<"\t loglik = "<<logPnew<<endl;
        cout<<"************************************************************************************************ delta logP = " <<smallnumber<<"\n";
        
		/* check convergence */
        if(fabs(smallnumber) < eps) //it has converged
        {
	  break;
        }
        else //not converged
        {
            logPold = logPnew;
        }
        
        /* compute backward probabilities */
        backward();
        
        /* UPDATES */
		//useful numbers
		for(iN=0; iN<this->N; iN++)
        {
			tempvector[iN]=this->logalpha[0][iN]+this->logbeta[0][iN];//is -inf if proba=0
			logSUMofGamma[iN]=logSumGamma(iN);
		}
		temp = Max(tempvector, N);
		probaDeno = 0.0;
        for(iN=0;iN<this->N;iN++)
        {
			probaDeno+=exp(this->logalpha[0][iN]+this->logbeta[0][iN]-temp);
		}
		probaDeno=temp+log(probaDeno);
		//updates
		for (iN=0; iN<this->N; iN++)
        {   //update logproba
			this->logproba[iN] = this->logalpha[0][iN]+this->logbeta[0][iN]-probaDeno;
			for(jN=0; jN<this->N; jN++)
            {   //update the A
				Anew[iN][jN]=exp(logSumEta(iN,jN)-logSUMofGamma[iN]);
			}
		}
		/* copy the new A to the model */
        for(iN=0; iN < this->N; iN++) {
            for(jN=0; jN < this->N; jN++) {
                this->A[iN][jN] = Anew[iN][jN];
            }
        }
        
    } /* main loop end */
    
    //Print the last results
    cout<<endl<<endl<<"---------------------------------------------------------"<<endl;
    cout<<"| FINAL ESTIMATION RESULTS \t\t\t\t |"<<endl;
    for(iN=0;iN<N;iN++){
        if(iN==0) cout<<"| unmodified component: \t\t\t\t |"<<endl;
        if(iN==1) cout<<"| modified component: \t\t\t\t\t |"<<endl;
        cout<<"|";
        cout<<" proba["<<iN<<"]="<<exp(this->logproba[iN])<<"\t\t\t\t\t\t |"<<endl;
        cout<<"|";
        for(jN=0;jN<N;jN++){
            cout<<" A["<<iN<<"]["<<jN<<"]="<<this->A[iN][jN]<<"\t";
        };
        cout<<"\t |"<<endl;
        cout << densityFunctions[iN] << "\t\t\t\t\t |" << endl;
        cout << this->densityFunctions[iN];
		//<<"| p["<<iN<<"]="<<this->densityFunctions[iN]->p<<"\t\t\t\t\t |"<<endl;
        //cout<<"| r["<<iN<<"]="<<this->densityFunctions[iN]->r<<"\t\t\t\t\t\t |"<<endl;
    }
    cout<<"---------------------------------------------------------"<<endl;
    
    /* free memory */
    free(logSUMofGamma);
    freeDoubleMatrix(Anew, this->N);
    free(tempvector);
    freeDoubleMatrix(Tauti,this->N);	
    return iteration;
}


void LogHMM::computeDensities() {
  for(int t=0; t<this->T; t++)
  {
      for(int iN=0; iN<this->N; iN++)
      {
          this->logdensities[t][iN] = this->densityFunctions[iN]->logdensity(t);
          //this->logdensities[t][iN] = this->densityFunctions[iN]->logdensity(t*Nmod);
      }
  }
/*for(int t=0;t<10;t++){
for(int iN=0;iN<this->N;iN++){
cout<<"logd["<<t<<"]["<<iN<<"]="<<this->logdensities[t][iN]<<"\t";
};cout<<endl;};
for(int t=T-10;t<T;t++){
for(int iN=0;iN<this->N;iN++){
cout<<"logd["<<t<<"]["<<iN<<"]="<<this->logdensities[t][iN]<<"\t";
};cout<<endl;};*/
}

/* for an initialized model compute the posterior probs
 posterior is a matrix [N x T] */
void LogHMM::posterior(double** post, int recompute) {
	double sum, temp;
	double *tempVec = (double*) calloc(this->N, sizeof(double));
    if (recompute) {
       // compute forward and backward variables
      computeDensities();
      forward();
      backward();
    }
	for(int t=0;t<this->T;t++)
    {
		sum=0.0;
	    for(int jN=0;jN<this->N;jN++)
        {
			tempVec[jN]=this->logalpha[t][jN]+this->logbeta[t][jN];
        }
	    temp=Max(tempVec,this->N);
        for(int jN=0;jN<this->N;jN++)
        {
			sum+=exp(this->logalpha[t][jN]+this->logbeta[t][jN]-temp);
		}
		for(int iN=0;iN<this->N;iN++)
        {
			tempVec[iN]=exp(this->logalpha[t][iN]+this->logbeta[t][iN]-temp)/sum;
			post[iN][t] = exp(this->logalpha[t][iN]+this->logbeta[t][iN]-temp)/sum;
		}
	}
    free(tempVec);
}

/* for an initialized model compute the viterbi path */
void LogHMM::viterbi(int* path, int recompute) {
    if (recompute) {
        for(int t=0; t<this->T; t++)
        {
            for(int iN=0; iN<this->N; iN++)
            {
                this->logdensities[t][iN] = this->densityFunctions[iN]->logdensity(t);
            }
        }
        /* compute forward and backward variables */
        forward();
        backward();
    }
    
    /* do the actual viterbi */
    
}

/*
 INPUT:
 O: observation size [T]
 T: size of the chromosomes, length = 1
 N: number of states, length = 1
 interationMAX: maximum number of iterations for the HMM, length = 1
 eps: tolerant number
 post: posterior, transposed matrix of matrix [T*N]
 OUTPUT:
 r: r parameters of the univariate distributions for every modifications, length = Nmod * 2
 p: p parameters of the univariate distributions for every modifications, length = Nmod * 2
 w: w parameters of the univariate distributions for every modifications, length = Nmod * 2
 A: transition probability matrix [N*N]
 proba: initial probability, length = N
 loglik: loglikelihood
 */
extern "C" {
  void R_univariate_hmm(double* O, int* T, int* N, int* densityNames, double* r, double* p, double* w, int* iterationMAX, double* eps, double* post, double* A, double* proba, double* loglik) {
  //void R_univariate_hmm(double* O, int* T, int* N, double* r, double* p, double* w, int* iterationMAX, double* eps, double* post, double* A, double* proba, double* loglik) {
    printf("seqlen %i, nstates %i, maxit %i, eps %f\n", *T, *N, *iterationMAX, *eps);

    /* initialize */
    printf("init..\n");
    LogHMM* model = new LogHMM(O, *T, *N);
  //  model->Nmod = 1;
    printf("model->T = %i\nmodel->N = %i\n", model->T, model->N);
      
  //  int i, j, t;
            
    if (densityNames[0] == NB && densityNames[1] == NB) //both are Negative Binomial
    {
        printf("using for both states the Negative Binomial distribution\n");
        for (int i=0; i<model->N; i++)
        {
            Density *d = new NegativeBinomial(O,(double)i+0.5,0.01);
            model->densityFunctions.push_back(d);            
        }
        /* estimate the parameters */
        printf("estimate parameters..\n");
        model->baumWelch(*iterationMAX, *eps);
        printf("OK\n");
    }//end_if
    else if (densityNames[0] == Z && densityNames[1] == NB) //ZiNB and NB
    {
        printf("using ZiNB and Negative Binomial.\n We do a grid search for w - the weight of the ZiNB distribution.\n");
        int maxNW = 10;
        double wList[maxNW];
        /*for (int iN=0; iN<5; iN++)
            wList[iN] = (iN+1) * 0.0002;
        for (int iN=5; iN<maxNW; iN++)
            wList[iN] = (iN-5+1) * 0.05;*/
	for(int iN=0;iN<maxNW;iN++)
	    wList[iN] = (iN) * 0.05;
        // double loglik[maxNW], rList[maxNW][N], pList[maxNW][N], Alist[maxNW][N][N], probaList[maxNW][N];
        double *loglikZN = (double*) calloc(maxNW, sizeof(double));
        int *lastIteListZN = (int*) calloc(maxNW, sizeof(int));
        double **rListZN = allocDoubleMatrix(maxNW, model->N);
        double **pListZN = allocDoubleMatrix(maxNW, model->N);
        double ***AlistZN = alloc3Ddouble(maxNW, model->N, model->N);
        double **probaListZN = allocDoubleMatrix(maxNW, model->N);
        for (int wID = 0; wID < maxNW; wID++) //for_wList
        {
	    printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  weight = %g \n", wList[wID]); 
            model->densityFunctions.clear();
            int i;
            for (i = 0; i < model->N; i++) //for_N_1
            {
                if (i == 0)
                {
                    ZiNB *d;
                    if(wID>0)
                    {
                        if(isnan(loglikZN[wID-1]))
                            d = new ZiNB(model->O,1.5,0.15,wList[wID]);
                        else
                            d = new ZiNB(model->O,rListZN[wID-1][i],pListZN[wID-1][i],wList[wID]);
                    }
                    else
                        d = new ZiNB(model->O,1.5,0.15,wList[wID]);
                    
                    model->densityFunctions.push_back(d);
                }
                if (i == 1)
                {
                    NegativeBinomial *d;
                    if(wID>0)
                    {
                        if(isnan(loglikZN[wID-1]))
                            d = new NegativeBinomial(model->O,3.0,0.01);
                        else
                            d = new NegativeBinomial(model->O,rListZN[wID-1][i],pListZN[wID-1][i]);
                    }
                    else
                        d = new NegativeBinomial(model->O,3.0,0.01);
                    
                    model->densityFunctions.push_back(d);
                }
            }//end_for_N_1
            /* estimate the parameters */
            printf("estimate parameters..\n");
            model->baumWelch(*iterationMAX, *eps);
            printf("OK\n");
            loglikZN[wID] = model->logLikelihood();
            // lastIteListZN[wID] = lastIteration;
            for( int iN=0; iN < model->N; iN++) //for_N_2
            {
             //   rListZN[wID][iN] = model->densityFunctions[iN]->getFirstParam();
              //  pListZN[wID][iN] = model->densityFunctions[iN]->getSecondParam();
                if (model->densityFunctions[iN]->getType() == NB)
                {
                    NegativeBinomial* temp = (NegativeBinomial*) model->densityFunctions[iN];
                    rListZN[wID][iN] = temp->getR();
                    pListZN[wID][iN] = temp->getP();
                }
                if (model->densityFunctions[iN]->getType() == Z)
                {
                    ZiNB* temp = (ZiNB*) model->densityFunctions[iN];
                    rListZN[wID][iN] = temp->getR();
                    pListZN[wID][iN] = temp->getP();
                }
                probaListZN[wID][iN] = model->logproba[iN];
                for (int jN = 0; jN < model->N; jN++)
                    AlistZN[wID][iN][jN] = model->A[iN][jN];
            }//end_for_N_2
        }//end_for_wList
        //set the optimal parameters to the model
        int curArgMaxZN = argMax(loglikZN,maxNW);
        for (int iN=0; iN < model->N; iN++)
        {
            if (model->densityFunctions[iN]->getType() == NB)
            {
                NegativeBinomial* tempDens = (NegativeBinomial*) model->densityFunctions[iN];
                tempDens->setR(rListZN[curArgMaxZN][iN]);
                tempDens->setP(pListZN[curArgMaxZN][iN]);
            }
            if (model->densityFunctions[iN]->getType() == Z)
            {
                ZiNB* tempDens = (ZiNB*) model->densityFunctions[iN];
                tempDens->setR(rListZN[curArgMaxZN][iN]);
                tempDens->setP(pListZN[curArgMaxZN][iN]);
                tempDens->setW(wList[curArgMaxZN]);
            }            
          
            model->logproba[iN] = probaListZN[curArgMaxZN][iN];
            for (int jN = 0; jN < model->N; jN++)
                model->A[iN][jN] = AlistZN[curArgMaxZN][iN][jN];
        }
        //model->posterior(post, postMax,true);
        free(loglikZN);
        free(lastIteListZN);
        freeDoubleMatrix(rListZN, maxNW);
        freeDoubleMatrix(pListZN, maxNW);
        freeDoubleMatrix(probaListZN, maxNW);
        free3Ddouble(AlistZN,maxNW, model->N);
    }//end_else
    if ( (densityNames[0] == Z && densityNames[1] == Z) || model->zn2zz == true ) //do ZiNB and ZiNB
    {
        printf("using ZiNB and ZiNB.\n We do a grid search for w - the weight of the ZiNB distribution.\n");        
        int maxNW = 10;
        double wList[maxNW];
        for (int iN=0; iN<5; iN++)
            wList[iN] = (iN+1) * 0.0002;
        for (int iN=5; iN<maxNW; iN++)
            wList[iN] = (iN-5+1) * 0.05;
        double *loglikZZ = (double*) calloc(maxNW*maxNW,sizeof(double));
        int *lastIteListZZ = (int*) calloc(maxNW*maxNW,sizeof(int));
        double **rListZZ = allocDoubleMatrix(maxNW*maxNW,model->N);
        double **pListZZ = allocDoubleMatrix(maxNW*maxNW,model->N);
        double ***AlistZZ = alloc3Ddouble(maxNW*maxNW,model->N,model->N);
        double **probaListZZ = allocDoubleMatrix(maxNW*maxNW,model->N);
        for (int wID = 0; wID < maxNW; wID++)
        {
            for (int w2 = 0; w2 < maxNW; w2++)
            {
		printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  weights = (%g , %g )", wList[wID],wList[w2]);
                model->densityFunctions.clear();
                int i;
                for (i = 0; i < model->N; i++)
                {
                    if (i == 0)
                    {
                        ZiNB *d;
                        if(wID*maxNW+w2>0)
                        {
                            if(isnan(loglikZZ[wID*maxNW+w2-1]))
                                d = new ZiNB(model->O,1.5,0.15,wList[wID]);
                            else
                                d = new ZiNB(model->O,rListZZ[wID*maxNW+w2-1][i],pListZZ[wID*maxNW+w2-1][i],wList[wID]);
                        }
                        else
                            d = new ZiNB(model->O,1.5,0.15,wList[wID]);
                        model->densityFunctions.push_back(d);
                    }
                    if (i == 1)
                    {
                        ZiNB *d;
                        if(wID*maxNW+w2>0)
                        {
                            if(isnan(loglikZZ[wID*maxNW+w2-1]))
                                d = new ZiNB(model->O,3.0,0.01,wList[w2]);
                            else
                                d = new ZiNB(model->O,rListZZ[wID*maxNW+w2-1][i],pListZZ[wID*maxNW+w2-1][i],wList[w2]);
                        }
                        else
                            d = new ZiNB(model->O,3.0,0.01,wList[w2]);
                        model->densityFunctions.push_back(d);
                    }
                }
                /* estimate the parameters */
                printf("estimate parameters..\n");
                model->baumWelch(*iterationMAX, *eps);
                printf("OK\n");
                loglikZZ[wID*maxNW+w2] = model->logLikelihood();
                for( int iN=0; iN < model->N; iN++)
                {
                    ZiNB *tempDens = (ZiNB*) model->densityFunctions[iN];
                    rListZZ[wID*maxNW+w2][iN] = tempDens->getR();
                    pListZZ[wID*maxNW+w2][iN] = tempDens->getP();
                    probaListZZ[wID*maxNW+w2][iN] = model->logproba[iN];
                    for (int jN = 0; jN < model->N; jN++)
                        AlistZZ[wID*maxNW+w2][iN][jN] = model->A[iN][jN];
                }
            }//end_for_wList (2)
        }//end_for_wList (1)
        int curArgMaxZZ = argMax(loglikZZ,maxNW*maxNW);
        for (int iN=0; iN < model->N; iN++)
        {
            ZiNB *tempDens = (ZiNB*) model->densityFunctions[iN];
            tempDens->setR(rListZZ[curArgMaxZZ][iN]);
            tempDens->setP(pListZZ[curArgMaxZZ][iN]);
            tempDens->setW(wList[curArgMaxZZ]);
            model->logproba[iN] = probaListZZ[curArgMaxZZ][iN];
            for (int jN = 0; jN < model->N; jN++)
                model->A[iN][jN] = AlistZZ[curArgMaxZZ][iN][jN];
        }
        free(loglikZZ);
        free(lastIteListZZ);
        freeDoubleMatrix(rListZZ, maxNW*maxNW);
        freeDoubleMatrix(pListZZ, maxNW*maxNW);
        freeDoubleMatrix(probaListZZ, maxNW*maxNW);
        free3Ddouble(AlistZZ,maxNW*maxNW, model->N);
    }//end_if
      
    /* initialize the distributions */
    
/*    for (i=0; i<model->N; i++)
    {
        Density* d;
        //if (densityNames[i] == Z)
        if (*w != 0.0)
        {
            printf("using zinba\n");
            // this is some stupid initialization just guessing..
            d = new ZiNB(O, (double)i + 0.5, 0.01, *w);
        }
        else
        {
            printf("using negbinom\n");
            // this is some stupid initialization just guessing..
            d = new NegativeBinomial(O, (double)i + 0.5, 0.01);
        }
        model->densityFunctions.push_back(d);
    }    
    printf("OK\n");      
    // estimate the parameters 
    printf("estimate parameters..\n");
    model->baumWelch(*iterationMAX, *eps);
    printf("OK\n");*/

    /* compute the posteriors and save results directly to the R pointer */
    /* we do need to recompute the forward and backward variables especially for case Z_NB and Z_Z*/
    printf("compute posterior.."); 
    double** post_matrix = allocDoubleMatrix(model->N, model->T);
    model->posterior(post_matrix, true);
    /* recode into column representation */
      int i, j, t;
    for (t=0; t<model->T; t++)
    {
      for (i=0; i<model->N; i++)
      {
          post[t + i * model->T] = post_matrix[i][t];
      }
    }
cout<<"after convergence, the posteriors are:"<<endl;
for(t=0;t<10;t++){
cout<<"post["<<t<<"][0]="<<post_matrix[0][t]<<", "<<post_matrix[1][t]<<endl;
};
    freeDoubleMatrix(post_matrix, model->N);

    /* also return the estimated transition matrix and the initial probs */
    for (i=0; i<model->N; i++)
    {
      proba[i] = exp(model->logproba[i]);
      for (j=0; j<model->N; j++)
      {
          A[i + j * model->N] = model->A[i][j];
      }
    }

    /* copy the estimated distribution params */
    for (i=0; i<model->N; i++)
    {
      if (model->densityFunctions[i]->getType() == Z)
        {
            printf("reading zinba parameters\n");
            ZiNB* d = (ZiNB*)(model->densityFunctions[i]);
            r[i] = d->getR();
            p[i] = d->getP();
            w[i] = d->getW();
        }
      if (model->densityFunctions[i]->getType() == NB) 
        {
            printf("reading negative binomial parameters\n");
            NegativeBinomial* d = (NegativeBinomial*)(model->densityFunctions[i]);
            r[i] = d->getR();
            p[i] = d->getP();
            w[i] = 0;
        }
    }
    *loglik = model->logLikelihood();
    printf("OK\n");
//cout to debug
for (int iN = 0; iN < model->N; iN++)
	cout << "State(r,p,w) = (" << r[iN] <<", " << p[iN] << ", " << w[iN] << ")" << endl;
  }
}

/*
 INPUT:
 O: observation in a transposed matrix of matrix [T*Nmod], length = Nmod * T
 T: size of the chromosomes, length = 1
 N: number of states, length = 1
 Nmod: number of modifications, length = 1
 r: r parameters of the univariate distributions for every modifications, length = Nmod * 2
 p: p parameters of the univariate distributions for every modifications, length = Nmod * 2
 w: w parameters of the univariate distributions for every modifications, length = Nmod * 2, in case of Negative Binomial, w = 0
 interationMAX: maximum number of iterations for the HMM, length = 1
 eps: tolerant number
 states: the combinatorial states obtained from univariate analyses, length = N
 OUTPUT:
 post: posterior, transposed matrix of matrix [T*N]
 A: transition probability matrix [N*N]
 proba: initial probability, length = N
 loglik: loglikelihood
 */
extern "C" {
    void R_multivariate_hmm(double* O, int* T, int* N, int *Nmod, int* states, double* r, double* p, double* w, double* cor_matrix_inv, double* det, int* iterationMAX, double* eps, double* post, double* A, double* proba, double* loglik){
        printf("seqlen %i, nstates %i, maxit %i, eps %f\n", *T, *N, *iterationMAX, *eps);

        /* initialize */
        printf("init..\n");
        LogHMM* model = new LogHMM(O, *T, *N, *Nmod);
       // model->Nmod = *Nmod;
        printf("model->T = %i\nmodel->N = %i\nmodel->Nmod = %i\n", model->T, model->N, model->Nmod);
        
        //Prepare the enrich (univariate) vector: enrich[N][Nmod], e.g., enrich[iN][imod] tells me at state iN, modification imod is non-enriched (0) or enirched (1)
        bool **enrich = allocBoolMatrix(model->N, model->Nmod);
        
        for(int iN=0; iN < model->N; iN++) //for each comb state considered
        {
            for(int imod=0; imod < model->Nmod; imod++) //for each modification of this comb state
            {
                enrich[iN][imod] = states[iN]&(int)pow(2,model->Nmod-imod-1);//if =0->this hidden state has modification imod non enriched; if !=0->enriched
                if (enrich[iN][imod] !=0 )
                    enrich[iN][imod] = 1;
            }
        }
        

        //Prepare the observation vector for each modification
        /*double **Ot = allocDoubleMatrix(model->Nmod, model->T);
        for (int imod = 0; imod < model->Nmod; imod ++)
        {
            for (int t=0; t < model->T; t++)
            {
                Ot[imod][t] = O[imod*T+t]
            }
        }*/
        

        /* initialize the distributions */
        int i, j, t;
        for (int iN=0; iN<model->N; iN++)
        {
            vector <Density*> tempMarginals;            
            for (int imod=0; imod < model->Nmod; imod++)
            {
                Density *d;
                if (enrich[iN][imod]) //construct the density function for modification imod being enriched
                {
                    if (w[2*imod+1] == 0) //NegativeBinomial
                        d = new NegativeBinomial(NULL, r[2*imod+1], p[2*imod+1]);
                    else //ZiNB
                        d = new ZiNB(NULL, r[2*imod+1], p[2*imod+1], w[2*imod+1]);
                }
                else //construct the density function for modification imod being non-enriched
                {
                    if (w[2*imod] == 0)
                        d = new NegativeBinomial(NULL, r[2*imod], p[2*imod]);
                    else
                        d = new ZiNB(NULL, r[2*imod], p[2*imod], w[2*imod]);
                }
                tempMarginals.push_back(d);
            }
            //MVCopulaApproximation *tempMVdens = new MVCopulaApproximation(O, tempMarginals, &(cor_matrix_inv[iN*Nmod*Nmod]), det[iN]);
            MVCopulaApproximation *tempMVdens = new MVCopulaApproximation(O, tempMarginals, &(cor_matrix_inv[iN*model->Nmod*model->Nmod]), det[iN]);
            model->densityFunctions.push_back(tempMVdens);
        }
        printf("OK\n");
        
        /* estimate the parameters */
        printf("estimate parameters...\n");
        model->estimateTransitions(*iterationMAX, *eps);
        printf("OK\n");
        
        /* compute the posteriors and save results directly to the R pointer */
        /* we do not need to recompute the forward and backward variables since
         they have been computed in the estimation already */
        printf("compute posterior...");
        double** post_matrix = allocDoubleMatrix(model->N, model->T);
        model->posterior(post_matrix, 0);
        /* recode into column representation */
        for (t=0; t<model->T; t++)
        {
            for (i=0; i<model->N; i++)
            {
                post[t + i * model->T] = post_matrix[i][t];
            }
        }
        freeDoubleMatrix(post_matrix, model->N);
        
        /* also return the estimated transition matrix and the initial probs */
        for (i=0; i<model->N; i++)
        {
            proba[i] = exp(model->logproba[i]);
            for (j=0; j<model->N; j++)
            {
                A[i + j * model->N] = model->A[i][j];
            }
        }
        *loglik = model->logLikelihood();
        printf("OK\n");
    }
}
