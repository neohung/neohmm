#include<stdio.h>
#include<stdlib.h>
#include"neohmm.h"
void readhmm(FILE *fp, HMM *phmm)
{
	int i,j;
	double **a,**b,**p;
	fscanf(fp, "M= %d\n", &(phmm->M));
	fscanf(fp, "N= %d\n", &(phmm->N));
	fscanf(fp,"A:\n");
	a = (double**) calloc((unsigned) phmm->N, sizeof(double*));
	for(i=0; i< phmm->N; i++){
		a[i] = (double*) calloc((unsigned) phmm->N, sizeof(double));
	}
	phmm->A = a;
	
	for (i=0;i<phmm->N;i++){
		for(j=0;j<phmm->N;j++){
			fscanf(fp, "%lf", &(phmm->A[i][j]));
		}
		fscanf(fp,"\n");
	}
	
	fscanf(fp,"B:\n");
	b = (double**) calloc((unsigned) phmm->N, sizeof(double*));
        for(i=0; i< phmm->N; i++){
                b[i] = (double*) calloc((unsigned) phmm->M, sizeof(double));
        }
        phmm->B = b;
	for (i=0;i<phmm->N;i++){
                for(j=0;j<phmm->M;j++){
                        fscanf(fp, "%lf", &(phmm->B[i][j]));
                }
		fscanf(fp,"\n");
        }
    fscanf(fp, "pi:\n");
	phmm->pi = (double *) calloc((unsigned) phmm->N,  sizeof(double*));
	for (i = 0; i < phmm->N; i++) {
		fscanf(fp, "%lf", &(phmm->pi[i])); 
	}
	//printHMM(phmm);	
}
void printHMM(HMM *phmm)
{
	int i,j;
       printf("M=%d\n",phmm->M);
       printf("N=%d\n",phmm->N);
       printf("A:\n");
        for (i=0;i<phmm->N;i++){
                for(j=0;j<phmm->N;j++){
                        printf("%lf ",phmm->A[i][j]);
                }
                        printf("\n");
        }
  	printf("B:\n");
        for (i=0;i<phmm->N;i++){
                for(j=0;j<phmm->M;j++){
                        printf("%lf ",phmm->B[i][j]);
                }
                        printf("\n");
        }
     printf("pi:\n");
      for (i=0;i<phmm->N;i++){
      	 printf("%lf ",phmm->pi[i]);
      }
}
void readsequence(FILE *fp, SEQUENCE *pseq)
{
	int i;
	int *O;
	fscanf(fp,"T= %d\n", &pseq->T);
	O = (int *)calloc((unsigned)pseq->T, sizeof(int));
	pseq->O = O;
        for (i=0;i<pseq->T;i++){
             fscanf(fp, "%d", &pseq->O[i]);
	}
	//printsequence(pseq);
}
void printsequence(SEQUENCE *pseq)
{
	int i;
	printf("T=%d\n",pseq->T);
	for (i=0;i<pseq->T;i++){
		printf("%d ",pseq->O[i] );
	}
	printf("\n");
}
void Forward(HMM *phmm, SEQUENCE *pseq, double *pprob)
{
	double **alpha, sum; 
	int i,j,t;
	alpha = (double**) calloc((unsigned) pseq->T, sizeof(double*));
	for(t=0; t< pseq->T; t++){
		alpha[t] = (double*) calloc((unsigned) phmm->N, sizeof(double));
	}
	//initial
    for (i = 0; i < phmm->N; i++){
            alpha[0][i] = phmm->pi[i]* phmm->B[i][pseq->O[0]];
    }
    //cal alpha[t][i], mean the prob of "see O[t] and the state is i"
    for (t = 0; t < pseq->T - 1; t++) {
                for (j = 0; j < phmm->N; j++) {
                        sum = 0.0;
                        for (i = 0; i < phmm->N; i++){
                               sum += alpha[t][i] * (phmm->A[i][j]);
 						}
                        alpha[t+1][j] = sum*(phmm->B[j][pseq->O[t+1]]);
                        //printf("alpha[%d][%d]=%lf\n",t+1,j,alpha[t+1][j]);
                }
    }
    *pprob = 0.0;
    for (i = 0; i < phmm->N; i++){
    	*pprob += alpha[pseq->T - 1][i];
    } 
    /*
	printf("alpha:\n");
	for (t = 0; t < pseq->T; t++){
		for (i = 0; i < phmm->N; i++){
			printf("%lf ", alpha[t][i]);
		}
		printf("\n");
	}
	printf("prob=%lf\n",*pprob);
	*/
}
void Backward(HMM *phmm, SEQUENCE *pseq, double *pprob)
{
	double **beta, sum; 
	int i,j,t;
	beta = (double**) calloc((unsigned) pseq->T, sizeof(double*));
	for(t=0; t< pseq->T; t++){
		beta[t] = (double*) calloc((unsigned) phmm->N, sizeof(double));
	}
	//initial
	for (i = 0; i < phmm->N; i++){
                beta[pseq->T-1][i] = 1.0;
    }
    for (t = pseq->T-2; t >= 0; t--) {
                for (i = 0; i < phmm->N; i++) {
                        sum = 0.0;
                        for (j = 0; j < phmm->N; j++){
                                sum += phmm->A[i][j] * (phmm->B[j][pseq->O[t+1]])*beta[t+1][j];
                        }
                        beta[t][i] = sum;
                }
    }
    *pprob = 0.0;
    for (i = 0; i < phmm->N; i++){
         *pprob +=(phmm->pi[i] * phmm->B[i][pseq->O[0]]* beta[0][i]);       
    }
    /*
    printf("beta:\n");
	for (t = pseq->T-1;t>=0; t--){
		for (i = 0; i < phmm->N; i++){
			printf("%lf ", beta[t][i]);
		}
		printf("\n");
	}
	printf("prob=%lf\n",*pprob);
	*/
}
int main()
{
	HMM hmm;
	SEQUENCE seq;
	char* hmmfilename="a.hmm";
	char* seqfilename="a.seq";
	FILE *fp = fopen(hmmfilename,"r");
	if (fp == NULL) {printf("fail to open %s\n", hmmfilename);exit(1);}
	readhmm(fp,&hmm);
	fclose(fp);
	fp = fopen(seqfilename,"r");
	if (fp == NULL) {printf("fail to open %s\n", seqfilename);exit(1);}
	readsequence(fp,&seq);
	fclose(fp);
	double prob_f,prob_b;
	Forward(&hmm, &seq, &prob_f);
	printf("test [%lf]\n",prob_f);
	Backward(&hmm, &seq, &prob_b);
	printf("test [%lf]\n",prob_b);

	exit(0);	
}
