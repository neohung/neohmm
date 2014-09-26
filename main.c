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
}
void readsequence(FILE *fp, SEQUENCE *pseq)
{
	int i;
	double *O;
	fscanf(fp,"T= %d\n", &pseq->T);
	O = (double *)calloc((unsigned)pseq->T, sizeof(double));
	pseq->O = O;
        for (i=0;i<pseq->T;i++){
             fscanf(fp, "%lf", &pseq->O[i]);
	}
	//printsequence(pseq);
}
void printsequence(SEQUENCE *pseq)
{
	int i;
	printf("T=%d\n",pseq->T);
	for (i=0;i<pseq->T;i++){
		printf("%lf ",pseq->O[i] );
	}
	printf("\n");
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
	printf("test\n");
	exit(0);	
}
