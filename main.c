#include<stdio.h>
#include<stdlib.h>
#include"neohmm.h"
void readhmm(FILE *fp, HMM *phmm)
{
	int i,j;
	double **a;
	fscanf(fp, "M= %d\n", &(phmm->M));
	fscanf(fp, "N= %d\n", &(phmm->N));
	fscanf(fp,"A:\n");
	a = (double**) calloc((unsigned) phmm->N, sizeof(double*));
	for(i=0; i< phmm->N; i++){
		a[i] = (double*) calloc((unsigned) phmm->N, sizeof(double));
		//a[i] = 0;
	}
	phmm->A = a;
	
	for (i=0;i<phmm->N;i++){
		for(j=0;j<phmm->N;j++){
			fscanf(fp, "%lf", &(phmm->A[i][j]));
		}
	}
	
	for (i=0;i<phmm->N;i++){
		for(j=0;j<phmm->N;j++){
			printf("%lf ",phmm->A[i][j]);
		}
			printf("\n");
	}
	

}
int main()
{
	HMM hmm;
	char* hmmfilename="a.hmm";
	FILE *fp = fopen(hmmfilename,"r");
	if (fp == NULL) {printf("fail to open %s\n", hmmfilename);exit(1);}
	readhmm(fp,&hmm);
	fclose(fp);
	printf("test\n");
	exit(0);	
}
