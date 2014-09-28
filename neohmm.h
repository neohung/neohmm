#ifndef NEOHMM_H
#define NEOHMM_H
typedef struct{
 int N;		//number of states
 int M;		//number of observation symbols
 double **A;  // transition matrix A[1..N][1..N]
 double **B;  // confususe matrix B[1..N][1..M]
 double *pi; // initial state distribution
}HMM;
typedef struct{
 int T;
 int *O;
}SEQUENCE;

void ReadHMM(FILE *fp, HMM *phmm);
void SaveHMM(FILE *fp,HMM *phmm);
void PrintHMM(HMM *phmm);
void ReadSequence(FILE *fp, SEQUENCE *pseq);
void PrintSequence(SEQUENCE *pseq);
double** Forward(HMM *phmm, SEQUENCE *pseq, double *pprob);
double** Backward(HMM *phmm, SEQUENCE *pseq,double *pprob);
void Viterbi(HMM *phmm, SEQUENCE *pseq, double *pprob, SEQUENCE *p_predict_seq);
HMM BaumWelch(HMM *phmm, SEQUENCE *pseq);
#endif
