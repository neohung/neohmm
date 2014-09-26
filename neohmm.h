typedef struct{
 int N;		//number of states
 int M;		//number of observation symbols
 double **A;  // transition matrix A[1..N][1..N]
 double **B;  // confususe matrix B[1..N][1..M]
 double **pi; // initial state distribution
}HMM;
typedef struct{
 int T;
  double *O;
}SEQUENCE;
