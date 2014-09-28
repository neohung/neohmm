#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "neohmm.h"
void do_training(char *trainedhmmfilename,char *readhmmfilename,char *readsequencefilename);
void do_viterbi(char *predictedfilename,char *readhmmfilename,char *readsequencefilename);

void Usage(char* processname)
{
	printf("\tTraining  : %s -T <trained.hmm> <initial.hmm> <data.seq> \n", processname);
	printf("\tPrediction: %s -V <predict.state> <trainid.hmm> <data.seq> \n", processname);
}
int main(int argc, char** argv)
{
	int c,errflg = 0, Tflg = 0, Vflg = 0;
	char *trainedhmmfilename,*predictedfilename, *readhmmfilename, *readsequencefilename;
	extern char *optarg;
    extern int optind, opterr, optopt;
    while((c=getopt(argc, argv, "hT:V:")) != -1)
	{
    	switch(c)
    	{
        	case 'h':
        		Usage(argv[0]);
				exit(1);
        		break;
        	case 'T':
        		if (Tflg){
                	errflg++;
                }else {
                    Tflg++;
                    //sscanf(optarg, "%s", &hmmfilename);
                    trainedhmmfilename=optarg;
                }
        		break;
        	case 'V':
        		if (Vflg){
                	errflg++;
                }else {
                    Vflg++;
                    //sscanf(optarg, "%s", &hmmfilename);
                    predictedfilename=optarg;
                }
        		break;
        	case ':':
        		printf("oops\n");
        		break;
        	case '?':
        		errflg++;
        		printf("Usage Error:\n");
        		break;
    	}
	}
	if (errflg){
		Usage(argv[0]);
		exit(1);
	}
	if (Tflg && !errflg){
		readhmmfilename = argv[optind];
		optind++;
		readsequencefilename = argv[optind];
		do_training(trainedhmmfilename,readhmmfilename,readsequencefilename);
	}
	if (Vflg && !errflg){
		readhmmfilename = argv[optind];
		optind++;
		readsequencefilename = argv[optind];
		do_viterbi(predictedfilename,readhmmfilename,readsequencefilename);
	}
	/*
	int t;
	HMM hmm;
	SEQUENCE seq, predict_seq,predict_seq_log;
	char* hmmfilename="a.hmm";
	char* seqfilename="a.seq";
	FILE *fp = fopen(hmmfilename,"r");
	if (fp == NULL) {printf("fail to open %s\n", hmmfilename);exit(1);}
	ReadHMM(fp,&hmm);
	fclose(fp);
	fp = fopen(seqfilename,"r");
	if (fp == NULL) {printf("fail to open %s\n", seqfilename);exit(1);}
	ReadSequence(fp,&seq);
	fclose(fp);
	double prob_f,prob_b,prob_v,prob_vl;
	double **alpha, **beta;
	alpha = Forward(&hmm, &seq, &prob_f);
	printf("prob_f [%lf]\n",prob_f);
	beta = Backward(&hmm, &seq,&prob_b);
	printf("prob_b [%lf]\n",prob_b);

	Viterbi(&hmm, &seq, &prob_v, &predict_seq);
	printf("before training: Viterbi=%lf\n",prob_v);
	for (t=0;t< predict_seq.T; t++){
		printf("%d->",predict_seq.O[t]);
	}
	printf("END\n");
//
	hmm = BaumWelch(&hmm, &seq);
//
	alpha = Forward(&hmm, &seq, &prob_f);
	printf("After training: prob_f [%lf]\n",prob_f);

    Viterbi(&hmm, &seq, &prob_v, &predict_seq);
	printf("After training: Viterbi=%lf\n",prob_v);
	for (t=0;t< predict_seq.T; t++){
		printf("%d->",predict_seq.O[t]);
	}
	printf("END\n");

	FILE *newfp = fopen("new.hmm","w");
	if (newfp == NULL) {printf("fail to open %s\n", "new.hmm");exit(1);}
	SaveHMM(newfp,&hmm);
	fclose(newfp);
	newfp = fopen("new.seq","w");
	if (newfp == NULL) {printf("fail to open %s\n", "new.seq");exit(1);}
	SaveSequence(newfp,&predict_seq);
	fclose(newfp);
*/	
	exit(0);	
}
void do_training(char *trainedhmmfilename,char *readhmmfilename,char *readsequencefilename)
{
	HMM hmm;
	SEQUENCE seq;
	FILE *fp = fopen(readhmmfilename,"r");
	if (fp == NULL) {printf("Fail to open <initial.hmm> :%s\n", readhmmfilename);exit(1);}
	ReadHMM(fp,&hmm);
	fclose(fp);
	fp = fopen(readsequencefilename,"r");
	if (fp == NULL) {printf("Fail to open <data.seq> :%s\n", readsequencefilename);exit(1);}
	ReadSequence(fp,&seq);
	fclose(fp);
	hmm = BaumWelch(&hmm, &seq);
	FILE *newfp = fopen(trainedhmmfilename,"w");
	if (newfp == NULL) {printf("Fail to open <trained.hmm> :%s\n", trainedhmmfilename);exit(1);}
	SaveHMM(newfp,&hmm);
	fclose(newfp);
	printf("save to %s, finish training\n",trainedhmmfilename);
}
void do_viterbi(char *predictedfilename,char *readhmmfilename,char *readsequencefilename)
{
	HMM hmm;
	SEQUENCE seq, predict_seq;
	double prob_v;
	FILE *fp = fopen(readhmmfilename,"r");
	if (fp == NULL) {printf("Fail to open <trainid.hmm> :%s\n", readhmmfilename);exit(1);}
	ReadHMM(fp,&hmm);
	fclose(fp);
	fp = fopen(readsequencefilename,"r");
	if (fp == NULL) {printf("Fail to open <data.seq> :%s\n", readsequencefilename);exit(1);}
	ReadSequence(fp,&seq);
	fclose(fp);
    Viterbi(&hmm, &seq, &prob_v, &predict_seq);
    printf("prob_v=%lf\n",prob_v);
	FILE *newfp = fopen(predictedfilename,"w");
	if (newfp == NULL) {printf("Fail to open <predict.state> :%s\n", predictedfilename);exit(1);}
	SaveSequence(newfp,&predict_seq);
	fclose(newfp);
	printf("save to %s, finish prediction\n",predictedfilename);
}
