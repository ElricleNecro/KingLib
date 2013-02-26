#include <king/utils.h>

unsigned int* unsignedint1d(unsigned int n)
{
	unsigned int * ptf=NULL;
	if (n >0) {
		ptf = (unsigned int *) calloc(n, sizeof(unsigned int)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation unsigned int1d\n");
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void unsignedint1d_libere(unsigned int *ptf)
{
	/* liberation d'un tableau 1D de unsigned ints */
	free(ptf);
	return ;
}

int* int1d(int n)
{
	int * ptf=NULL;
	if (n >0) {
		ptf = (int *) calloc(n, sizeof(int)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation int1d\n");
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void int1d_libere(int *ptf)
{
	/* liberation d'un tableau 1D de ints */
	free(ptf);
	return ;
}

float* float1d(int n)
{
	float * ptf=NULL;
	if (n >0) {
		ptf = (float *) calloc(n, sizeof(float)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation float1d\n");
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void float1d_libere(float *ptf)
{
	/* liberation d'un tableau 1D de floats */
	free(ptf);
	return ;
}

double* double1d(int n)
{
	double * ptf=NULL;
	if (n >0) {
		ptf = (double *) calloc(n, sizeof(double)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation double1d\n");
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void double1d_libere(double *ptf)
{
	/* liberation d'un tableau 1D de doubles */
	free(ptf);
	return ;
}

int maxlocdouble1d(double * x, int n)
{
	double tmp;
	int imax;
	int i;
	tmp = x[0];
	imax = 0 ;
	for (i=0; i<n; i++){
		if( x[i] > tmp) {
			tmp = x[i];
			imax = i ;
		}
	}
	return imax;
}

void float2d_libere(float ** mat) {
	free(mat[0]); /* liberation du pointeur de la matrice */
	free(mat) ; /* liberation du vecteur des pointeurs de debut de ligne */
}

float ** float2d(int nlignes, int mcolonnes){
	/* allocation d'un tableau 2D de nlignes x mcolonnes floats */
	float * ptf=NULL, ** pt=NULL;
	int k;
	if (nlignes >0 && mcolonnes >0) { /* Allocation globale de la matrice */
		ptf = (float *) calloc(nlignes*mcolonnes, sizeof(float)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation float2d \n");
		exit(EXIT_FAILURE);
	}
	/* allocation des n pointeurs de debut de ligne */
	pt = (float **) calloc(nlignes, sizeof(float*));
	if (pt==NULL) {
		fprintf(stderr, "erreur allocation float2d \n");
		exit(EXIT_FAILURE);
	}
	/* affectation des pointeurs de debut de ligne */
	for (k=0; k < nlignes; k++) {
		pt[k] = &(ptf[k*mcolonnes]);
	}
	return pt;
}

void double2d_libere(double ** mat) {
	free(mat[0]); /* liberation du pointeur de la matrice */
	free(mat) ; /* liberation du vecteur des pointeurs de debut de ligne */
}

double ** double2d(int nlignes, int mcolonnes){
	/* allocation d'un tableau 2D de nlignes x mcolonnes doubles */
	double * ptf=NULL, ** pt=NULL;
	int k;
	if (nlignes >0 && mcolonnes >0) { /* Allocation globale de la matrice */
		ptf = (double *) calloc(nlignes*mcolonnes, sizeof(double)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation double2d \n");
		exit(EXIT_FAILURE);
	}
	/* allocation des n pointeurs de debut de ligne */
	pt = (double **) calloc(nlignes, sizeof(double*));
	if (pt==NULL) {
		fprintf(stderr, "erreur allocation double2d \n");
		exit(EXIT_FAILURE);
	}
	/* affectation des pointeurs de debut de ligne */
	for (k=0; k < nlignes; k++) {
		pt[k] = &(ptf[k*mcolonnes]);
	}
	return pt;
}
