#ifndef __UTILS_H_GUI__28082012_170000
#define __UTILS_H_GUI__28082012_170000

#include <stdio.h>
#include <stdlib.h>

/****************************************************************************************\
 * 	Code récupéré de la librairie mnitab de Mr lefrére, enseignant du MNI du M1	*
 * 	Physique Générale et du M1 SDUEE de l'université Pierre et Marie Curie Paris 6	*
\****************************************************************************************/
unsigned int* unsignedint1d(unsigned int n);
void          unsignedint1d_libere(unsigned int *ptf);

int*          int1d(int n);
void          int1d_libere(int *ptf);

float*        float1d(int n);
void          float1d_libere(float *ptf);

double*       double1d(int n);
void          double1d_libere(double *ptf);
int           maxlocdouble1d(double * x, int n);

float**       float2d(int nlignes, int mcolonnes);
void          float2d_libere(float ** mat);

double**      double2d(int nlignes, int mcolonnes);
void          double2d_libere(double ** mat);

#endif
