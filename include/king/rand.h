#ifndef RAND_PERSO__H_28082012_170000
#define RAND_PERSO__H_28082012_170000

/*********************************************************\
 * Générateurs aléatoire tiré des numérical recipes in C *
\*********************************************************/

#include <stdlib.h> //Change to math.h in K&R C.
#if defined __STDC__ && defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#include <tgmath.h>
#else
#include <math.h>
#endif

/*********************************************************\
 * 	Quelques define nécessaires aux algorithme	 *
\*********************************************************/
/*
#ifdef RAN0
#define    IA 16807
#define    IM 2147483647
#define    AM (1.0/IM)
#define    IQ 127773
#define    IR 2836
#define    MASK 123459876
/ * “Minimal” random number generator of Park and Miller. Returns a uniform random deviate
 * between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
 * to initialize the sequence; idum must not be altered between calls for successive deviates in
 * a sequence.
 * /
float ran0(long *idum);

#endif

#ifdef RAN1
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/ * “Minimal” random number generator of Park and Miller with Bays-Durham shuﬄe and added
 * safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
 * successive deviates in a sequence. RNMX should approximate the largest ﬂoating value that is
 * less than 1.
 * /
float ran1(long *idum);
#endif
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


/*
 * Long period (> 2 × 1018 ) random number generator of L’Ecuyer with Bays-Durham shuffle
 * and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
 * the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
 * idum between successive deviates in a sequence. RNMX should approximate the largest ﬂoating
 * value that is less than 1.
*/
float ran2(long *idum);

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

/*
 * According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted
 * for the above values :
 * 		ran3 has one nice feature: if your machine is poor on integer arithmetic (i.e.,
 * 		is limited to 16-bit integers), you can declare mj, mk, and ma[] as float, deﬁne
 * 		mbig and mseed as 4000000 and 1618033, respectively, and the routine will be
 * 		rendered entirely ﬂoating-point.
 *
 *
 * Returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value to
 * initialize or reinitialize the sequence.
 */
float ran3(long *idum);

#endif
