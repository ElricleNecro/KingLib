#ifndef __CTE_PHYS_UTILES_28082012_170000
#define __CTE_PHYS_UTILES_28082012_170000

#if defined __STDC__ && defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#include <tgmath.h>
#else
#include <math.h>
#endif

/********************************************************************************\
 *	 	Quelques constantes de la physique (en SI)			*
\********************************************************************************/

//Constante de la gravitation :
extern double G_SI 6.67e-11;		// Nm^2kg^{-2}
//Vitesse de la lumiere :
extern double C_LUM 3.0e8;		// m/s
//Constante de Boltzmann :
extern double KB 1.380664e-23;		// J/K

//Valeur de pi :
#ifndef M_PI
#define PI acos(-1.0)
#endif

#endif
