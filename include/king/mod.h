#ifndef __MODULE_H_GUIGUI_28082012_170000
#define __MODULE_H_GUIGUI_28082012_170000

#include <stdio.h>
#include <stdlib.h>
#if defined __STDC__ && defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#include <tgmath.h>
#else
#include <math.h>
#endif

#ifdef USE_FENV
#include <fenv.h>
#endif

#include "rk4.h"
#include "cte_phys.h"

//#define _DEBUG_CALC_POT__
//#define _DEBUG__RK4_SPH_PERSO__

/**
 * \struct phys_king
 */
/**
 * \typedef Phys_King: Structure contenant toute les quantités physiques associées à un modèle de King.
 */
typedef struct phys_king {
	double m;		/**< Masse d'une particule de l'amas**/
	double El;		/**< Énergie de libération**/
	double rc;		/**< Rayon de coeur**/
	double sigma2;		/**< Ecart-type de l'énergie.**/
	double rho0;		/**< Densité centrale de la sphère.**/
	double W0;		/**< Condition initiale sur (El - m\phi)/\sigma^2**/
	double rho_ori;		/**< Valeur de la densité adimensionnee en x=0**/
	double eps;		/**< Précision pour l'équation autour de zero.**/
	double G;		/**< Constante de la gravitation.**/
        double rmax;            /**< Taille spatiale de l'objet.**/
        double vmax;            /**< Taille en vitesse de l'objet.**/
        double Mtot;            /**< Masse totale de l'objet.**/
}Phys_King;

double get_pas(void);
void set_pas(double ndt);
int Resol_King(Phys_King *obj, double **res, double *tmax, int *taille);

// Densité adimensionnee du modéle de King en fonction du potentiel :
double rhoking(const double phi);

/*
 * Le systéme d'équation de la sphere isotherme en boîte :
 * u*(3-u-v)/x
 * v*(u-1)/x
 * u <=> v[0]
 * v <=> v[1]
 *
 * \param x temps
 * \param *v valeurs au temps t+dt
 * \param *y valeurs au temps t
 * \param N taille du tableau de valeurs
 * \param *param paramètre de l'utilisateur
 */
void f_sphereiso(const double x, const double *v, double *y, const int N, void *param);

/**
 * Le système d'équation pour le modéle de King :
 * \f{align}
 * \frac{d \phi}{dr}
 * -4*pi*G*rho(r)
 *  \f{align}
 *
 * \param x temps
 * \param *v valeurs au temps t+dt
 * \param *y valeurs au temps t
 * \param N taille du tableau de valeurs
 * \param *param paramètre de l'utilisateur
 */
void f_sphereking(const double x, const double *v, double *y, const int N, void *param);

/**
 * Le système d'équation pour le modéle de King :
 * \f{align}
 * \frac{d \phi}{dr}
 * -4*pi*G*rho(r)
 * \f{align}
 *
 *  Cette version rajoute une option de précision pour décider de l'intervalle sur lequel utiliser le développement limité de la solution plutôt que l'équation
 *
 * \param x temps
 * \param *v valeurs au temps t+dt
 * \param *y valeurs au temps t
 * \param N taille du tableau de valeurs
 * \param *param paramètre de l'utilisateur
 */
void f_sphereking_eps(const double x, const double *v, double *y, const int N, void *param);

#endif
