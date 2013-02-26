#ifndef __RK4_H__GUIGUI_28082012_170000
#define __RK4_H__GUIUGI_28082012_170000

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>

#include"utils.h"

//#define __DEBUG_RK4_VEC

/*
 *              !!!!! Fonction d'initialisation !!!!!
 */
int init_rk4v(const int N);
void free_rk4v(void);

/*
 *              !!!!! VERSION SCALAIRE !!!!!
 * Fonction de résolution d'EDO au 1er ordre utilisant Runge Kutta 4
 * x : paramétre de "temps" (y(x))
 * *y : la fonction recherché
 * dx : le pas de variation
 * double (*func)(double,double,void*) : fonction second membre
 * *param : paramétre à passer à la fonction second membre
 */
double rk4(double x, double y, double dx, double(*func)(double, double, void*), void *param);

/*
 *              !!!!! VERSION VECTORIELLE !!!!!
 * Fonction de résolution d'EDO au 1er ordre utilisant Runge Kutta 4
 * x : paramétre de "temps" (y(x))
 * *y : la fonction recherché
 * dx : le pas de variation
 * double (*func)(double,double*,double*,const int,void*) : fonction second membre
 * *param : paramétre à passer à la fonction second membre
 */
void rk4v(double x, double *y, double dx, const int N, void(*func)(const double, const double*, double*, const int, void*), void *param);

/*
 *              !!!!! VERSION VECTORIELLE A PAS VARIABLE !!!!!
 * Fonction de résolution d'EDO au 1er ordre utilisant Runge Kutta 4
 * x : paramétre de "temps" (y(x))
 * *y : la fonction recherché
 * dx : le pas de variation
 * double (*func)(double,double,void*) : fonction second membre
 * eps : précision recherché
 * *param : paramétre à passer à la fonction second membre
 */
void pas_var(double *t, double *y, double dx, const int N, void(*func)(const double, const double*, double*, const int, void*), const double eps, void *param);

#endif
