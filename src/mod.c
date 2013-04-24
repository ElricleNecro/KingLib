#include <king/mod.h>

/*
Phys initPhys(void)
{
	Phys param;
	param.G  = 1.;
	param.c  = 1.;
	param.M  = 1.;
	param.bc = sqrt(27.);
	return param;
}
*/

double dt = 1e-4;

double get_pas(void)
{
	return dt;
}

void set_pas(const double ndt)
{
	dt = ndt;
}

int Resol_King(Phys_King *obj, double **res, double *tmax, int *taille)
{
	double  t     = 0.0,				/**<-Conditions initiales*/
		t_max = *tmax,				/**<-Pas et temps max*/
		prec  = 1e-8;				/**<-Précision de la résolution*/
	double *u2    = NULL;
	int     N     = 2;				/**<-Taille du système à résoudre*/
	int     i_max = (int)((t_max-t)/dt);


	/*********************************************************\
	 *    Creation du tableaux des positions et vitesses	 *
	\*********************************************************/
	if( (u2       = (double*)calloc(N, sizeof(double))) == NULL )
	{
		fprintf(stderr, "\033[31mErreur d'allocation du tableau u2\n\033[00m");
		return -1;
	}

	*taille       = i_max; //(int)*tmax + 1;

	/**************************************************************************\
	 *				Initialisation				  *
	\**************************************************************************/
	int  err      = 1;
	if( (err      = init_rk4v(N)) != 1 )
	{
		fprintf(stderr, "\033[31mErreur à l'initialisation des routines RK4 :: erreur::%d\033[00m\n", err);
		free(u2);
		return err;
	}

	u2[0]         = obj->W0;
	u2[1]         = 0.0;

/*
	printf("\033[32mParamètre :\033[00m\n\033[33m\t- u(0) : %g v(0) = %g\n\t- xi : %g x_max : %g\n\t- dx : %g, precision : %g\n\t-> Nombre d'itération : %d\n\033[00m",
			u2[0], u2[1], t, t_max, dt, prec, i_max);
*/

#ifdef KING_LOG
	FILE *fich = NULL;
	if( (fich  = fopen("King.log", "w")) == NULL)
	{
		fprintf(stderr, "\033[31mImpossible d'écrire dans le fichier '%s'\n\033[00m", "King.log");
		free(u2);
		double2d_libere(res);
		free_rk4v();
		return -30;
	}
#endif

	/***************************************************************************\
	 * 			Physique de la sphére				   *
	\***************************************************************************/
	obj->rho_ori = rhoking(u2[0]);

	int     j = 0;
	double rp = u2[0];

	res[j][0] = t;
       	res[j][1] = u2[0];
	res[j][2] = u2[1];
	res[j][3] = rhoking(u2[0]);

#ifdef KING_LOG
	fprintf(fich, "%g\t%g\t%g\t%g\n", t, u2[0], u2[1], rhoking(u2[0]));
#endif

#ifdef _DEBUG_CALC_POT__
	fprintf(stderr, "%g\t%g\t%g\t%g\n", t, u2[0], u2[1], rhoking(u2[0]));
#endif

	for(int i = 1; i<i_max; i++)
	{
		j++;
		pas_var(&t, u2, dt, N, f_sphereking, prec, &(*obj));

#ifdef USE_FENV
		if( fetestexcept(FE_DIVBYZERO|FE_INEXACT|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW) )//isnan(rhoking(u2[0])/obj->rho_ori) != 0 )
#else
		if( isnan(rhoking(u2[0])/obj->rho_ori) != 0 )
#endif
		{
			t  -= dt;
			//dt /= 2.0;
			break;
		}
		//else if( fabs(rp - u2[0])/rp <= 0.01 )
		//	break;
		else
	                rp  = u2[0];

#ifdef KING_LOG
		fprintf(fich, "%g\t%g\t%g\t%g\n", t, u2[0], u2[1], rhoking(u2[0]));
#endif

#ifdef _DEBUG_CALC_POT__
		fprintf(stderr, "%g\t%g\t%g\t%g\n", t, u2[0], u2[1], rhoking(u2[0]));
#endif

		res[j][0] = t;
	       	res[j][1] = u2[0];
		res[j][2] = u2[1];
		res[j][3] = rhoking(u2[0]);
	}

#ifdef _DEBUG_CALC_POT__
	fprintf(stderr, "%g\n", res[0][0]);
#endif

	*taille = j;

#ifdef _DEBUG_CALC_POT__
	fprintf(stderr, "%i\n\n", *taille);
#endif

	*tmax   = t;

	(void)rp;

#ifdef KING_LOG
	fclose(fich);
#endif
	free_rk4v();
	free(u2);

	return 1;
}

// Densité du modéle de King en fonction du potentiel :
double rhoking(const double phi)
{
	return (exp(phi)*erf(sqrt(phi)) - sqrt(4.0*phi/(PI))*(1.0+2.0*phi/(3.0)));
}

// Le systéme d'équation de la sphere :
// Le prototype est imposé par les routines rk4?.
void f_sphereiso(const double x, const double *v, double *y, const int N, void *param)
{
	if( x/2.0 >= /*DBL_EPSILON*/1.0e-3)
	{
		y[0] = v[0]*(3.0 - v[0] - v[1])/x;	// u(x)
		y[1] = v[1]*(v[0] - 1.0)/x;		// v(x)
	}
	else
	{
		y[0] = -2.0/5.0*x + 19.0*4.0/1050.0*x*x*x - 823.0*6.0/189000.0*x*x*x*x*x;
		y[1] = 2.0/3.0*x - 4.0/30.0*x*x*x + 6.0/315.0*x*x*x*x*x;
	}

	(void)N;
	(void)param;
}

// Le systéme d'équation de King :
// Le prototype est imposé par les routines rk4?.
void f_sphereking(const double x, const double *v, double *y, const int N, void *param)
{
	if(x > 1e-6)
	{
		y[0]           = v[1];
		y[1]           = -rhoking(v[0]) - 2.0*v[1]/x;
	}
	else
	{
		Phys_King *par = (Phys_King*)param;
#ifdef _DEBUG__RK4_SPH_PERSO__
		fprintf(stderr, "%g ", par->W0);
#endif
		double g0      = par->W0,
		       g2      = -rhoking(g0)/3.0,
		       g4      = 2.0 * rhoking(g0) * (sqrt(g0/PI) - exp(g0) * erf(sqrt(g0))/2.0)/5.0;
		// Les termes g1, g3, g5 sont nuls.

		//gamma           = g0 + g2*x*x/2.0 + g4*x*x*x*x/24.0;
		y[0]           = g2*x + g4*x*x*x/6.0;
		y[1]           = g2   + g4*x*x/2.0;
	}

	(void)N;
}

void f_sphereking_eps(const double x, const double *v, double *y, const int N, void *param)
{
	Phys_King *par = (Phys_King*)param;
	y[0]           = v[1];
	y[1]           = -rhoking(v[0]/*, par*/) - 2.0*v[1]/fabs(x-par->eps);

	(void)N;
	(void)param;
}
