#include <king/rk4.h>

double *k1       = NULL,
       *k2       = NULL,
       *k3       = NULL,
       *k4       = NULL,
       *k5       = NULL,
       *k6       = NULL,
       *tmp      = NULL;
double *pas1     = NULL,
       *pas2     = NULL;
FILE   *log_fich = NULL;

int init_rk4v(const int N)
{
	if( k1 != NULL )
		free_rk4v();

	if( (k1  = (double*)calloc(N, sizeof(double))) == NULL)
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}
	if( (k2  = (double*)calloc(N, sizeof(double))) == NULL)
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}
	if( (k3  = (double*)calloc(N, sizeof(double))) == NULL)
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}
	if( (k4  = (double*)calloc(N, sizeof(double))) == NULL)
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}
	if( (k5  = (double*)calloc(N, sizeof(double))) == NULL)
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}
	if( (k6  = (double*)calloc(N, sizeof(double))) == NULL)
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}

	if( (tmp  = (double*)calloc(N, sizeof(double))) == NULL)
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}

	if( (pas1 = (double*)calloc(N, sizeof(double))) == NULL )
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}
	if( (pas2 = (double*)calloc(N, sizeof(double))) == NULL )
	{
		fprintf(stderr,
			"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",
			__FILE__,__LINE__-1);
		return -1;
	}
/*	if( (log_fich = fopen("rk4_pas.log", "w")) == NULL )
	{
		fprintf(log_fich, "\033[31mImpossible d'ouvrir le fichier log : %s <- %s:%d\033[00m\n", "rk4_pas.log", __FILE__, __LINE__);
		return -1;
	}
*/	return 1;
}

void free_rk4v(void)
{
	free(k1),k1=NULL;
	free(k2),k2=NULL;
	free(k3),k3=NULL;
	free(k4),k4=NULL;
	free(k5),k5=NULL;
	free(k6),k6=NULL;
	free(tmp),tmp=NULL;
	free(pas1),pas1=NULL;
	free(pas2),pas2=NULL;

//	fclose(log_fich);
}

double rk4(double x, double y, double dx, double(*func)(double t, double x, void *param), void *param)
{
	double sk1, sk2, sk3, sk4;

	sk1=func(x, y, param);
	sk2=func(x + dx/2., y + sk1*dx/2., param);
	sk3=func(x + dx/2., y + sk2*dx/2., param);
	sk4=func(x + dx, y + sk3*dx, param);

	return y + (sk1 + 2.*sk2 + 2.*sk3 + sk4)*dx/6.;
}

void rk4v(double x, double *y, double dx, const int N, void(*func)(const double, const double*, double*, const int, void*), void *param)
{
	for(int i=0; i<N; i++)
		tmp[i] = y[i];

	func(x, y, k1, N, param);
	for(int i=0; i<N; i++)
	{
		tmp[i] = y[i] + k1[i]*dx/2.;
	}

	func(x + dx/2., tmp, k2, N, param);
	for(int i=0; i<N; i++)
	{
		tmp[i] = y[i] + k2[i]*dx/2.;
	}

	func(x + dx/2., tmp, k3, N, param);
	for(int i=0; i<N; i++)
	{
		tmp[i] = y[i] + k3[i]*dx;
	}

	func(x + dx, tmp, k4, N, param);

	for(int i=0; i<N; i++)
		y[i] = y[i] + (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i])*dx/6.;

}

void pas_var(double *t, double *y, double dx, const int N, void(*func)(const double, const double*, double*, const int, void*), const double eps, void *param)
{
	int    nbite = 0 ;
	double err ;
	double err0 ;
	double *Et  = NULL ;
	double tb   = *t ;
	double hmax = dx / 3. ;
	double h    = hmax ;
	double hn ;
	double to, tout ;
	double *anc = NULL;
	int    ii ;

	Et    = double1d(N);
	anc   = double1d(N);
	to    = tb ;
	tout  = *t + dx ;

	nbite = 0 ;
	while (tb < tout && nbite < 100) {
		nbite++ ;
		for(int i=0; i<N; i++) {
			pas1[i] = pas2[i] = y[i] ;
		}
		rk4v(tb, pas1, h, N, func, param) ;
		rk4v(tb, pas2, h*0.5, N, func, param) ;
		rk4v(tb + h*0.5, pas2, h*0.5, N, func, param) ;

		func(tb + h, pas2, tmp, N, param) ;

		for (int i=0; i<N; i++) {
			Et[i]  = fabs(pas2[i] - pas1[i]) ;
		}
		ii = maxlocdouble1d(Et, N) ;
		err0 = eps * (fabs(pas2[ii]) + h * fabs(tmp[ii])) ;
		err = Et[ii] ;
		if (err > DBL_EPSILON ) {
			hn  = 0.95 * h * pow(err0/err, 1./5.) ;
		} else {
			hn  = 10. * h ;
		}
		//	printf (" %d %g %g %g %g %g %g\n", nbite, tb, err0, err, h, hn, hmax) ;
		if (hn > hmax) {
			hn = hmax ;
		}

		if (err > err0) {
			h = hn ;
		} else {
			for(int i=0; i<N; i++) {
				anc[i] = y[i] ;
				y[i] = (16. * pas2[i] - pas1[i]) / 15. ;
			}
			to = tb ;
			tb += h ;
			h = hn ;
		}
	}
	if (tb > tout) {
		for (int i=0; i<N; i++) {
			pas1[i] = anc[i] ;
			pas2[i] = anc[i] ;
		}
		h = tout - to ;
		rk4v(to, pas1, h, N, func, param) ;
		rk4v(to, pas2, h*0.5, N, func, param) ;
		rk4v(to + h*0.5, pas2, h*0.5, N, func, param) ;
		for (int i=0; i<N; i++) {
			y[i] = (16. * pas2[i] - pas1[i]) / 15. ;
		}
		tb = to + h ;
	}
	*t = tb ;

	free(Et) ;
	free(anc) ;
}

/*
void rkfe(double x, double *y, double dx, const int N, void(*func)(double, const double*, double*, const int, void*), void *param)
{
	double *k1  = NULL,
	       *k2  = NULL,
	       *k3  = NULL,
	       *k4  = NULL,
	       *k5  = NULL,
	       *k6  = NULL,
	       *tmp = NULL;

	if( (k1  = (double*)calloc(N, sizeof(double))) == NULL)
		fprintf(stderr,"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",__FILE__,__LINE__),exit(EXIT_FAILURE);
	if( (k2  = (double*)calloc(N, sizeof(double))) == NULL)
		fprintf(stderr,"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",__FILE__,__LINE__),exit(EXIT_FAILURE);
	if( (k3  = (double*)calloc(N, sizeof(double))) == NULL)
		fprintf(stderr,"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",__FILE__,__LINE__),exit(EXIT_FAILURE);
	if( (k4  = (double*)calloc(N, sizeof(double))) == NULL)
		fprintf(stderr,"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",__FILE__,__LINE__),exit(EXIT_FAILURE);
	if( (k5  = (double*)calloc(N, sizeof(double))) == NULL)
		fprintf(stderr,"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",__FILE__,__LINE__),exit(EXIT_FAILURE);
	if( (k6  = (double*)calloc(N, sizeof(double))) == NULL)
		fprintf(stderr,"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",__FILE__,__LINE__),exit(EXIT_FAILURE);
	if( (tmp  = (double*)calloc(N, sizeof(double))) == NULL)
		fprintf(stderr,"\033[31mImpossible d'allouer la memoire : %s:%d\033[00m\n",__FILE__,__LINE__),exit(EXIT_FAILURE);

	do
	{
		func(x, y, k1, N, param);
		for(int i=0; i<N; i++)
		{
			k1[i]  = dx*k1[i];
			tmp[i] = y[i] + k1[i]*1./4.;
		}

		func(x + dx/4., tmp, k2, N, param);
		for(int i=0; i<N; i++)
		{
			k2[i]  = dx*k2[i];
			tmp[i] = y[i] + k2[i]*9./32. + k1[i]*9./32.;
		}

		func(x + dx*3./8., tmp, k3, N, param);
		for(int i=0; i<N; i++)
		{
			k3[i]  = dx*k3[i];
			tmp[i] = y[i] + k3[i]*7296./2197 + k2[i]*7200./2197. + k1[i]*1932./2197;
		}

		func(x + dx*12./13., tmp, k4, N, param);
		for(int i=0; i<N; i++)
		{
			k4[i]  = dx*k4[i];
			tmp[i] = y[i] + 439./216.*k1[i] - 8.*k2[i] + 3680./513. - 845./4104.*k4[i];
		}

		func(x + dx, tmp, k5, N, param);
		for(int i=0; i<N; i++)
		{
			k5[i]  = dx*k5[i];
			tmp[i] = y[i] - 8./27.*k1[i] + 2.*k2[i] - 3544./2565.*k3[i] + 1859./4104.*k4[i] - 11./40.*k5[i];
		}

		func(x + dx*0.5, tmp, k6, N, param);
		for(int i=0; i<N; i++)
		{
			k6[i]  = dx*k6[i];
		}

		for(int i=0; i<N; i++)
			y[i] = y[i] + (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i])*dx/6.;

	}while(true);

	free(k1);
	free(k2);
	free(k3);
	free(k4);
	free(k5);
	free(k6);
	free(tmp);
}
*/

