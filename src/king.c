#include <king/king.h>

/*extern*/ int K_r    = 0;
/*extern*/ int K_pot  = 1;
/*extern*/ int K_dpot = 2;
/*extern*/ int K_rho  = 3;
/*extern*/ int K_mu   = 4;
/*extern*/ int K_temp = 5;

static double G = G_SI;

double King_get_PhysVal(const King *obj, int i, int j)
{
        return obj->don[i][j];
}

void King_SetG(const double G_new)
{
	G = G_new;
}

King King_New(const double W0, const double rc, const double sig_v)
{
	King Amas;

	Amas.amas.W0     = W0;
	Amas.amas.rc     = rc;
	Amas.amas.sigma2 = sig_v*sig_v;
	Amas.amas.G      = G;
	Amas.amas.rho0   = Amas.amas.sigma2 / ( 8.0 * M_PI * Amas.amas.G * Amas.amas.rc * Amas.amas.rc );
	Amas.don         = NULL;

	return Amas;
}

int read_utile(King *Amas, const char const *str)
{
	FILE *fich = NULL;

	// Coefficient pour la relation ajusté T_c <=> W_O (cf Rapprt Chapitre 4) :
	double a         = -10.0698;
	double b         = 0.220152;
	double c         = -1.63409;
	double d         = -2.3341;
	double e         = 16.913;
        double G         = 0.0;

	//Fraction de El (cf cahiers + rapport) :
	//        double frac      = 2.0;
	double sig_v     = 2.9,                                         /**<- En km/s */
	       rc        = 1.35,                                        /**<- En arc minutes */
//               m         = 2e30,                                      /**<- En kg */
	       R_sun     = 8.9,                                         /**<- En kpc */
	       Tc        = pow(10.0, 9.32);                             /**<- En années*/

	if( (fich = fopen(str, "r")) == NULL )
	{
		perror("Erreur : ");
		return -1;
	}

	if( fscanf(fich, "%lf", &rc) != 1)
	{
                fputs("\033[31mImpossible de lire la premiére ligne du fichier (rc)\033[00m", stderr);
		fclose(fich);
		return -2;
	}
	if( fscanf(fich, "%lf", &sig_v) != 1)
	{
                fputs("\033[31mImpossible de lire la seconde ligne du fichier (sig_v)\033[00m", stderr);
		fclose(fich);
		return -3;
	}
	if( fscanf(fich, "%lf %lf", &Tc, &R_sun) != 2)
	{
                fputs("\033[31mImpossible de lire la troisiéme ligne du fichier (Tc et R_sun)\033[00m", stderr);
		fclose(fich);
		return -4;
	}
	if( fscanf(fich, "%lf", &G) != 1)
	{
                fputs("\033[31mImpossible de lire la derniére ligne du fichier (G)\033[00m", stderr);
		fclose(fich);
		return -5;
	}

	fclose(fich);

/*
	fprintf(stderr,"Tc : %g\n",Tc);
	fprintf(stderr,"a : %g\n",a);
	fprintf(stderr,"b : %g\n",b);
	fprintf(stderr,"c : %g\n",c);
	fprintf(stderr,"d : %g\n",d);
	fprintf(stderr,"e : %g\n",e);
	fprintf(stderr,"1/b : %g\n",1/b);
	fprintf(stderr,"(d*Tc + e - c)/a : %g\n", (d*Tc + e - c)/a);
	fprintf(stderr,"log( (d*Tc + e - c)/a ) : %g\n",log( (d*Tc + e - c)/a ));
*/

	Amas->frac        = 2.0;
	Amas->amas.rc     = R_sun * tan(rc/60.0 /*degrés*/ * PI/180.0/*radians*/) /*kpc*/ * 1e3 /*pc*/ * 3.086e13 /*km*/ * 1e3 /*m*/;
	Amas->amas.W0     = fabs((1.0/b) * log( (d*Tc + e - c)/a ));
	//21.0; //9.0; //1.0/0.220152 log( -(-2.3341 * log10(pow(10.0, 9.32)) + 16.913 + 1.63409)/10.0698 );
	//	Amas->amas.m      = m /*kg*/;
	Amas->amas.sigma2 = (sig_v /*km/s*/ * 1e3 /*m/s*/)*(sig_v /*km/s*/ * 1e3 /*m/s*/); // Vitesse au carré
	Amas->amas.rho0   = Amas->amas.sigma2 / (8.0 * PI * (G) * Amas->amas.rc * Amas->amas.rc);
	//Amas.amas.El     = Amas.amas.W0 * Amas.amas.sigma2 / (1-frac);
	//	Amas->amas.El     = -frac*Amas->amas.W0 * Amas->amas.sigma2 / (frac-1);
	//

	Amas->amas.G      = G;

	return 1;
}

double King_Temp(const King *Amas, int i)
{
	double msig = Amas->amas.m * Amas->amas.sigma2, 							/** $m \sigma^2$ */
	       cte  = sqrt(2.0 * msig/acos(-1.0)) / pow(Amas->amas.m, 2.0), 					/** $\sqrt{\dfrac{2 m \sigma^2}{\pi}}\dfrac{1}{m^2}$ */
	       cte2 = 3.0/2.0 * sqrt(2.0 * PI * msig), 								/** $\dfrac{3 \sqrt{2m\pi\sigma^2}}{2}$ */
	       phi  = Amas->amas.El - Amas->amas.m * Amas->don[i][1],						/** $E_l - m\psi(r)$ */
	       gam  = phi/Amas->amas.sigma2, 									/** $\dfrac{\phi(r)}{\sigma^2}$ */
	       fcts = pow(2.0 * Amas->amas.m * phi, 5.0/2.0) * pow(msig, -2.0) / (5.0 * rhoking(gam));		/** $\dfrac{2m\phi(r)^{5/2} (m\sigma^2)^{-2}}{5\rho(\gamma(r))}$ */

#ifdef DBG_Temp_FUNC
	fprintf(stderr, "\033[31m%s::debug\n\tCTE  :: %g\n\tCTE2 :: %g\n\tFCTS :: %g\n\tdon[%i][1] :: %g\n\tGamma :: %g\n\trhoking :: %g\t (2.0*(El-m*psi))**5/2 :: %g\033[00m\n",
			__func__,
			cte,
			cte2,
			fcts,
			i,
			Amas->don[i][1],
			gam,
			rhoking(gam),
			pow(2.0*(Amas->amas.El - Amas->amas.m*Amas->amas.m*Amas->don[i][1]), 5.0/2.0)
	       );
	fprintf(stderr, "\033[31m\tEl :: %g\n\tm :: %g\n\tsigma2 :: %g\033[00m\n",
			Amas->amas.El,
			Amas->amas.m,
			Amas->amas.sigma2
	       );
#endif

#ifdef KB_UTIL
	return Amas->amas.m * cte * (cte2 - fcts) / (3.0 * KB);
#else
	return cte * (cte2 - fcts);
#endif
}

double King_TempMoy(const King *Amas/*, double Mtot*/)
{
	/***********************************************\
	 * 	Calcul de la température moyenne       *
	\***********************************************/
	double sum = 0.0,
	       s   = 0.0,
	       dr  = Amas->don[1][0] - Amas->don[0][0];
	double a   = Amas->don[0][0],
	       b   = Amas->don[Amas->lig-1][0];
	int    i,
	       it  = Amas->lig - 2;
#ifdef DBG_Temp_FUNC
	fprintf(stderr, "\033[31m%s:%s:%d :: a = %g, b = %g, dr = %g, it = %d\n\tAmas.lig = %d, %p\033[00m\n",
			__FILE__,
			__func__,
			__LINE__,
			a,
			b,
			dr,
			it,
			Amas->lig,
			&(Amas->lig)
	       );
#endif
/*#ifdef DBG_Temp_FUNC
	it = 5;
	fprintf(stderr, "\033[32ma :: %g\tb :: %g\tdr :: %g\033[00m\n", a, b, dr);
	fprintf(stderr, "\033[32m%s:%s:%d: Sum :: %g\033[00m\n", __FILE__, __func__, __LINE__, sum);
	for(sum=0.0,i=1;i<=it;i++)
	{
#	ifdef USE_NR
		sum += Amas->don[i][0]*Amas->don[i][0] * pow( (Amas->amas.El - Amas->amas.m * Amas->don[i][1]), 5.0/2.0);
#	else
		sum += dr * Amas->don[i][0]*Amas->don[i][0] * pow( (Amas->amas.El - Amas->amas.m * Amas->don[i][1]), 5.0/2.0) ;
		fprintf(stderr, "\033[32m%s:%s:%d: Sum :: %g (%g, %g, %g)\033[00m\n",
				__FILE__,
				__func__,
				__LINE__,
				sum,
				Amas->don[i][0]*Amas->don[i][0],
				pow( (Amas->amas.El - Amas->amas.m * Amas->don[i][1]), 5.0/2.0),
				Amas->don[i][0]*Amas->don[i][0] * pow( (Amas->amas.El - Amas->amas.m * Amas->don[i][1]), 5.0/2.0)
		       );
#	endif
	}
#else*/
	for(i=1;i<=it;i++)
	{
#	ifdef USE_NR
		sum += Amas->don[i][0]*Amas->don[i][0] * pow( (Amas->amas.El - Amas->amas.m * Amas->don[i][1]), 5.0/2.0);
#	else
		sum += dr * Amas->don[i][0]*Amas->don[i][0] * pow( (Amas->amas.El - Amas->amas.m * Amas->don[i][1]), 5.0/2.0) ;
#	endif
	}
//#endif

#ifdef USE_NR
	s = 0.5*(s+(b-a)*sum/it);
#else
	s = sum + dr*( a*a * pow( (Amas->amas.El - Amas->amas.m * Amas->don[0][1]), 5.0/2.0)
		       +  b*b * pow( (Amas->amas.El - Amas->amas.m * Amas->don[Amas->lig-1][1]), 5.0/2.0) )/2.0;
#endif

#ifdef DBG_Temp_FUNC
	fprintf(stderr, "sum = %g, s = %g\n\t%g\n\t%g",
			sum,
			s,
			a*a*pow(Amas->amas.El - Amas->amas.m * Amas->don[0][1] , 5.0/2.0),//Amas->don[0][0]*Amas->don[0][0] * pow( (Amas->amas.El - Amas->amas.m * Amas->don[0][1]), 5.0/2.0),
		       	b*b*pow(Amas->amas.El - Amas->amas.m * Amas->don[Amas->lig-1][1], 5.0/2.0)//Amas->don[Amas->lig-2][0]*Amas->don[Amas->lig-2][0] * pow( (Amas->amas.El - Amas->amas.m * Amas->don[Amas->lig-2][1]), 5.0/2.0)
		);
#endif

	return  3.0*Amas->amas.sigma2/Amas->amas.m - 32.0 * PI * PI * Amas->amas.rho0 / (5.0*Amas->amas.m * Amas->amas.Mtot * pow(PI*Amas->amas.sigma2, 3.0/2.0)) * s;
}

double King_mu(const King *Amas, const int a, const int b)
{
	double sum = 0.0,
	       s   = 0.0,
	       xa  = Amas->don[a][0],
	       xb  = Amas->don[b][0],
	       dr  = Amas->don[1][0] - Amas->don[0][0];

	for(int i = a + 1; i <= b - 1; i++)
	{
#	ifdef USE_NR
		sum += 4*PI * Amas->don[i][0]*Amas->don[i][0] * Amas->don[i][3];
#	else
		sum += 4*PI * dr * Amas->don[i][0]*Amas->don[i][0] * Amas->don[i][3];
#	endif
	}

#ifdef USE_NR
	s = 0.5*(s+(b-a)*sum/it);
#else
	s = sum + dr*( xa*xa * 4*PI * Amas->don[a][3] + xb*xb * 4*PI * Amas->don[b][3])/2.0;
#endif

	return s;
}

double King_Ec(const King *Amas, const int a, const int b)
{
	double sum = 0.0,
	       s   = 0.0,
	       xa  = Amas->don[a][0],
	       xb  = Amas->don[b][0],
	       dr  = Amas->don[1][0] - Amas->don[0][0];

	for(int i = a + 1; i <= b - 1; i++)
	{
#	ifdef USE_NR
		sum += 4*PI * Amas->don[i][0]*Amas->don[i][0] * King_Temp(Amas, i);
#	else
		sum += 4*PI * dr * Amas->don[i][0]*Amas->don[i][0] * King_Temp(Amas, i);
#	endif
	}

#ifdef USE_NR
	s = 0.5*(s+(b-a)*sum/it);
#else
	s = sum + dr*( xa*xa * 4*PI * King_Temp(Amas, a) + xb*xb * 4*PI * King_Temp(Amas, b))/2.0;
#endif

	return Amas->amas.m * s / 2.0;
}

double King_Ep(const King *Amas, const int a, const int b)
{
	double sum = 0.0,
	       s   = 0.0,
	       xa  = Amas->don[a][0],
	       xb  = Amas->don[b][0],
	       dr  = Amas->don[1][0] - Amas->don[0][0];

	for(int i = a + 1; i <= b - 1; i++)
	{
#	ifdef USE_NR
		sum += 4*PI * Amas->don[i][0] * Amas->don[i][3] * Amas->don[i][K_mu];//King_mu(Amas, a, i);
#	else
		sum += 4*PI * dr * Amas->don[i][0] * Amas->don[i][3] * Amas->don[i][K_mu];//King_mu(Amas, a, i);
#	endif
	}

#ifdef USE_NR
	s = 0.5*(s+(b-a)*sum/it);
#else
	s = sum + dr*( xa * 4*PI * Amas->don[a][3] * Amas->don[a][K_mu]/*King_mu(Amas, a, b)*/ + xb * 4*PI * Amas->don[b][3] * Amas->don[b][K_mu]/* King_mu(Amas, a, b)*/)/2.0;
//	s = sum + dr*( xa * 4*PI * Amas->don[a][3] * King_mu(Amas, a, b) + xb * 4*PI * Amas->don[b][3] * King_mu(Amas, a, b))/2.0;
#endif

	return -4.0 * PI * Amas->amas.G * s;
}

//#ifdef SWIG_RB
//void freeKing(King *obj){
//	double2d_libere(obj->don);
//}

//#else
void King_free(King *obj)
{
	double2d_libere(obj->don), obj->don=NULL;
}

//#endif

//void King_ugbs(King *obj, double *rmax, double *vmax, double *Mtot, const int NbPart, const double G)
void King_ugbs(King *obj)
{
	int taille, err;
	obj->amas.rmax    = 100.0;

	obj->col = 5;

	if( (obj->don = double2d((int)((obj->amas.rmax)/get_pas()), obj->col)) == NULL)
	{
		fprintf(stderr, "\033[31mErreur d'allocation du tableau res\n\033[00m");
		exit(EXIT_FAILURE);
	}

	if( (err = Resol_King(&(obj->amas), obj->don, &obj->amas.rmax, &taille)) != 1)
	{
		fprintf(stderr, "\033[31mErreur lors du calcul du potentiel :: %d\033[00m\n", err);
		King_free(obj);
		exit(err);
	}
	obj->lig = taille;
}

//void get_box_size(King *obj, double *rmax, double *vmax, double *Mtot, const int NbPart, const double G)
void King_gbs(King *obj, const int NbPart)
{
	int taille, err;
	obj->amas.rmax    = 100.0;

	obj->col = 5;

	if( (obj->don = double2d((int)((obj->amas.rmax)/get_pas()), obj->col)) == NULL)
	{
		fprintf(stderr, "\033[31mErreur d'allocation du tableau res\n\033[00m");
		exit(EXIT_FAILURE);
	}

	if( (err = Resol_King(&(obj->amas), obj->don, &(obj->amas.rmax), &taille)) != 1)
	{
		fprintf(stderr, "\033[31mErreur lors du calcul du potentiel :: %d\033[00m\n", err);
		King_free(obj);
		exit(err);
	}

#ifdef KING_DIM_LOG
	FILE *fich = NULL;
	if( (fich  = fopen("King_dim.log", "w")) == NULL)
	{
		fprintf(stderr, "\033[31mImpossible d'écrire dans le fichier '%s'\n\033[00m", "King.log");
		double2d_libere(obj->don);
		exit(EXIT_FAILURE);
	}
#endif
	obj->amas.rmax    = obj->amas.rmax * obj->amas.rc;
	obj->lig = taille;

	for(int j = 0; j < taille; j++)
	{
		obj->don[j][0] = obj->don[j][0] * obj->amas.rc;
		obj->don[j][3] = obj->don[j][3] * obj->amas.rho0;
	}

	double dr               = obj->don[1][0] - obj->don[0][0];
	obj->amas.Mtot          = 0.0;
	for(int j = 1; j < taille; j++)
		obj->amas.Mtot += 4.0*PI * obj->don[j][0] * obj->don[j][0] * obj->don[j][3] /*f(a + i*dr)*/ * dr;

	obj->amas.Mtot         += ( 4.0*PI * obj->don[0][0] * obj->don[0][0] * obj->don[0][3] + 4.0*PI * obj->don[taille - 1][0] * obj->don[taille - 1][0] * obj->don[taille - 1][3] ) * dr/2.0;
			//(obj->don[j][0] - obj->don[j-1][0]) * 4.0*PI*(obj->don[j][0] * obj->don[j][0] * obj->don[j][3] + obj->don[j-1][0] * obj->don[j-1][0] * obj->don[j-1][3])/2.0;		/** <- Intégration par la méthode des trapéze */
//	*Mtot                 *= (obj->don[j][0] - obj->don[j-1][0]);
//	*Mtot                 += (obj->don[j][0] - obj->don[j-1][0])/2.0 * 4.0*PI*(obj->don[taille-1][0] * obj->don[taille-1][3] + obj->don[0][0] * obj->don[0][3])

	obj->amas.m            = obj->amas.Mtot/NbPart;

        obj->amas.sigma2       = obj->amas.m * obj->amas.sigma2 / 2.0 /*Joules*/;
        //obj->amas.El           = Amas.amas.W0 * Amas.amas.sigma2 / (1-frac);
        obj->amas.El           = - obj->amas.G * obj->amas.m * obj->amas.Mtot / obj->amas.rmax;//-obj->frac*obj->amas.W0 * obj->amas.sigma2 / (obj->frac-1);


	// On redimensionne le probléme :
	for(int j = 0; j < taille; j++)
	{
	       	obj->don[j][1] = (obj->amas.El - obj->don[j][1] * obj->amas.sigma2)/obj->amas.m;
		obj->don[j][2] = -obj->amas.sigma2 * obj->don[j][2]/obj->amas.m;
		obj->don[j][4] = obj->don[j][0] * obj->don[j][0] * obj->don[j][2] / obj->amas.G; //King_mu(obj, 0, j);
#ifdef KING_DIM_LOG
		fprintf(fich, "%g\t%g\t%g\t%g\n", obj->don[j][0], obj->don[j][1], obj->don[j][2], obj->don[j][3]);
#endif
	}

#ifdef KING_DIM_LOG
	fclose(fich);
#endif

	obj->amas.vmax    = sqrt(2.0 * (obj->amas.El/obj->amas.m - obj->don[0][1]));

#ifdef DBG_Temp_FUNC
	fprintf(stderr, "%s::%s::%d :: Nb de ligne : %d, %p\n", __FILE__, __func__, __LINE__, obj->lig, &obj->lig);
#endif
}

void King_gbs_cb(King *obj, int *NbPart, gbs_cb f /*double (*f)(King*, const int, void*)*/, void *data)
{
	int taille, err;
	obj->amas.rmax    = 100.0;

	obj->col = 5;

	if( (obj->don = double2d((int)((obj->amas.rmax)/get_pas()), obj->col)) == NULL)
	{
		fprintf(stderr, "\033[31mErreur d'allocation du tableau res\n\033[00m");
		exit(EXIT_FAILURE);
	}

	if( (err = Resol_King(&(obj->amas), obj->don, &(obj->amas.rmax), &taille)) != 1)
	{
		fprintf(stderr, "\033[31mErreur lors du calcul du potentiel :: %d\033[00m\n", err);
		King_free(obj);
		exit(err);
	}

#ifdef KING_DIM_LOG
	FILE *fich = NULL;
	if( (fich  = fopen("King_dim.log", "w")) == NULL)
	{
		fprintf(stderr, "\033[31mImpossible d'écrire dans le fichier '%s'\n\033[00m", "King.log");
		double2d_libere(obj->don);
		exit(EXIT_FAILURE);
	}
#endif
	obj->amas.rmax    = obj->amas.rmax * obj->amas.rc;
	obj->lig = taille;

	for(int j = 0; j < taille; j++)
	{
		obj->don[j][0] = obj->don[j][0] * obj->amas.rc;
		obj->don[j][3] = obj->don[j][3] * obj->amas.rho0;
	}

	double dr               = obj->don[1][0] - obj->don[0][0];
	obj->amas.Mtot          = 0.0;
	for(int j = 1; j < taille; j++)
		obj->amas.Mtot += 4.0*PI * obj->don[j][0] * obj->don[j][0] * obj->don[j][3] /*f(a + i*dr)*/ * dr;

	obj->amas.Mtot         += ( 4.0*PI * obj->don[0][0] * obj->don[0][0] * obj->don[0][3] + 4.0*PI * obj->don[taille - 1][0] * obj->don[taille - 1][0] * obj->don[taille - 1][3] ) * dr/2.0;
			//(obj->don[j][0] - obj->don[j-1][0]) * 4.0*PI*(obj->don[j][0] * obj->don[j][0] * obj->don[j][3] + obj->don[j-1][0] * obj->don[j-1][0] * obj->don[j-1][3])/2.0;		/** <- Intégration par la méthode des trapéze */
//	*Mtot                 *= (obj->don[j][0] - obj->don[j-1][0]);
//	*Mtot                 += (obj->don[j][0] - obj->don[j-1][0])/2.0 * 4.0*PI*(obj->don[taille-1][0] * obj->don[taille-1][3] + obj->don[0][0] * obj->don[0][3])

	*NbPart                = f(obj, data);
	obj->amas.m            = obj->amas.Mtot/(*NbPart);

        obj->amas.sigma2       = obj->amas.m * obj->amas.sigma2 / 2.0 /*Joules*/;
        //obj->amas.El           = Amas.amas.W0 * Amas.amas.sigma2 / (1-frac);
        obj->amas.El           = - obj->amas.G * obj->amas.m * obj->amas.Mtot / obj->amas.rmax;//-obj->frac*obj->amas.W0 * obj->amas.sigma2 / (obj->frac-1);


	// On redimensionne le probléme :
	for(int j = 0; j < taille; j++)
	{
		obj->don[j][1] = (obj->amas.El - obj->don[j][1] * obj->amas.sigma2)/obj->amas.m;
		obj->don[j][2] = -obj->amas.sigma2 * obj->don[j][2]/obj->amas.m;
		obj->don[j][4] = obj->don[j][0] * obj->don[j][0] * obj->don[j][2] / obj->amas.G; //King_mu(obj, 0, j);
#ifdef KING_DIM_LOG
		fprintf(fich, "%g\t%g\t%g\t%g\n", obj->don[j][0], obj->don[j][1], obj->don[j][2], obj->don[j][3]);
#endif
	}

#ifdef KING_DIM_LOG
	fclose(fich);
#endif

	obj->amas.vmax    = sqrt(2.0 * (obj->amas.El/obj->amas.m - obj->don[0][1]));

#ifdef DBG_Temp_FUNC
	fprintf(stderr, "%s::%s::%d :: Nb de ligne : %d, %p\n", __FILE__, __func__, __LINE__, obj->lig, &obj->lig);
#endif
}

void King_gud(King *obj)
{
	King_ugbs(obj);
	obj->amas.rmax    = obj->amas.rmax * obj->amas.rc;
}

void King_CalcPot(King *obj)
{
	for(int j = 0; j < obj->lig; j++)
		obj->don[j][1] = (obj->amas.El - obj->don[j][1] * obj->amas.sigma2)/obj->amas.m;
}

void King_CalcVMax(King *obj)
{
	obj->amas.vmax    = sqrt(2.0 * (obj->amas.El/obj->amas.m - obj->don[0][1]));
}

void King_CalcRMax(King *obj)
{
	obj->amas.rmax    = obj->amas.rmax * obj->amas.rc;
}

void King_CalcDPot(King *obj)
{
	for(int j = 0; j < obj->lig; j++)
		obj->don[j][2] = -obj->amas.sigma2 * obj->don[j][2]/obj->amas.m;
}

void King_SetM(King *obj, int NbPart)
{
	obj->amas.m = obj->amas.Mtot / NbPart;
}

void King_CalcMu(King *obj)
{
	for(int j = 0; j < obj->lig; j++)
		obj->don[j][4] = obj->don[j][0] * obj->don[j][0] * obj->don[j][2] / obj->amas.G; //King_mu(obj, 0, j);
}

void King_CalcEl(King *obj)
{
        obj->amas.El           = - obj->amas.G * obj->amas.m * obj->amas.Mtot / obj->amas.rmax;//-obj->frac*obj->amas.W0 * obj->amas.sigma2 / (obj->frac-1);
}

void King_CalcSig2(King *obj)
{
        obj->amas.sigma2       = obj->amas.m * obj->amas.sigma2 / 2.0 /*Joules*/;
}

void King_CalcMtot(King *obj)
{
	double dr               = obj->don[1][0] - obj->don[0][0];
	obj->amas.Mtot          = 0.0;
	for(int j = 1; j < obj->lig; j++)
		obj->amas.Mtot += 4.0*M_PI * obj->don[j][0] * obj->don[j][0] * obj->don[j][3] /*f(a + i*dr)*/ * dr;

	obj->amas.Mtot         += ( 4.0*M_PI * obj->don[0][0] * obj->don[0][0] * obj->don[0][3] + 4.0*M_PI * obj->don[obj->lig - 1][0] * obj->don[obj->lig - 1][0] * obj->don[obj->lig - 1][3] ) * dr/2.0;
}

void King_CalcR(King *obj)
{
	for(int i = 0; i < obj->lig; i++)
		obj->don[i][0] = obj->don[i][0] * obj->amas.rc;
}

void King_CalcRho(King *obj)
{
	for(int i = 0; i < obj->lig; i++)
		obj->don[i][3] = obj->don[i][3] * obj->amas.rho0;
}

double King_distrib(const King *obj, double E)
{
	if( E > obj->amas.El )
		return 0.0;

	return obj->amas.rho0 * pow(2.0 * M_PI * obj->amas.m * obj->amas.sigma2, -1.5)
			     * (exp((obj->amas.El - E)/obj->amas.sigma2) -1.0);
}

double King_don_pot(const King *obj, double r)
{
	if( r == 0.0 )
		return (obj->amas.El - obj->amas.W0*obj->amas.sigma2) / obj->amas.m;
	// Pour obtenir la valeur du potentiel en r, nous utilisons une interpolation linéaire :
	int r1 = floor(r/(get_pas() * obj->amas.rc)),
	    r2 = ceil(r/(get_pas() * obj->amas.rc));
	if( r2 == obj->lig )
	{
		r1--;
		r2--;
	}
	else if( r2 == 0 )
		r2++;
	if( r1 == r2 && r1 != 0 )
		r1--;
	else if( r1 == r2 && r1 == 0 )
		r2++;

	return obj->don[r1][1] + (r - obj->don[r1][0])*( obj->don[r2][1] - obj->don[r1][1] )
				/(obj->don[r2][0] - obj->don[r1][0]);
	// // potentiel[r1] + (posvit[i][R] - r1*dr)*( potentiel[r2] - potentiel[r1] ) /(r2*dr - r1*dr)
	//return obj.don[(int)round(r/(obj.amas.rc * get_pas()))][1];
}

