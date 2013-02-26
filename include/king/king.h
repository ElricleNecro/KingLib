#ifndef __KING_MODELE_POT__H_28082012_170000
#define __KING_MODELE_POT__H_28082012_170000

#include <stdlib.h>
#if defined __STDC__ && defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#include <tgmath.h>
#else
#include <math.h>
#endif
#include <stdio.h>

#include "mod.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

/**
 * \file king.h
 *
 * Ce fichier contient toutes les fonctions permettant de calculer des quantités physiques
 * associées au modèle de King.
 *
 * \author Guillaume Plum
 * \date Mardi 6 décembre 2011
 * \version 1.0
 *
 * \todo Rajouter les quantités telle la masse totale ou les vitesses et rayons maximums
 * \todo Fusionner avec mod.h
 */

/**
 * \struct king
 */
/**
 * \typedef King: Structure contenant toute les quantités physiques associées à un modèle de King.
 * \sa king
 */
typedef struct king
{
	Phys_King amas;		/*! Contient les paramètres physiques d'un King */
	double    **don;	/*! Contient les rayon, potentiel, dérivée du potentiel et densité volumique */
	double    frac;		/*! \deprecated Servait pour calculer l'énergie de libération */
	int       lig;		/*! Nombre de ligne du tableau don */
	int       col;		/*! Nombre de colonne du tableau don */
} King;

/**
 * \typedef gbs_cb: prototype de la fonction callback pour la fonction King_gbs_cb.
 */
typedef int(*gbs_cb)(King*,void*);

/*
extern int K_r    = 0;
extern int K_pot  = 1;
extern int K_dpot = 2;
extern int K_rho  = 3;
extern int K_mu   = 4;
*/

/**
 * Permet de récuperer la valeur du tableau de donnée de la structure King
 * pour la coordonée \f$(i,j)\f$.
 *
 * \param *obj Structure de type King
 * \param i Donnée physique (K_r, K_pot, K_dpot, K_rho, K_mu)
 * \param j Numero de ligne voulu
 *
 * \return La valeur voulu.
 *
 * \sa King
 */
double King_get_PhysVal(const King *obj, int i, int j);

/**
 * Crée une nouvelle structure de type King en calculant tous les paramètres nécessaire à partir de :
 * @param W0 La profondeur du puits de potentiel (sans dimension),
 * @param rc La longueur caractéristique (rayon de cœur) du modèle,
 * @param sig_v la dispersion de vitesse.
 */
King NewKing(const double W0, const double rc, const double sig_v);

/**
 * Calcul la température cinétique au rayon Amas.don[i][0]
 *
 * \param *Amas Structure de type King
 * \param i Indice pour obtenir le rayon à l'indice i
 *
 * \return la valeur de la température correspondante
 *
 * \sa King
 */
double King_Temp(const King *Amas, int i);

/**
 * Calcul de la température moyenne du Modèle de King.
 *
 * \param *Amas Structure de type King
 *
 * \return la température moyenne
 *
 * \sa King
 */
// * \param Mtot Masse totale de l'amas.
double King_TempMoy(const King *Amas/*, double Mtot*/);

/**
 * Calcul de \f$\mu(r)\f$, la masse contenue dans le rayon r.
 *
 * \param *Amas Structure de type King
 * \param a Indice de la borne inférieur (O généralement)
 * \param b Indice de la borne supérieur
 *
 * \return La masse contenue dans la sphére de rayon r
 *
 * \sa King
 */
double King_mu(const King *Amas, const int a, const int b);

/**
 * Calcul de \f$E_c(r)\f$
 *
 * \param *Amas Structure de type King
 * \param a Indice de la borne inférieur (O généralement)
 * \param b Indice de la borne supérieur
 *
 * \return L'énergie cinétique pour les particules contenues dans la sphére de rayon r
 *
 * \sa King
 */
double King_Ec(const King *Amas, const int a, const int b);

/**
 * Calcul de \f$E_p(r)\f$
 *
 * \param *Amas Structure de type King
 * \param a Indice de la borne inférieur (O généralement)
 * \param b Indice de la borne supérieur
 *
 * \return L'énergie cinétique pour les particules contenues dans la sphére de rayon r
 *
 * \sa King
 */
double King_Ep(const King *Amas, const int a, const int b);

/**
 * Résout les équations du modèles de King et donne les solutions à la solution obtenue. Il calcule aussi les bornes en distances et vitesses
 * maximum ainsi que la masse totale du système.
 *
 * \param[inout] *obj Structure de type King
 * \param[in]  NbPart	Nombre de particule du système
 *
 * \sa King
 */
void King_gbs(King *obj, const int NbPart);

//void King_gbs(King *obj, double *rmax, double *vmax, double *Mtot, const int NbPart, const double G);
//void get_box_size(King *obj, double *rmax, double *vmax);

/**
 * Résout les équations du modèles de King et donne les solutions à la solution obtenue. Il calcule aussi les bornes en distances et vitesses
 * maximum ainsi que la masse totale du système.
 *
 * \param[inout] *obj Structure de type King.
 * \param[inout] *NbPart Nombre de particule du système.
 * \param[in] f Fonction callback pour calculer le nombre de particule du King.
 * \param[in] *data Paramètre a passer en plus à la fonction callback.
 *
 * \sa King
 */
void King_gbs_cb(King *obj, int *NbPart, gbs_cb f /*double (*f)(King*, const int, void*)*/, void *data);

/**
 * Résout les équations du modèles de King et donne à la solution obtenue.
 *
 * \param[inout] *obj Structure de type King
 *
 * \sa King
 */
void King_ugbs(King *obj);

/**
 * Fonction de Distribution du modèle de King, adimensionnée :
 * \f{align}{
 *	\rho_0 \left(2m\pi\sigma^2\right)^{-3/2} \left(\mathrm{e}^{\dfrac{E_l - E}{\sigma^2}} - 1\right)
 * \f}
 *
 * \param[in] obj Structure de type King
 * \param[in] E Energie
 *
 * \return $f_K(E)$ si $E < E_l$, $0$ sinon
 *
 * \sa King
 */
double King_distrib(const King *obj, double E);

/**
 * Interpole le potentiel pour une valeur r du rayon.
 *
 * \param[in] obj Structure de type King
 * \param[in] E Energie
 *
 * \return Le potentiel en r
 *
 * \sa King
 */
double King_don_pot(const King *obj, double r);

/*
 * Cette fonction lit dans un fichier les paramètres du modèle de King à simuler
 * ainsi que la valeur de G, pour rester cohérent avec les unités donné par les autres variables dans le fichier.
 *
 * \param *Amas Structure de type King
 * \param[in]  *str	Nom du fichier à lire
 *
 * \return 1 ou un nombre négatif correspondant à l'erreur
 *
 * \sa King
 */
int read_utile(King *Amas, const char const *str);
//int read_utile(const char const *str, King *Amas, double *G);

/**
 * Résout les équations du modèles de King et donne les solutions à la solution obtenue.
 *
 * \param[inout] *obj Structure de type King
 * \param[in]  NbPart	Nombre de particule du système
 *
 * \sa King
 */
void King_gud(King *obj);

/**
 * Cette fonction calcule la densité volumique de masse du profil en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcRho(King *obj);

/**
 * Cette fonction redimensionne l'axe \f$\vec{r}\f$ en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcR(King *obj);

/**
 * Cette fonction calcule la Masse totale de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcMtot(King *obj);

/**
 * Cette fonction calcule la dispersion d'énergie de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcSig2(King *obj);

/**
 * Cette fonction calcule l'énergie de libération de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcEl(King *obj);

/**
 * Cette fonction calcule le profil de masse de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcMu(King *obj);

/**
 * Cette fonction calcule la dérivée du potentiel de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcDPot(King *obj);

/**
 * Cette fonction calcule le potentiel de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcPot(King *obj);

/**
 * Cette fonction calcule la vitesse maximale des particules de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcVMax(King *obj);

/**
 * Cette fonction calcule le rayon de l'objet en se servant des résultat de Resol_King.
 * \param *obj Structure de type King.
 */
void King_CalcRMax(King *obj);

/**
 * Cette fonction calcule la masse d'une particule à l'aide du nombre d'étoiles souhaité dans l'objet.
 * \param *obj Structure de type King.
 * \param[in] NbPart Nombre de Particule de l'objet.
 */
void King_SetM(King *obj, const int NbPart);

/**
 * S'occupe de désallouer la mémoire utilisée pour le modèle de King
 *
 * \param[in] *obj Structure de type King
 *
 * \sa King
 */
void King_free(King *obj);

#endif
