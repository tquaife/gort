#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gortt.h>


/*
n.b. the relative azimuth is passed in rather than taken from g: gortt_kc
needs Kg on the principal plane (raa = 0 and 180 deg) as well as at the
actual geometry, and it used to get them by writing to g->raa and putting
it back afterwards.
*/
double gortt_kg( const gortt_parameters *p, const gortt_geometry *g, double raa )
{

	double Kg, overlap ;

	overlap = gortt_overlap( p, g, raa );

	Kg = exp( -( p->lambda*pow( p->r, 2 )*M_PI*( SEC(g->sza_prime)+SEC(g->vza_prime)-overlap ) ) );

	return( Kg );

}


/*
O-bar: the overlap between the illumination and view shadows of a single
crown, in units of pi r^2, so that Kg = exp(-lambda pi r^2 (sec i' + sec v' - O)).
All angles are primed.

For a crown centred at height h (Li & Strahler 1992 eqs 5-6):

	O(h)   = ( t - sin(t) cos(t) ) secs / pi
	cos(t) = (h/b) sep / secs
	secs   = sec i' + sec v'
	sep    = sqrt( D^2 + (tan i' tan v' sin raa)^2 )
	D^2    = tan^2 i' + tan^2 v' - 2 tan i' tan v' cos raa

D is the separation of the two shadow centres per unit height, as in the
linear kernels, and sep is the effective separation the overlap depends on.
Its sin(raa) term (the Ambrals/LiSparse form) matters only off the principal
plane; in the principal plane sep = D = |tan i' - tan v' cos raa|, which is
L&S eq. 6.  O(h) falls to zero when cos(t) reaches 1, where the shadows stop
touching, and at the hotspot (sep = 0) it is secs/2 whatever the height.

L&S eq. 3 wants the overlap averaged over the centre heights, which gortt
takes as uniform in [h1,h2].  The average is a Gauss-Legendre sum over the
NH_OVERLAP heights p->h_overlap and weights p->w_overlap, set up in
gortt_init_params; only h varies inside the loop.

Previously a single h = 2b + h2 - h1 (the canopy depth) was used, and two
alternative forms of sep and secs sat here behind if(TRUE)/if(FALSE).  All
three were scored against a Monte Carlo and are worse: KC-INVESTIGATION.md,
section 4, and figures/goms_kernels.py, where they survive as switches.
*/
double gortt_overlap( const gortt_parameters *p, const gortt_geometry *g, double raa )
{

	double overlap, t, D, d ;
	double sep, secs, cos_t, ratio ;
	int i ;


	/* As in the linear kernels*/
	d = pow( tan( g->sza_prime ),2 ) + pow( tan( g->vza_prime ),2 ) \
		- 2.0 * tan( g->sza_prime ) * tan( g->vza_prime ) * cos( raa ) ;

	D = sqrt( MAX( 0.0, d ) );

	sep  = sqrt( D*D + pow( ( tan(g->sza_prime)*tan(g->vza_prime)*sin(raa) ), 2) );
	secs = SEC(g->sza_prime) + SEC(g->vza_prime) ;

	/*cos(t) = h * ratio: everything but the height, hoisted out of the loop*/
	ratio = sep / ( secs * p->b ) ;

	overlap = 0.0 ;
	for( i=0; i<NH_OVERLAP; i++ ){

		cos_t = p->h_overlap[ i ] * ratio ;

		cos_t = MAX( -1.0, cos_t );
		cos_t = MIN(  1.0, cos_t );

		t = acos( cos_t );

		overlap += p->w_overlap[ i ] * MAX( 0.0, ( t - sin( t )*cos_t )*secs/M_PI ) ;
	}

	return( overlap );

}

double gortt_kc_ambrals( const gortt_parameters *p, const gortt_geometry *g )
{

	double Kc, cos_xi_p ;

	cos_xi_p = cos( g->sza_prime )*cos( g->vza_prime )+sin( g->sza_prime )*sin( g->vza_prime )*cos( g->raa );

	Kc = ( 1.0 - exp( -p->lambda*M_PI*p->rr*SEC( g->vza_prime ) ) ) * 0.5 * ( 1.0 + cos_xi_p );

	return( Kc );

}


/*
This version interpolates between values of
f calulated on the principal plane to produce
the actual Kc value. Based on the method used
by Wenge in her code, which is Li & Strahler
(1992) IEEE TGRS 30(2), sec. H, eq. (25).
*/
double gortt_kc( const gortt_parameters *p, const gortt_geometry *g, double Kg )
{

	double Kc ;
	double beta, junk ;
	double f, f0, f180 ;
	double F, F0, F180 ;
	double Kg0, Kg180 ;
	double frac ;

	/*at the actual geometry: only F and beta are used, f is discarded*/
	gortt_kc_fFbeta( p, g, g->raa, Kg, &f, &F, &beta );

	/*and on the two halves of the principal plane*/
	Kg0 = gortt_kg( p, g, DTOR( 0. ) ) ;
	gortt_kc_fFbeta( p, g, DTOR( 0. ), Kg0, &f0, &F0, &junk );

	Kg180 = gortt_kg( p, g, DTOR( 180. ) ) ;
	gortt_kc_fFbeta( p, g, DTOR( 180. ), Kg180, &f180, &F180, &junk );

	frac = g->raa / M_PI ;
	if( frac > 1.0 ) frac = 2.0 - frac ;

	/*slot in user defined beta*/
	if( p->use_user_beta == TRUE ) beta = p->beta ;

	/*
	Li & Strahler (1992) sec. H, eq. (25): f is computed at raa = 0, pi
	and interpolated to the actual azimuth. gortt_kc_fFbeta already
	returns f with its factor of F in it (eq. 17), so it is not applied
	again here. This used to read f0*F0 and f180*F180, i.e. F squared.
	*/
	f = (1.-frac)*f0 + frac*f180 ;

	f = beta*f+(1.0-beta)*F ;

	Kc = f * ( 1.0 - Kg );

	return( Kc ) ;
}

/*
Only calculates f, F and beta of the kc function
this is done so we can interpolate bewteen various
f's to get the true value of Kc for OPP cases.
*/
void gortt_kc_fFbeta( const gortt_parameters *p, const gortt_geometry *g, double raa, double Kg, double *f, double *F, double *beta )
{

	double phase_prime, Gamma, Gamma_v, Gamma_c, Gamma_i ;
	double overlap, Mv, Mi, M ;
	double PvMv, PiMi, Po ;
	double theta_Mi, D ;

	/*double f, F, beta ; */

	overlap = gortt_overlap( p, g, raa );
	phase_prime = cos( g->vza_prime )*cos( g->sza_prime )+sin( g->vza_prime )*sin( g->sza_prime )*cos( raa ) ;

	Mi = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )  ;
	Mv = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )  ;

	Gamma   = M_PI * p->rr * ( SEC(g->sza_prime)+SEC(g->vza_prime)-overlap ) ;
	Gamma_c = M_PI * p->rr *   SEC(g->vza_prime) * 0.5 * ( 1.0 + phase_prime ) ;
	Gamma_v = M_PI * p->rr *   SEC(g->vza_prime) ;

	*F = Gamma_c / Gamma ;

	M = 1.0 - ( 1.0 - Kg )/( p->lambda * Gamma ) ;

	theta_Mi = acos( 1.0 - 2.0*Mi );

	/*IGARSS 92 paper: PvMv is its eq. 4, PiMi its eq. 5*/

	/*
	Gamma_i for beta: the shadow area a crown casts in the illumination
	direction. IGARSS '92 eq. 7 defines it as the shadowing area on a level,
	and Abuelgasim & Strahler (1994) call lambda*Gamma_i the "coverage and sun
	position". It used to be Gamma_v, which made beta depend on the view
	zenith. KC-INVESTIGATION.md section 8; OPEN-QUESTIONS #9.
	*/
	Gamma_i = M_PI * p->rr * SEC(g->sza_prime) ;

	PiMi = ( 1 - cos( theta_Mi * ( 1 - ( g->sza_prime - g->vza_prime * cos( raa ) )/M_PI ) ) ) / 2.0 ;
	PvMv = Mv -( 1.0 - cos( g->vza_prime * cos( raa ) - g->sza_prime ) ) / 2.0 ;


	/*
	Po: the part of the crown shaded from both directions. The smaller of the
	two mutually shadowed sets lies (almost) inside the larger, which is what
	Li & Strahler's cases 2(b) and 2(c) assume and what a Monte Carlo of the
	scene shows everywhere (KC-INVESTIGATION.md, section 2).

	It used to be PvMv throughout forward scatter (L&S p. 282). That cancels
	the view-direction shadowing, and once PvMv overtakes PiMi at large
	forward zeniths Kc grows like sec(vza') -- past 1 - Kg - Kz, where
	gortt_rsurf clamps Kt. With the minimum, Gamma_c/Gamma_v - PvMv =
	1 - Mv makes Kc <= 1 - Kg - Kz exact in the principal plane at beta = 1.
	*/
	Po = MIN( PvMv, PiMi ) ;

	if( g->sza_prime < 0.000000001 ){

		*beta = 0.0 ;

	}else{

		D = p->r * ( 1.0 / tan ( g->sza_prime/2.0 ) );
		*beta = (p->lambda*Gamma_i)/(p->lambda*Gamma_i+(p->h2-p->h1)/D)*(1.0-exp(-p->lambda*Gamma_i-(p->h2-p->h1)/D))/(1.0-exp(-p->lambda*Gamma_i)) ;

	}

	*f = *F * ( 1.0 - Gamma_v*( PvMv + PiMi - Po )/Gamma_c ) / ( 1.0 - M ) ;


	return ;
}


/*
Kuusk's bidirectional within-crown gap probability, as used in Ni et al.
(1999) eqs 3-7:

	P(n>0 | h=0, i, v) = P(n>0 | h=0, i) P(n>0 | h=0, v) H

	H         = exp( tau(i) * path_mean * hotspot )
	path_mean = sqrt( lsza * lvza )
	hotspot   = ( 1 - exp(-path_sep/r) ) / ( path_sep/r )
	path_sep  = sqrt( lsza^2 + lvza^2 - 2 lsza lvza cos(xi) )

tau(i) = k(i) FAVD.  lsza and lvza are the within-crown path lengths along the
sun and view directions, recovered by inverting Beer's law on the one-way
within-crown gap probabilities: l = -ln(P) / (k FAVD) (eqs 5-6).  path_sep is
the distance between the ends of the two paths (eq. 7), and the crown radius r
plays the role of Kuusk's correlation length (k_l = 1/r).

At the hotspot path_sep = 0 and hotspot = 1, so H = exp(tau path_mean), which
undoes the double counting of a shared path -- exactly, while k is the same for
the sun and view directions (true as long as the leaf angle distribution is
hardwired; OPEN-QUESTIONS #4-5).  As the directions separate, hotspot falls
towards 0 and H towards 1: the two gaps become independent.

n.b. cos(xi) uses the unprimed angles; a primed version is left commented out.
*/
double gortt_kuusk( const gortt_parameters *p, const gortt_geometry *g )
{

	double kuusk, lsza, lvza, path_sep, cos_xi, H ;
	double path_mean, hotspot ;

	cos_xi = cos( g->sza )*cos( g->vza )+sin( g->sza )*sin( g->vza )*cos( g->raa );
	/*05/06/06*/
	/*cos_xi = cos( g->sza_prime )*cos( g->vza_prime )+sin( g->sza_prime )*sin( g->vza_prime )*cos( g->raa );*/

	lsza = -log( g->p_ngt0_heq0_sza ) / ( p->k * p->favd ) ;
	lvza = -log( g->p_ngt0_heq0_vza ) / ( g->k_vza * p->favd ) ;


	/*sqrt trap for -ve; zero separation is the hotspot, where the limit is 1*/
	if( ( lsza*lsza + lvza*lvza - 2.*lsza*lvza*cos_xi ) > 0.0 ){
		path_sep = sqrt( lsza*lsza + lvza*lvza - 2.*lsza*lvza*cos_xi );
		hotspot = ( 1.0 - exp( -path_sep/p->r ) ) / ( path_sep/p->r ) ;
	}else{
		hotspot = 1.0 ;
	}

	if( ( lsza * lvza ) > 0.0 ){
		path_mean = sqrt( lsza * lvza ) ;
	}else{
		path_mean = 0.0 ;
	}

	H = exp( p->k * p->favd * path_mean * hotspot ) ;

	kuusk = g->p_ngt0_heq0_sza * g->p_ngt0_heq0_vza * H ;

	/*fprintf( stderr, "### %f\n", kuusk );*/

	return( kuusk );

}

/*
Fill in the radiative transfer terms for one geometry and one waveband.

These used to be a call tree of small functions, each recomputing its
dependencies from scratch every time it was asked for: one pass of the
waveband loop in gortt_rsurf evaluated the semi-infinite transmittance
T_inf_ff thirty-nine times for a value that does not change within it.
They are filled here instead, once, in dependency order -- which is the
order Ni et al. (1999) derive them in, so the block still reads against
the paper.
*/
void gortt_rt_init( const gortt_parameters *p, const gortt_geometry *g,
                    const gortt_spectra *s, int iw, gortt_rt *rt )
{

	double k, k_open ;


	/*
	leaf single scattering albedo, and the two quantities
	derived from it that are used throughout
	*/

	rt->w = *( s->rleaf + iw ) + *( s->tleaf + iw ) ;
	rt->gamma = sqrt( 1 - rt->w ) ;
	rt->g_asym = -(4.0/9.0)*( *(s->rleaf+iw) - *(s->tleaf+iw) )/( rt->w ) ;

	k = gortt_leaf_angle_distribution( p, g->sza );


	/*direct unscattered transmittance*/

	/*rt->t_0 = exp( -( k*p->favd*SEC(g->sza) ) ) ;*/
	/*05/06/06*/
	rt->t_0 = exp( -( k*p->elai*SEC( g->sza_prime ) ) ) ;


	/*
	direct-unscattered transmittance for a discontinuous
	plant canopy of height h at given geometry
	*/

	rt->t_prime_0 = g->p_neq0_heq0_sza + g->p_ngt0_heq0_sza ;


	/*
	------------------------------
	semi-infinite homogeneous layer
	------------------------------
	*/

	/*hemispherical-hemispherical reflectance*/

	rt->R_inf_ff = ( 1.0 - rt->gamma )/( 1.0 + rt->gamma ) ;

	/*directional-hemispherical reflectance*/

	/*rt->R_inf_df = ( 1.0 - rt->gamma )/( 1.0 + 2.0 *cos( g->sza )*rt->gamma ) ;*/
	/*05/06/06*/
	rt->R_inf_df = ( 1.0 - rt->gamma )/( 1.0 + 2.0 *cos( g->sza_prime )*rt->gamma ) ;

	/*hemispherical-hemispherical transmittance*/

	/*rt->T_inf_ff = exp(-( 2.0*rt->gamma*k*p->favd )) ;*/
	/*05/06/06*/
	rt->T_inf_ff = exp(-( 2.0*rt->gamma*k*p->elai )) ;

	/*directional-hemispherical transmittance*/

	rt->T_inf_df  = ( rt->w / 2.0 ) ;
	/*rt->T_inf_df *= ( 1. + 2. * cos( g->sza ) ) / ( 1. - pow( ( 2.*rt->gamma*cos( g->sza ) ), 2 ) )  ;*/
	/*05/06/06*/
	rt->T_inf_df *= ( 1. + 2. * cos( g->sza_prime ) ) / ( 1. - pow( ( 2.*rt->gamma*cos( g->sza_prime ) ), 2 ) )  ;
	rt->T_inf_df *= ( rt->T_inf_ff - rt->t_0 ) ;


	/*
	------------------------
	finite homogeneous layer
	------------------------
	*/

	/*hemispherical hemispherical transmittance*/

	rt->t_ff  = rt->T_inf_ff ;
	rt->t_ff *= ( 1. - pow( rt->R_inf_ff, 2 ) ) ;
	rt->t_ff /= ( 1. - pow( rt->R_inf_ff*rt->T_inf_ff, 2 ) );

	/*hemispherical-hemispherical reflectance*/

	rt->p_ff  = rt->R_inf_ff ;
	rt->p_ff *= ( 1. - pow( rt->T_inf_ff, 2 ) ) ;
	rt->p_ff /= ( 1. - pow( rt->T_inf_ff*rt->R_inf_ff, 2 ) )  ;

	/*directional hemispherical transmittance*/

	rt->t_df = rt->T_inf_df - rt->p_ff \
			* ( rt->t_0 * rt->R_inf_df + rt->T_inf_df * rt->R_inf_ff ) ;

	/*direction hemispherical reflectance*/

	rt->p_df = rt->R_inf_df - rt->t_ff * \
		 	(  rt->t_0     * rt->R_inf_df + \
				 rt->T_inf_df * rt->R_inf_ff ) ;


	/*
	--------------------
	discontinuous canopy
	--------------------
	*/

	/*
	directional hemispherical path transmittance
	for a discontinuous plant canopy of height h
	at given geometry
	*/

	rt->t_prime_df = rt->t_df * ( 1 - rt->t_prime_0 ) ;

	/*
	hemispherical-hemispherical path transmittance
	for a discontinuous plant canopy of height h
	at given geometry
	*/

	k_open = p->k_open[ 0 ] + p->k_openep[ 0 ] ;
	rt->t_prime_ff = rt->t_ff * ( 1.0 - k_open ) + k_open ;

	/*
	the path reflectances of a discontinuous canopy are taken
	to be those of the equivalent homogeneous one
	*/

	rt->p_prime_df = rt->p_df ;
	rt->p_prime_ff = rt->p_ff ;


	rt->kuusk = gortt_kuusk( p, g ) ;

	return ;

}


