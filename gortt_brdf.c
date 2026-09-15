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


double gortt_overlap( const gortt_parameters *p, const gortt_geometry *g, double raa )
{

	double overlap, t, D, d ;
	double t1, t2, cos_t;


	/* As in the linear kernels*/
	d = pow( tan( g->sza_prime ),2 ) + pow( tan( g->vza_prime ),2 ) \
		- 2.0 * tan( g->sza_prime ) * tan( g->vza_prime ) * cos( raa ) ;

	D = sqrt( MAX( 0.0, d ) );

/*
	printf( "+++ %f %f %f\n", D, d, cos( raa ) );
*/

	if( TRUE ){ /*ambrals style*/

		t2 = sqrt( D*D + pow( ( tan(g->sza_prime)*tan(g->vza_prime)*sin(raa) ), 2) );

	}else{ /*Li&Strahler '92*/

		/*Principal plane only??*/
		t2 = fabs( tan( g->sza_prime) - tan(g->vza_prime)*cos(raa) ) ;

	}


	if( FALSE ){ /*from Wenge's code*/

		/*
		Might be required for off-principal plane...
		*/

		double tmp ;

		tmp = pow( tan( g->sza_prime ) * sin( raa ) / D, 2 ) \
			  + pow((tan( g->sza_prime ) * cos( raa ) - tan( g->vza_prime ))/D, 2 ) \
				* pow( SEC( g->sza_prime ), 2 );

		if( tmp < 0.0 ) tmp = 0.0 ;

		t1 = sqrt( tmp ) ;

		tmp = pow( tan( g->vza_prime ) * sin( raa ) / D, 2 ) \
			  + pow((tan( g->vza_prime ) * cos( raa ) - tan( g->sza_prime ))/D, 2 ) \
				* pow( SEC( g->vza_prime ), 2 );

		if( tmp < 0.0 ) tmp = 0.0 ;

		t1 += sqrt( tmp ) ;

	}else{ /*Li&Strhaler'92*/

		t1 = ( SEC(g->sza_prime) + SEC(g->vza_prime) ) ;

	}


	cos_t = (p->h/p->b) * t2 / t1  ;

	cos_t = MAX( -1.0, cos_t );
	cos_t = MIN(  1.0, cos_t );

	t = acos( cos_t );

	/*
	overlap = MAX( 0.0, ( t - sin( t )*cos( t ) )*( SEC(g->sza_prime)+SEC(g->vza_prime) )/M_PI ) ;
	*/

	overlap = MAX( 0.0, ( t - sin( t )*cos_t )*( SEC(g->sza_prime)+SEC(g->vza_prime) )/M_PI ) ;

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

	/*IGARSS 92 paper*/

	Gamma_i = Gamma_v ;

	PiMi = ( 1 - cos( theta_Mi * ( 1 - ( g->sza_prime - g->vza_prime * cos( raa ) )/M_PI ) ) ) / 2.0 ;
	PvMv = Mv -( 1.0 - cos( g->vza_prime * cos( raa ) - g->sza_prime ) ) / 2.0 ;


	/*set Po*/

	if( ( raa < DTOR( 270. ) ) && ( raa > DTOR( 90. ) ) ) Po = PvMv ;
	else if( fabs( g->vza ) > fabs( g->sza ) ) Po = PiMi  ;
	else Po = PvMv  ;

	if( g->sza_prime < 0.000000001 ){

		*beta = 0.0 ;

	}else{

		D = p->r * ( 1.0 / tan ( g->sza_prime/2.0 ) );
		*beta = (p->lambda*Gamma_i)/(p->lambda*Gamma_i+(p->h2-p->h1)/D)*(1.0-exp(-p->lambda*Gamma_i-(p->h2-p->h1)/D))/(1.0-exp(-p->lambda*Gamma_i)) ;

	}

	*f = *F * ( 1.0 - Gamma_v*( PvMv + PiMi - Po )/Gamma_c ) / ( 1.0 - M ) ;


	return ;
}


double gortt_kc_PP_only( const gortt_parameters *p, const gortt_geometry *g, double Kg )
{

	double Kc ;

	double phase_prime, Gamma, Gamma_v, Gamma_c ;
	double overlap, Mv, Mi, M, f, F, f1, f2  ;
	double Pv, PvMv, Pi, PiMi, Po ;
	double theta_Mi, theta_Mv, beta, D ;

	overlap = gortt_overlap( p, g, g->raa );
	phase_prime = cos( g->vza_prime )*cos( g->sza_prime )+sin( g->vza_prime )*sin( g->sza_prime )*cos( g->raa ) ;

	Mi = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )  ;
	Mv = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )  ;

	Gamma   = M_PI * p->rr * ( SEC(g->sza_prime)+SEC(g->vza_prime)-overlap ) ;
	Gamma_c = M_PI * p->rr *   SEC(g->vza_prime) * 0.5 * ( 1.0 + phase_prime ) ;
	Gamma_v = M_PI * p->rr *   SEC(g->vza_prime) ;

	F = Gamma_c / Gamma ;

	M = 1.0 - ( 1.0 - Kg )/( p->lambda * Gamma ) ;

	theta_Mi = acos( 1.0 - 2.0*Mi );
	theta_Mv = acos( 1.0 - 2.0*Mv );

	if( FALSE ){

		/*TGARSS '92 paper*/

		/*
		This is really not correct...!
		*/

		Pi = ( 1.0 - cos( theta_Mi - g->sza_prime + g->vza_prime * cos( g->raa ) ) ) / ( 1.0 - cos( theta_Mi ) ) ;

		if( ( theta_Mv - g->sza_prime + g->vza_prime * cos( g->raa ) ) >= 0 ){
			Pv = ( 1.0 - cos( theta_Mv - g->vza_prime + g->vza_prime * cos( g->raa ) ) ) / ( 1.0 - cos( theta_Mv ) ) ;
		}else{
			Pv = 0 ;
		}

		PiMi = Pi*Mi ;
		PvMv = Pv*Mv ;

		f = F * ( 1.0 - Gamma_v*( PvMv + PiMi - Po )/Gamma_c ) / ( 1.0 - M ) ;

	}else{

		/*IGARSS 92 paper*/

		double Gamma_i ;
		Gamma_i = Gamma_v ;

		PiMi = ( 1 - cos( theta_Mi * ( 1 - ( g->sza_prime - g->vza_prime * cos( g->raa ) )/M_PI ) ) ) / 2.0 ;
		PvMv = Mv -( 1.0 - cos( g->vza_prime * cos( g->raa ) - g->sza_prime ) ) / 2.0 ;


		/*
		Angles now redefined, so this is incorrect
		if( ( g->vza > 0 && g->sza < 0 )||( g->vza < 0 && g->sza > 0 ) ) Po = PvMv ;
		else if( fabs( g->vza ) > fabs( g->sza ) ) Po = PiMi ;
		else Po = PvMv  ;
		*/


		/*set Po*/

		if( ( g->raa < DTOR( 270 ) ) && ( g->raa > DTOR( 90 ) ) ) Po = PvMv ;
		else if( fabs( g->vza ) > fabs( g->sza ) ) Po = PiMi ;
		else Po = PvMv  ;


		D = p->r * ( 1.0 / tan ( g->sza_prime/2.0 ) );


		if( TRUE ){ /*IGARSS paper*/

			beta = (p->lambda*Gamma_i)/(p->lambda*Gamma_i+(p->h2-p->h1)/D)*(1.0-exp(-p->lambda*Gamma_i-(p->h2-p->h1)/D))/(1.0-exp(-p->lambda*Gamma_i)) ;

		}else{ /*RSE '94*/
			/*totally wrong as far as I can make out...*/
			beta = (p->lambda*Gamma_i)/(p->lambda*Gamma_i+(p->h2-p->h1)/D)*(1.0-exp(-p->lambda*Gamma_i+(p->h2-p->h1)/D))/(1.0-exp(p->lambda*Gamma_i)) ;
		}

		f1 = F * ( 1.0 - Gamma_v*( PvMv + PiMi - Po )/Gamma_c ) / ( 1.0 - M ) ;
		f2 = F ;

		f = beta*f1+(1.0-beta)*f2 ;

	}


	/*fprintf( stderr, "%f %f\n", Pv, Pi );*/

	Kc = f * ( 1.0 - Kg );

	return Kc ;

}


/*
Kuusk's canopy bi-directional gap probability
*/
double gortt_kuusk( const gortt_parameters *p, const gortt_geometry *g )
{

	double kuusk, lsza, lvza, lsv, cos_xi, H ;
	double t1, t2 ;

	cos_xi = cos( g->sza )*cos( g->vza )+sin( g->sza )*sin( g->vza )*cos( g->raa );
	/*05/06/06*/
	/*cos_xi = cos( g->sza_prime )*cos( g->vza_prime )+sin( g->sza_prime )*sin( g->vza_prime )*cos( g->raa );*/

	/*
	lsza = -log( g->p_ngt0_heq0_sza ) / ( p->k * p->favd ) ;
	lvza = -log( g->p_ngt0_heq0_vza ) / ( g->k_vza * p->favd ) ;
	*/

	lsza = -log( g->p_ngt0_heq0_sza ) / ( p->k * p->favd ) ;
	lvza = -log( g->p_ngt0_heq0_vza ) / ( g->k_vza * p->favd ) ;


	/*sqrt trap for -ve*/
	if( ( lsza*lsza + lvza*lvza - 2.*lsza*lvza*cos_xi ) > 0.0 ){
		lsv = sqrt( lsza*lsza + lvza*lvza - 2.*lsza*lvza*cos_xi );
		t2 = ( 1.0 - exp( -lsv/p->r ) ) / ( lsv/p->r ) ;
	}else{
		t2 = 1.0 ;
	}

	if( ( lsza * lvza ) > 0.0 ){
		t1 = sqrt( lsza * lvza ) ;
	}else{
		t1 = 0.0 ;
	}

	/*
	H = exp( p->k * p->favd * t1 * t2 ) ;
	*/


	H = exp( p->k * p->favd * t1 * t2 ) ;


	/*
	H = exp( p->k * p->favd * sqrt( li * lv ) * ( 1.0 - exp( -liv/p->r ) ) / ( liv/p->r ) );
	*/
	/*
	if( liv > 0.0 )
		H = exp( p->k * p->favd * sqrt( li * lv ) * ( 1.0 - exp( -liv/p->r ) ) / ( liv/p->r ) );
	else
		H = exp( p->k * p->favd * sqrt( li * lv ) );
	*/


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


