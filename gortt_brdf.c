#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gortt.h>


double gortt_kg( p, g )
gortt_parameters *p ;
gortt_geometry *g ;
{

	double Kg, overlap ;

	overlap = gortt_overlap( p, g );
		
	Kg = exp( -( p->lambda*pow( p->r, 2 )*M_PI*( SEC(g->sza_prime)+SEC(g->vza_prime)-overlap ) ) );
	
	return( Kg );

}


double gortt_overlap( p, g )
gortt_parameters *p ;
gortt_geometry *g ;
{

	double overlap, t, D, d ;
	double t1, t2, cos_t;	
	
	
	/* As in the linear kernels*/
	d = pow( tan( g->sza_prime ),2 ) + pow( tan( g->vza_prime ),2 ) \
		- 2.0 * tan( g->sza_prime ) * tan( g->vza_prime ) * cos( g->raa ) ;
			
	D = sqrt( MAX( 0.0, d ) );

/*
	printf( "+++ %f %f %f\n", D, d, cos( g->raa ) );
*/

	if( TRUE ){ /*ambrals style*/		
	
		t2 = sqrt( D*D + pow( ( tan(g->sza_prime)*tan(g->vza_prime)*sin(g->raa) ), 2) );
		
	}else{ /*Li&Strahler '92*/
	
		/*Principal plane only??*/	
		t2 = fabs( tan( g->sza_prime) - tan(g->vza_prime)*cos(g->raa) ) ;
	
	}
	
	
	if( FALSE ){ /*from Wenge's code*/
	
		/*
		Might be required for off-principal plane...
		*/
	
		double tmp ;
	
		tmp = pow( tan( g->sza_prime ) * sin( g->raa ) / D, 2 ) \
			  + pow((tan( g->sza_prime ) * cos( g->raa ) - tan( g->vza_prime ))/D, 2 ) \
				* pow( SEC( g->sza_prime ), 2 );
	
		if( tmp < 0.0 ) tmp = 0.0 ;
		
		t1 = sqrt( tmp ) ;
	
		tmp = pow( tan( g->vza_prime ) * sin( g->raa ) / D, 2 ) \
			  + pow((tan( g->vza_prime ) * cos( g->raa ) - tan( g->sza_prime ))/D, 2 ) \
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

double gortt_kc_ambrals( p, g )
gortt_parameters *p ;
gortt_geometry *g ;
{

	double Kc, cos_xi_p ;

	cos_xi_p = cos( g->sza_prime )*cos( g->vza_prime )+sin( g->sza_prime )*sin( g->vza_prime )*cos( g->raa );

	Kc = ( 1.0 - exp( -p->lambda*M_PI*p->rr*SEC( g->vza_prime ) ) ) * 0.5 * ( 1.0 + cos_xi_p );

	return( Kc );

}


double gortt_kc( p, g, Kg )
/*
This version interpolates between values of
f calulated on the principal plane to produce
the actual Kc value. Based on the method used 
by Wenge in her code.
*/
gortt_parameters *p ;
gortt_geometry *g ;
double Kg ;
{
	
	double Kc ;
	double beta, junk ;
	double f, f0, f180 ;
	double F, F0, F180 ;
	double Kg0, Kg180 ;
	double raa_swap, frac ;
	
	
	raa_swap = g->raa ;
	
	
	gortt_kc_fFbeta( p, g, Kg, &f, &F, &beta );
	
	g->raa = DTOR( 0. );
	Kg0 = gortt_kg( p, g ) ;
	gortt_kc_fFbeta( p, g, Kg0, &f0, &F0, &junk );
	
	g->raa = DTOR( 180. );
	Kg180 = gortt_kg( p, g ) ;
	gortt_kc_fFbeta( p, g, Kg180, &f180, &F180, &junk );
	
	g->raa = raa_swap ;
	
	frac = g->raa / M_PI ;
	if( frac > 1.0 ) frac = 2.0 - frac ;
	
	/*slot in user defined beta*/
	if( p->use_user_beta == TRUE ) beta = p->beta ;	
	
	f = (1.-frac)*f0*F0 + frac*f180*F180 ;


	f = beta*f+(1.0-beta)*F ;

	Kc = f * ( 1.0 - Kg ); 



	return( Kc ) ;
}

void gortt_kc_fFbeta( p, g, Kg, f, F, beta )
/*
Only calculates f, F and beta of the kc function
this is done so we can interpolate bewteen various
f's to get the true value of Kc for OPP cases.
*/
gortt_parameters *p ;
gortt_geometry *g ;
double Kg ;
double *f ;
double *F ;
double *beta ;
{
	
	double phase_prime, Gamma, Gamma_v, Gamma_c, Gamma_i ;
	double overlap, Mv, Mi, M ;
	double PvMv, PiMi, Po ;
	double theta_Mi, theta_Mv, D ;
	
	/*double f, F, beta ; */
	
	overlap = gortt_overlap( p, g );
	phase_prime = cos( g->vza_prime )*cos( g->sza_prime )+sin( g->vza_prime )*sin( g->sza_prime )*cos( g->raa ) ;
	
	Mi = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )  ;
	Mv = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )  ;

	Gamma   = M_PI * p->rr * ( SEC(g->sza_prime)+SEC(g->vza_prime)-overlap ) ; 
	Gamma_c = M_PI * p->rr *   SEC(g->vza_prime) * 0.5 * ( 1.0 + phase_prime ) ;
	Gamma_v = M_PI * p->rr *   SEC(g->vza_prime) ;

	*F = Gamma_c / Gamma ;
	
	M = 1.0 - ( 1.0 - Kg )/( p->lambda * Gamma ) ;
	
	theta_Mi = acos( 1.0 - 2.0*Mi );
	theta_Mv = acos( 1.0 - 2.0*Mv );

	/*IGARSS 92 paper*/
	
	Gamma_i = Gamma_v ;
	
	PiMi = ( 1 - cos( theta_Mi * ( 1 - ( g->sza_prime - g->vza_prime * cos( g->raa ) )/M_PI ) ) ) / 2.0 ; 
	PvMv = Mv -( 1.0 - cos( g->vza_prime * cos( g->raa ) - g->sza_prime ) ) / 2.0 ;
		
		
	/*set Po*/
	
	if( ( g->raa < DTOR( 270. ) ) && ( g->raa > DTOR( 90. ) ) ) Po = PvMv ;
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


double gortt_kc_PP_only( p, g, Kg )
gortt_parameters *p ;
gortt_geometry *g ;
double Kg ;
{

	double Kc ;

	double phase_prime, Gamma, Gamma_v, Gamma_c ;
	double overlap, Mv, Mi, M, f, F, f1, f2  ;
	double Pv, PvMv, Pi, PiMi, Po ;
	double theta_Mi, theta_Mv, beta, D ;
	
	overlap = gortt_overlap( p, g );
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



double gortt_t_prime_df( p, g, s )
/*
directional hemispherical path transmittance
for a discontinuous plant canopy of height h
at given geometry
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double t_prime_df ;

	t_prime_df = gortt_t_df( p, g, s ) * ( 1 - gortt_t_prime_0( p, g, s ) ) ;

	return( t_prime_df );

}


double gortt_t_prime_ff( p, g, s )
/*
hemispherical-hemispherical path transmittance
for a discontinuous plant canopy of height h
at given geometry
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double t_prime_ff, k_open ;

	k_open = p->k_open[ 0 ] + p->k_openep[ 0 ] ;
	t_prime_ff = gortt_t_ff( p, g, s ) * ( 1.0 - k_open ) + k_open ;
	
	return( t_prime_ff );

}



double gortt_t_ff( p, g, s )
/*
hemispherical hemispherical transmittance
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double t_ff ;

	t_ff  = gortt_T_inf_ff( p, g, s ) ;
	t_ff *= ( 1. - pow( gortt_R_inf_ff( p, g, s ), 2 ) ) ;
	t_ff /= ( 1. - pow( gortt_R_inf_ff( p, g, s )*gortt_T_inf_ff( p, g, s ), 2 ) );

	return( t_ff );

}




double gortt_t_df( p, g, s )
/*
directional hemispherical transmittance
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double t_df ;

	t_df  = gortt_T_inf_df( p, g, s ) - gortt_p_ff( p, g, s ) \
				* ( gortt_t_0( p, g, s ) * gortt_R_inf_df( p, g, s ) + gortt_T_inf_df( p, g, s ) * gortt_R_inf_ff( p, g, s ) ) ;


	return( t_df );

}




double gortt_t_prime_0( p, g, s )
/*
direct-unscattered transmittance
for a discontinuous plant canopy of 
height h at given geometry
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double t_prime_0 ;
	
	t_prime_0 = p->p_neq0_heq0_sza + p->p_ngt0_heq0_sza ;
	
	return( t_prime_0 );

}


double gortt_T_inf_df( p, g, s )
/*
directional-hemispherical transmittance
(semi-infinite homogeneous layer)
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double T_inf_df ;


	T_inf_df  = ( s->leaf_sscat_albedo / 2.0 ) ;
	/*T_inf_df *= ( 1. + 2. * cos( g->sza ) ) / ( 1. - pow( ( 2.*s->leaf_gamma*cos( g->sza ) ), 2 ) )  ;*/
	/*05/06/06*/
	T_inf_df *= ( 1. + 2. * cos( g->sza_prime ) ) / ( 1. - pow( ( 2.*s->leaf_gamma*cos( g->sza_prime ) ), 2 ) )  ;
	T_inf_df *= ( gortt_T_inf_ff( p, g, s ) - gortt_t_0( p, g, s ) ) ;

	return( T_inf_df );

}

double gortt_T_inf_ff( p, g, s )
/*
hemispherical-hemispherical transmittance
(semi-infinite homogeneous layer)
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double k;

	k = gortt_leaf_angle_distribution( p, g->sza );
	/*return( exp(-( 2.0*s->leaf_gamma*k*p->favd )) );*/
	/*05/06/06*/
	return( exp(-( 2.0*s->leaf_gamma*k*p->elai )) );

}



double gortt_p_ff( p, g, s )
/*
hemispherical-hemispherical reflectance
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double p_ff ;


	p_ff  = gortt_R_inf_ff( p, g, s ) ;
	p_ff *= ( 1. - pow( gortt_T_inf_ff( p, g, s ), 2 ) ) ;
	p_ff /= ( 1. - pow( gortt_T_inf_ff( p, g, s )*gortt_R_inf_ff( p, g, s ), 2 ) )  ;	

	return( p_ff );

}



double gortt_t_0( p, g, s )
/*
direct unscattered transmittance
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double k ;

	k = gortt_leaf_angle_distribution( p, g->sza );
	/*return( exp( -( k*p->favd*SEC(g->sza) ) ) );*/
	/*05/06/06*/
	return( exp( -( k*p->elai*SEC( g->sza_prime ) ) ) );

}

double gortt_R_inf_df( p, g, s )
/*
directional-hemispherical reflectance 
(semi-infinite homogeneous layer)
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double R_inf_df ;

	/*R_inf_df = ( 1.0 - s->leaf_gamma )/( 1.0 + 2.0 *cos( g->sza )*s->leaf_gamma ) ;*/
	/*05/06/06*/
	R_inf_df = ( 1.0 - s->leaf_gamma )/( 1.0 + 2.0 *cos( g->sza_prime )*s->leaf_gamma ) ;

	return( R_inf_df );

}




double gortt_R_inf_ff( p, g, s )
/*
hemispherical-hemispherical reflectance 
(semi-infinite homogeneous layer)
*/

gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double R_inf_ff ;

	R_inf_ff = ( 1.0 - s->leaf_gamma )/( 1.0 + s->leaf_gamma ) ;

	return( R_inf_ff );

}




double gortt_phase_function_assymetry( p, g, s )
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double gfunc ;

	gfunc = -(4.0/9.0)*( *(s->rleaf+s->n) - *(s->tleaf+s->n) )/(s->leaf_sscat_albedo) ;

	return( gfunc );

}


double gortt_p_prime_df( p, g, s )
/*
Direction hemispherical path reflectance 
for a discontinuous canopy
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double p_p_df ;

	p_p_df = gortt_p_df(  p, g, s  );

	return( p_p_df );

}

double gortt_p_df( p, g, s )
/*
Direction hemispherical reflectance 
for a discontinuous canopy
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double p_df ;

	p_df = gortt_R_inf_df( p, g, s ) - gortt_t_ff( p, g, s ) * \
		 	(  gortt_t_0( p, g, s )      * gortt_R_inf_df( p, g, s ) + \
				 gortt_T_inf_df( p, g, s ) * gortt_R_inf_ff( p, g, s ) ) ;
	
	return( p_df );

}



double gortt_kuusk( p, g, s )
/*
Kuusk's canopy bi-directional gap probability
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double kuusk, lsza, lvza, lsv, cos_xi, H ;
	double t1, t2 ;
	
	cos_xi = cos( g->sza )*cos( g->vza )+sin( g->sza )*sin( g->vza )*cos( g->raa );
	/*05/06/06*/
	/*cos_xi = cos( g->sza_prime )*cos( g->vza_prime )+sin( g->sza_prime )*sin( g->vza_prime )*cos( g->raa );*/
	
	/*
	lsza = -log( p->p_ngt0_heq0_sza ) / ( p->k * p->favd ) ;
	lvza = -log( p->p_ngt0_heq0_vza ) / ( p->k_vza * p->favd ) ;
	*/
	
	lsza = -log( p->p_ngt0_heq0_sza ) / ( p->k * p->favd ) ;
	lvza = -log( p->p_ngt0_heq0_vza ) / ( p->k_vza * p->favd ) ;
	
	
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
	
	
	kuusk = p->p_ngt0_heq0_sza * p->p_ngt0_heq0_vza * H ;
	
	/*fprintf( stderr, "### %f\n", kuusk );*/
	
	return( kuusk );

}


double gortt_p_prime_ff( p, g, s )
/*
hemispherical hemispherical path reflectance 
for a discontinuous canopy
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double p_p_ff ;

	p_p_ff = gortt_p_ff(  p, g, s  );

	return( p_p_ff );

}
