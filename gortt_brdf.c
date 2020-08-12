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


double gortt_kc_wenge( ppp, ggg, Kg )
/*
Stripped from code by Wenge and Conghe
I think this is wrong!
*/
gortt_parameters *ppp ;
gortt_geometry *ggg ;
double Kg ;
{

  float x,y,sinpv,sinpi,cospi,cospv,cost,temp;
  float h, b, beta, hbratio;
  float sunzn, sunaz,secti,sinti,costi,tanti;
  float vaz, vzn, mv, mi, htdif, D;
  float cosphi,sinphi,overpp0,overpp180,iv0;
  float iv180, f0, f180, kg0, kg180;
  float dist0, dist180, absdist0, absdist180, t0, t180, cost0, cost180;
  float bigF,bigF0, bigF180, m0, m180, thetai, thetav, fhat0;
  float fhat180, gammac180, pm0,pm180;
  float pvmv0, pvmv180, pimi0, pimi180, gamma0, gamma180, gammac0;
  float sintv, costv, sectv, tantv, sumax, thetami, thetamv;
  float gammax, gammai, gammav;
  float over, kg;
  float f, iv ; 
  
  
	double lambda0, vazpi;
 
  /*for C, G,Z */
  float  H;
 

	double h1, h2, phip, brratio, tip ;
	double lambda, radius, Kc ;
	
	lambda = ppp->lambda ;
	radius = ppp->r ;
	
	h1=ppp->h1 ;
	h2=ppp->h2 ;

  h = (h1+h2)/2.0;
 /* b = brratio*radius; */
 	brratio = ppp->b / ppp->r ;
 
 	b = ppp->b ;
  hbratio = ppp->h/ppp->b ;
  
  htdif = (h2-h1)/b; 
  H = htdif*b;
  
  
	
	/*sunzn = DTOR(SZN);*/
	sunzn = ggg->sza ;
  /*phip =  DTOR(SAZ);*/
  phip =  ggg->saa;
  tanti = tan(sunzn);

  tanti *= brratio;
  if(tanti < 0) tanti = 0.;
  tip = atan(tanti);

  sunaz=phip;
  thetai = tip;
  
  sinti = sin(thetai);
  costi = cos(thetai);
  secti = 1/costi;
  
  
  /* This is for hot-spot quantit, C and beta*/
  gammai = M_PI* secti;
  gamma0 = gammai;
  
  lambda0 = lambda*radius*radius;
  kg0 = exp( -lambda0*gammai );
  
  /* Calculate the beta - Based on Li and Strahler IGARSS92 paper and 
     formulation by Abdelgadir Abuelgasim (RSE 93 paper)*/
  
  D= brratio*htdif*sin(thetai*.5) / cos(thetai*.5);
  temp = lambda0*gammai ;
  if(temp < 0.00001) beta = 1.0;
  else
    beta=(lambda0*gammai/(lambda0*gammai+D))*
      ((1.0-exp(-(lambda0*gammai+D)))/(1.0-exp(-lambda0*gammai)));
  
  /* This is M - eqn 14*/
  temp = lambda0*gamma0 ;
  if(temp < 0.00001) m0 = 0.0;
  else m0 = 1.0-((1.0-kg0)/ temp);
  mi = m0;
  
  /* This is theta sub Mi - eqn 18*/ 
  if(mi < 0.) mi = 0.;
  thetami = 2.0*asin(sqrt(mi));
  /*
    if(thetami < PI/2.0) thetami = atan(tan(thetami)*radius/b);
    else thetami = PI - atan(tan(PI-thetami)*radius/b);
    */
  

  /* This is the Ac -- the physical surface area of an illuminated crown*/
  vzn =  ggg->vza;
  tantv = tan(vzn);
  tantv *= brratio;
  if(tantv < 0) tantv = 0.;
  tip=vzn;
  phip =  ggg->vaa ;
  
  
  thetav = tip;
  vaz = phip - sunaz;
  /*
    if(vaz > 0.) vaz *= -1;
    */
  if(vaz < 0.) vaz *= -1;
  cosphi = cos(vaz);
  sinphi = sin(vaz);
  
  sintv = sin(thetav);
  costv = cos(thetav);
  sectv = 1/costv;
  /* if(Iter == 100) printf("%f %f %f\n ",vaz,cosphi,sinphi); */
  x = tanti*tanti + tantv*tantv - 2*tanti*tantv*cosphi;
  if (x > 0.00001) {
    x = sqrt(x);
    /* x: distance between centers of two projection */
    sinpv = tanti*sinphi/x;
    cospv = (tanti* cosphi - tantv ) / x;
    sinpi = tantv*sinphi/x;
    cospi = (tantv* cosphi - tanti ) / x;
  }
  else {
    x = 0.;
    sinpv = 0;
    cospv = 1.;
    sinpi = 0;
    cospi = 1.;
  }
  temp =sinpv*sinpv + cospv*cospv*sectv*sectv;
  if (temp < 0.) temp = 0.;
  y = sqrt(temp);
  temp = sinpi*sinpi + cospi*cospi*secti*secti;
  if (temp < 0.) temp = 0.;
  y += sqrt(temp);
  /* sum of two half axis */
  cost = hbratio*x/y;
  if(cost < 0.) cost = 0.;
  if (cost >= 1.0) temp = 0.0;
  else temp = acos(cost);
  over = (temp - sin(temp)*cos(temp)) * (secti+sectv)/M_PI;
  /*  accurate PP, very good approximation on PC, reasonably
      good off PP/PC unless b/r very large AND h/b very small */
  sumax = secti+ sectv;
  
  /*Calculate iv explicitly --eqn 9 -- the illuminated portion 
    of the area of a single spheriod */
  iv=costi*costv +sinti*sintv*cosphi;
  
  /*Calculate F (with gamma and gammac) explicitly--eqn 12*/
  
  /* This is gamma - pg 14, at given sunzn and vzn, vaz*/
  gammax = M_PI*(sumax-over);
  
  /*Use this gammax in the kg calculation*/
  kg=exp(-lambda0*gammax);
  
  /* Calculate f0, bigF0 on PP at give sunzn and vzn*/
  
  dist0 = tanti - tantv;
  
  /*Absdist is the second part of eqn 4 -- (cos phi = +1 along PP)*/
  
  absdist0 = dist0*hbratio;
  if (absdist0 < 0.0) {
    absdist0 *= -1.0;
  }
  
  /* eqn 6*/
  cost0 = absdist0 /sumax;
  if(cost0 < 0.) cost0 = 0.;
  if (cost0 >= 1.0) t0 =0.0;
  else  t0 = acos(cost0);
  /* eqn 5 -- the exact overlap for the PP*/
  overpp0 = (t0-sin(2.0*t0)/2.0)*sumax;
  
  /* Calculate f180, bigF180 */
       
  dist180 = tanti + tantv;
  
  /*Absdist is the second part of eqn 4 -- (cos phi = -1 along PP)*/
  
  absdist180 = dist180*hbratio;
  
  /* eqn 6*/
  cost180 = absdist180 /sumax;
  if(cost180 < 0.) cost180 = 0.;
  if (cost180 >= 1.0) t180 = 0.0;
  else t180 = acos(cost180);
  
  /* eqn 5 -- the exact overlap for the PP*/
  overpp180 = (t180-sin(2.0*t180)/2.0)*sumax;
  
  /*Calculate the mutual shadowing effect along the PP*/
  
  /* This is gamma - pg 14*/
  gamma0 = M_PI*sumax- overpp0;
  gamma180 = M_PI*sumax- overpp180;
  
  /* the kg computation -- eqn 3 (the overlap should be Over divided by PI)*/
  kg0 = exp( -lambda0*gamma0);
  kg180 = exp( -lambda0*gamma180);
  
  /* the illuminated portion of the area of a single spheriod - eqn 9*/
  iv0 =  costi*costv + sinti*sintv;
  iv180 =  costi*costv - sinti*sintv;
  /* This is gamma sub v - pg 14*/
  gammav = M_PI* sectv;
  /* This is gamma sub i - pg 14*/
  gammai = M_PI* secti;
  /* This is gamma sub c - top of eqn 12*/
  gammac0 = 0.5 * (1+iv0)*gammav;
  gammac180 = 0.5 * (1+iv180)*gammav;
  
  /* F - eqn 12*/
  bigF0 = gammac0/gamma0;
  bigF180 = gammac180/gamma180;
  
  
  /* This is M - eqn 14*/
  temp = lambda0*gamma0 ;
  if(temp < 0.00001) m0 = 0.0;
  else
    m0 = 1.0-((1.0-kg0)/(temp));
  temp = lambda0*gamma180 ;
  if(temp < 0.00001) m180 = 0.0;
  else
    m180 = 1.0-((1.0-kg180)/(temp));
  
  /* This is Mv - eqn 11*/
  temp = lambda0*gammav ;
  if(temp < 0.00001) mv = 0.0;
  else
    mv = 1.0 - ((1.0 - exp(-temp))/(temp));
  
  /* This is Mi - eqn 10*/
  temp = lambda0*gammai ;
  if(temp < 0.00001) mi = 0.0;
  else
    mi = 1.0 - ((1.0 - exp(-temp))/(temp));
  
  /* This is theta sub Mi - eqn 18*/ 
  if(mi < 0.) mi = 0.;
  thetami = 2.0*asin(sqrt(mi));
  
  /* This is theta sub Mv - eqn 18*/
  if(mv < 0.) mv = 0.;
  thetamv = 2.0*asin(sqrt(mv));
  
  /*pvmv and pimi are eqns 4 and 5 in the IGARSS 1992 paper*/
  
  pvmv0=mv-(1-cos(thetav-thetai))/2.0;
  pvmv180=mv-(1-cos(-thetav-thetai))/2.0;
  pimi0=(1-cos(thetami*(1 - (thetai - thetav)/M_PI)))/2.0;
  pimi180=(1-cos(thetami*(1 - (thetai + thetav)/M_PI)))/2.0;
  pm0 = pvmv0;
  if(pimi0 > pvmv0) pm0 = pimi0;
  pm180 = pvmv180;
  if(pimi180 > pm180) pm180 = pimi180;
  
  /* This f is from eqn 17 and the first part of eqn 20*/
  if (thetav > thetai) { 
    f0= (1-(gammav*pvmv0/gammac0))/(1-m0);
  }
  else  {
    f0= (1-(gammav*pm0/gammac0))/(1-m0);
         }
  
  f180= (1-(gammav*pm180/gammac180))/(1-m180);  
  fhat0=beta*f0*bigF0 + (1.0-beta)*bigF0;      
  fhat180=beta*f180*bigF180 + (1.0-beta)*bigF180;
  
  /* This is gamma sub c - top of eqn 12*/
  gammac0 = 0.5 * (1+iv)*gammav;
  
  bigF = gammac0/gammax;
  /* The following is for preparing interpretation of Mutual shadowing */
  
  /*Linearly interpolate f (the non F portion)*/
  vazpi = vaz/M_PI;
  if(vazpi > 1.0) vazpi = 2.0-vazpi;
  f = (1.0-vazpi)*f0*bigF0+ (vazpi)*f180*bigF180;
  fhat0 = beta*f + (1.0-beta)*bigF;
  Kc = fhat0*(1.0-Kg);
	

	
	return( Kc );

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
