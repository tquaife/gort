#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<ctype.h>
#include<gortt.h>
#include<soil_rho.h>

int main( int argc, char **argv )
{

	gortt_geometry geom ;
	gortt_parameters params ;
	gortt_spectra spect ;
	gortt_control cont ;
	
	int i, j, nw_check, na_check, na=0 ;
	
	double x1, x2 ;
	
	char line[ MAX_LINE_LEN ], out_head[ MAX_LINE_LEN ], tmp_str1[ MAX_LINE_LEN ];
	
	FILE *prob_fp ;
	
	double *soil_vector_1, *soil_vector_2, *soil_vector_3, *soil_vector_4 ;
	double out_vza, out_vaa, out_sza, out_saa ;
	
	/*
	Initialise control structure
	*/
	
    cont.print_comp_spectra=FALSE ;
    cont.print_comp_proport=FALSE ;
    cont.calc_integrals=FALSE ;
    cont.use_q08_pn_kopen=FALSE ;
    cont.lidar_mode=FALSE ;

	spect.rsl1=0.2 ;
	spect.rsl2=0.1 ;
	spect.rsl3=0.03726 ;
	spect.rsl4=-0.002426 ;

	spect.is_user_leaf=FALSE ;
	spect.is_user_soil=FALSE ;
	spect.is_file_soil=FALSE ;
	
	
	/*
	default prospect values taken from the main.f90
	file in the PROSPECT-D main repositpry
	*/
	
	spect.p_N= 1.2;    // structure coefficient
	spect.p_Cab= 30.;  // chlorophyll content (µg.cm-2) 
	spect.p_Car= 10.;  // carotenoid content (µg.cm-2) 
	spect.p_Anth= 1.0; // anthocyanin content (µg.cm-2) 
	spect.p_Cbrown = 0.0; //brown pigment content (arbitrary units)
	spect.p_Cw = 0.015;	// EWT (cm)
	spect.p_Cm = 0.009;	// leaf mass per unit area (g.cm-2)
			
		
	geom.vza = 0. ;
	geom.vaa = 0. ;
	geom.sza = 0. ;
	geom.saa = 0. ;
	
	params.lambda = 0.405 ;
	params.r = 0.76 ;
	params.b = 3.55263 * params.r ;
	params.h1 = 3.0 ;
	params.h2 = 8.5 ;
	params.favd = 0.858 ;

	params.dz  = 0.20 ;	
	params.ds  = 0.20 ;
	params.dth = DTOR( 1 ) ;

	params.nlayers = 15 ;

	/*direct/diffuse term*/
	/*params.fd = 1.0 ;*/
	/*params.fd = 0.903 ;*/

	/*params.lad = LAD_EXTREMOPHILE ;*/
	params.lad = LAD_05 ;

	/*the following prob. 
	don't need changing*/
	params.maxcrowns = 30 ;
	params.nh_es = 20 ;

  /*number integral points*/
	params.npoints=32 ;

  /*diffuse radiation*/
	params.use_user_fd=FALSE ;
	
	gortt_cl_parser( argc, argv, &params, &geom, &spect, &cont ) ;
	
	
	/*read soil spectra if requested*/
	
	if( spect.is_file_soil == TRUE ) gortt_read_soil_lut( &spect );
	
	
	/*initialise some GORT stuff*/
	
	gortt_init_params( &params, &geom ) ;
	
	/**
	Three options for dealing with photon scattering
	probabilities:
	**/
	
	/*1) calculate own probabilities if not reading from */
	if( params.read_prob_file == FALSE )
	    if( cont.use_q08_pn_kopen == TRUE )
		    gortt_gap_probabilities_Q08( &params, &geom ) ;
        else
		    gortt_gap_probabilities( &params, &geom ) ;

  /*2) write probability file for later use*/
	if( params.write_prob_file == TRUE ){
		for( j=0; j<90; j++ )
			printf( "%d %0.40f %0.40f\n", j, params.p_n0[0][j], params.epgap[0][j] );
		printf( "-1 %0.40f %0.40f\n",  params.k_open[0], params.k_openep[0]  );
		return( EXIT_SUCCESS );
	}
	
  /*3) read existing probability file*/	
	if( params.read_prob_file == TRUE ){
		if( ( prob_fp = fopen( params.prob_fn, "r") ) == NULL ){
			fprintf( stderr, "%s: error opening probability file: %s\n", *argv, params.prob_fn );
			exit( EXIT_FAILURE );
		}
		while( fscanf( prob_fp, "%d %lf %lf", &j, &x1, &x2 ) == 3 ){
			if( j>=0 ){
				params.p_n0[0][j]=x1;
				params.epgap[0][j]=x2;
			}else{
				params.k_open[0]=x1;
				params.k_openep[0] =x2;
			}
		}
		fclose( prob_fp );
	}

	
	/**
	Parse first line in the BRDF input file
	**/
	
	if( fgets( line, MAX_LINE_LEN, stdin ) == NULL ){
		fprintf( stderr, "%s: error reading data on stdin\n", *argv );
		exit( EXIT_FAILURE );
	}
	strcpy( out_head, line );
	
	/*read number of angles*/
	if( get_first_string_token( line, tmp_str1 ) == 0 ){
		fprintf( stderr, "%s: error reading number of angles from line 1\n", *argv );
		exit( EXIT_FAILURE );
	}
	na_check = atoi( tmp_str1 );
		
	/*read number of wavelengths*/		
	if( get_first_string_token( line, tmp_str1 ) == 0 ){
		fprintf( stderr, "%s: error reading number of wavebands from line 1\n", *argv );
		exit( EXIT_FAILURE );
	}
	nw_check = atoi( tmp_str1 );
	
  /*alloc to hold wavelengths*/
	spect.wavelength = (double *)malloc( nw_check * sizeof( double ) );
	
	/*read wavelengths*/
	spect.nw = 0;
	while( get_first_string_token( line, tmp_str1 ) )
		spect.wavelength[ spect.nw++ ] = atof( tmp_str1 );
	
	if( nw_check != spect.nw ){
		fprintf( stderr, "%s: expected number of wavelengths (%d) does not match with number found (%d)\n", *argv, nw_check, spect.nw );
		exit( EXIT_FAILURE );
	}
	
	
	/*memory allocations*/
	spect.rsurf = (double *)malloc( nw_check * sizeof( double ) );
	spect.rsoil = (double *)malloc( nw_check * sizeof( double ) );
	spect.rleaf = (double *)malloc( nw_check * sizeof( double ) );
	spect.tleaf = (double *)malloc( nw_check * sizeof( double ) ); 

	/* alloc for spectra of each component at each waveband*/
	spect.scomp = (double *)malloc( 4 * nw_check * sizeof( double ) ); 

  /*allocs for integral terms*/
	if( cont.calc_integrals ){
	  spect.albedo = (double *)malloc( nw_check * sizeof( double ) );
	  spect.favegt = (double *)malloc( nw_check * sizeof( double ) );
	  spect.fasoil = (double *)malloc( nw_check * sizeof( double ) );
	  params.abscissa = (double *)malloc( params.npoints * sizeof( double ) );
	  params.weights  = (double *)malloc( params.npoints * sizeof( double ) );
	}

	
	/*set up integration */

  if( cont.calc_integrals ){
    	gauleg( -1., 1., params.abscissa, params.weights, params.npoints ); 	
    	/* Test print....
    	for (i=0;i<params.npoints;i++)
    	  printf( "%d %f %f\n", i, *(params.abscissa+i), *(params.weights+i) );
    	exit(0);
    	*/
	}
	
	/*set up spectral stuff here*/
	soil_vector_1 = default_soil_vector_1 ;
	soil_vector_2 = default_soil_vector_2 ;
	soil_vector_3 = default_soil_vector_3 ;
	soil_vector_4 = default_soil_vector_4 ;

	/*soil reflectance*/
	gortt_price_soil( &spect, soil_vector_1, soil_vector_2, soil_vector_3, soil_vector_4 );	

	/*leaf reflectance & transmittance*/
	gortt_prospect_interface( &spect );
	
	
	printf( "%s", out_head );

	while( fgets( line, MAX_LINE_LEN, stdin ) != NULL ){
	
		if( sscanf( line, "%lf %lf %lf %lf", &out_vza, &out_vaa, &out_sza, &out_saa ) != 4 ){
			fprintf( stderr, "%s: error on input, line %d\n", *argv, na+1 );
			exit( EXIT_FAILURE );
		}
		na++;
		
		geom.vza = DTOR( out_vza ) ;
		geom.vaa = DTOR( out_vaa ) ;
		geom.sza = DTOR( out_sza ) ;
		geom.saa = DTOR( out_saa ) ;

		/*
		correct angles for suitability with code (see walthall.c)
		make sure raa is correct too...
		*/

		/*corrections to the angles*/
				
		/* zenith angles +ve */
		if ( geom.sza < 0.0 ) {
			geom.saa += M_PI ;
			geom.sza *= -1.0 ;
		}
		if ( geom.vza < 0.0 ) {
			geom.vaa += M_PI ;
			geom.vza *= -1.0 ;
		}
		
		/* azimuth angles between 0 & 360 degrees */
		
		while ( geom.saa > 2 * M_PI )
			geom.saa -= 2 * M_PI;
   
		while ( geom.vaa > 2 * M_PI )
			geom.vaa -= 2 * M_PI;
   
		while ( geom.saa < 0 )
			geom.saa += 2 * M_PI;
   
		while ( geom.vaa < 0 )
			geom.vaa += 2 * M_PI;
		
				
				
		geom.raa = geom.saa - geom.vaa ;
		geom.raa = fabs((geom.raa - 2 * M_PI * (int) (0.5 + geom.raa * M_1_PI * 0.5)));



		geom.vza_prime = gortt_prime_theta( &params, geom.vza );
		geom.sza_prime = gortt_prime_theta( &params, geom.sza );

		/*this is needed by the kuusk function*/
		params.k_vza = gortt_leaf_angle_distribution( &params, geom.vza ) ;

		/*From Ni et al. '99*/
		if( ! params.use_user_fd )
		  params.fd = cos( geom.sza ) / ( cos( geom.sza ) + 0.09 );

		/*gortt*/
		gortt_set_zenith_dependant_probabilities( &params, &geom );
		gortt_rsurf( &params, &geom, &spect );

	
	/**
	LIDAR calculations
	**/
	
	
	gortt_lidar( &params, &geom, &spect );
	
	
    /*BRF outputs 
    n.b. the integral routines write over the BRFs so
    must print anything from previous calculations now
    */
		printf( "%f %f %f %f ", out_vza, out_vaa, out_sza, out_saa );
		for( i=0; i<spect.nw; i++ ){
			printf( "%f ", *( spect.rsurf + i ) );
      if( cont.print_comp_spectra==TRUE )
			  printf( "{ %f %f %f %f } ",*(spect.scomp+i*4+0),*(spect.scomp+i*4+1),*(spect.scomp+i*4+2),*(spect.scomp+i*4+3));
		}
    if( cont.print_comp_proport==TRUE )
		  printf( "[ %f %f %f %f ]", geom.Kc, geom.Kg, geom.Kt, geom.Kz );


		/*integral terms*/
		if( cont.calc_integrals ){
		  gortt_energy( &params, &geom, &spect );
	    for( i=0; i<spect.nw; i++ )
			  printf( "%f %f %f", *( spect.albedo + i ), *( spect.favegt + i ), *( spect.fasoil + i ) );
	  }
	
	  printf( "\n" );
	
	}

	if( na_check != na ){
		fprintf( stderr, "%s: expected number of angles (%d) does not match with number found (%d)\n", *argv, na_check, na );
		exit( EXIT_FAILURE );
	}	
	
	
	/*deallocations*/
	
	gortt_dealloc_1d( spect.wavelength );
	gortt_dealloc_1d( spect.rsurf );
	gortt_dealloc_1d( spect.rsoil );
	gortt_dealloc_1d( spect.rleaf );
	gortt_dealloc_1d( spect.tleaf );
	gortt_dealloc_1d( spect.scomp );

  if( cont.calc_integrals ){
	  gortt_dealloc_1d( spect.albedo );
	  gortt_dealloc_1d( spect.favegt );
	  gortt_dealloc_1d( spect.fasoil );
	  gortt_dealloc_1d( params.abscissa );
	  gortt_dealloc_1d( params.weights );

	}

	if( params.height    != NULL )gortt_dealloc_1d( params.height ); 
	if( params.height_p  != NULL )gortt_dealloc_1d( params.height_p );  
	if( params.theta     != NULL )gortt_dealloc_1d( params.theta );  
	if( params.theta_p   != NULL )gortt_dealloc_1d( params.theta_p );  

	if( params.vb        != NULL )gortt_dealloc_1d( params.vb ); 
	if( params.k_open    != NULL )gortt_dealloc_1d( params.k_open );  
	if( params.dk_open   != NULL )gortt_dealloc_1d( params.dk_open );  
	if( params.lk_up     != NULL )gortt_dealloc_1d( params.lk_up );  
	if( params.lk_down   != NULL )gortt_dealloc_1d( params.lk_down ); 
	if( params.k_openep  != NULL )gortt_dealloc_1d( params.k_openep ); 
	if( params.es        != NULL )gortt_dealloc_1d( params.es );  
	if( params.factorial != NULL )gortt_dealloc_1d( params.factorial );  

	if( params.fb        != NULL ) gortt_dealloc_2d( params.fb, params.nlayers );  
	if( params.t_open    != NULL ) gortt_dealloc_2d( params.t_open, params.nlayers ); 
	if( params.dt_open   != NULL ) gortt_dealloc_2d( params.dt_open, params.nlayers );  
	if( params.s_p       != NULL ) gortt_dealloc_2d( params.s_p, params.nlayers );  
	if( params.v_g       != NULL ) gortt_dealloc_2d( params.v_g, params.nlayers );  
	if( params.p_s0      != NULL ) gortt_dealloc_2d( params.p_s0, params.nlayers );  
	if( params.p_n0      != NULL ) gortt_dealloc_2d( params.p_n0, params.nlayers );  
	if( params.epgap     != NULL ) gortt_dealloc_2d( params.epgap, params.nlayers );  

	if( params.pd_s    	 != NULL ) gortt_dealloc_3d( params.pd_s, params.nlayers, params.nth );
	
	return( EXIT_SUCCESS );

}


void gortt_rsurf( p, g, s )
/*
Top level GORT function.

Kc = areal proportion of sunlit and viewed crown
Kg = areal proportion of the sunlit and viewed background
Kt = areal proportion of shaded and viewed crown
Kz = areal proportion of shaded and viewed background

C = signature of sunlit crown
T = signature of shaded crown
G = signature of sunlit soil
Z = signature of shaded soil

*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

	double Kc, Kg, Kt, Kz ;
	double C, G, T, Z;
	double Cd, Gd, Td, Zd;
	double Cf, Gf, Tf, Zf;
	double CdC, CdG, CdCG ;
	double CfC, CfG, CfCG ;

	double Kprime_g, Kprime_z ;

	int i ;
	

	
	/*
	======================
	Geometric "kernels"
	======================
	*/
	
	g->vza_prime = gortt_prime_theta( p, g->vza );
	g->sza_prime = gortt_prime_theta( p, g->sza );
	
	
	/*Sunlit & viewed ground*/
	Kg = gortt_kg( p, g ) ;
	
	
	/*sunlit & viewed crown*/
	Kc = gortt_kc( p, g, Kg ) ;
	/*Kc = gortt_kc_wenge( p, g, Kg ) ;*/
	/*Kc = gortt_kc_ambrals( p, g  ) ;*/
	
	
	/*shaded and viewd ground*/
	Kz = exp( -( p->lambda * M_PI * pow( p->r, 2 ) ) / cos( g->vza_prime ) ) - Kg ;
	
	
	/*shaded and viewd crown*/
	Kt = 1.0 - Kc - Kz - Kg ;
	Kt = MAX( 0.0, Kt );
	

	/*Kg and Kz normallised by crown shape*/
	Kprime_g = exp( -( p->lambda*M_PI*p->rr )/cos( g->sza_prime ) ) - Kg ;
	Kprime_z = 1.0 - exp( -( p->lambda*M_PI*p->rr )/cos( g->vza_prime ) ) - Kprime_g ;


	
	/*
	========================
	Wavelength dependent 
	bits go inside this loop
	========================
	*/
	
	for( i=0; i<s->nw; i++ ){
	
		s->n = i ;
	
		/*
		set up some spectral values for
		the efficiency (gamma used a lot).
		*/
	
		s->leaf_sscat_albedo = *( s->rleaf + i ) + *( s->tleaf + i ) ;
		s->leaf_gamma = sqrt( 1 - s->leaf_sscat_albedo ) ;
	
	
		/*sunlit ground
		 Could include brdf model here*.
		 
		 *However - note that gort_energy() assumes a
		 Lambertian background and so this would have 
		 to change to accomadate this correctly.
		*/
		
		Gd = *( s->rsoil + i ) ;
		Gf = *( s->rsoil + i ) ;
	
		G = p->fd*Gd + (1-p->fd)*Gf ;

	
		
		/*shaded ground*/
	
		/*Zd = ( gortt_t_prime_df( p, g, s ) + p->p_neq0_heq0_sza ) * *( s->rsoil + i ) ;*/
		Zd = ( gortt_t_prime_df( p, g, s ) + p->p_ngt0_heq0_sza ) * *( s->rsoil + i ) ;
		Zf = ( gortt_t_prime_ff( p, g, s ) - p->k_openep[ 0 ] ) * *( s->rsoil + i ) ;
			
		Z = p->fd*Zd + (1-p->fd)*Zf ;

	
		/*Sunlit Crown*/
		
		/*CdC = gortt_p_prime_df( p, g, s ) + ( ( 1.0-s->leaf_sscat_albedo ) \
			  * gortt_kuusk( p, g, s ) * s->leaf_sscat_albedo \
				* ( 1.0 - gortt_phase_function_assymetry( p, g, s ) ) ) \
				/ ( 2.0 * cos( g->sza ) * cos( g->vza ) ) ;*/
				
		CdC = gortt_p_prime_df( p, g, s ) + ( ( 1.0-s->leaf_sscat_albedo ) \
			  * gortt_kuusk( p, g, s ) * s->leaf_sscat_albedo \
				* ( 1.0 - gortt_phase_function_assymetry( p, g, s ) ) ) \
				/ ( 2.0 * cos( g->sza_prime ) * cos( g->vza_prime ) ) ;
		
						
		
		
		CfC = gortt_p_prime_ff( p, g, s ) ;
		
		CdG = ( Z * Kprime_z + G * Kprime_g ) * p->k_openep[ 0 ] ;
		
		CfG = ( ( p->k_openep[ 0 ]+p->k_open[ 0 ] )*G + \
					( 1-(p->k_openep[ 0 ]+p->k_open[ 0 ] ) )*Z ) * p->k_openep[ 0 ] ;
		
		CdCG = ( gortt_t_prime_df( p, g, s ) + gortt_t_prime_0( p, g, s ) ) \
				 * ( *( s->rsoil + i ) / ( 1.0 - *( s->rsoil + i ) * gortt_p_prime_ff( p, g, s ) ) ) \
				 * ( gortt_t_prime_ff( p, g, s ) - p->k_open[ 0 ] );
		
		CfCG = gortt_t_prime_ff( p, g, s ) \
				 * ( *( s->rsoil + i ) / ( 1.0 - *( s->rsoil + i ) * gortt_p_prime_ff( p, g, s ) ) ) \
				 * ( gortt_t_prime_ff( p, g, s ) - p->k_open[ 0 ] );
		
		
		Cd = CdC + CdG + CdCG ;
		Cf = CfC + CfG + CfCG ;
	
		C = p->fd*Cd + (1-p->fd)*Cf ;
		
		
		
		/*Shaded crown*/
		
		/*
		n.b. - this are the same as CdCG and CfCG
		*/
		
		Td = ( gortt_t_prime_df( p, g, s ) + gortt_t_prime_0( p, g, s ) ) \
			 * ( *( s->rsoil + i ) / ( 1.0 - *( s->rsoil + i ) * gortt_p_prime_ff( p, g, s ) ) ) \
			 * ( gortt_t_prime_ff( p, g, s ) - p->k_open[ 0 ] );
		
		Tf = gortt_t_prime_ff( p, g, s ) \
			 * ( *( s->rsoil + i ) / ( 1.0 - *( s->rsoil + i ) * gortt_p_prime_ff( p, g, s ) ) ) \
			 * ( gortt_t_prime_ff( p, g, s ) - p->k_open[ 0 ] );
 
		
		T = p->fd*Td + (1-p->fd)*Tf ;
	
	
		/*
		put the answers in the appropriate place...
		*/
	
		*( s->rsurf + i ) = Kc*C + Kg*G + Kt*T + Kz*Z ;
		

		/*fprintf( stderr, "%lf %lf %lf %lf\n", C, G, T, Z);*/
		
		*( s->scomp + i*4 + 0 ) = C ;
		*( s->scomp + i*4 + 1 ) = G ;
		*( s->scomp + i*4 + 2 ) = T ;
		*( s->scomp + i*4 + 3 ) = Z ;
	
	}


	g->Kc = Kc ;
	g->Kg = Kg ;
	g->Kt = Kt ;
	g->Kz = Kz ;


	return ;

}


double gortt_prime_theta( p, za )
gortt_parameters *p ;
double za ;
{

	return( atan( ( p->b/p->r ) * tan( za ) ) );

}



double gortt_leaf_angle_distribution( p, za )
gortt_parameters *p ;
double za ;
{

	double ag, bg, cg, dg ;
	
	ag = bg = cg = dg = 1 ;
	
	if( p->lad != LAD_05 ){
		fprintf( stderr, "gortt currently assumes a G function of 0.5\n" ) ;
		return( 0.5 );
	}
	
	if( p->lad == LAD_PLANOPHILE  ){
		dg = 0.0 ;
	}else if( p->lad == LAD_ERECTOPHILE  ){
		bg = -1.0 ;
		dg = 0.0 ;
	}else if( p->lad == LAD_PLAGIOPHILE  ){
		bg = -1.0 ;
		cg = 2.0 ;
		dg = 0.0 ;
	}else if( p->lad == LAD_EXTREMOPHILE ){
		cg = 2.0 ;
		dg = 0.0 ;
	}else if( p->lad == LAD_UNIFORM  	 ){
		ag = 0.0 ;
		bg = 0.0 ;
		cg = 0.0 ;
	}else if( p->lad == LAD_05  	 ){
		return( 0.5 );	
	}
	
	
	return( 2.0 / M_PI*( ag+bg*cos( 2.0*cg*za ) ) + dg*sin( za ) );

}


void gortt_init_params( p, g ) 
gortt_parameters *p ;
gortt_geometry *g ;

{

	int i, j, s_max ;


	p->ellipticity = p->b / p->r ;
	p->rr = p->r * p->r ;
	p->rrr = p->rr * p->r;
	p->h =  2.0 * p->r * p->ellipticity  + p->h2 - p->h1;

	/*
	the above = 2. * p->b + p->h2 - p->h1
	*/


	/*
	* figure out tau.
	*/ 
	
	p->k = gortt_leaf_angle_distribution( p, g->sza ) ;
	/*p->elai = p->favd * ( ( 4.0/3.0 ) * p->lambda * M_PI * p->ellipticity * p->rrr ) ; */
	p->elai = p->favd * ( ( 1.333333 ) * p->lambda * M_PI * p->ellipticity * p->rrr ) ;
	p->tau  =  p->k * p->favd ;


	/*
	* printf("g.FAVD=%7.4f g.ELAI=%7.4f",g.FAVD, g.ELAI) ; 
	* While h1 and h2 define the volume in which we find crown centers, we
	* need to integrate from h1-r to h2+r in order to consider all heights at
	* which scattering by foliage may occur.  These will be referred to as z1
	* and z2.
	*/

	p->z1 = p->h1 - p->r * p->ellipticity;
	p->z2 = p->h2 + p->r * p->ellipticity;
	p->lv = p->lambda / ( p->h2 - p->h1 );

	/*Corresponded variables in transformed space */

	p->favd_p = p->favd * p->ellipticity;
	p->tau_p  = p->k * p->favd_p;
	p->lv_p = p->lv * p->ellipticity;

	p->z1_p = p->z1 / p->ellipticity;
	p->z2_p = p->z2 / p->ellipticity;
	p->h1_p = p->h1 / p->ellipticity;
	p->h2_p = p->h2 / p->ellipticity;

	/*
	would the following be better done with a ceil()?
	(for the case where there is no remainder?)
	*/

	/*
	new version - fix nlayers instead of dz and ds
	*/
	/*p->nlayers = (int) ((p->z2 - p->z1) / p->dz) + 1;*/
 	
	
 	p->dz =  (double)(p->z2 - p->z1) / ((double) p->nlayers - 1.0 )  ;
	p->ds = p->dz ;
 	p->dz_p = p->dz / p->ellipticity ; 

 
  p->nh1 =   (int) ((p->h1 - p->z1) / p->dz) + 1 ;
	p->nh2 =   (int) ((p->h2 - p->z1) / p->dz) + 1 ;
   
	/*
	* The following variable relates to the discrete approximation to P(s) at
	* a given height.  The variable 'ds' is used in the conversions between
	* real values of 's' and array indices.
	*/

	/*
	* Figure out the number of angles theta (nth) for use in
	* various numerical integrations. 
	*/
   
	p->nth = (int) ( DTOR(90.0) / p->dth + 0.5 ) + 1;


	/*init all array pointers to NULL*/

	p->height = NULL;
	p->height_p = NULL ;
	p->theta = NULL ;
	p->theta_p = NULL ;
	p->vb = NULL ;
	p->k_open = NULL ;
	p->k_openep = NULL ;
	p->dk_open = NULL ;
	p->lk_up = NULL ;
	p->lk_down = NULL ;
	p->es = NULL ;
	p->factorial = NULL ;
	p->fb = NULL ;
	p->t_open = NULL ;
	p->dt_open = NULL ;
	p->s_p = NULL ;
	p->v_g = NULL;
	p->p_s0 = NULL ;
	p->p_n0 = NULL ;
	p->epgap = NULL ;
	p->pd_s = NULL ;


	/*
	 Do the factorial array
	*/

	/* p->factorial = ( double *)malloc( ( p->maxcrowns+1 ) * sizeof( double ) );   */
	if( !( p->factorial = ( double *)calloc( ( p->maxcrowns+1 ) , sizeof( double ) ) ) ){
		fprintf( stderr, "unable to alloc p->factorial\n" );
		exit( EXIT_FAILURE );
	}

	 p->factorial[0] = 1;
   for ( i = 1; i<=p->maxcrowns; i++) 
      p->factorial[i] = p->factorial[i - 1]*(double) i;
			

	/*
	* Initialize height and theta arrays.
	*/

	if( !( p->height = (double *)calloc( p->nlayers, sizeof( double ) ) ) ){
		fprintf( stderr, "unable to alloc p->height\n" );
		exit( EXIT_FAILURE );
	}
	if( !( p->height_p = (double *)calloc( p->nlayers, sizeof( double ) ) ) ){
		fprintf( stderr, "unable to alloc p->height_p\n" );
		exit( EXIT_FAILURE );
	}
	if( !( p->theta = (double *)calloc( p->nth, sizeof( double ) ) ) ){
		fprintf( stderr, "unable to alloc p->theta\n" );
		exit( EXIT_FAILURE );
	}
	if( !( p->theta_p = (double *)calloc( p->nth, sizeof( double ) ) ) ){
		fprintf( stderr, "unable to alloc p->theta_p\n" );
		exit( EXIT_FAILURE );
	}

	for ( i=p->nlayers-1; i>=0; i--) {
		p->height[i] = p->z2 - p->dz * (double)(p->nlayers-1-i);
		p->height_p[i] = p->height[i]/p->ellipticity ; 
	}

	for (i = 0; i < p->nth; i++) {
		p->theta[i] = p->dth * (double) i;
		if( p->theta[i] >= M_PI/2.0 ) { 
			p->theta[i]  = M_PI/2.0  - 1.0*M_PI/180.0;
		}

		/*
		* Get theta prime -- original angles adjusted by ellipticity.
		*/
		p->theta_p[i] = atan( tan(p->theta[i]) * p->ellipticity );
		if( p->theta_p[i] >= M_PI/2.0 ) { 
			p->theta_p[i] = M_PI/2.0 - 1.0*M_PI/180.0;
		}
	
	}


	/*
	Other allocations
	*/

	 p->vb       = gortt_alloc_1d( p->nlayers ); 
	 p->k_open   = gortt_alloc_1d( p->nlayers );  
	 p->dk_open  = gortt_alloc_1d( p->nlayers );  
	 p->lk_up    = gortt_alloc_1d( p->nlayers );  
	 p->lk_down  = gortt_alloc_1d( p->nlayers ); 
	 p->k_openep = gortt_alloc_1d( p->nlayers ); 
	 p->es       = gortt_alloc_1d( p->nlayers );  
	 p->fb       = gortt_alloc_2d( p->nlayers , p->nth );  
	 p->t_open   = gortt_alloc_2d( p->nlayers , p->nlayers ); 
	 p->dt_open  = gortt_alloc_2d( p->nlayers , p->nlayers );  
	 p->s_p      = gortt_alloc_2d( p->nlayers , p->nth );  
	 p->v_g      = gortt_alloc_2d( p->nlayers , p->nth );  
	 p->p_s0     = gortt_alloc_2d( p->nlayers , p->nth );  
	 p->p_n0     = gortt_alloc_2d( p->nlayers , p->nth );  
	 p->epgap    = gortt_alloc_2d( p->nlayers , p->nth );  


	/*
	* Allocate the probability distribution table.
	*/
	if( !(  p->pd_s = (double ***)calloc( p->nlayers, sizeof(double ** ) ) ) ){;
		fprintf(stderr, "unable to alloc p->pd_s\n");
		exit( EXIT_FAILURE );
	}
	 
	for( i=0; i<p->nlayers; i++ ){
      
		if( !(  p->pd_s[i] = (double **)calloc( p->nth, sizeof( double * ) ) ) ){;
			fprintf(stderr, "unable to alloc p->pd_s[%d]\n", i );
			exit( EXIT_FAILURE );
		}


		for ( j=0; j<p->nth; j++ ){

			/*
			* figure out maximum path length at this height and this angle.
			* Allocate enough memory to hold all possible s values.	 
			*/
 
			/*s_max = gortt_s_to_index( p ,(p->z2 - p->height[i]) / cos( p->theta_p[j] ) );*/

			s_max = gortt_s_to_index( p ,(p->z2_p - p->height_p[i]) / cos( p->theta_p[j] ) );

			/*
			fprintf(stderr,"s_max %d\n",s_max);
			fflush( stderr );
			*/
			
			p->pd_s[i][j] = gortt_alloc_1d( s_max + PD_S_BUFF ) ;
			
			/*
			if( !(  p->pd_s[i][j] = ( double *)calloc( ( s_max+PD_S_BUFF ) , sizeof( double ) ) ) ){;
				fprintf(stderr, "unable to alloc p->pd_s[%d]\n", i );
				exit( EXIT_FAILURE );
			}	
			*/

		}
   
	}


	return;
}



void gortt_set_zenith_dependant_probabilities( p, g ) 
/*
linearly interoplate the gap probabilities to get the
probs. at the required zenith angle.
Requires that p->p_n0 and p->epgap have been populated.
*/
gortt_parameters *p ;
gortt_geometry *g ;
{

	int cindex ;
	int findex ;
	double d, pos ;
		
		
	/*solar zenith*/	
	
	pos = fabs( g->sza ) / p->dth ;
	
	cindex = ceil( pos );
	findex = floor( pos );
	
	d = pos - findex ;
	
	p->p_neq0_heq0_sza = d * p->p_n0[ 0 ][ cindex ] + ( 1.0 - d ) * p->p_n0[ 0 ][ findex ] ;
	p->p_ngt0_heq0_sza = d * p->epgap[ 0 ][ cindex ] + ( 1.0 - d ) * p->epgap[ 0 ][ findex ] ;


	/*view zenith*/

	pos = fabs( g->vza ) / p->dth ;
	
	cindex = ceil( pos );
	findex = floor( pos );
	
	d = pos - findex ;
	
	p->p_neq0_heq0_vza = d * p->p_n0[ 0 ][ cindex ] + ( 1.0 - d ) * p->p_n0[ 0 ][ findex ] ;
	p->p_ngt0_heq0_vza = d * p->epgap[ 0 ][ cindex ] + ( 1.0 - d ) * p->epgap[ 0 ][ findex ] ;


	return ;

}




double * gortt_alloc_1d( n )
int n ;
{

	double *r;
	
	if (! ( r = (double *) calloc( n, sizeof(double) ) ) ) {
		fprintf(stderr, "Memory allocation failed (gortt_alloc_1d)\n");
		exit( EXIT_FAILURE );
	}

	return ( r );

}

void gortt_dealloc_1d( d )
double *d ;
{

	free( d );
	
	return;

}

double ** gortt_alloc_2d( n, m )
int n, m ;
{


	int q;
	double **r;

   
	if (! ( r = (double **) calloc( (unsigned) n, sizeof(double) ) ) ) {
		fprintf(stderr, "Memory allocation failed (gortt_alloc_2d)\n");
		exit( EXIT_FAILURE );
	}
	 
	for( q=0; q<n; q++)
		r[q] = gortt_alloc_1d( m );
	

	return( r );

}

void gortt_dealloc_2d( d, n )
double **d ;
int n;
{

	int q ;

	for( q=0; q<n; q++) gortt_dealloc_1d( d[ q ] );
	free( d );
	
	return;

}

void gortt_dealloc_3d( d, n, m )
double ***d ;
int n, m;
{

	int q, r ;

	for( q=0; q<n; q++){ 
		for( r=0; r<m; r++){ 
			free( d[ q ][ r ] );
		}
		free( d[ q ] );
	}
	
	free( d );
	
	return;

}



void gortt_cl_parser( argc, argv, p, g, s, c )
int argc;
char **argv;
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
gortt_control *c ;
{


	int i, use_true_p=FALSE, use_lai=FALSE;
	float hb=2.0, br=1.0, pcc=0.5, lai=2.0;
	
	p->read_prob_file=FALSE ;
	p->write_prob_file=FALSE;
	p->use_user_beta=FALSE ;
	
	
	
	for( i=1; i<argc; i++ ){

		if( *argv[ i ] == '-' ){
		
			/**/ if( !strncasecmp( argv[ i ], "-favd", 5 ) )   p->favd     = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-h1", 3 ) )     p->h1       = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-h2", 3 ) )     p->h2       = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-lambda", 7 ) ) p->lambda   = atof( argv[ ++i ] ); 


			else if( !strncmp( argv[ i ], "-HB", 3 ) )   { use_true_p=TRUE ;  hb = atof( argv[ ++i ] ); } 
			else if( !strncmp( argv[ i ], "-BR", 3 ) )   { use_true_p=TRUE ;  br = atof( argv[ ++i ] ); }
			else if( !strncmp( argv[ i ], "-PCC", 7 ) )  { use_true_p=TRUE ;  pcc = atof( argv[ ++i ] ); }
			
			else if( !strncmp( argv[ i ], "-LAI", 7 ) )  { use_lai=TRUE ;  lai = atof( argv[ ++i ] ); }

			
			else if( !strncasecmp( argv[ i ], "-beta", 5 ) ){ p->use_user_beta=TRUE; p->beta=atof( argv[ ++i ] );} 

			else if( !strncasecmp( argv[ i ], "-diffuse", 5 ) ){ p->use_user_fd=TRUE; p->fd=1.0-atof( argv[ ++i ] );} 
			
			else if( !strncmp( argv[ i ], "-alb_leaf", 9 ) ){
				s->is_user_leaf = TRUE ;
				s->user_r_leaf = atof( argv[ ++i ] ) ;
				
			}
			
			else if( !strncmp( argv[ i ], "-alb_soil", 9 ) ){
				s->is_user_soil = TRUE ;
				s->is_file_soil = FALSE ;
				s->user_r_soil = atof( argv[ ++i ] ) ;
				
			}

			else if( !strncmp( argv[ i ], "-soil_spectra", 10 ) ){
				s->is_user_soil = FALSE ;
				s->is_file_soil = TRUE ;
				strcpy( s->soil_file, argv[ ++i ] );
				
			}
			

			else if( !strncmp( argv[ i ], "-prnspec", 7 ) ) c->print_comp_spectra=TRUE;
			else if( !strncmp( argv[ i ], "-prnprop", 7 ) ) c->print_comp_proport=TRUE;			
			else if( !strncmp( argv[ i ], "-energy", 7 ) ) c->calc_integrals=TRUE;
			else if( !strncmp( argv[ i ], "-q08_pn_kopen", 7 ) ) c->use_q08_pn_kopen=TRUE;			
			else if( !strncmp( argv[ i ], "-lidar", 6 ) ) c->lidar_mode=TRUE;			
			
			else if( !strncmp( argv[ i ], "-P", 2 ) ){ p->read_prob_file=TRUE; strcpy( p->prob_fn, argv[ ++i ] );}
			else if( !strncmp( argv[ i ], "-W", 2 ) ) p->write_prob_file=TRUE;
			
			/*prospect*/
			else if( !strncasecmp( argv[ i ], "-N", 2 ) )   s->p_N     = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-cab", 4 ) ) s->p_Cab     = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-car", 4 ) ) s->p_Car     = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-canth", 3 ) )  s->p_Anth      = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-cbrown", 3 ) )  s->p_Cbrown      = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-cw", 3 ) )  s->p_Cw      = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-cm", 3 ) )  s->p_Cm      = atof( argv[ ++i ] ); 

			/*price*/
			else if( !strncasecmp( argv[ i ], "-rsl1", 5 ) )  s->rsl1    = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-rsl2", 5 ) )  s->rsl2    = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-rsl3", 5 ) )  s->rsl3    = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-rsl4", 5 ) )  s->rsl4    = atof( argv[ ++i ] ); 

			/*leave these at the bottom! (so as not to catch other opts)*/
			else if( !strncasecmp( argv[ i ], "-b", 2 ) )      p->b        = atof( argv[ ++i ] ); 
			else if( !strncasecmp( argv[ i ], "-r", 2 ) )      p->r        = atof( argv[ ++i ] ); 





			
			else if( !strncasecmp( argv[ i ], "-u", 2 ) ){  gortt_usage( argv[ 0 ] ); exit( EXIT_SUCCESS ) ; }
			
			else{
			
				fprintf( stderr, "%s: unknown option on command line: %s\n", argv[ 0 ], argv[ i ] );
				fprintf( stderr, "(use the option -u to see brief usage instructions)\n" );
				exit( EXIT_FAILURE );

			}

		}else{
		
			fprintf( stderr, "%s: unknown argument on command line: %s\n", argv[ 0 ], argv[ 1 ] );
			fprintf( stderr, "(use the option -u to see brief usage instructions)\n" );
			exit( EXIT_FAILURE );
		
		}

	}

	if( use_true_p ){
	
		p->r=10. ;
		p->b=br*p->r ;
		p->h1=p->b*2. ;
		p->h2=hb*p->b+p->h1 ;
		p->lambda=pcc/(p->r*p->r*M_PI) ;
	
	}

	if( use_lai ){
	
		p->favd=lai*3./(p->lambda*p->r*p->r*M_PI*p->b*4.0) ;
	
	}


	return ;

}



void gortt_usage( char *bin_name )
{


	fprintf( stderr, "usage: %s [options] < angles.dat\n\n", bin_name );
	fprintf( stderr, "The first line of the input data reads:\nN M W_1 W_2 [...] W_M\n"  );
	fprintf( stderr, "where N is the number of view--illumination geometries\nM is the number of wavelengths and\n"  );
	fprintf( stderr, "W_i (i=1,M) are the wavelengths at which to predict the canopy reflectance\n"  );
	fprintf( stderr, "The rest of the input data is four columns of ascii:\nview_zenith view_azimuth solar_zenith solar azimuth\n\n"  );
	

	fprintf( stderr, "The command line options are:\n"  );

	
  fprintf( stderr, "\n============ Crown geometry options:\n" );

	fprintf( stderr, "-beta arg\tforce the proportion of mutual shadowing to arg\n"  );
	fprintf( stderr, "         \t[n.b. if beta is not set gortt uses the model by Li and Strahler (IGARSS'92) to determine mutual shadowing]\n"  );

  fprintf( stderr, "\n------------ EITHER (old style):\n" );

	fprintf( stderr, "-h1 arg    \tset the lower boundary of the crown centres (m) to arg\n"  );
	fprintf( stderr, "-h2 arg    \tset the upper boundary of the crown centres (m) to arg\n"  );
	fprintf( stderr, "-b arg     \tset the vertical crown radius (m) to arg\n"  );
	fprintf( stderr, "-r arg     \tset the horizontal crown radius (m) to arg\n"  );
	fprintf( stderr, "-lambda arg\tset the tree stem density (1/m2) to arg\n"  );


  fprintf( stderr, "\n------------ OR (new style):\n" );
	fprintf( stderr, "-HB  arg\tset the ratio of the centroid height range to the vertical crown radius to arg\n"  );
	fprintf( stderr, "-BR  arg\tset the ratio of the vertical to horizontal crown radius to arg\n"  );
	fprintf( stderr, "-PCC arg\tset the projected crown cover (at nadir) to arg\n"  );


	fprintf( stderr, "\nThe above old style options refer to the original GORT paper and the \n"  );
	fprintf( stderr, "new style options are as they are expressed in the Quaife et al. (2008) \n"  );	
	fprintf( stderr, "DALEC paper. As soon as a new style option is specified all old style \n"  );
	fprintf( stderr, "options are ignored. Note that they are same thing - the code simply uses\n"  );
	fprintf( stderr, "the new style options to calculate the old ones but the new ones are \n"  );
	fprintf( stderr, "preferred because it reduces equifinality. \n"  );

  fprintf( stderr, "\n============ Amount of leaf material:\n" );
//  fprintf( stderr, "Two ways of specifying the amount of leaf material. Do not mix. \n" );
  fprintf( stderr, "------------ EITHER:\n" );
	fprintf( stderr, "-favd arg\tset the foliage volume area density (1/m2) within crown to arg\n"  );
  fprintf( stderr, "\n------------ OR:\n" );
	fprintf( stderr, "-LAI  arg\tset the leaf area index (m2/m2) for the scene to arg\n"  );

	fprintf( stderr, "\n============ Prospect leaf options:\n" );
	
	fprintf( stderr, "-N arg  \tset the leaf structure variable to arg\n"  );
	fprintf( stderr, "-Cab arg\tset the leaf chlorophyl content (µg.cm-2) to arg\n"  );
	fprintf( stderr, "-Cw arg \tset the equivelant leaf water thickness (cm) to arg\n"  );
	fprintf( stderr, "-Car arg \tset the carotenoid content (µg.cm-2) to arg\n"  );
	fprintf( stderr, "-Anth arg \tset the anthocyanin content (µg.cm-2) to arg\n"  );
	fprintf( stderr, "-Cbrown arg \tset the /brown pigment content (arbitrary units) to arg\n"  );
	fprintf( stderr, "-Cm arg \tset the leaf mass per unit area (g.cm-2) to arg\n"  );
	
	
	
	fprintf( stderr, "\n============ Price soil spectra options:\n" );
	fprintf( stderr, "-rsl1 arg \tset the weight of the first soil vector to arg\n"  );
	fprintf( stderr, "-rsl2 arg \tset the weight of the second soil vector to arg\n"  );
	fprintf( stderr, "-rsl3 arg \tset the weight of the third soil vector to arg\n"  );
	fprintf( stderr, "-rsl4 arg \tset the weight of the fourth soil vector to arg\n"  );
	
	fprintf( stderr, "\n============ User override for spectral properties:\n" );
	fprintf( stderr, "-alb_leaf     arg \tset leaf albedo to arg (turns prospect off)**\n"  );
	fprintf( stderr, "-alb_soil     arg \tset soil albedo to arg (turns price off)**\n"  );
  fprintf( stderr, "-soil_spectra arg \tread soil spectra from file arg (turns price off)\n"  );
	fprintf( stderr, "                  \t[** n.b. use only one of the above soil options]\n"  );



	fprintf( stderr, "\n============ Read/write gap probabilities:\n" );
	fprintf( stderr, "-W       \tgo as far as calculating the gap probabilities and write these to the stdout and exit\n"  );
	fprintf( stderr, "-P file  \tread gap probabilities from file that have been written using the -W option\n"  );
	fprintf( stderr, "         \tn.b. the above two options are included to allow fast BRF calculation based on\n"  );
	fprintf( stderr, "         \tpre-computed gap probabilities. The model must be run using the same crown and canopy\n"  );
	fprintf( stderr, "         \tgeometry options for the read and write. Spectral options (i.e. Prospect and Price options)\n"  );
	fprintf( stderr, "         \tmay be varied whist running using a give probability file.\n"  );

	fprintf( stderr, "\n============ Input/output options:\n" );


	fprintf( stderr, "-prnspec\tprint the scene component spectra for each wavelength inside {}\n" );
	fprintf( stderr, "-prnprop\tprint the viewed proportions of scene components inside []\n" );
	fprintf( stderr, "-energy \tprint the spectral albedo, absoption by veg and absorption by soil for each wavelength after other outputs.\n" );
	fprintf( stderr, "-u      \tprint this message and exit\n"  );

  fprintf( stderr, "\n" );

	return ;

}


char	get_first_string_token( line, token )
char	*line;
char	*token;
{

	int	i=0, j=0;

	/* Check that we are not already */
	/* at the end of the line        */

	if( *line == '\0' ) return( 0 );
	
	/* Skip leading white space and */
	/* return zero if end of string */
	/* is encountered               */
	
	while( isspace( *(line+i) ) ) 
		if( *(line+( ++i )) == '\0' ) return( 0 );
	
	
	/* Put the first block of non white 
	   space characters into element */
	
	while( ! isspace( *(line+i) ) ){
		if( ( token[j++] = *(line+(i++) ) ) == '\0' ){
		
			*line = '\0';
			return( 1 );
		
		}
	}
	
	
	/*put a NULL terminator in at j - (*) */
	
	token[ j ] = '\0' ;
	
	/* Copy remainder of line into the beggining */
	/* of itself */
	
	for( j=0; ( *(line+j)=*(line+i) ) != '\0' ; ++j, ++i );
	
	
	
	return( 1 );

}


void	gortt_price_soil( s, sv1, sv2, sv3, sv4 )
gortt_spectra *s ;
double *sv1, *sv2, *sv3, *sv4 ;
{


	int upper, lower, i ;
	double fraction ;
	
	double rs_lower, rs_upper ;
	
	for( i=0; i<s->nw; i++ ){
	
		if( *( s->wavelength +i ) < 400 || *( s->wavelength +i ) > 2500 ){
			fprintf( stderr, "gortt_price_soil: wavlength out of range (400-2500)\n" );
			exit( EXIT_FAILURE );
		}

	
		if( s->is_user_soil ){
		
			*(s->rsoil + i) = s->user_r_soil ;
		
		}else{
		
			upper = 1. + (*( s->wavelength +i )-400) / 5.0 ;
			lower = (*( s->wavelength +i )-400) / 5.0 ;
		
			fraction = (double) ( *( s->wavelength +i ) - 400. ) / 5.0 - lower ;
	
	
			rs_lower = s->rsl1 * *( sv1 + lower ) + s->rsl2 * *( sv2 + lower ) +  s->rsl3 * *( sv3 + lower ) + s->rsl4 * *( sv4 + lower ) ;
			rs_upper = s->rsl1 * *( sv1 + upper ) + s->rsl2 * *( sv2 + upper ) +  s->rsl3 * *( sv3 + upper ) + s->rsl4 * *( sv4 + upper ) ;
	
			*(s->rsoil + i) = rs_lower * ( 1 - fraction ) + rs_upper * fraction ;
		
		}

	}

	return ;

}


void gortt_prospect_interface( s )
gortt_spectra *s ;
{

	double  prospect_out[PROSPECT_NWBANDS*2];
	
	int upper, lower, i ;
	float fraction ;

    if( s->is_user_leaf == FALSE )
            prospect_DB_(&(s->p_N),&(s->p_Cab),&(s->p_Car),&(s->p_Anth),\
                         &(s->p_Cbrown),&(s->p_Cw),&(s->p_Cm), prospect_out);

	/*for( i=0; i<2100; i++ )
	        fprintf(stderr, "%d %f %f\n",400+i,prospect_out[i],prospect_out[i+PROSPECT_NWBANDS]);
	*/
	
	
	for( i=0; i<s->nw; i++ ){
		if( *( s->wavelength +i ) < 400 || *( s->wavelength +i ) > 2500 ){
			fprintf( stderr, "gortt_prospect_interface: wavlength out of range (400-2500)\n" );
			exit( EXIT_FAILURE );
		}

		if( s->is_user_leaf ){
			*(s->rleaf + i) = s->user_r_leaf / 2.0 ;
			*(s->tleaf + i) = s->user_r_leaf / 2.0 ;		
		
		}else{
			/*get the leaf optics for the desired wavelength*/	
			upper = 1 + (*( s->wavelength +i )-PROSPECT_LOWER_WL) / PROSPECT_SPECTRAL_RESLN ;
			lower = (*( s->wavelength +i )-PROSPECT_LOWER_WL) / PROSPECT_SPECTRAL_RESLN ;
		
			fraction = (float) ( *( s->wavelength +i )-PROSPECT_LOWER_WL) / PROSPECT_SPECTRAL_RESLN - lower ;
	
			*(s->rleaf + i) = prospect_out[lower] * (1-fraction) + prospect_out[upper] * fraction ;
			*(s->tleaf + i) = prospect_out[lower+PROSPECT_NWBANDS] * (1-fraction)\
			                + prospect_out[upper+PROSPECT_NWBANDS] * fraction ;
		}
	
	}
	
	return ;
}


/**
gortt_read_soil_lut:

read the soil spectra from the file specified on the command line.
First column should be wavelenght in nm and the second column should
be the soil albedo.

File must be in ascending order and cover the range 400<->2500. It
may have values outside this range and the sampling interval is arbitrary
and need not be constant.
**/
void	gortt_read_soil_lut( s )
gortt_spectra *s;
{


	FILE *fp ;
	char line[ MAX_LINE_LEN ];
	int n=0, i, index ;
	double this_wl, this_rs ;
	double last_wl, last_rs ;
	
	if( ( fp = fopen( s->soil_file, "r" ) ) == NULL ){
		fprintf( stderr, "gortt: cannot open file: %s\n", s->soil_file );
		exit( EXIT_FAILURE );
	}

	while( fgets( line, MAX_LINE_LEN, fp ) != NULL ){

		n++;
		if( sscanf( line, "%lf %lf", &this_wl, &this_rs ) != 2 ){
			fprintf( stderr, "gortt: error in soil file (%s), line %d\n", s->soil_file, n+1 );
			exit( EXIT_FAILURE );
		}
		
		if( n==1 && this_wl > 400 ){
			fprintf( stderr, "gortt: error in soil file (%s), first wavelength (%lf) should be <=400\n", s->soil_file, this_wl );
			exit( EXIT_FAILURE );
			
		}
		
		if( n>1 ){
		
			for( i=ceil( last_wl ); i<=floor( this_wl ); i++ ){
			
				index=i-400;
				
				if( (index>=0) && (index<=2100) )
					s->soil_spectra[ index ] = last_rs+(i-last_wl)/(this_wl-last_wl)*(this_rs-last_rs) ; 
			}
		
		}
		
		last_wl=this_wl;
		last_rs=this_rs;

	}

	if( last_wl < 2500 ){
		fprintf( stderr, "gortt: error in soil file (%s), last wavelength (%lf) should be >=2500\n", s->soil_file, last_wl );
		exit( EXIT_FAILURE );		
	}

	
	for( i=0;i<=2100;i++ ) printf( "%d %lf\n", i+400, s->soil_spectra[ i ] ) ;
	exit( EXIT_FAILURE );
	


	fclose( fp );

	return ;


}


