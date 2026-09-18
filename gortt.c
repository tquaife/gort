#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<ctype.h>
#include<gortt.h>
#include<soil_rho.h>
/*
Driver-internal helpers. These are the decomposition of what used to be a
single ~370 line main(). They are the program's plumbing rather than part
of the model, so they are static here rather than declared in gortt.h.
*/

static void gortt_set_output_format( char *, int, char * );
static void gortt_set_defaults( gortt_parameters *, gortt_geometry *, gortt_spectra *, gortt_control *,
                                gortt_lidar_params * );
static void gortt_compute_probabilities( gortt_parameters *, gortt_control * );
static void gortt_init_squeeze_lai( gortt_parameters * );
static double gortt_bg_coupling( const gortt_parameters *, const gortt_rt *, double );
static void gortt_write_prob_file( gortt_parameters * );
static void gortt_read_prob_file( gortt_parameters *, char * );
static void gortt_read_header( FILE *, gortt_spectra *, int *, char *, char * );
static void gortt_lidar_spectra( gortt_spectra *, const gortt_lidar_params * );
static void gortt_alloc_outputs( gortt_parameters *, gortt_spectra *, gortt_control * );
static void gortt_free_outputs( gortt_parameters *, gortt_spectra *, gortt_control * );
static void gortt_free_params( gortt_parameters * );
static void gortt_set_geometry( gortt_geometry *, double, double, double, double );
static void gortt_run_geometry( gortt_parameters *, gortt_geometry *, gortt_spectra * );
static void gortt_print_result( gortt_parameters *, gortt_geometry *, gortt_spectra *,
                                gortt_control *, char *, double, double, double, double );


int main( int argc, char **argv )
{

	gortt_geometry geom ;
	gortt_parameters params ;
	gortt_spectra spect ;
	gortt_control cont ;
	gortt_lidar_params lid ;

	char line[ MAX_LINE_LEN ], out_head[ MAX_LINE_LEN ] ;
	char fmt[ 16 ] ;

	int na_check, na=0 ;
	double out_vza, out_vaa, out_sza, out_saa ;


	gortt_set_output_format( fmt, sizeof( fmt ), *argv );
	gortt_set_defaults( &params, &geom, &spect, &cont, &lid );

	gortt_cl_parser( argc, argv, &params, &spect, &cont, &lid ) ;


	/*read soil spectra if requested*/

	if( spect.is_file_soil == TRUE ) gortt_read_soil_lut( &spect );


	/*initialise some GORT stuff*/

	gortt_init_params( &params, &geom ) ;

	/**
	Three options for dealing with photon scattering
	probabilities:
	**/

	/*1) calculate own probabilities if not reading from file*/
	if( params.read_prob_file == FALSE )
		gortt_compute_probabilities( &params, &cont );

	/*2) write probability file for later use.
	n.b. this returns before stdin is read at all*/
	if( params.write_prob_file == TRUE ){
		gortt_write_prob_file( &params );
		gortt_free_params( &params );
		return( EXIT_SUCCESS );
	}

	/*3) read existing probability file*/
	if( params.read_prob_file == TRUE )
		gortt_read_prob_file( &params, *argv );

	/*
	Both of those settle K_open, which is the last thing -rt_crown_path's
	RT_PATH_SQUEEZE depth needs. It sits here rather than in either branch
	so that -P and a fresh run cannot disagree about it.
	*/
	gortt_init_squeeze_lai( &params ) ;


	/**
	Lidar mode short-circuits the rest of main() the same way -W does,
	and for the same kind of reason: the device is always at nadir, so
	no view/sun geometry on stdin would mean anything and stdin is never
	read at all. The spectral grid is the single lidar wavelength, which
	is all gortt_price_soil and gortt_prospect_interface need to give the
	waveform its leaf transmittance and soil reflectance.
	**/

	if( cont.lidar_mode == TRUE ){

		gortt_lidar_spectra( &spect, &lid );
		gortt_alloc_outputs( &params, &spect, &cont );

		gortt_price_soil( &spect, default_soil_vector_1, default_soil_vector_2,
		                          default_soil_vector_3, default_soil_vector_4 );
		gortt_prospect_interface( &spect );

		gortt_lidar( &params, &spect, &lid, fmt );

		gortt_free_outputs( &params, &spect, &cont );
		gortt_free_params( &params );
		return( EXIT_SUCCESS );

	}


	/**
	Parse first line in the BRDF input file
	**/

	gortt_read_header( stdin, &spect, &na_check, out_head, *argv );

	gortt_alloc_outputs( &params, &spect, &cont );

	/*
	Soil and leaf spectra do not depend on the viewing geometry,
	so are calculated once here rather than inside the loop.
	*/

	/*soil reflectance*/
	gortt_price_soil( &spect, default_soil_vector_1, default_soil_vector_2,
	                          default_soil_vector_3, default_soil_vector_4 );

	/*leaf reflectance & transmittance*/
	gortt_prospect_interface( &spect );


	printf( "%s", out_head );

	while( fgets( line, MAX_LINE_LEN, stdin ) != NULL ){

		if( sscanf( line, "%lf %lf %lf %lf", &out_vza, &out_vaa, &out_sza, &out_saa ) != 4 ){
			fprintf( stderr, "%s: error on input, line %d\n", *argv, na+1 );
			exit( EXIT_FAILURE );
		}
		na++;

		gortt_set_geometry( &geom, out_vza, out_vaa, out_sza, out_saa );

		/*gortt*/
		gortt_run_geometry( &params, &geom, &spect );

		gortt_print_result( &params, &geom, &spect, &cont, fmt,
		                    out_vza, out_vaa, out_sza, out_saa );

	}

	if( na_check != na ){
		fprintf( stderr, "%s: expected number of angles (%d) does not match with number found (%d)\n", *argv, na_check, na );
		exit( EXIT_FAILURE );
	}

	gortt_free_outputs( &params, &spect, &cont );
	gortt_free_params( &params );

	return( EXIT_SUCCESS );

}


/*
Output precision for the reflectance/albedo columns. The default "%f"
gives 4 significant figures, which is too coarse to see a small change
in the numerics. Setting GORTT_PREC to a number of significant digits
(17 round-trips a double exactly) raises it. Intended for regression
testing - see tests/run.sh - hence an environment variable rather than
a command line option.
*/
static void gortt_set_output_format( char *fmt, int len, char *bin_name )
{

	char *prec_env ;
	int prec ;

	strcpy( fmt, "%f " );

	if( ( prec_env = getenv( "GORTT_PREC" ) ) == NULL ) return ;

	prec = atoi( prec_env );

	if( prec < 1 || prec > 17 ){
		fprintf( stderr, "%s: GORTT_PREC must be between 1 and 17\n", bin_name );
		exit( EXIT_FAILURE );
	}

	/*n.b. built from an integer, never from the string itself*/
	snprintf( fmt, len, "%%.%dg ", prec );

	return ;

}


/*
Model defaults, applied before the command line is parsed so that any
option overrides them.
*/
static void gortt_set_defaults( gortt_parameters *p, gortt_geometry *g, gortt_spectra *s, gortt_control *c,
                                gortt_lidar_params *l )
{

	c->print_comp_spectra=FALSE ;
	c->print_comp_proport=FALSE ;
	c->calc_integrals=FALSE ;
	c->use_q08_pn_kopen=FALSE ;
	c->lidar_mode=FALSE ;

	s->rsl1=0.2 ;
	s->rsl2=0.1 ;
	s->rsl3=0.03726 ;
	s->rsl4=-0.002426 ;

	s->is_user_leaf=FALSE ;
	s->user_refl_leaf=0.0 ;
	s->user_trans_leaf=0.0 ;
	s->is_user_soil=FALSE ;
	s->is_file_soil=FALSE ;


	/*
	default prospect values taken from the main.f90
	file in the PROSPECT-D main repositpry
	*/

	s->p_N= 1.2;    /* structure coefficient */
	s->p_Cab= 30.;  /* chlorophyll content (µg.cm-2) */
	s->p_Car= 10.;  /* carotenoid content (µg.cm-2) */
	s->p_Anth= 1.0; /* anthocyanin content (µg.cm-2) */
	s->p_Cbrown = 0.0; /* brown pigment content (arbitrary units) */
	s->p_Cw = 0.015;	/* EWT (cm) */
	s->p_Cm = 0.009;	/* leaf mass per unit area (g.cm-2) */


	g->vza = 0. ;
	g->vaa = 0. ;
	g->sza = 0. ;
	g->saa = 0. ;

	p->lambda = 0.405 ;
	p->r = 0.76 ;
	p->b = 3.55263 * p->r ;
	p->h1 = 3.0 ;
	p->h2 = 8.5 ;
	p->favd = 0.858 ;

	/*n.b. dz and ds are not set here: gortt_init_params derives both
	from nlayers. They would need a default again if the commented-out
	"fix dz, derive nlayers" alternative there were ever reinstated.*/
	p->dth = DTOR( 1 ) ;

	p->nlayers = 15 ;

	/*the BRF only needs the gap probabilities at the ground; -lidar
	raises this to nlayers*/
	p->nh_probs = 1 ;

	/*direct/diffuse term*/
	/*p->user_fd = 1.0 ;*/
	/*p->user_fd = 0.903 ;*/

	/*p->lad = LAD_EXTREMOPHILE ;*/
	p->lad = LAD_05 ;

	/*the following prob.
	don't need changing*/
	p->maxcrowns = 30 ;

	/*multiplier on the NQ_* quadrature orders. 1 is the converged
	setting; -nquad exists to be able to show that it is.*/
	p->nquad = 1 ;

	/*number integral points*/
	p->npoints=32 ;

	/*diffuse radiation: without -diffuse the direct fraction is derived
	per geometry from the solar zenith, so user_fd is never read*/
	p->use_user_fd=FALSE ;
	p->user_fd=1.0 ;

	/*the experimental scattering formulations are all off: with none of
	the -rt_* options given the model is Ni et al. (1999) as published*/
	p->rt_crown_path=FALSE ;
	p->rt_crown_path_mode=RT_PATH_CROWN ;
	p->rt_recollision=FALSE ;
	p->rt_bg_coupling=FALSE ;
	p->rt_crown_crown=FALSE ;
	p->rt_fd_paper=FALSE ;
	p->rt_crown_brdf=FALSE ;
	p->rt_upscatter=FALSE ;

	/*the lidar instrument. 1064nm is the VCL/LVIS/SLICER line, and the
	one both Ni-Meister et al. (2001) and (2018) evaluate their rho_v at.
	J0 is a pure scale factor on both returns, so the default of one makes
	the output a fraction of the emitted pulse energy. binsize <= 0 means
	"use p->dz", which is not known until gortt_init_params has run.*/
	l->wl=1064.0 ;
	l->J0=1.0 ;
	l->pulse=0.0 ;
	l->binsize=-1.0 ;
	l->binned_output=FALSE ;

	return ;

}


static void gortt_compute_probabilities( gortt_parameters *p, gortt_control *c )
{

	if( c->use_q08_pn_kopen == TRUE )
		gortt_gap_probabilities_Q08( p ) ;
	else
		gortt_gap_probabilities( p ) ;

	return ;

}


/*
The third of the optical depths -rt_crown_path can pick (RT_PATH_SQUEEZE):
all the foliage spread over the part of the hemisphere that is not open,
rather than over everything (elai) or over the crowns' nadir footprint
(crown_lai). It needs K_open, so it cannot be done in gortt_init_params
with the other derived constants.

K_open -> 1 is an empty canopy and the quotient runs away; there is no
foliage to put anywhere in that limit, so fall back on elai.
*/
static void gortt_init_squeeze_lai( gortt_parameters *p )
{

	double k_open ;


	k_open = p->k_open[ 0 ] + p->k_openep[ 0 ] ;

	if( k_open < 1.0 - 1.0e-9 )
		p->squeeze_lai = p->elai / ( 1.0 - k_open ) ;
	else
		p->squeeze_lai = p->elai ;

	return ;

}


/*
Write the gap probabilities to stdout so that a later run can pick them
up with -P and skip recomputing them.

Two things about the format matter, and both used to be wrong:

The loop runs to nth, not to 90. The tables cover 0..90 degrees
inclusive, so there are nth = 91 rows; stopping at 90 left the last one
unwritten and a -P run read it back as the calloc zero.

The values are written as %.17g, not %0.40f. The gap probabilities fall
off superexponentially towards grazing -- by 88 degrees P(n=0) is around
7e-33 -- and %0.40f is *fixed* point, so everything below 1e-40 was
written as exactly zero. -P then took the log of it in gortt_kuusk and
returned -nan. %.17g round-trips a double exactly whatever its exponent.
*/
static void gortt_write_prob_file( gortt_parameters *p )
{

	int j ;

	for( j=0; j<p->nth; j++ )
		printf( "%d %.17g %.17g\n", j, p->p_n0[0][j], p->epgap[0][j] );

	printf( "-1 %.17g %.17g\n",  p->k_open[0], p->k_openep[0]  );

	return ;

}


static void gortt_read_prob_file( gortt_parameters *p, char *bin_name )
{

	FILE *fp ;
	int j ;
	double x1, x2 ;

	if( ( fp = fopen( p->prob_fn, "r") ) == NULL ){
		fprintf( stderr, "%s: error opening probability file: %s\n", bin_name, p->prob_fn );
		exit( EXIT_FAILURE );
	}

	while( fscanf( fp, "%d %lf %lf", &j, &x1, &x2 ) == 3 ){
		if( j>=0 ){
			p->p_n0[0][j]=x1;
			p->epgap[0][j]=x2;
		}else{
			p->k_open[0]=x1;
			p->k_openep[0] =x2;
		}
	}

	fclose( fp );

	return ;

}


/*
Read the first line of the input: the number of angles that should follow,
the number of wavebands, and then the wavelengths themselves. The line is
copied to out_head so that it can be echoed to the output.
*/
static void gortt_read_header( FILE *fp, gortt_spectra *s, int *na_check, char *out_head, char *bin_name )
{

	char line[ MAX_LINE_LEN ], tmp_str1[ MAX_LINE_LEN ];
	int nw_check ;

	if( fgets( line, MAX_LINE_LEN, fp ) == NULL ){
		fprintf( stderr, "%s: error reading data on stdin\n", bin_name );
		exit( EXIT_FAILURE );
	}
	strcpy( out_head, line );

	/*read number of angles*/
	if( get_first_string_token( line, tmp_str1 ) == 0 ){
		fprintf( stderr, "%s: error reading number of angles from line 1\n", bin_name );
		exit( EXIT_FAILURE );
	}
	*na_check = atoi( tmp_str1 );

	/*read number of wavelengths*/
	if( get_first_string_token( line, tmp_str1 ) == 0 ){
		fprintf( stderr, "%s: error reading number of wavebands from line 1\n", bin_name );
		exit( EXIT_FAILURE );
	}
	nw_check = atoi( tmp_str1 );

	/*alloc to hold wavelengths*/
	if( nw_check < 1 ){
		fprintf( stderr, "%s: number of wavebands on line 1 must be at least 1 (got %d)\n", bin_name, nw_check );
		exit( EXIT_FAILURE );
	}
	s->wavelength = (double *)malloc( nw_check * sizeof( double ) );

	/*read wavelengths.
	n.b. keep counting past nw_check, so that a header line listing more
	wavebands than it declares is still caught below, but do not store
	them - the array only has room for nw_check of them.*/
	s->nw = 0;
	while( get_first_string_token( line, tmp_str1 ) ){
		if( s->nw < nw_check )
			s->wavelength[ s->nw ] = atof( tmp_str1 );
		s->nw++ ;
	}

	if( nw_check != s->nw ){
		fprintf( stderr, "%s: expected number of wavelengths (%d) does not match with number found (%d)\n", bin_name, nw_check, s->nw );
		exit( EXIT_FAILURE );
	}

	return ;

}


/*
The lidar equivalent of gortt_read_header: lidar mode never reads stdin,
so the spectral grid is built from -lidar_wl instead of from a header
line. One waveband is all the waveform needs -- rho_v and rho_g are
evaluated at the laser line and nowhere else -- and setting it up this
way means gortt_alloc_outputs, gortt_price_soil, gortt_prospect_interface
and gortt_free_outputs are all reused untouched.
*/
static void gortt_lidar_spectra( gortt_spectra *s, const gortt_lidar_params *l )
{

	s->nw = 1 ;
	s->wavelength = (double *)malloc( s->nw * sizeof( double ) );
	s->wavelength[ 0 ] = l->wl ;

	return ;

}


static void gortt_alloc_outputs( gortt_parameters *p, gortt_spectra *s, gortt_control *c )
{

	/*memory allocations*/
	s->rsurf = (double *)malloc( s->nw * sizeof( double ) );
	s->rsoil = (double *)malloc( s->nw * sizeof( double ) );
	s->rleaf = (double *)malloc( s->nw * sizeof( double ) );
	s->tleaf = (double *)malloc( s->nw * sizeof( double ) );

	/* alloc for spectra of each component at each waveband*/
	s->scomp = (double *)malloc( 4 * s->nw * sizeof( double ) );

	/*written by gortt_rsurf whether or not an integral was asked for*/
	s->ztrans = (double *)malloc( s->nw * sizeof( double ) );

	/*allocs for integral terms, and set up the integration*/
	if( c->calc_integrals ){
		s->albedo = (double *)malloc( s->nw * sizeof( double ) );
		s->favegt = (double *)malloc( s->nw * sizeof( double ) );
		s->fasoil = (double *)malloc( s->nw * sizeof( double ) );
		p->abscissa = (double *)malloc( p->npoints * sizeof( double ) );
		p->weights  = (double *)malloc( p->npoints * sizeof( double ) );

		gauleg( -1., 1., p->abscissa, p->weights, p->npoints );
	}

	return ;

}


static void gortt_free_outputs( gortt_parameters *p, gortt_spectra *s, gortt_control *c )
{

	gortt_dealloc_1d( s->wavelength );
	gortt_dealloc_1d( s->rsurf );
	gortt_dealloc_1d( s->rsoil );
	gortt_dealloc_1d( s->rleaf );
	gortt_dealloc_1d( s->tleaf );
	gortt_dealloc_1d( s->scomp );
	gortt_dealloc_1d( s->ztrans );

	if( c->calc_integrals ){
		gortt_dealloc_1d( s->albedo );
		gortt_dealloc_1d( s->favegt );
		gortt_dealloc_1d( s->fasoil );
		gortt_dealloc_1d( p->abscissa );
		gortt_dealloc_1d( p->weights );
	}

	return ;

}


/*
The canopy arrays, which are allocated by gortt_init_params rather than
by gortt_alloc_outputs. Kept separate because the -W path exits before
stdin is read at all, so it has these but not the spectral outputs.
*/
static void gortt_free_params( gortt_parameters *p )
{

	if( p->height    != NULL )gortt_dealloc_1d( p->height );
	if( p->height_p  != NULL )gortt_dealloc_1d( p->height_p );
	if( p->theta     != NULL )gortt_dealloc_1d( p->theta );
	if( p->theta_p   != NULL )gortt_dealloc_1d( p->theta_p );

	if( p->zq        != NULL )gortt_dealloc_1d( p->zq );
	if( p->xq_gl     != NULL )gortt_dealloc_1d( p->xq_gl );
	if( p->wq_gl     != NULL )gortt_dealloc_1d( p->wq_gl );
	if( p->k_open    != NULL )gortt_dealloc_1d( p->k_open );
	if( p->k_openep  != NULL )gortt_dealloc_1d( p->k_openep );
	if( p->factorial != NULL )gortt_dealloc_1d( p->factorial );

	if( p->pn0q      != NULL ) gortt_dealloc_2d( p->pn0q, p->nq_ct + 1 );
	if( p->p_n0      != NULL ) gortt_dealloc_2d( p->p_n0, p->nlayers );
	if( p->epgap     != NULL ) gortt_dealloc_2d( p->epgap, p->nlayers );


	return ;

}


/*
Take one set of view and solar angles in degrees and put them into the
geometry structure in radians, with the sign and range conventions the
rest of the code assumes: zenith angles positive, azimuths in [0,2pi),
and the relative azimuth folded into [0,pi].

(see walthall.c)
*/
static void gortt_set_geometry( gortt_geometry *g, double vza, double vaa, double sza, double saa )
{

	g->vza = DTOR( vza ) ;
	g->vaa = DTOR( vaa ) ;
	g->sza = DTOR( sza ) ;
	g->saa = DTOR( saa ) ;

	/* zenith angles +ve */
	if ( g->sza < 0.0 ) {
		g->saa += M_PI ;
		g->sza *= -1.0 ;
	}
	if ( g->vza < 0.0 ) {
		g->vaa += M_PI ;
		g->vza *= -1.0 ;
	}

	/* azimuth angles between 0 & 360 degrees */

	while ( g->saa > 2 * M_PI )
		g->saa -= 2 * M_PI;

	while ( g->vaa > 2 * M_PI )
		g->vaa -= 2 * M_PI;

	while ( g->saa < 0 )
		g->saa += 2 * M_PI;

	while ( g->vaa < 0 )
		g->vaa += 2 * M_PI;


	g->raa = g->saa - g->vaa ;
	g->raa = fabs((g->raa - 2 * M_PI * (int) (0.5 + g->raa * M_1_PI * 0.5)));

	return ;

}


/*
Evaluate the model for the geometry currently in g, leaving the
reflectances in s.
*/
static void gortt_run_geometry( gortt_parameters *p, gortt_geometry *g, gortt_spectra *s )
{

	g->vza_prime = gortt_prime_theta( p, g->vza );
	g->sza_prime = gortt_prime_theta( p, g->sza );

	/*this is needed by the kuusk function*/
	g->k_vza = gortt_leaf_angle_distribution( p, g->vza ) ;

	/*From Ni et al. '99*/
	if( p->use_user_fd )
		g->fd = p->user_fd ;
	else if( p->rt_fd_paper )
		g->fd = cos( g->sza ) / ( cos( g->sza ) + 0.05 );
	else
		g->fd = cos( g->sza ) / ( cos( g->sza ) + 0.09 );

	gortt_set_zenith_dependant_probabilities( p, g );
	gortt_rsurf( p, g, s );

	return ;

}


/*
One output line: the angles as they were given, then the BRF per waveband,
then whatever the -prnspec/-prnprop/-energy options ask for.

n.b. the integral routines write over the BRFs, so everything from the
previous calculation has to be printed before gortt_energy is called.
*/
static void gortt_print_result( gortt_parameters *p, gortt_geometry *g, gortt_spectra *s,
                                gortt_control *c, char *fmt,
                                double vza, double vaa, double sza, double saa )
{

	int i ;

	/*the angles are echoes of the input, so are left at %f*/
	printf( "%f %f %f %f ", vza, vaa, sza, saa );

	for( i=0; i<s->nw; i++ ){
		printf( fmt, *( s->rsurf + i ) );
		if( c->print_comp_spectra==TRUE ){
			printf( "{ " );
			printf( fmt, *(s->scomp+i*4+0) );
			printf( fmt, *(s->scomp+i*4+1) );
			printf( fmt, *(s->scomp+i*4+2) );
			printf( fmt, *(s->scomp+i*4+3) );
			printf( "} " );
		}
	}

	if( c->print_comp_proport==TRUE ){
		printf( "[ " );
		printf( fmt, g->Kc );
		printf( fmt, g->Kg );
		printf( fmt, g->Kt );
		printf( fmt, g->Kz );
		printf( "] " );
	}

	/*integral terms*/
	if( c->calc_integrals ){
		gortt_energy( p, g, s );
		for( i=0; i<s->nw; i++ ){
			printf( fmt, *( s->albedo + i ) );
			printf( fmt, *( s->favegt + i ) );
			printf( fmt, *( s->fasoil + i ) );
		}
	}

	printf( "\n" );

	return ;

}


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
/*
The background-canopy coupling factor 1/(1 - rho_s rho'_ff), or one when
-rt_bg_coupling is off. See the comment beside G in gortt_rsurf.

rho_s rho'_ff is a product of two reflectances and so is below one for any
physical pair, but clamp anyway: a user-supplied soil albedo is not checked
against one anywhere, and a runaway here would be silent.
*/
static double gortt_bg_coupling( const gortt_parameters *p, const gortt_rt *rt, double rsoil )
{

	double x ;


	if( p->rt_bg_coupling == FALSE ) return( 1.0 );

	x = rsoil * rt->p_prime_ff ;

	return( ( x < 1.0 ) ? 1.0 / ( 1.0 - x ) : 1.0 );

}


void gortt_rsurf( const gortt_parameters *p, gortt_geometry *g, gortt_spectra *s )
{

	double Kc, Kg, Kt, Kz ;
	double C, G, T, Z;
	double Cd, Gd, Zd;
	double Cf, Gf, Zf;
	double CdC, CdG, CdCG ;
	double CfC, CfG, CfCG ;

	double Kprime_g, Kprime_z ;
	double bg, cc ;

	gortt_rt rt ;

	int i ;



	/*
	======================
	Geometric "kernels"
	======================
	*/

	g->vza_prime = gortt_prime_theta( p, g->vza );
	g->sza_prime = gortt_prime_theta( p, g->sza );


	/*Sunlit & viewed ground*/
	Kg = gortt_kg( p, g, g->raa ) ;


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

		/*
		everything the radiative transfer half of the model
		contributes, for this waveband and this geometry
		*/

		gortt_rt_init( p, g, s, i, &rt );


		/*sunlit ground
		 Could include brdf model here*.

		 *However - note that gort_energy() assumes a
		 Lambertian background and so this would have
		 to change to accomadate this correctly.
		*/

		/*
		coupling: the background does not just reflect what reaches it,
		it reflects it back and forth against the canopy above. The
		series sums to 1/(1 - rho_s rho'_ff), which is exactly the factor
		eqs (36)-(37) already give the crown-background term -- the
		background components are the only ones without it.

		Ni et al. (1999) section III-B says so: the background is "simply
		assumed to be Lambertian", which is "not a problem for denser
		vegetation", with anything better left to future studies. That
		is what -rt_bg_coupling turns on. It matters most where the
		background dominates, i.e. exactly the sparse canopies the paper
		says it has not validated.
		*/
		bg = gortt_bg_coupling( p, &rt, *( s->rsoil + i ) ) ;

		Gd = *( s->rsoil + i ) * bg ;
		Gf = *( s->rsoil + i ) * bg ;

		G = g->fd*Gd + (1-g->fd)*Gf ;



		/*shaded ground*/

		/*Zd = ( rt.t_prime_df + g->p_neq0_heq0_sza ) * *( s->rsoil + i ) ;*/
		Zd = ( rt.t_prime_df + g->p_ngt0_heq0_sza ) * *( s->rsoil + i ) * bg ;
		Zf = ( rt.t_prime_ff - p->k_openep[ 0 ] ) * *( s->rsoil + i ) * bg ;

		Z = g->fd*Zd + (1-g->fd)*Zf ;

		/*
		the same mixture without the soil reflectance in it. Z is
		linear in rsoil, so this is Z/rsoil -- but it stays finite for
		a black soil, where that division is 0/0. gortt_energy uses it.
		*/
		*( s->ztrans + i ) = ( g->fd * ( rt.t_prime_df + g->p_ngt0_heq0_sza )
		                     + (1-g->fd) * ( rt.t_prime_ff - p->k_openep[ 0 ] ) ) * bg ;


		/*Sunlit Crown*/

		/*CdC = rt.p_prime_df + ( ( 1.0-rt.w ) \
			  * rt.kuusk * rt.w \
				* ( 1.0 - rt.g_asym ) ) \
				/ ( 2.0 * cos( g->sza ) * cos( g->vza ) ) ;*/

		CdC = rt.p_prime_df + ( ( 1.0-rt.w ) \
			  * rt.kuusk * rt.w \
				* ( 1.0 - rt.g_asym ) ) \
				/ ( 2.0 * cos( g->sza_prime ) * cos( g->vza_prime ) ) ;




		CfC = rt.p_prime_ff ;

		CdG = ( Z * Kprime_z + G * Kprime_g ) * p->k_openep[ 0 ] ;

		CfG = ( ( p->k_openep[ 0 ]+p->k_open[ 0 ] )*G + \
					( 1-(p->k_openep[ 0 ]+p->k_open[ 0 ] ) )*Z ) * p->k_openep[ 0 ] ;

		CdCG = ( rt.t_prime_df + rt.t_prime_0 ) \
				 * ( *( s->rsoil + i ) / ( 1.0 - *( s->rsoil + i ) * rt.p_prime_ff ) ) \
				 * ( rt.t_prime_ff - p->k_open[ 0 ] );

		CfCG = rt.t_prime_ff \
				 * ( *( s->rsoil + i ) / ( 1.0 - *( s->rsoil + i ) * rt.p_prime_ff ) ) \
				 * ( rt.t_prime_ff - p->k_open[ 0 ] );


		Cd = CdC + CdG + CdCG ;
		Cf = CfC + CfG + CfCG ;

		C = g->fd*Cd + (1-g->fd)*Cf ;



		/*
		Shaded crown.

		Td and Tf used to be written out again in full here; they were
		term for term identical to CdCG and CfCG, which the comment
		alongside them already noted. The shaded crown signature is
		exactly the crown-ground interaction without the crown itself.
		*/

		T = g->fd*CdCG + (1-g->fd)*CfCG ;


		/*
		Neighbour crowns. Ni et al. (1999), just under eq. (39): "Notice
		that the component terms C and T are also affected by the
		neighboring crowns. This factor is usually small and very
		complicated to model, so we ignore it here."

		-rt_crown_crown puts a first-order version back. A crown radiates
		rho'_ff of what falls on it; a fraction (1 - K_open) of that meets
		another crown rather than the sky or the ground, and the exchange
		repeats, so the geometric series adds

			rho'_ff (1-K_open) / ( 1 - rho'_ff (1-K_open) )

		of the crown's own brightness. It goes on both crown components:
		a shaded crown is lit by its neighbours as much as a sunlit one
		is, which is why T -- black in the near infrared today -- is
		where it shows up most. This is the most speculative of the -rt_*
		options and the one to drop first if it does not earn its place.
		*/
		if( p->rt_crown_crown == TRUE ){
			cc = rt.p_prime_ff * ( 1.0 - ( p->k_open[ 0 ] + p->k_openep[ 0 ] ) ) ;
			cc = ( cc < 1.0 ) ? cc / ( 1.0 - cc ) : 0.0 ;
			T += cc * rt.p_prime_ff ;
			C += cc * rt.p_prime_ff ;
		}


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


double gortt_prime_theta( const gortt_parameters *p, double za )
{

	return( atan( ( p->b/p->r ) * tan( za ) ) );

}



double gortt_leaf_angle_distribution( const gortt_parameters *p, double za )
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


void gortt_init_params( gortt_parameters *p, gortt_geometry *g )
{

	int i ;


	p->ellipticity = p->b / p->r ;
	p->rr = p->r * p->r ;
	p->rrr = p->rr * p->r;

	/*
	Heights for the overlap average.  Li & Strahler (1992) eqs 4-6 use the
	height of the crown centre, and eq. 3 the overlap averaged over the
	distribution of centre heights, which gortt takes as uniform in [h1,h2].
	(This used to be a single h = 2b + h2 - h1, the canopy depth, which
	narrowed the hotspot by about half for the SOBS stand.)  The rule is set
	up on [-1,1] and mapped, so that h1 == h2 gives every point the same
	height rather than zero weights.
	*/
	gauleg( -1.0, 1.0, p->h_overlap, p->w_overlap, NH_OVERLAP ) ;
	for( i=0; i<NH_OVERLAP; i++ ){
		p->h_overlap[ i ] = 0.5*( p->h1 + p->h2 ) + 0.5*( p->h2 - p->h1 )*p->h_overlap[ i ] ;
		p->w_overlap[ i ] = 0.5*p->w_overlap[ i ] ;
	}



	/*
	* figure out tau.
	*/

	/*
	n.b. this is evaluated once, at whatever sza g happens to hold, which
	is the default of zero: gortt_cl_parser does not set the geometry, the
	per-angle loop does, and that runs later. Harmless while the leaf angle
	distribution is hardwired to 0.5, but p->k would have to move into the
	per-geometry path along with g->k_vza before any real LAD is added.
	*/
	p->k = gortt_leaf_angle_distribution( p, g->sza ) ;
	/*p->elai = p->favd * ( ( 4.0/3.0 ) * p->lambda * M_PI * p->ellipticity * p->rrr ) ; */
	p->elai = p->favd * ( ( 1.333333 ) * p->lambda * M_PI * p->ellipticity * p->rrr ) ;
	p->tau  =  p->k * p->favd ;

	/*
	The crown's own optical depth, for -rt_crown_path. 4b/3 is the mean
	vertical chord of a spheroid of semi-axes (r,r,b): the chord at impact
	parameter s is 2b sqrt(1-(s/r)^2) and the area weight over the projected
	disc is 2s ds/r^2, which integrates to 4b/3. For a sphere that is the
	4r/3 gortt_get_es returns at h=0, so the two agree by construction.

	Equivalently this is elai/PCC with PCC = lambda pi r^2, i.e. the same
	foliage spread over the crowns' footprint instead of over everything.
	squeeze_lai is filled later, in gortt_compute_probabilities, because it
	needs K_open.
	*/
	p->crown_lai   = p->favd * ( 4.0 / 3.0 ) * p->b ;
	p->squeeze_lai = p->elai ;


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
	p->xq_gl = NULL ;
	p->wq_gl = NULL ;
	p->k_open = NULL ;
	p->k_openep = NULL ;
	p->factorial = NULL ;
	p->zq = NULL ;
	p->pn0q = NULL ;
	p->p_n0 = NULL ;
	p->epgap = NULL ;

	/*
	The rule gortt_crown_proj_volume and gortt_get_es use for the one
	part of their integrands that is not integrable in closed form. It is
	set up on [-1,1] once here rather than per call, because between them
	they are evaluated tens of thousands of times a run and gauleg finds
	its nodes by Newton iteration.
	*/
	p->nq_gl = NQ_GL * p->nquad ;
	p->xq_gl = gortt_alloc_1d( p->nq_gl ) ;
	p->wq_gl = gortt_alloc_1d( p->nq_gl ) ;
	gauleg( -1.0, 1.0, p->xq_gl, p->wq_gl, p->nq_gl ) ;

	/*
	The height grid that gortt_calc_epgap_ht partitions the first-contact
	probability over, and P(n=0) on it. Uniform over the full canopy
	[z1_p, z2_p] rather than relative to the height being asked about, so
	that the table can be filled once for every h, and so that z1_p is a
	node: the ground layer then needs no partial panel and the BRF comes
	out exactly independent of nlayers.
	*/
	p->nq_ct = NQ_CONTACT * p->nquad ;
	p->zq    = gortt_alloc_1d( p->nq_ct + 1 ) ;
	p->pn0q  = gortt_alloc_2d( p->nq_ct + 1, p->nth ) ;
	for( i=0; i<=p->nq_ct; i++ )
		p->zq[ i ] = p->z1_p
		           + ( p->z2_p - p->z1_p ) * (double) i / (double) p->nq_ct ;


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

	 p->k_open   = gortt_alloc_1d( p->nlayers );
	 p->k_openep = gortt_alloc_1d( p->nlayers );
	 p->p_n0     = gortt_alloc_2d( p->nlayers , p->nth );
	 p->epgap    = gortt_alloc_2d( p->nlayers , p->nth );


	return;
}



/*
linearly interoplate the gap probabilities to get the
probs. at the required zenith angle.
Requires that p->p_n0 and p->epgap have been populated.
*/
void gortt_set_zenith_dependant_probabilities( const gortt_parameters *p, gortt_geometry *g )
{

	int cindex ;
	int findex ;
	double d, pos ;


	/*solar zenith*/

	pos = fabs( g->sza ) / p->dth ;

	cindex = ceil( pos );
	findex = floor( pos );

	d = pos - findex ;

	g->p_neq0_heq0_sza = d * p->p_n0[ 0 ][ cindex ] + ( 1.0 - d ) * p->p_n0[ 0 ][ findex ] ;
	g->p_ngt0_heq0_sza = d * p->epgap[ 0 ][ cindex ] + ( 1.0 - d ) * p->epgap[ 0 ][ findex ] ;


	/*view zenith*/

	pos = fabs( g->vza ) / p->dth ;

	cindex = ceil( pos );
	findex = floor( pos );

	d = pos - findex ;

	g->p_neq0_heq0_vza = d * p->p_n0[ 0 ][ cindex ] + ( 1.0 - d ) * p->p_n0[ 0 ][ findex ] ;
	g->p_ngt0_heq0_vza = d * p->epgap[ 0 ][ cindex ] + ( 1.0 - d ) * p->epgap[ 0 ][ findex ] ;


	return ;

}




double * gortt_alloc_1d( int n )
{

	double *r;

	if (! ( r = (double *) calloc( n, sizeof(double) ) ) ) {
		fprintf(stderr, "Memory allocation failed (gortt_alloc_1d)\n");
		exit( EXIT_FAILURE );
	}

	return ( r );

}

void gortt_dealloc_1d( double *d )
{

	free( d );

	return;

}

double ** gortt_alloc_2d( int n, int m )
{


	int q;
	double **r;


	/*n.b. sizeof(double *), not sizeof(double) - this is a row of pointers*/
	if (! ( r = (double **) calloc( (unsigned) n, sizeof(double *) ) ) ) {
		fprintf(stderr, "Memory allocation failed (gortt_alloc_2d)\n");
		exit( EXIT_FAILURE );
	}

	for( q=0; q<n; q++)
		r[q] = gortt_alloc_1d( m );


	return( r );

}

void gortt_dealloc_2d( double **d, int n )
{

	int q ;

	for( q=0; q<n; q++) gortt_dealloc_1d( d[ q ] );
	free( d );

	return;

}

void gortt_dealloc_3d( double ***d, int n, int m )
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



static char *gortt_next_arg( int argc, char **argv, int *i )
/*
Return the value following the option at argv[*i], advancing *i onto it.

Without the bounds check, an option that takes a value and is given as the
last word on the command line ("gortt -Cab") reads argv[argc], which is
NULL, and atof/strcpy then segfault.
*/
{
	if( *i + 1 >= argc ){
		fprintf( stderr, "%s: option %s requires a value\n", argv[ 0 ], argv[ *i ] );
		exit( EXIT_FAILURE );
	}

	return( argv[ ++(*i) ] );
}


void gortt_cl_parser( int argc, char **argv, gortt_parameters *p, gortt_spectra *s, gortt_control *c,
                      gortt_lidar_params *l )
{


	int i, use_true_p=FALSE, use_lai=FALSE;
	int got_refl_leaf=FALSE, got_trans_leaf=FALSE;
	int got_lidar_opt=FALSE, got_rt_path=FALSE;
	float hb=2.0, br=1.0, pcc=0.5, lai=2.0;
	char *arg;

	p->read_prob_file=FALSE ;
	p->write_prob_file=FALSE;
	p->use_user_beta=FALSE ;



	for( i=1; i<argc; i++ ){

		if( *argv[ i ] == '-' ){

			/*n.b. -nlayers and -nquad must stay at the TOP of this chain:
			-N below is matched case-insensitively on two characters, and
			would otherwise swallow both of them as the PROSPECT structure
			parameter.*/
			/**/ if( !strncmp( argv[ i ], "-nlayers", 8 ) )    p->nlayers  = atoi( gortt_next_arg( argc, argv, &i ) );
			else if( !strncmp( argv[ i ], "-nquad", 6 ) )      p->nquad    = atoi( gortt_next_arg( argc, argv, &i ) );

			else if( !strncasecmp( argv[ i ], "-favd", 5 ) )   p->favd     = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-h1", 3 ) )     p->h1       = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-h2", 3 ) )     p->h2       = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-lambda", 7 ) ) p->lambda   = atof( gortt_next_arg( argc, argv, &i ) );


			else if( !strncmp( argv[ i ], "-HB", 3 ) )   { use_true_p=TRUE ;  hb = atof( gortt_next_arg( argc, argv, &i ) ); }
			else if( !strncmp( argv[ i ], "-BR", 3 ) )   { use_true_p=TRUE ;  br = atof( gortt_next_arg( argc, argv, &i ) ); }
			else if( !strncmp( argv[ i ], "-PCC", 7 ) )  { use_true_p=TRUE ;  pcc = atof( gortt_next_arg( argc, argv, &i ) ); }

			else if( !strncmp( argv[ i ], "-LAI", 7 ) )  { use_lai=TRUE ;  lai = atof( gortt_next_arg( argc, argv, &i ) ); }


			else if( !strncasecmp( argv[ i ], "-beta", 5 ) ){ p->use_user_beta=TRUE; p->beta=atof( gortt_next_arg( argc, argv, &i ) );}

			else if( !strncasecmp( argv[ i ], "-diffuse", 5 ) ){ p->use_user_fd=TRUE; p->user_fd=1.0-atof( gortt_next_arg( argc, argv, &i ) );}

			/*n.b. both of these must stay above the "-r" test at the
			bottom, which would otherwise swallow -refl_leaf.*/
			else if( !strncmp( argv[ i ], "-refl_leaf", 10 ) ){
				s->is_user_leaf = TRUE ;
				got_refl_leaf = TRUE ;
				s->user_refl_leaf = atof( gortt_next_arg( argc, argv, &i ) ) ;

			}

			else if( !strncmp( argv[ i ], "-trans_leaf", 11 ) ){
				s->is_user_leaf = TRUE ;
				got_trans_leaf = TRUE ;
				s->user_trans_leaf = atof( gortt_next_arg( argc, argv, &i ) ) ;

			}

			else if( !strncmp( argv[ i ], "-alb_soil", 9 ) ){
				s->is_user_soil = TRUE ;
				s->is_file_soil = FALSE ;
				s->user_r_soil = atof( gortt_next_arg( argc, argv, &i ) ) ;

			}

			else if( !strncmp( argv[ i ], "-soil_spectra", 10 ) ){
				s->is_user_soil = FALSE ;
				s->is_file_soil = TRUE ;
				strcpy( s->soil_file, gortt_next_arg( argc, argv, &i ) );

			}


			else if( !strncmp( argv[ i ], "-prnspec", 7 ) ) c->print_comp_spectra=TRUE;
			else if( !strncmp( argv[ i ], "-prnprop", 7 ) ) c->print_comp_proport=TRUE;
			else if( !strncmp( argv[ i ], "-energy", 7 ) ) c->calc_integrals=TRUE;
			else if( !strncmp( argv[ i ], "-q08_pn_kopen", 7 ) ) c->use_q08_pn_kopen=TRUE;

			/*
			The experimental scattering formulations. n.b. all of these
			must stay above the "-r" test at the bottom, which is matched
			on two characters and would otherwise swallow every one of
			them AND eat the following argv word as a crown radius --
			the same trap -refl_leaf is commented for above. They are
			matched on strlen+1, i.e. including the terminator, so that
			-rt_crown_pathologically is an error rather than a synonym.
			*/
			else if( !strncmp( argv[ i ], "-rt_crown_crown", 16 ) ) p->rt_crown_crown=TRUE;
			else if( !strncmp( argv[ i ], "-rt_crown_path", 15 ) )  p->rt_crown_path=TRUE;
			else if( !strncmp( argv[ i ], "-rt_crown_brdf", 15 ) )  p->rt_crown_brdf=TRUE;
			else if( !strncmp( argv[ i ], "-rt_recollision", 16 ) ) p->rt_recollision=TRUE;
			else if( !strncmp( argv[ i ], "-rt_upscatter", 14 ) )   p->rt_upscatter=TRUE;
			else if( !strncmp( argv[ i ], "-rt_bg_coupling", 16 ) ) p->rt_bg_coupling=TRUE;
			else if( !strncmp( argv[ i ], "-rt_fd_paper", 13 ) )    p->rt_fd_paper=TRUE;

			else if( !strncmp( argv[ i ], "-rt_path", 9 ) ){

				got_rt_path = TRUE ;
				arg = gortt_next_arg( argc, argv, &i ) ;

				if( !strcasecmp( arg, "scene" ) )        p->rt_crown_path_mode = RT_PATH_SCENE ;
				else if( !strcasecmp( arg, "squeeze" ) ) p->rt_crown_path_mode = RT_PATH_SQUEEZE ;
				else if( !strcasecmp( arg, "crown" ) )   p->rt_crown_path_mode = RT_PATH_CROWN ;
				else{
					fprintf( stderr, "%s: -rt_path must be scene, squeeze or crown (got %s)\n",
					         argv[ 0 ], arg );
					exit( EXIT_FAILURE );
				}

			}

			else if( !strncmp( argv[ i ], "-rt_all", 8 ) ){
				p->rt_crown_path=TRUE;
				p->rt_recollision=TRUE;
				p->rt_bg_coupling=TRUE;
				p->rt_crown_crown=TRUE;
				p->rt_fd_paper=TRUE;
				p->rt_crown_brdf=TRUE;
				p->rt_upscatter=TRUE;
			}

			/*the lidar options. n.b. every one of these must stay ABOVE
			the "-lidar" test below, which is matched on six characters
			and would otherwise swallow all of them. None of them is a
			prefix of another: -lidar_wl and -lidar_waveform diverge at
			index 8, so their relative order does not matter.
			n.b. p->nh_probs is NOT raised here -- see the post-loop
			block, which has to wait for -nlayers.*/
			else if( !strncmp( argv[ i ], "-lidar_waveform", 15 ) ){ got_lidar_opt=TRUE; l->binned_output=TRUE; }
			else if( !strncmp( argv[ i ], "-lidar_wl", 9 ) )      { got_lidar_opt=TRUE; l->wl      = atof( gortt_next_arg( argc, argv, &i ) ); }
			else if( !strncmp( argv[ i ], "-lidar_J0", 9 ) )      { got_lidar_opt=TRUE; l->J0      = atof( gortt_next_arg( argc, argv, &i ) ); }
			else if( !strncmp( argv[ i ], "-lidar_pulse", 12 ) )  { got_lidar_opt=TRUE; l->pulse   = atof( gortt_next_arg( argc, argv, &i ) ); }
			else if( !strncmp( argv[ i ], "-lidar_binsize", 14 ) ){ got_lidar_opt=TRUE; l->binsize = atof( gortt_next_arg( argc, argv, &i ) ); }

			else if( !strncmp( argv[ i ], "-lidar", 6 ) ) c->lidar_mode=TRUE;

			else if( !strncmp( argv[ i ], "-P", 2 ) ){ p->read_prob_file=TRUE; strcpy( p->prob_fn, gortt_next_arg( argc, argv, &i ) );}
			else if( !strncmp( argv[ i ], "-W", 2 ) ) p->write_prob_file=TRUE;

			/*prospect*/
			else if( !strncasecmp( argv[ i ], "-N", 2 ) )   s->p_N     = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-cab", 4 ) ) s->p_Cab     = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-car", 4 ) ) s->p_Car     = atof( gortt_next_arg( argc, argv, &i ) );
			/*PROSPECT-D calls this one Anth, without the C ("content") prefix
			the other pigments carry, and so does gortt_usage. -Canth is kept
			as an alias because it was the only spelling the parser accepted
			for a long time. n.b. -Canth is matched on "-ca" alone, so it has
			to stay below -cab and -car.*/
			else if( !strncasecmp( argv[ i ], "-anth", 5 ) )   s->p_Anth      = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-canth", 3 ) )  s->p_Anth      = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-cbrown", 3 ) )  s->p_Cbrown      = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-cw", 3 ) )  s->p_Cw      = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-cm", 3 ) )  s->p_Cm      = atof( gortt_next_arg( argc, argv, &i ) );

			/*price*/
			else if( !strncasecmp( argv[ i ], "-rsl1", 5 ) )  s->rsl1    = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-rsl2", 5 ) )  s->rsl2    = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-rsl3", 5 ) )  s->rsl3    = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-rsl4", 5 ) )  s->rsl4    = atof( gortt_next_arg( argc, argv, &i ) );

			/*leave these at the bottom! (so as not to catch other opts)*/
			else if( !strncasecmp( argv[ i ], "-b", 2 ) )      p->b        = atof( gortt_next_arg( argc, argv, &i ) );
			else if( !strncasecmp( argv[ i ], "-r", 2 ) )      p->r        = atof( gortt_next_arg( argc, argv, &i ) );






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

	/*
	Refuse a half-specified leaf. Defaulting the missing one to zero would
	silently change the single-scattering albedo -- which is what the whole
	scattering calculation is built on -- and the mistake is easy to make
	when converting a script from the old single-albedo option.
	*/
	if( got_refl_leaf != got_trans_leaf ){
		fprintf( stderr, "%s: -refl_leaf and -trans_leaf must be given together\n", argv[ 0 ] );
		exit( EXIT_FAILURE );
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

	/*
	A -lidar_* option on its own does nothing, and "nothing" here means the
	run falls through to the BRF path and blocks on stdin with no output and
	no message -- which looks exactly like a hang. None of these options has
	any meaning outside lidar mode, so refuse instead.

	n.b. this is checked after the loop, so the order on the command line
	does not matter: -lidar_wl 1550 -lidar is as good as -lidar -lidar_wl 1550.
	*/
	if( got_lidar_opt == TRUE && c->lidar_mode == FALSE ){
		fprintf( stderr, "%s: the -lidar_* options require -lidar\n", argv[ 0 ] );
		exit( EXIT_FAILURE );
	}

	/*
	-rt_path only says which optical depth -rt_crown_path should use, so on
	its own it does nothing at all. Refuse it rather than let it be inert:
	the run would otherwise produce the published numbers while looking
	like it had been asked for something else.
	*/
	if( got_rt_path == TRUE && p->rt_crown_path == FALSE ){
		fprintf( stderr, "%s: -rt_path requires -rt_crown_path\n", argv[ 0 ] );
		exit( EXIT_FAILURE );
	}

	/*
	The recollision probability describes a crown, so it can only rescale
	a crown: applied to the stand-sized layer of eqs (15)-(16) it divides
	a crown quantity by a scene one and the result means nothing. Refuse
	the combination rather than let it look like it worked.
	*/
	if( p->rt_recollision == TRUE && p->rt_crown_path == FALSE ){
		fprintf( stderr, "%s: -rt_recollision requires -rt_crown_path\n", argv[ 0 ] );
		exit( EXIT_FAILURE );
	}

	/*
	Same argument for the crown BRDF: it replaces rho_df on a chord
	through a crown, and without -rt_crown_path there is no chord -- the
	layer is the whole stand and there is nothing for a phase function to
	be a property of.
	*/
	if( p->rt_crown_brdf == TRUE && p->rt_crown_path == FALSE ){
		fprintf( stderr, "%s: -rt_crown_brdf requires -rt_crown_path\n", argv[ 0 ] );
		exit( EXIT_FAILURE );
	}

	/*
	gortt_init_params divides by nlayers-1 to get dz, and the probability
	tables are indexed h and h+1.
	*/
	if( p->nlayers < 2 ){
		fprintf( stderr, "%s: -nlayers must be at least 2 (got %d)\n", argv[ 0 ], p->nlayers );
		exit( EXIT_FAILURE );
	}

	/*
	-nquad scales the internal quadrature orders, which are otherwise
	independent of nlayers. 1 is already converged; larger values are for
	demonstrating that.
	*/
	if( p->nquad < 1 ){
		fprintf( stderr, "%s: -nquad must be at least 1 (got %d)\n", argv[ 0 ], p->nquad );
		exit( EXIT_FAILURE );
	}

	/*
	A waveform is a function of height, so lidar mode needs the gap
	probabilities at every layer rather than only at the ground.

	n.b. this has to happen HERE rather than in the -lidar branch of the
	loop above, because -nlayers may be given after -lidar on the command
	line. Doing it in the branch would leave nh_probs at whatever nlayers
	happened to be at the time, and gortt_lidar would silently read the
	calloc zeros above that row.
	*/
	if( c->lidar_mode == TRUE ) p->nh_probs = p->nlayers ;

	/*
	Lidar mode is incompatible with three things, all for the same reason:
	they only ever deal with the gap probabilities at the ground.

	  -P / -W   gortt_write_prob_file and gortt_read_prob_file serialise
	            row 0 alone, so a -P run has zeros above the ground. -W
	            also short-circuits main(), as lidar mode does.
	  -q08      gortt_gap_probabilities_Q08 is only coded for h=0 and
	            says so itself.

	Silently returning zeros for the whole canopy is a much worse outcome
	than refusing, so refuse.
	*/
	if( c->lidar_mode == TRUE ){

		if( p->read_prob_file == TRUE || p->write_prob_file == TRUE ){
			fprintf( stderr, "%s: -lidar cannot be combined with -P or -W: the gap probability\n", argv[ 0 ] );
			fprintf( stderr, "file format holds the ground layer only, and a waveform needs every layer\n" );
			exit( EXIT_FAILURE );
		}

		if( c->use_q08_pn_kopen == TRUE ){
			fprintf( stderr, "%s: -lidar cannot be combined with -q08_pn_kopen: the Q08\n", argv[ 0 ] );
			fprintf( stderr, "approximation is only valid at the ground (h=0)\n" );
			exit( EXIT_FAILURE );
		}

		/*PROSPECT-D and the Price soil EOFs both run 400-2500nm*/
		if( l->wl < PROSPECT_LOWER_WL ||
		    l->wl > PROSPECT_LOWER_WL + ( PROSPECT_NWBANDS - 1 ) * PROSPECT_SPECTRAL_RESLN ){
			fprintf( stderr, "%s: -lidar_wl must be between %g and %gnm (got %g)\n", argv[ 0 ],
			         PROSPECT_LOWER_WL,
			         PROSPECT_LOWER_WL + ( PROSPECT_NWBANDS - 1 ) * PROSPECT_SPECTRAL_RESLN, l->wl );
			exit( EXIT_FAILURE );
		}

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
	fprintf( stderr, "-favd arg\tset the foliage area volume density (1/m) within crown to arg\n"  );
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
	fprintf( stderr, "-refl_leaf    arg \tset leaf reflectance to arg (turns prospect off)*\n"  );
	fprintf( stderr, "-trans_leaf   arg \tset leaf transmittance to arg (turns prospect off)*\n"  );
	fprintf( stderr, "                  \t[* n.b. these two must be given together]\n"  );
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

	fprintf( stderr, "\n============ Numerical resolution:\n" );
	fprintf( stderr, "-nlayers arg\tset the number of canopy height layers to arg (default 15, minimum 2)\n"  );
	fprintf( stderr, "            \tthis is the REPORTING grid: it sets how finely the lidar profile is written\n"  );
	fprintf( stderr, "            \tout, and nothing else. No quadrature steps with it, so the BRF does not\n"  );
	fprintf( stderr, "            \tmove when it changes. Only lidar mode gets slower as it rises.\n"  );
	fprintf( stderr, "-nquad arg  \tmultiply the internal quadrature orders by arg (default 1, minimum 1)\n"  );
	fprintf( stderr, "            \tthe default is converged; raise it to check that for a given canopy\n"  );

	fprintf( stderr, "\n============ Lidar options:\n" );
	fprintf( stderr, "-lidar           \tpredict the lidar waveform instead of the BRF, and exit\n"  );
	fprintf( stderr, "                 \tn.b. lidar mode does not read stdin at all: the device is assumed to be\n"  );
	fprintf( stderr, "                 \tat nadir, so no view/illumination geometry is required. Normal output is\n"  );
	fprintf( stderr, "                 \tsuppressed. Cannot be combined with -W, -P or -q08_pn_kopen, all of which\n"  );
	fprintf( stderr, "                 \tonly deal with the gap probabilities at the ground.\n"  );
	fprintf( stderr, "-lidar_wl arg    \tset the lidar wavelength (nm) to arg (default 1064, must be 400-2500)\n"  );
	fprintf( stderr, "-lidar_J0 arg    \tset the emitted pulse energy to arg (default 1, a pure scale factor)\n"  );
	fprintf( stderr, "-lidar_waveform  \temit the binned waveform rather than the per-layer profile\n"  );
	fprintf( stderr, "-lidar_pulse arg \tset the Gaussian pulse width (sigma, m) to arg (default 0, no convolution)\n"  );
	fprintf( stderr, "-lidar_binsize arg\tset the range bin size (m) to arg (default: the canopy layer thickness)\n"  );
	fprintf( stderr, "                 \tn.b. the last two apply to -lidar_waveform only\n"  );
	fprintf( stderr, "                 \tn.b. all of the -lidar_* options require -lidar, and are refused without it\n"  );

	fprintf( stderr, "\n============ Experimental scattering options:\n" );

  fprintf( stderr, "  Ni et al. (1999) section V says the model is validated over dense forests, and\n"  );
  fprintf( stderr, "  defers three things that matter as a canopy opens up; most of these switch one\n"  );
  fprintf( stderr, "  of them on. -rt_crown_brdf and -rt_upscatter do not: the first comes from\n"  );
  fprintf( stderr, "  measuring a single crown against a Monte Carlo, the second from the two-stream\n"  );
  fprintf( stderr, "  closure underneath eqs (8)-(11). All are off by default: with none of them\n"  );
  fprintf( stderr, "  given the model is the paper as published.\n"  );

	fprintf( stderr, "-rt_crown_path \ttake the crown path reflectances from a crown rather than from a layer\n"  );
	fprintf( stderr, "               \tholding the whole scene's foliage. Eqs (15)-(16) use the latter, so a\n"  );
	fprintf( stderr, "               \tstand of widely spaced crowns comes out with black crowns.\n"  );
	fprintf( stderr, "-rt_path arg   \twhich layer -rt_crown_path means: scene (as published, a control),\n"  );
	fprintf( stderr, "               \tsqueeze (elai/(1-Kopen)) or crown (the default, chord-averaged)\n"  );
	fprintf( stderr, "-rt_recollision\trescale the crown path reflectances by the crown's recollision\n"  );
	fprintf( stderr, "               \tprobability, so that they carry every scattering order. Requires\n"  );
	fprintf( stderr, "               \t-rt_crown_path: it describes a crown, so it can only rescale a crown.\n"  );
	fprintf( stderr, "-rt_crown_brdf \tgive rho'_df a direction: single scattering exactly, multiple\n"  );
	fprintf( stderr, "               \tscattering isotropically, normalised so the hemispherical integral\n"  );
	fprintf( stderr, "               \tis still rho_df. Requires -rt_crown_path.\n"  );
	fprintf( stderr, "-rt_bg_coupling\tlet the sunlit and shaded background exchange radiation with the canopy,\n"  );
	fprintf( stderr, "               \tthe 1/(1-rho_s rho'_ff) that eqs (36)-(37) already give the crown\n"  );
	fprintf( stderr, "-rt_crown_crown\tadd the neighbour-crown term the paper drops just below eq. (39)\n"  );
	fprintf( stderr, "-rt_upscatter  \ttwo-stream upscatter beta = 1/2 + (rl-tl)/(8w) instead of the 1/2\n"  );
	fprintf( stderr, "               \teqs (8)-(11) assume. Inert where the leaf reflects and transmits\n"  );
	fprintf( stderr, "               \tequally, so it is a red-band option. n.b. NOT -beta, which is the\n"  );
	fprintf( stderr, "               \tgeometric-optical mutual shadowing weight.\n"  );
	fprintf( stderr, "-rt_fd_paper   \tdirect beam fraction mu/(mu+0.05), which is what section IV gives;\n"  );
	fprintf( stderr, "               \twithout it gortt uses mu/(mu+0.09). n.b. -diffuse overrides both.\n"  );
	fprintf( stderr, "-rt_all        \tall seven of the above\n"  );

	fprintf( stderr, "\n============ Input/output options:\n" );


	fprintf( stderr, "-prnspec\tprint the scene component spectra for each wavelength inside {}\n" );
	fprintf( stderr, "-prnprop\tprint the viewed proportions of scene components inside []\n" );
	fprintf( stderr, "-energy \tprint the spectral albedo, absoption by veg and absorption by soil for each wavelength after other outputs.\n" );
	fprintf( stderr, "-u      \tprint this message and exit\n"  );

  fprintf( stderr, "\n" );

	return ;

}


char get_first_string_token( char *line, char *token )
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


void gortt_price_soil( gortt_spectra *s, double *sv1, double *sv2, double *sv3, double *sv4 )
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

		}else if( s->is_file_soil ){

			/*
			-soil_spectra: gortt_read_soil_lut has already resampled
			the file onto a 1nm grid, so this is the same linear
			interpolation as the PROSPECT one, on that grid instead
			of the Price 5nm one.
			*/

			upper = 1 + ( *( s->wavelength +i ) - SOIL_LUT_LOWER_WL ) / SOIL_LUT_SPECTRAL_RESLN ;
			lower = ( *( s->wavelength +i ) - SOIL_LUT_LOWER_WL ) / SOIL_LUT_SPECTRAL_RESLN ;

			if( upper > SOIL_LUT_NBANDS - 1 ) upper = SOIL_LUT_NBANDS - 1 ;

			fraction = ( *( s->wavelength +i ) - SOIL_LUT_LOWER_WL ) / SOIL_LUT_SPECTRAL_RESLN - lower ;

			*(s->rsoil + i) = s->soil_spectra[ lower ] * ( 1 - fraction )
			                + s->soil_spectra[ upper ] * fraction ;

		}else{

			upper = 1. + (*( s->wavelength +i )-400) / 5.0 ;
			lower = (*( s->wavelength +i )-400) / 5.0 ;

			/*
			At exactly 2500nm upper lands one past the end of the
			Price vectors. fraction is zero there, so the out of
			range element was always multiplied by zero and never
			affected the answer - but it was still read.
			*/
			if( upper > PRICE_NBANDS - 1 ) upper = PRICE_NBANDS - 1 ;

			fraction = (double) ( *( s->wavelength +i ) - 400. ) / 5.0 - lower ;


			rs_lower = s->rsl1 * *( sv1 + lower ) + s->rsl2 * *( sv2 + lower ) +  s->rsl3 * *( sv3 + lower ) + s->rsl4 * *( sv4 + lower ) ;
			rs_upper = s->rsl1 * *( sv1 + upper ) + s->rsl2 * *( sv2 + upper ) +  s->rsl3 * *( sv3 + upper ) + s->rsl4 * *( sv4 + upper ) ;

			*(s->rsoil + i) = rs_lower * ( 1 - fraction ) + rs_upper * fraction ;

		}

	}

	return ;

}


void gortt_prospect_interface( gortt_spectra *s )
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
			*(s->rleaf + i) = s->user_refl_leaf ;
			*(s->tleaf + i) = s->user_trans_leaf ;

		}else{
			/*get the leaf optics for the desired wavelength*/
			upper = 1 + (*( s->wavelength +i )-PROSPECT_LOWER_WL) / PROSPECT_SPECTRAL_RESLN ;
			lower = (*( s->wavelength +i )-PROSPECT_LOWER_WL) / PROSPECT_SPECTRAL_RESLN ;

			/*
			As in gortt_price_soil: at exactly 2500nm upper lands one
			past the end of the reflectance half of prospect_out, and
			upper+PROSPECT_NWBANDS one past the end of the array
			itself. fraction is zero there so the answer was right,
			but the read was not.
			*/
			if( upper > PROSPECT_NWBANDS - 1 ) upper = PROSPECT_NWBANDS - 1 ;

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
void gortt_read_soil_lut( gortt_spectra *s )
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

				if( (index>=0) && (index < SOIL_LUT_NBANDS) )
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




	fclose( fp );

	return ;


}


