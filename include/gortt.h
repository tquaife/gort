#ifndef GORTT
#define GORTT

#define PD_S_BUFF 3

#define DTOR(x)  ((x)*M_PI/180.0)
#define RTOD(x)  ((x)*180.0/M_PI)
#define SEC(x)   ((1.0/cos((x))))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#define LAD_PLANOPHILE 1
#define LAD_ERECTOPHILE 2
#define LAD_PLAGIOPHILE 3
#define LAD_EXTREMOPHILE 4
#define LAD_UNIFORM 5

#define LAD_05 100

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#define MAX_LINE_LEN 1000

/*prospect*/
#define PROSPECT_NWBANDS 2101
#define PROSPECT_LOWER_WL 400.0
#define PROSPECT_SPECTRAL_RESLN 1.0


/* new 09/10/12:
A structure to hold control variables,
e.g. to decide what ouput to print etc
*/

typedef struct  {

    unsigned short print_comp_spectra ;
    unsigned short print_comp_proport ;
    unsigned short calc_integrals ;
    unsigned short use_q08_pn_kopen ;
    unsigned short lidar_mode ;

} gortt_control ;



/*
One view/sun geometry, plus everything derived from it. Every field below
is rewritten by gortt_set_geometry / gortt_run_geometry for each angle in
the input, so nothing here survives from one output line to the next.
*/
typedef struct  {
		double vza ;
		double vaa ;
		double sza ;
		double saa ;

		double raa ;
		double vza_prime ;
		double sza_prime ;

		/*viewed areal proportions: sunlit crown, sunlit
		background, shaded crown, shaded background*/
		double Kc ;
		double Kg ;
		double Kt ;
		double Kz ;

		/*leaf projection G at the view zenith (the solar
		zenith equivalent is p->k, which does not vary)*/
		double k_vza ;

		/*fraction of the incoming radiation that is direct*/
		double fd ;

		/*gap probabilities interpolated to this geometry's
		zenith angles by gortt_set_zenith_dependant_probabilities*/
		double p_neq0_heq0_sza ;
		double p_ngt0_heq0_sza ;
		double p_neq0_heq0_vza ;
		double p_ngt0_heq0_vza ;

} gortt_geometry ;


typedef struct  {
		double *wavelength ;
		double *rsoil ;
		double *rleaf ;
		double *tleaf ;
		double *rsurf ;
		double *scomp ;

		/*intgral terms*/
		double *albedo ;
		double *favegt ;
		double *fasoil ;

		/*number of wavebands*/
		int nw ;
		/*counter!*/
		int n ;

		/*functions of the above*/
		double leaf_sscat_albedo ;
		double leaf_gamma ;

		/*leaf & soil properties*/

		double p_N ;
		double p_Cab;
		double p_Car;
		double p_Anth;
		double p_Cbrown;
		double p_Cw;
		double p_Cm;

		double rsl1;
		double rsl2;
		double rsl3;
		double rsl4;

		/*user specified rho leaf & soil*/

		int is_user_leaf ;
		int is_user_soil ;
		double user_r_leaf ;
		double user_r_soil ;

		int is_file_soil ;
		char soil_file[ MAX_LINE_LEN ] ;
		double soil_spectra[2101] ;

		} gortt_spectra ;


/*
Canopy and crown description, and everything precomputed from it. The
fields are grouped by lifetime, and the whole struct is read-only once
gortt_init_params and the gap-probability routines have run -- which is
why the BRF, albedo and lidar code can take it as const.

Names follow Ni et al. (1999): r and b are the crown radii, lambda the
crown number density, favd the foliage area volume density, and so on.
Anything suffixed _p lives in the ellipsoid->sphere transformed space
(see gortt_prime_theta); mixing primed and unprimed quantities in one
expression is almost always a bug.
*/
typedef struct  {

		/*
		---- 1. user inputs -------------------------------------------
		Written by gortt_set_defaults and then gortt_cl_parser.
		Read-only from gortt_init_params onwards.
		*/

		/*crown and canopy geometry*/
		double r ;			/*horizontal crown radius*/
		double b ;			/*vertical crown radius*/
		double h1 ;			/*lower bound of crown centres*/
		double h2 ;			/*upper bound of crown centres*/
		double lambda ;		/*crown number density*/
		double favd ;		/*foliage area volume density*/
		int lad ;			/*leaf angle distribution, LAD_*/

		/*numerical resolution*/
		double dth ;		/*zenith step of the probability tables*/
		int nlayers ;		/*number of height layers*/
		int maxcrowns ;		/*truncation of the Poisson crown count*/
		int nh_es ;
		int npoints ;		/*Gauss-Legendre points for the albedo*/

		/*diffuse fraction: if the user supplied one it is held here and
		copied into g->fd for every geometry, otherwise g->fd is derived
		from the solar zenith*/
		unsigned short use_user_fd ;
		double user_fd ;

		/*user override of the hotspot/overlap parameter*/
		int use_user_beta ;
		double beta ;

		/*gap probability file i/o (-W and -P)*/
		int read_prob_file ;
		int write_prob_file ;
		char prob_fn[ MAX_LINE_LEN ] ;


		/*
		---- 2. derived constants -------------------------------------
		Written only by gortt_init_params, from group 1 above.
		Read-only afterwards.
		*/

		double ellipticity ;	/*b/r: the vertical scaling*/
		double rr ;				/*r^2*/
		double rrr ;			/*r^3*/
		double h ;				/*total canopy depth*/
		double k ;				/*leaf projection G at the solar zenith*/
		double elai ;			/*effective LAI*/
		double tau ;
		double tau_p ;
		double z1 ;				/*lowest height with any foliage*/
		double z2 ;				/*highest height with any foliage*/
		double lv ;
		double favd_p ;
		double lv_p ;
		double z1_p ;
		double z2_p ;
		double h1_p ;
		double h2_p ;

		/*discretisation, all fixed by nlayers and dth*/
		double dz ;
		double dz_p ;
		double ds ;
		int nh1 ;
		int nh2 ;
		int nth ;				/*number of zenith angles, 0..90 deg*/

		double *height ;
		double *height_p ;
		double *theta ;
		double *theta_p ;
		double *factorial ;

		/*Gauss-Legendre abscissa and weights, allocated only when an
		integral is actually asked for*/
		double *abscissa ;
		double *weights ;


		/*
		---- 3. gap probabilities -------------------------------------
		Written once, either by gortt_gap_probabilities(_Q08) or by
		reading a -P file. The expensive part of a run, and the reason
		it is worth caching. Read-only afterwards.
		*/

		double **p_s0 ;
		double **p_n0 ;
		double **epgap ;
		double *k_open ;
		double *k_openep ;
		double ***pd_s ;

		} gortt_parameters ;



/*TQ's implimentation of the GORT paper*/
void   gortt_rsurf( gortt_parameters *, gortt_geometry *, gortt_spectra * );
void   gortt_gap_probabilities( gortt_parameters *, gortt_geometry * );
void   gortt_gap_probabilities_Q08( gortt_parameters *, gortt_geometry * );
void   gortt_set_zenith_dependant_probabilities( gortt_parameters *, gortt_geometry * );
double gortt_prime_theta( gortt_parameters *, double  );
double gortt_kg( gortt_parameters *, gortt_geometry *, double );
double gortt_kc( gortt_parameters *, gortt_geometry *, double );
void   gortt_kc_fFbeta( gortt_parameters *, gortt_geometry *, double, double, double *, double *, double * );
double gortt_kc_PP_only( gortt_parameters *, gortt_geometry *, double );
double gortt_kc_ambrals( gortt_parameters *, gortt_geometry * 	);
double gortt_overlap( gortt_parameters *, gortt_geometry *, double );
double gortt_t_prime_df( gortt_parameters *, gortt_geometry *, gortt_spectra * );
double gortt_t_prime_ff( gortt_parameters *, gortt_geometry *, gortt_spectra * );
double gortt_t_prime_0( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_t_df( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_t_ff( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_t_0( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_T_inf_df( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_T_inf_ff( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_R_inf_df( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_R_inf_ff( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_p_ff( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_p_prime_df( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_p_df( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_p_prime_ff( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_phase_function_assymetry( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_kuusk( gortt_parameters *, gortt_geometry *, gortt_spectra *  );
double gortt_leaf_angle_distribution( gortt_parameters *, double );

void  gortt_lidar( gortt_parameters *, gortt_geometry *, gortt_spectra * );


/*TQ's integral routines*/

void   gortt_energy( gortt_parameters *, gortt_geometry *, gortt_spectra * );
void   gortt_albedo( gortt_parameters *, gortt_geometry *, gortt_spectra * );
void   gauleg( double, double, double *, double *, int );

/*Wenge's functions for within crown probs*/
void   gortt_init_params( gortt_parameters *, gortt_geometry * );
int    gortt_s_to_index( gortt_parameters *, double );
double gortt_index_to_s( gortt_parameters *, int );

double gortt_crown_proj_volume( gortt_parameters *, double, double );
double gortt_crown_proj_cross_section( gortt_parameters *, double, double, double );
double gortt_weird_cross_section( gortt_parameters *, double, double, double );
double gortt_left_circle_area( double, double );
double gortt_right_ellipse_area( double, double, double );
void   gortt_get_pd_s( gortt_parameters *, int, int );

double gortt_get_es( gortt_parameters *, int, int );
double gortt_get_s( gortt_parameters *, int, double, int );
double gortt_get_pcc( gortt_parameters *, double  );
double gortt_vol( gortt_parameters *, int, int, int, double ) ;

void   gortt_calc_epgap( gortt_parameters * );
double gortt_calc_pgap( gortt_parameters *, double, double );
void   gortt_calc_kopen( gortt_parameters * );


/*memory alloc/dealloc*/
double *  gortt_alloc_1d( int ) ;
void      gortt_dealloc_1d( double * ) ;
double ** gortt_alloc_2d( int, int ) ;
void      gortt_dealloc_2d( double **, int ) ;
void      gortt_dealloc_3d( double ***, int, int ) ;


/*price and prospect*/
void   gortt_price_soil( gortt_spectra *, double *, double *, double *, double * );
void   gortt_prospect_interface( gortt_spectra * );
//void   leaftwo_( double *, double *, double *, double *, double *, double *, double * );
void   prospect_DB_( double *, double *, double *, double *, double *, double *, double *, double * );

/*read soil spectra*/
void	gortt_read_soil_lut( gortt_spectra * );


/*parsing and usage*/

void   gortt_cl_parser( int, char **, gortt_parameters *, gortt_geometry *, gortt_spectra *, gortt_control * );
void   gortt_usage( char * );
char	 get_first_string_token( char *, char * );


#endif /*GORTT*/





