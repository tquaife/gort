#ifndef GORTT
#define GORTT

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

/*Gauss-Legendre points for averaging the shadow overlap over crown-centre
heights in [h1,h2] (gortt_overlap).  24 puts Kg within ~1e-5 of the limit;
convergence is only algebraic, because the overlap has a kink where it
reaches zero.*/
#define NH_OVERLAP 24

/*
Quadrature resolutions for the gap probabilities. These are deliberately
NOT derived from nlayers: nlayers is the output grid, and the accuracy of
the model must not depend on how finely the user asks for it to be
reported. All three are scaled by p->nquad (the -nquad option), which
exists so that convergence can be demonstrated from the command line.
*/
#define NQ_VGAMMA  24	/*GL points across the non-smooth band of Vgamma*/

/*prospect*/
#define PROSPECT_NWBANDS 2101
#define PROSPECT_LOWER_WL 400.0
#define PROSPECT_SPECTRAL_RESLN 1.0

/*the grid -soil_spectra is resampled onto: 400-2500nm at 1nm*/
#define SOIL_LUT_NBANDS 2101
#define SOIL_LUT_LOWER_WL 400.0
#define SOIL_LUT_SPECTRAL_RESLN 1.0


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
The lidar instrument: everything the waveform needs that is neither the
canopy nor the scene spectra. Written by gortt_cl_parser, read by
gortt_lidar.

Unlike p/g/s/c this is not threaded through the model. Lidar mode
short-circuits main() -- the device is always at nadir, so there is no
per-angle loop for it to be threaded through -- and gortt_lidar is the
only function that ever sees it.
*/

typedef struct  {

    double wl ;                     /*wavelength, nm*/
    double J0 ;                     /*emitted pulse energy*/
    double pulse ;                  /*Gaussian pulse sigma, m; 0 = no convolution*/
    double binsize ;                /*range bin size, m; <=0 means use p->dz*/
    unsigned short binned_output ;  /*emit the binned waveform, not the profile*/

} gortt_lidar_params ;



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

		/*
		the shaded-background signature Z divided by the soil
		reflectance, i.e. the fraction of the incident radiation that
		reaches the background. gortt_rsurf builds Z as this times
		rsoil, so recovering it by division is exact for any soil but
		0/0 for a black one, which is why it is carried separately.
		gortt_energy needs it; nothing else does.
		*/
		double *ztrans ;

		/*number of wavebands*/
		int nw ;

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
		/*leaf reflectance and transmittance are held separately: the
		scattering asymmetry in gortt_rt_init is built from their
		difference, so a single albedo cannot express a real leaf.*/
		double user_refl_leaf ;
		double user_trans_leaf ;
		double user_r_soil ;

		int is_file_soil ;
		char soil_file[ MAX_LINE_LEN ] ;
		double soil_spectra[ SOIL_LUT_NBANDS ] ;

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
		int nquad ;			/*multiplier on the NQ_* quadrature orders*/

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
		/*crown-centre heights and weights (summing to one) at which
		gortt_overlap evaluates and averages the overlap: uniform in [h1,h2]*/
		double h_overlap[ NH_OVERLAP ] ;
		double w_overlap[ NH_OVERLAP ] ;
		/*Gauss-Legendre rule on [-1,1], NQ_VGAMMA*nquad points, mapped
		onto whatever interval gortt_crown_proj_volume needs it for*/
		int nq_vg ;
		double *xq_vg ;
		double *wq_vg ;
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

		Only rows h < nh_probs of the height-indexed tables are filled;
		the rest are left at their calloc zero. The BRF only ever needs
		h = 0, so nh_probs defaults to 1, but a lidar waveform is a
		function of height and -lidar raises it to nlayers.
		*/

		int nh_probs ;

		double **p_s0 ;
		double **p_n0 ;
		double **epgap ;
		double *k_open ;
		double *k_openep ;

		} gortt_parameters ;


/*
The radiative transfer half of the model: every term Ni et al. (1999)
derive between the leaf optics and the four scene signatures, for one
geometry and one waveband. Filled by gortt_rt_init.

Reading the names: _df is direct-to-diffuse, _ff diffuse-to-diffuse and
_0 direct unscattered; t/T is a transmittance, p/R a reflectance; _inf
is the semi-infinite homogeneous layer; and _prime is the discontinuous
canopy as opposed to the equivalent homogeneous one.
*/
typedef struct  {

		double w ;			/*leaf single scattering albedo*/
		double gamma ;		/*sqrt( 1 - w )*/
		double g_asym ;		/*phase function asymmetry*/

		/*semi-infinite homogeneous layer*/
		double R_inf_df ;
		double R_inf_ff ;
		double T_inf_df ;
		double T_inf_ff ;

		/*finite homogeneous layer*/
		double t_0 ;
		double t_df ;
		double t_ff ;
		double p_df ;
		double p_ff ;

		/*discontinuous canopy*/
		double t_prime_0 ;
		double t_prime_df ;
		double t_prime_ff ;
		double p_prime_df ;
		double p_prime_ff ;

		/*Kuusk's bi-directional gap probability*/
		double kuusk ;

		} gortt_rt ;



/*
TQ's implimentation of the GORT paper.

Everything below takes gortt_parameters as const: once gortt_init_params
and the gap-probability routines have run, the canopy description does
not change for the rest of the run. gortt_geometry is const wherever the
function only evaluates at the geometry it is given, i.e. everywhere
except gortt_rsurf, which fills in the four viewed proportions.
*/
void   gortt_rsurf( const gortt_parameters *, gortt_geometry *, gortt_spectra * );
void   gortt_gap_probabilities( gortt_parameters * );
void   gortt_gap_probabilities_Q08( gortt_parameters * );
void   gortt_set_zenith_dependant_probabilities( const gortt_parameters *, gortt_geometry * );
double gortt_prime_theta( const gortt_parameters *, double  );
double gortt_kg( const gortt_parameters *, const gortt_geometry *, double );
double gortt_kc( const gortt_parameters *, const gortt_geometry *, double );
void   gortt_kc_fFbeta( const gortt_parameters *, const gortt_geometry *, double, double, double *, double *, double * );
double gortt_kc_ambrals( const gortt_parameters *, const gortt_geometry * 	);
double gortt_overlap( const gortt_parameters *, const gortt_geometry *, double );
void   gortt_rt_init( const gortt_parameters *, const gortt_geometry *, const gortt_spectra *, int, gortt_rt * );
double gortt_kuusk( const gortt_parameters *, const gortt_geometry * );
double gortt_leaf_angle_distribution( const gortt_parameters *, double );

/*
The lidar waveform. No gortt_geometry: the device is at nadir, which is
theta[0] == 0 exactly, so nothing has to be interpolated. The char * is
the GORTT_PREC output format, so the waveform prints as sharply as the
BRF does.
*/
void  gortt_lidar( const gortt_parameters *, const gortt_spectra *,
                   const gortt_lidar_params *, char * );


/*TQ's integral routines*/

void   gortt_energy( const gortt_parameters *, const gortt_geometry *, gortt_spectra * );
void   gortt_albedo( const gortt_parameters *, const gortt_geometry *, gortt_spectra * );
void   gauleg( double, double, double *, double *, int );

/*Wenge's functions for within crown probs*/
void   gortt_init_params( gortt_parameters *, gortt_geometry * );

double gortt_crown_proj_volume( gortt_parameters *, double, double );
double gortt_crown_proj_cross_section( gortt_parameters *, double, double, double );
double gortt_weird_cross_section( gortt_parameters *, double, double, double );
double gortt_left_circle_area( double, double );
double gortt_right_ellipse_area( double, double, double );
double gortt_calc_epgap_ht( gortt_parameters *, int, int );

double gortt_get_es( gortt_parameters *, double, int );
double gortt_get_s( gortt_parameters *, double, double, int );
double gortt_get_pcc( gortt_parameters *, double  );
double gortt_vol( gortt_parameters *, double, double, int, double ) ;

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

void   gortt_cl_parser( int, char **, gortt_parameters *, gortt_spectra *, gortt_control *,
                        gortt_lidar_params * );
void   gortt_usage( char * );
char	 get_first_string_token( char *, char * );


#endif /*GORTT*/





