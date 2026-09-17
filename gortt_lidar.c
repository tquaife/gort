#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gortt.h>


/*
The canopy lidar waveform of Ni-Meister, Jupp and Dubayah (2001), IEEE
TGRS 39(9), 1943-1958, built on the gap probabilities GORT has already
computed.

The whole model is three equations. With P(z) the probability of a gap
above height z at nadir, J0 the emitted pulse energy, rho_v the volume
backscattering coefficient of the foliage and rho_g that of the ground:

	-dRv(z)/dz = J0 rho_v dP(z)/dz          the canopy return   (eq 1)
	Rg         = J0 rho_g P(0)              the ground return   (eq 2)
	Fapp(z)    = dlogP(z)/dz                apparent foliage    (eq 4)

This is a one-way attenuation model, and that is only right at the
hotspot -- every shadow in the direction of the outgoing pulse is
occluded, so a scattering surface the pulse can reach is also one the
receiver can see. Nadir viewing guarantees it, which is why gortt_lidar
does not take a gortt_geometry: the device points straight down, nadir
is theta[0] == 0 exactly, and nothing has to be interpolated.

The two backscattering coefficients come from the phase function at zero
scattering angle. For randomly oriented Lambertian leaves (eq 9), which
is the G = 0.5 case the rest of gortt is hardwired to, p(xi) at xi = 0
is 8t/(3 omega) (eq 12), so from eq 6

	rho_v = omega p(0) / 4 = 2t/3

with t the leaf transmittance. The paper quotes rho_v = omega/5, which
is this under its own assumption that t = 0.3 omega at 1.064um; taking
t straight from PROSPECT-D instead is the same expression without the
assumption. For a Lambertian ground p(xi) = 2 in backscatter (eq 14),
so rho_g = omega_g/2 (eqs 15, 16).
*/


static int    gortt_lidar_nbins( const gortt_parameters *, const gortt_lidar_params *, double, double * );
static void   gortt_lidar_deposit( double *, int, double, double, double, double );
static void   gortt_lidar_convolve( double *, int, double, double );
static void   gortt_lidar_print_profile( const gortt_parameters *, const gortt_lidar_params *,
                                         const double *, const double *, double, double, double, char * );
static void   gortt_lidar_print_waveform( const gortt_parameters *, const gortt_lidar_params *,
                                          const double *, double, double, double, char * );


void gortt_lidar( const gortt_parameters *p, const gortt_spectra *s,
                  const gortt_lidar_params *l, char *fmt )
{

	int h ;

	double rho_v, rho_g, Rg ;
	double *P, *dPdz ;


	/*
	The within-crown term is only filled for rows h < nh_probs, and the
	rest of the table is left at its calloc zero. gortt_cl_parser raises
	nh_probs to nlayers for -lidar; if that ever stops happening the
	waveform would come out quietly wrong rather than obviously wrong,
	so check it here rather than trusting it.
	*/
	if( p->nh_probs < p->nlayers ){
		fprintf( stderr, "gortt_lidar: gap probabilities are only available for %d of %d layers\n",
		         p->nh_probs, p->nlayers );
		exit( EXIT_FAILURE );
	}

	P    = gortt_alloc_1d( p->nlayers );
	dPdz = gortt_alloc_1d( p->nlayers );

	/*
	Total gap probability at nadir: the between-crown term P(n=0) plus
	the within-crown term P(n>0), which is what eq 3 and the paper's
	Fig. 2 mean by P(z). theta[0] is zero exactly.
	*/
	for( h=0; h<p->nlayers; h++ )
		P[ h ] = p->p_n0[ h ][ 0 ] + p->epgap[ h ][ 0 ] ;

	/*centred difference inside, one-sided at both ends*/
	for( h=1; h<p->nlayers-1; h++ )
		dPdz[ h ] = ( P[ h+1 ] - P[ h-1 ] ) / ( 2.0 * p->dz ) ;

	dPdz[ 0 ]            = ( P[ 1 ] - P[ 0 ] ) / p->dz ;
	dPdz[ p->nlayers-1 ] = ( P[ p->nlayers-1 ] - P[ p->nlayers-2 ] ) / p->dz ;

	rho_v = 2.0 * s->tleaf[ 0 ] / 3.0 ;
	rho_g = 0.5 * s->rsoil[ 0 ] ;

	/*
	P[0] is the gap probability at z1 = h1 - b, the lowest height at
	which there is any foliage at all. Nothing below it attenuates, so
	it is also the gap probability at the ground: the P(0) of eqs 2, 21.
	*/
	Rg = l->J0 * rho_g * P[ 0 ] ;

	if( l->binned_output == TRUE )
		gortt_lidar_print_waveform( p, l, P, rho_v, rho_g, Rg, fmt );
	else
		gortt_lidar_print_profile( p, l, P, dPdz, rho_v, rho_g, Rg, fmt );

	gortt_dealloc_1d( P );
	gortt_dealloc_1d( dPdz );

	return;

}


/*
The raw per-layer profile: the quantities Ni-Meister et al. plot in their
Figs. 4 to 6, on gortt's own height grid and with no instrument model
applied.

Rows run top-down, h = nlayers-1 to 0, which is the order in which the
return arrives. Rg goes on the header line rather than in a data row
because it is an energy while Rv is an energy density -- putting it in
the Rv column would mix units within a column.
*/
static void gortt_lidar_print_profile( const gortt_parameters *p, const gortt_lidar_params *l,
                                       const double *P, const double *dPdz,
                                       double rho_v, double rho_g, double Rg, char *fmt )
{

	int h ;

	printf( "%d ", p->nlayers );
	printf( fmt, l->wl );
	printf( fmt, rho_v );
	printf( fmt, rho_g );
	printf( fmt, l->J0 );
	printf( fmt, Rg );
	printf( "\n" );

	for( h=p->nlayers-1; h>=0; h-- ){

		printf( fmt, p->height[ h ] );
		printf( fmt, P[ h ] );
		printf( fmt, dPdz[ h ] );

		/*the apparent foliage profile, eq 4. P is a probability and
		cannot be negative, but a black canopy drives it to zero and
		the ratio with it.*/
		printf( fmt, P[ h ] > 0.0 ? dPdz[ h ] / P[ h ] : 0.0 );

		/*the canopy return density, eq 1*/
		printf( fmt, l->J0 * rho_v * dPdz[ h ] );

		printf( "\n" );

	}

	return ;

}


/*
The waveform as an instrument would record it: canopy and ground on one
uniform range grid, convolved with the transmitted pulse (the paper's
Fig. 7, where a finite pulse width is what lets the ground delta and the
canopy density be plotted together).

The canopy energy is deposited from the layer-to-layer increments of P
rather than by differentiating and then re-integrating dP/dz. That is
both more direct and exactly conservative: the canopy bins sum to
J0 rho_v ( P(top) - P(0) ), which is eq 22's Rv(0) = J0 rho_v (1-P(0))
whenever P reaches one at the canopy top.
*/
static void gortt_lidar_print_waveform( const gortt_parameters *p, const gortt_lidar_params *l,
                                        const double *P, double rho_v, double rho_g, double Rg,
                                        char *fmt )
{

	int h, k, nbins ;

	double db, ztop, z ;
	double *w ;


	db = ( l->binsize > 0.0 ) ? l->binsize : p->dz ;

	nbins = gortt_lidar_nbins( p, l, db, &ztop );
	w = gortt_alloc_1d( nbins );

	/*
	Canopy: the beam intercepted between layers h and h+1 is
	P[h+1] - P[h], and it is returned from somewhere in between, so put
	it at the midpoint.
	*/
	for( h=0; h<p->nlayers-1; h++ ){

		z = 0.5 * ( p->height[ h ] + p->height[ h+1 ] ) ;

		gortt_lidar_deposit( w, nbins, db, ztop, z,
		                     l->J0 * rho_v * ( P[ h+1 ] - P[ h ] ) );

	}

	/*ground: a delta at z = 0 until the pulse smears it*/
	gortt_lidar_deposit( w, nbins, db, ztop, 0.0, Rg );

	if( l->pulse > 0.0 ) gortt_lidar_convolve( w, nbins, db, l->pulse );

	printf( "%d ", nbins );
	printf( fmt, l->wl );
	printf( fmt, rho_v );
	printf( fmt, rho_g );
	printf( fmt, l->J0 );
	printf( fmt, db );
	printf( fmt, l->pulse );
	printf( "\n" );

	/*
	Printed as a density, w/db, so that it is directly comparable with
	the Rv column of the profile output; multiplying back by db and
	summing recovers the total returned energy.
	*/
	for( k=0; k<nbins; k++ ){
		printf( fmt, ztop - (double)k * db );
		printf( fmt, w[ k ] / db );
		printf( "\n" );
	}

	gortt_dealloc_1d( w );

	return ;

}


/*
Put a point source of energy e at height z onto the range grid, split
between the two bins it falls between in proportion to how close it is
to each (linear, or "cloud in cell", deposition). Bin k is centred on
ztop - k db, so a larger z is a smaller k.

Dropping the whole of e into the nearest bin instead is a sampling of
the waveform rather than a binning of it, and it aliases: the canopy
sources are the nlayers-1 layer midpoints, evenly spaced by dz, so when
dz > db the bins alternate between catching one source and catching
none and the waveform comes out as a staircase. It bites whenever
-lidar_binsize is set finer than the layer spacing, and gortt spreads
its layers over z1..z2, so dz depends on how deep the canopy is and
tall canopies get there first: over 120 layers, b = 5 m crowns give
dz = 0.12 m, already coarser than a 0.1 m bin, where b = 1 m crowns
give 0.05 m. Splitting linearly removes it: the result
is exact for a source at a bin centre and conserves both the total
energy and its first moment everywhere else, at the cost of convolving
the waveform with a db-wide triangle, which adds db^2/6 to the variance
of a feature and is well inside the pulse width whenever the binsize is.

The bounds are checked per bin rather than per source so that a source
in the outermost half-bin still deposits the share that lands on the
grid. gortt_lidar_nbins pads by three sigma at both ends, so in practice
nothing is near an edge.
*/
static void gortt_lidar_deposit( double *w, int nbins, double db, double ztop,
                                 double z, double e )
{

	int k ;

	double x, f ;


	x = ( ztop - z ) / db ;
	k = (int)floor( x ) ;
	f = x - (double)k ;

	if( k   >= 0 && k   < nbins ) w[ k   ] += ( 1.0 - f ) * e ;
	if( k+1 >= 0 && k+1 < nbins ) w[ k+1 ] +=         f   * e ;

	return ;

}


/*
Size the range grid and return the height of its first (highest) bin.
It has to reach beyond the canopy at the top and below the ground at
the bottom, because the pulse moves energy both ways; three sigma is
where the Gaussian used below is truncated, and one bin is the minimum
that keeps the canopy top and the ground off the very edge.
*/
static int gortt_lidar_nbins( const gortt_parameters *p, const gortt_lidar_params *l,
                              double db, double *ztop )
{

	double pad ;

	pad = MAX( 3.0 * l->pulse, db );

	*ztop = p->z2 + pad ;

	/*
	n.b. round rather than truncate, and keep a spare bin at each end.
	The pad is exactly three sigma, so the kernel of a source sitting in
	the very last bin reaches exactly the end of the grid -- and a
	truncating divide here is then one bin short of that whenever the
	division lands just under an integer, which silently drops energy
	off the bottom. Two spare bins cost nothing and make the convolution
	conservative for any binsize.
	*/
	return (int)( ( *ztop + pad ) / db + 0.5 ) + 2 ;

}


/*
Convolve the binned waveform with a normalised Gaussian of width sigma,
truncated at three sigma. Energy that lands outside the grid is dropped
rather than folded back, which is why gortt_lidar_nbins pads by three
sigma at both ends: with that padding there is nothing left to drop.
*/
static void gortt_lidar_convolve( double *w, int nbins, double db, double sigma )
{

	int i, j, k, half ;

	double x, sum ;
	double *kern, *out ;


	half = (int)( 3.0 * sigma / db + 0.5 ) ;
	if( half < 1 ) return ;

	kern = gortt_alloc_1d( 2*half + 1 );
	out  = gortt_alloc_1d( nbins );

	sum = 0.0 ;
	for( i=-half; i<=half; i++ ){
		x = (double)i * db / sigma ;
		kern[ i+half ] = exp( -0.5 * x * x ) ;
		sum += kern[ i+half ] ;
	}

	/*normalise on the truncated kernel, not analytically, so that the
	discrete convolution conserves energy exactly*/
	for( i=0; i<2*half+1; i++ ) kern[ i ] /= sum ;

	for( k=0; k<nbins; k++ ){
		for( i=-half; i<=half; i++ ){
			j = k + i ;
			if( j >= 0 && j < nbins ) out[ j ] += w[ k ] * kern[ i+half ] ;
		}
	}

	for( k=0; k<nbins; k++ ) w[ k ] = out[ k ] ;

	gortt_dealloc_1d( kern );
	gortt_dealloc_1d( out );

	return ;

}
