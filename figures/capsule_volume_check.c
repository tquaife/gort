/*
 * gortt_vol against an independent Monte Carlo, and the capsule cross
 * section against a pixel count.
 *
 *   gcc -Wall -g -I../include -o capsule_volume_check capsule_volume_check.c \
 *       ../gortt_pn_kopen.o -lm
 *   ./capsule_volume_check
 *
 * gortt_vol(z_h, z_s, t, h_b) is the volume of the region of crown centres
 * whose crown meets the beam segment from z_h up to z_s, truncated below the
 * plane h_b.  gortt_calc_epgap_ht reads it twice, as a difference across
 * [h1_p, h2_p], to get the expected number of crowns a beam penetrates.  The
 * region is a capsule: the Minkowski sum of the segment with a ball of radius
 * r, which is what a crown is in the primed frame.
 *
 * This program checks two things.
 *
 * 1. gortt_capsule_cross_section, the horizontal cross section of that region,
 *    against a pixel count of the set {x : dist(x, segment) <= r} at the same
 *    height.  This is the integrand gortt_vol integrates, and it is what the
 *    model uses on its own when h2 == h1 collapses the crown centre layer
 *    onto a plane.  It agrees to the resolution of the count.
 *
 * 2. gortt_vol against a Monte Carlo of the same volume.  It agrees to the
 *    noise of the Monte Carlo when h_b is above the capsule -- where it
 *    returns its own V_0 = pi r^2 s + 4/3 pi r^3, the capsule volume -- and
 *    it does NOT agree when h_b cuts the capsule.  See OPEN-QUESTIONS.md #14.
 *
 * The two are independent of each other: the first tests a formula against a
 * pixel count, the second tests gortt's closed form against random points.
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gortt.h>

static gortt_parameters P;

/*
 * distance from (x,y,z) to the segment from (0,0,z_h) to (dx,0,z_s), the beam
 * lying in the x-z plane at zenith t so that dx = (z_s-z_h) tan t
 */
static double dseg( double x, double y, double z, double dx, double z_h, double z_s )
{
	double ux = dx, uz = z_s - z_h, L2 = ux*ux + uz*uz, q ;
	double px = x, pz = z - z_h ;

	q = ( px*ux + pz*uz ) / L2 ;
	if( q < 0.0 ) q = 0.0 ;
	if( q > 1.0 ) q = 1.0 ;
	px -= q*ux ; pz -= q*uz ;

	return sqrt( px*px + y*y + pz*pz ) ;
}

/* pixel count of the cross section at height z */
static double area_count( double t, double z_h, double z_s, double z )
{
	int    i, j, n = 3000 ;
	double dx = ( z_s - z_h )*tan( t ) ;
	double xlo = -P.r - 1.0, xhi = dx + P.r + 1.0 ;
	double ylo = -P.r - 1.0, yhi = P.r + 1.0 ;
	double hx = ( xhi - xlo )/n, hy = ( yhi - ylo )/n, area = 0.0 ;

	for( i=0; i<n; i++ )
		for( j=0; j<n; j++ )
			if( dseg( xlo + (i+0.5)*hx, ylo + (j+0.5)*hy, z, dx, z_h, z_s ) <= P.r )
				area += hx*hy ;

	return area ;
}

/* Monte Carlo of the volume below h_b */
static double vol_mc( double t, double z_h, double z_s, double h_b, long n )
{
	double dx = ( z_s - z_h )*tan( t ) ;
	double xlo = -P.r, xhi = dx + P.r ;
	double ylo = -P.r, yhi = P.r ;
	double zlo = z_h - P.r, zhi = MIN( h_b, z_s + P.r ) ;
	double box, hit = 0.0 ;
	long   i ;

	if( zhi <= zlo ) return 0.0 ;
	box = ( xhi - xlo )*( yhi - ylo )*( zhi - zlo ) ;

	srand48( 12345 ) ;
	for( i=0; i<n; i++ ){
		double x = xlo + drand48()*( xhi - xlo ) ;
		double y = ylo + drand48()*( yhi - ylo ) ;
		double z = zlo + drand48()*( zhi - zlo ) ;
		if( dseg( x, y, z, dx, z_h, z_s ) <= P.r ) hit += 1.0 ;
	}

	return box*hit/(double) n ;
}

int main( void )
{
	double tdeg[] = { 0.0, 20.0, 45.0, 70.0 } ;
	double z_h = 0.0, z_s = 2.5, t, a, c, worst_a = 0.0, worst_v = 0.0 ;
	int    it, k ;

	P.r = 1.7 ; P.rr = P.r*P.r ; P.rrr = P.rr*P.r ;
	P.theta_p = calloc( 8, sizeof( double ) ) ;

	printf( "1. gortt_capsule_cross_section against a pixel count\n\n" ) ;
	printf( "   %5s %7s %12s %12s %10s\n", "t", "z", "formula", "count", "diff" ) ;

	for( it=0; it<4; it++ ){
		t = tdeg[it]*M_PI/180.0 ;
		for( k=-2; k<=5; k++ ){
			double z = z_h - P.r + k*( z_s - z_h + 2.0*P.r )/6.0 ;
			a = gortt_capsule_cross_section( &P, t, z_h, z_s, z ) ;
			c = area_count( t, z_h, z_s, z ) ;
			if( fabs( a - c ) > worst_a ) worst_a = fabs( a - c ) ;
			printf( "   %5.1f %7.2f %12.6f %12.6f %10.2e\n", tdeg[it], z, a, c, a - c ) ;
		}
	}
	printf( "\n   worst = %.2e, against a pixel of %.2e\n\n", worst_a, ( 2.0*P.r + 2.0 )/3000.0 ) ;

	printf( "2. gortt_vol against a Monte Carlo of the same volume\n\n" ) ;
	printf( "   %5s %7s %12s %12s %10s %8s\n", "t", "h_b", "gortt_vol", "MC", "diff", "%" ) ;

	for( it=0; it<4; it++ ){
		double hbs[] = { -0.26, 0.5, 1.11, 1.5, 4.3 } ;
		t = tdeg[it]*M_PI/180.0 ;
		P.theta_p[0] = t ;
		for( k=0; k<5; k++ ){
			a = gortt_vol( &P, z_h, z_s, 0, hbs[k] ) ;
			c = vol_mc( t, z_h, z_s, hbs[k], 40000000L ) ;
			if( c > 0.0 && fabs( a - c )/c > worst_v ) worst_v = fabs( a - c )/c ;
			printf( "   %5.1f %7.2f %12.6f %12.6f %10.2e %7.2f%%\n",
			        tdeg[it], hbs[k], a, c, a - c, c > 0.0 ? 100.0*( a - c )/c : 0.0 ) ;
		}
	}
	printf( "\n   worst = %.2f%%.  h_b = 4.30 is above the capsule; the rest cut it.\n", 100.0*worst_v ) ;

	return 0 ;
}
