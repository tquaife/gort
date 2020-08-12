#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gortt.h>


void gortt_energy( p, g, s )
/*
Calculate the SW energy balance in the canopy
The variables Fd1, Fd2, Fu1, Fd2 are the flux 
terms normalised to the incoming (i.e. Fd1=1.),
where:

d=down-welling
u=up-welling
1=above canopy
2=below canopy

Fu1 is the albedo.

Note that this routine assumes a Lambertian 
background reflectance at the moment. This 
would have to change if either Z or G are
modified to have a BRDF.
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{


  double Fd1=1., Fd2, Fu1, Fu2, G, Z, Pn0, rsoil ; 
	int k;

  gortt_albedo( p, g, s );

  Pn0=p->p_neq0_heq0_sza ;
  
  for( k=0;k<s->nw;k++ ){
  
    Fu1 = *(s->albedo+k) ;

    rsoil= *(s->rsoil+k) ;

    G = *( s->scomp + k*4+1 ) ;
    Z = *( s->scomp + k*4+3 ) ;

    Fu2 = G*Pn0+Z*(1.-Pn0);
    Fd2 = Pn0+Z*(1.-Pn0)/rsoil;

    *(s->favegt+k)=Fd1-Fu1-Fd2+Fu2;
    *(s->fasoil+k)=Fd2-Fu2 ;

    fprintf( stderr, "++ %f %f %f %f %f %f %f \n", Fd1, Fd2, Fu1, Fu2, G, Z, Pn0 );

  }

}

void gortt_albedo( p, g, s )
/*
Hemispherical integral of GORT reflectances.
Calculated for each spectral channel in s and
the solar geometry in g (n.b. we are integrating
over the viewing hemisphere so only the solar
zenith angle matters).
*/
gortt_parameters *p ;
gortt_geometry *g ;
gortt_spectra *s ;
{

  double xm, xr, ym, yr, x ,y ;
  double *sum_y, *sum_x ;
  int i, j, k ;

	sum_y = (double *)malloc( p->npoints * sizeof( double ) );
	sum_x = (double *)malloc( p->npoints * sizeof( double ) );

  xm=0.5*(1.-1.);
  xr=0.5*(1.+1.);
  ym=0.5*(2.*M_PI-0.);
  yr=0.5*(2.*M_PI+0.);


  for( k=0;k<s->nw;k++ ) *(sum_y+k) = 0.;
  for( i=0;i<p->npoints;i++ ){
    
    y = ym + yr* *(p->abscissa+i);
    
    g->vaa = y;
    
    /*corrections...*/
    while ( g->vaa > 2 * M_PI ) g->vaa -= 2 * M_PI;
    g->raa = g->saa - g->vaa ;
		g->raa = fabs((g->raa - 2 * M_PI * (int) (0.5 + g->raa * M_1_PI * 0.5)));
    
    for( k=0;k<s->nw;k++ ) *(sum_x+k)=0.;
    for( j=p->npoints/2.;j<p->npoints;j++ ){ 
    
      x=xm + xr* *(p->abscissa+j) ;

      g->vza = acos( x );
      
      /*corrections... un-needed!
      - but exercising caution...*/
      if ( g->vza < 0.0 ) {
			  g->vaa += M_PI ;
			  g->vza *= -1.0 ;
		  }
		
      g->vza_prime = gortt_prime_theta( p, g->vza );
		  g->sza_prime = gortt_prime_theta( p, g->sza );

  		/*this is needed by the kuusk function*/
  		p->k_vza = gortt_leaf_angle_distribution( p, g->vza ) ;

		  gortt_set_zenith_dependant_probabilities( p, g );
		 
		  /* Checking probability terms...
		  fprintf( stderr, "+++ %f\n", p->p_neq0_heq0_sza );
		  */
		 
		  gortt_rsurf( p, g, s );

      for( k=0;k<s->nw;k++ )
        *(sum_x+k)=*(sum_x+k)+*(s->rsurf+k)* *(p->weights+j)*fabs(x)*xr ;
    
    }
    for( k=0;k<s->nw;k++ )
      *(sum_y+k) = *(sum_y+k) + *(sum_x+k)* *(p->weights+i)*yr ;
  }     
  
  for( k=0;k<s->nw;k++ ) *(s->albedo+k) = *(sum_y+k) / M_PI ;

}


#define EPS 3.0e-11
void gauleg(x1,x2,x,w,n)
/**
Gauss-Legendre ascissa and weights using the numerical
recepies algorithm. This code taken from:
http://www.codeforge.com/read/33081/GAULEG.C__html

x1 = lower limit of integral
x2 = upper limit of integral
x  = array of abscissa
w  = array of weights
n  = number of points

NB - this version is modifed to use zero indexed arrays

**/
double x1,x2,*x,*w;
int n;
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;


	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
/*
	for (i=1;i<=m;i++)  {
*/
	for (i=0;i<m;i++)  {
/*
		z=cos(3.141592654*(i-0.25)/(n+0.5));
*/
		z=cos(3.141592654*(i+0.75)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		/*
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
		*/
		x[i]=xm-xl*z;
		x[n-1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}
#undef EPS

