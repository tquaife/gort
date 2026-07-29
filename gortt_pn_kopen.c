#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gortt.h>
#include<float.h>

void gortt_gap_probabilities( gortt_parameters *p )
{

	int i ;

	int h;	/* 'height' index */
	int t;	/* 'theta' index */

	int s_max_p ;

	double v_g ;	/* Vgamma at this height and theta */

	/*
	Find probabilities of photon not hitting a
	crown, given height and theta
	*/

	for (t = 0; t < p->nth; t++) {
		for (h = 0; h < p->nlayers; h++) {
			/*
			* Find P(n=0), using Vgamma.
			*/
			v_g = gortt_crown_proj_volume( p, p->theta_p[t], p->height_p[h]  );
			p->p_n0[h][t] = exp(-1.0 * p->lv_p * v_g );
		}
	}


	/*
	* Get P(s=0) for all heights.  It's just the differential of P(n=0) for
	* decreasing height.
	*/

	for (t=0; t<p->nth; t++){
		p->p_s0[p->nlayers - 1][t] = 0.0;
		for ( h=p->nlayers-2; h>=0; h-- ){
			p->p_s0[h][t] = p->p_n0[h + 1][t] - p->p_n0[h][t];
		}
	}

	/*
	* Get probability distributions for s for each possible h and theta.
	*
	* n.b. this is the expensive loop, and it is bounded by nh_probs
	* rather than nlayers: nothing downstream reads pd_s above h=0
	* unless a height profile has been asked for. P(s=0) above is still
	* filled for every height, because gortt_get_pd_s at h=0 integrates
	* over the p_s0 values at all the heights above it.
	*/
	for (t = 0; t < p->nth; t++) {
		for (h = 0; h < p->nh_probs; h++) {

			/*
			* figure out maximum path length at this height and this angle.
			*/

			s_max_p = gortt_s_to_index( p, (p->z2_p - p->height_p[h]) / cos( p->theta_p[t]) );

			/*
			fprintf(stderr,"s_max_p %d\n",s_max_p);
			fflush( stderr );
			*/

			for ( i=0; i<(s_max_p+PD_S_BUFF); i++ ){
				p->pd_s[h][t][i] = 0.0;
			}

			gortt_get_pd_s( p, h, t );

		}
	}



	gortt_calc_epgap( p );
	gortt_calc_kopen( p );

	return ;
}




int gortt_s_to_index( gortt_parameters *p, double s )
{
	return ( (int) ( s / p->ds + 0.5 ) );
}

double gortt_index_to_s( gortt_parameters *p, int index )
{
	return ( (double) index * p->ds );
}


double gortt_crown_proj_volume( gortt_parameters *p, double t, double h )
{

	/*
	* Calculate crown projection volume for a crown at the given height, by
	* numerical integration, using the midpoint rule.
	*/

	double	vol = 0.0;
	double	z;

	for ( z=p->h1_p+p->dz_p/2.0; z<=p->h2_p; z+=p->dz_p ) {
		vol += gortt_crown_proj_cross_section( p, t, h, z ) * (p->dz_p);
	}

	return( vol );
}


double gortt_crown_proj_cross_section( gortt_parameters *p, double t, double h, double z )
{

	/* all the calculation is in transformed dimension */
	/*
	* Calculate and return the cross-sectional area of the crown projection
	* for height h at height z.  The shape is very strange.  For a given
	* height h, the horizontal cross section at z < h - r sin (th) is
	* circular. For z > h + r sin (th), the cross section is elliptical.
	* In between, the shape is part circle and part ellipse.
	*/

	double	h_low;
	double	h_high;
	double	csa;		/* cross-sectional area */

	if (z < h-p->r ) return( 0.0 );

	h_low = h - p->r * sin( t );
	h_high = h + p->r * sin( t );


	if (z <= h_low) {

		/*
		* Cross-section is circular.
		*/

		double r_p, a;

		a = p->rr - (h - z) * (h - z) ;

		if( a <= 0 ) r_p = 0 ;
		else r_p = sqrt( a );

		csa = M_PI * r_p * r_p;

	}else if (z > h_low && z < h_high) {

		/*
		* Cross-section is neither a circle nor an ellipse.
		*/

		csa = gortt_weird_cross_section( p, t, h, z );

	}else{

		/*
		* Cross-section is an ellipse.
		*/

		csa = M_PI * p->rr * SEC( t );
   }


   return( csa );

}



double gortt_weird_cross_section( gortt_parameters *p, double t, double h, double z )
{

	/* r' - radius of the circular part of the cross section at the height 'z' */
	double	r_p;

	/* x-coord of the circle center, relative to the ellipse center */
	double	x_cc;

	/* x-coord of the line separating the circle
	part and the ellipse part of the cross section. */
	double	x_p;

	/* area of the circle part of the cross-section. */
	double	a_cp;

	/* area of the ellipse part of the cross section. */
	double	a_ep;

	/* total area */
	double	a;

	double	zdiff;

	zdiff = h - z;

	r_p = sqrt( p->rr - zdiff * zdiff);
	x_cc = zdiff * tan(t);
	x_p = x_cc / (1.0 - cos(t) * cos(t));

	/*
	* So now we know the cut-off line for each partial area. We need to find:
	*
	* 1. The area of the circle of radius r', centered at zero, and cut off by
	* the line x = x_p - x_cc.
	*
	* 2. The area of the ellipse to the right of the line x = x_p.
	*
	* Of course, the terms 'left' and 'right' here are only relative.
	*/

	a_cp = gortt_left_circle_area( r_p, x_p - x_cc );
	a_ep = gortt_right_ellipse_area( p->r, p->r * SEC( t ), x_p );
	a = a_cp + a_ep;

	return( a );

}


double gortt_left_circle_area( double r, double x_cut )
{

	double            area_tot;
	double            ang_sector;
	double            area_sector;
	double            area_triangle;
	double            a;

	area_tot = M_PI * r * r;
	ang_sector = acos(fabs(x_cut) / r) * 2.0;
	area_sector = area_tot * ang_sector / ( 2.0 * M_PI );
	area_triangle = fabs(x_cut) * sqrt(r * r - x_cut * x_cut);
	if (x_cut > 0.0)
		a = area_tot - (area_sector - area_triangle);
	else
		a = area_sector - area_triangle;

	return( a );
}



double gortt_right_ellipse_area( double r, double b, double x_cut )
{

	double	x_cut_p;
	double	a_p;

	x_cut_p = x_cut / (b / r);
	a_p = M_PI * r * r;

	a_p -= gortt_left_circle_area( r, x_cut_p );

	return( a_p * (b / r) );

}





/*
* Get k_open and delta(K_open).  This involves integrating P(n=0) and P(s=0)
* over theta.  We re-transform theta back to its original space -- that is,
* not adjusted by crown ellipticity.  Use the trapezoidal rule.
*/
void gortt_calc_kopen( gortt_parameters *p )
{

	/*
	these are gobal in Wenge's code...
	*/

	double tmp1, tmp1_last;
	double tmp2, tmp2_last;
	int h, t;


	for (h = 0; h < p->nh_probs; h++) {

		p->k_open[h] = 0.0 ;
		p->k_openep[h] = 0.0;

		tmp1_last = p->p_n0[h][0] * sin(2.0*p->theta[0]);
		tmp2_last = p->epgap[h][0] * sin(2.0*p->theta[0]);

		for (t = 1; t < p->nth; t++) {

			tmp1 = p->p_n0[h][t] * sin(2.0*p->theta[t]);
			p->k_open[h] += (tmp1 + tmp1_last) / 2.0 * p->dth;
			tmp1_last = tmp1;

			tmp2 = p->epgap[h][t] * sin(2.0*p->theta[t]);
			p->k_openep[h] += (tmp2 + tmp2_last) / 2.0 * p->dth;
			tmp2_last = tmp2;

		}

		/*
    if(h==0) printf("\n%8.4f %8.4f",p->k_open[h],p->k_openep[h]);
   	*/

	}

	return ;

}

static double gortt_trisec( double, double, double, double );
static double gortt_cylind( double, double, double, double );
static double gortt_triang_fcn( double, double, double, double );
static double gortt_cylind_fcn( double, double );
static double gortt_sector( double, double, double );
static double gortt_triang( double, double, double, int );

void gortt_get_pd_s( gortt_parameters *p, int h, int t )
{

	/*
	* Get a probability distribution for the within-crown path length 's' at
	* the given height and zenith angle. This routine fills the array
	* PD_s[h][t][...].
	*/

	/*static int	first_time = 0; */

	/* index to s' array */
	int	sp_i;

	/* index to P_s0 array */
	/*int	ps0_i;*/

	/*int	i, maxs;*/

	/* number-of-crowns counter */
	int	n ;

	/* probability associated with a given s'  value at a given height. */
	double	P_s_p;

	/* probability of penetrating 'n' crowns at a given s' value */
	double	P_n;

	double	temp1;
	/*double	temp2;*/
	double	s;

	/*double	tmp_Lv; */

	/*
	* es is the average distance that a beam must pass through a single
	* crown in order to reach this height. s_p is the after-entering-the-
	* crown pathlength for the current sp_i. Both used to be arrays on
	* gortt_parameters, but neither is read outside this function and
	* neither was indexed by everything it varies with.
	*/

	double	es;
	double	s_p;


	es = gortt_get_es( p, h, t );

	/*
	* Iterate over our possible s' (s_p[]) values.  The s' values are
	* after-entering-the-crown pathlengths, which are not equivalent to
	* within-crown pathlengths since they don't account for exiting from
	* crowns.
	*/

	for( sp_i = p->nlayers-1; sp_i>h; sp_i--) {

		/*
		* The s' array (s_p[][]) holds after-entering-the-crown pathlengths
		* associated with a given zenith angle and within-canopy height.
		*/

		s_p = (double)(p->height_p[sp_i] - p->height_p[h]) / cos(p->theta_p[t]);

		if(sp_i == ( p->nlayers-1 ) ) {
			/*
			* Path length is zero.  Calculations wouldn't work, but we know
			* that the appropriate index is zero and we know what the
			* probability is.
			*/
			p->pd_s[h][t][0] += p->p_s0[sp_i][t];
			continue;

		}

		/*
		* Get the probability associated with the current value of s' at this
		* height.
		*/

		P_s_p = p->p_s0[sp_i][t];

		/*
		* Iterate over number of crowns penetrated by a beam given the
		* pathlength s'.
		*/

		for( n=1; n<=p->maxcrowns; n++ ) {

			/*
			* Get the probability that 'n' crowns are penetrated, given the
			* current after-entering-crown pathlength s_p[i].
			*/

			temp1 = gortt_vol( p, h, sp_i, t, p->h2_p) - gortt_vol( p, h, sp_i, t, p->h1_p);
			temp1 *= p->lv_p;



			P_n = (pow(temp1, (double) n) * exp(-temp1)) /
				(p->factorial[n] * (1.0 - exp(-temp1)));

			/*
			* Get the mean within-crown pathlength.
			*/

			s = s_p * (1.0 - exp(-1.0 * (double) n * es / s_p));

			/* convert back to original elliptical dimension
			s *= g.ellipticity * cos(theta_p[t])/cos(theta[t]); */

			/*
			* Accumulate the appropriate value.  We use the calculated path
			* length 's' as an index into the array that we are filling, and
			* increment its value by the product of the probability that 'n'
			* crowns are penetrated and the probability of finding this
			* particular after-entering-the-crown pathlength.
			*/

			p->pd_s[h][t][gortt_s_to_index( p, s )] += P_n * P_s_p;



		}
  }

	return;

}


double gortt_get_es( gortt_parameters *p, int z, int t )
{


	/*
	* Return the expected value of S(z, h) -- on average, the distance that a
	* beam passes through a single crown in order to reach the height z. This
	* could be seriously optimized.  But anyway, we integrate over heights at
	* which such a sphere could be centered.
	*/


	double h;
	double dh;
	double ES = 0.0;



	dh = (p->h2_p - p->h1_p) / (double) p->nh_es;


	for (h = p->h1_p + dh / 2.0; h <= p->h2_p; h += dh)
		ES += gortt_get_s( p, z, h, t) * ( gortt_get_pcc( p, h ) * dh );



   return (ES);
}


double gortt_get_s( gortt_parameters *p, int z, double h, int th )
{
	double /*th_p,*/ S;



	/*
	* Return the average distance a beam must pass through a single crown
	* centered at height h to reach height z.  We calculate and return the
	* volume of a sphere centered at h above the height z, divided by its
	* projected area.
	*/

	if (p->height_p[z] > h + p->r - 0.0001){
		/*
		* Plane is above the sphere.  Average distance is zero.
		*/
		S = 0.0;


	}else if(p->height_p[z] < h - p->r + 0.0001){
      /*
       * Plane is below the sphere.  Average distance is the volume of the
       * sphere divided by its projected area.
       */
      S = 4.0 * p->r/ 3.0;


	}else{

		double proj_area = 0.0;
		double V_sphere = 0.0;
		double V_slice = 0.0;
		double V_tot = 0.0;
		double ht = 0.0;
		double zdiff = 0.0;

		V_sphere = 4.0 * M_PI * p->rrr / 3.0;

		zdiff = fabs(h - p->height_p[z]);

		ht = p->r - zdiff;
		V_slice = M_PI * ht * ht / 3.0 * (3.0 * p->r - ht);

		if (p->height_p[z] > h)
			V_tot = V_slice;
		else
			V_tot = V_sphere - V_slice;


		/* Adjust V_tot by the solar zenith angle */
		V_tot /= cos(p->theta_p[th]);


		/*
		* Get projected area. The 'crown_proj_cross_section' function assumes
		* that we are projecting _towards_ the sun, as is the case for the
		* calculation of V_gamma and other things.  But here we are projecting
		* _away_ from the sun, so we flip things around.
		*/

		if( h < p->height_p[z] )
			proj_area = gortt_crown_proj_cross_section( p, p->theta_p[th], h, (h - zdiff) );
		else
	 		proj_area = gortt_crown_proj_cross_section( p, p->theta_p[th], h, (h + zdiff) );

      S = V_tot / proj_area ;


	}



	return (S);
}


double gortt_get_pcc( gortt_parameters *p, double h )
{



	/*
	* return the probability of finding a crown center at the height 'h'
	*
	* n.b. h is deliberately unused: crown centres are uniformly
	* distributed between h1 and h2, so the density does not depend on
	* height. The argument is kept because the caller is integrating
	* over h and any non-uniform distribution would need it.
	*/
	return (1.0 / (p->h2_p - p->h1_p));
}





double gortt_vol( gortt_parameters *p, int h, int h_s, int t, double h_b )
{

	double   V, V_0, V_sp1, V_sp2, V_cyln;
	double   tmp_s;
	double   h_t, h_tt;

	tmp_s = (p->height_p[h_s] - p->height_p[h]) / cos(p->theta_p[t]) ;
	V_0 =  M_PI * p->rr * tmp_s;
	V_0 += (4.0/3.0) * M_PI * p->rrr;

	if( (p->height_p[h] - p->r) >= h_b ){

		V = 0.0;

	}else if( (p->height_p[h] - p->r*sin(p->theta_p[t])) >= h_b ) {

		h_t = p->r - ( p->height_p[h] - h_b );
		V = (M_PI/3.0) * h_t * h_t * (3.0*p->r - h_t);



	}else if( (p->height_p[h] + p->r*sin(p->theta_p[t])) >= h_b ) {
		V_sp1 = (2.0/3.0)*M_PI*p->rrr ;
		V_sp1 -= gortt_trisec(p->height_p[h], h_b, p->theta_p[t], p->r);
		h_tt = (h_b-(p->height_p[h]-p->r*sin(p->theta_p[t])))/cos(p->theta_p[t]);



		if( p->height_p[h_s]-p->r*sin(p->theta_p[t]) >= h_b ) {

			double hh1, hh2, hh;

			hh1 = (p->height_p[h]-h_b)/sin(p->theta_p[t]);
			hh2 = p->r;
			hh = h_tt;
			V_cyln = gortt_cylind(p->r, hh1, hh2, hh);
			V_sp2 = 0.0;



	  }else{

			double hh1, hh2, hh;

			V_sp2 = 0.0;
			hh1 = (p->height_p[h]-h_b)/sin(p->theta_p[t]);
			hh2 =(p->height_p[h_s]-h_b)/sin(p->theta_p[t]);
			hh = (p->height_p[h_s]-p->height_p[h])/cos(p->theta_p[t]) ;
			V_cyln = gortt_cylind(p->r, hh1, hh2, hh);
			V_sp2 = gortt_trisec( h_b,p->height_p[h_s], p->theta_p[t], p->r);


		}

	V = V_sp1 + V_cyln + V_sp2;

	}else if( p->height_p[h_s] - p->r*sin(p->theta_p[t]) >= h_b ) {

		double tmp_h;

		tmp_h = (h_b-p->height_p[h])/cos(p->theta_p[t]);
		V_cyln = M_PI*p->r*p->r*tmp_h;
		V_sp1 = (2.0/3.0)*M_PI*p->rrr ;
		V = V_sp1 + V_cyln;



	}else if( p->height_p[h_s] + p->r*sin(p->theta_p[t]) >= h_b ) {

		double hh1, hh2, hh;
		double tmp_h;

		h_tt = (p->height_p[h_s]+p->r*sin(p->theta_p[t]) - h_b)/cos(p->theta_p[t]);
		hh1 = (h_b-p->height_p[h_s])/sin(p->theta_p[t]);
		hh2 = p->r;
		hh = h_tt;
		tmp_h = (p->height_p[h_s]-p->height_p[h])/cos(p->theta_p[t]);
		V_cyln = M_PI*p->r*p->r*tmp_h - gortt_cylind(p->r, hh1, hh2, hh);
		V_sp2 = gortt_trisec( h_b, p->height_p[h_s], p->theta_p[t], p->r);
		V_sp1 = (2.0/3.0)*M_PI*p->rrr ;
		V = V_cyln + V_sp2 + V_sp1;


	}else if( p->height_p[h_s] + p->r  >= h_b ){

		h_t = p->r - (h_b - p->height_p[h_s]);
		V_sp1 = (M_PI/3.0) * h_t * h_t * ( 3.0*p->r - h_t);
		V = V_0 - V_sp1;


	}else {      /* height[h_s] + g.R < g.h2_p */

		V = V_0;


	}

  return (V);

}


static double gortt_trisec( double hh, double hh_b, double th, double r )
{
    double x, tmp;
    double b, a1, a2;
/*    int noint = 20; */
    int noint = 20;

    tmp = (hh - hh_b);

/*    printf("\n trisec: %6.2f", r*r-tmp*tmp); */
    x   = -1.0 * tmp*sin(th) + sqrt(r*r- tmp*tmp)*cos(th);

    b = - tmp/sin(th);    /* input parm for triangle volume calculation  */
/*    h_0 = tmp*cos(th) + sqrt(r*r-tmp*tmp)*sin(th);
      b   = x - h_0/tan(th);   alternative input parm for triangle volume  */
    a1 = x ;              /* input parm for sector   volume calculation  */
    a2 = r;               /* input parm for sector   volume calculation  */

    return( gortt_triang(b,r,th,noint)+gortt_sector(a1,a2,r));

}


/*     a1  and  a2  lie in the range (-r,r)  */
static double gortt_sector( double a1, double a2, double r )
{
 double b1, b2,volume;

      b1=r*r*a1-(a1*a1*a1)/3.0;
      b2=r*r*a2-(a2*a2*a2)/3.0;
      volume=M_PI*(b2-b1)/2.0;

      return(volume);
}

/*    b  lies in the range (-r,r)
      theta  in (0,pi/2)
 */
static double gortt_triang( double b, double r, double the, int noint )
{

	double sint, cost, a1, x0, h, volume;
	int m, i;
	double sum1, sum2;

	sint=sin(the);
	cost=cos(the);

	a1 = r*r - b*b*sint*sint;
	x0=b*(sint*sint) + sqrt(a1)*cost;

	/*composite Simpson's rule */
	m=noint;
	h=.50*(x0-b)/(float) m;

	sum1=0.0;
	for(i=0; i<m; i++)
		sum1+=gortt_triang_fcn(b+(float)(2*i+1)*h, b, r,the);

	volume=4.0*sum1;

	sum2=0.0;
	for(i=0; i<m-1; i++)
		sum2 += gortt_triang_fcn(b+(float)(2*(i+1))*h, b, r, the);

	volume += 2.0*sum2;

	/*
	note  that   fcn(x0,b,r,theta)=0  by definition of x0
	and   fcn(b,b,r,theta)=(r**2-x**2)*pi/2
	volume += (r*r - b*b)*PI*.50;
	*/

	volume +=  gortt_triang_fcn(x0,b,r,the);
	volume +=  gortt_triang_fcn(b,b,r,the);
	volume *=h/3.0;

/*     end of composite Simpson's rule */
	return(volume);
}



static double gortt_triang_fcn( double x, double b, double r, double the )
{

	double a1,a2,a3,func;

	a1=tan(the)*(x-b);
	a2=r*r-x*x ;
	a3=a2-a1*a1;
	if(fabs(a3)<0.0000000001) a3 = 0.0;
	func = 2.0* a1*sqrt(a3);

	return(func);

}



/*     integrand of  sqrt(r**2-x**2) */
static double gortt_cylind_fcn( double x, double r )
{

	double func;
	func =.50*x*sqrt(r*r-x*x)+.50*r*r*asin(x/r);

	return(func);

}




/*
      plane cuts both ends of the cylinder
      h1        is the bottom x-intercept, lies in (-r,r)
      h2        is the bottom x-intercept, lies in (-r,r)
      h2 > h1
      height    is the height of the cylinder
 */
static double gortt_cylind( double r, double h1, double h2, double h )
{

	double slope, volume;
	double tmp1, tmp2;

	slope=h/(h2-h1);
	tmp1 = sqrt(r*r-h1*h1);
	tmp2 = sqrt(r*r-h2*h2);

	volume= tmp1*tmp1*tmp1 - tmp2*tmp2*tmp2;

	volume /= 3.0;
	volume -= h1*(gortt_cylind_fcn(h2,r)-gortt_cylind_fcn(h1,r));
	volume *= 2.0*slope;

	if(h2 < r) {
		double phi, s1,s2;
		phi = acos(h2/r);
		s1 = r*r*phi;
		s2 = r*sin(phi)*h2;
		volume += (s1-s2)*h;
	}

	return(volume);
}





/*
* Fill the EPgap array.  each element of the table holds the expected gap
* probability for a given height and angle.
*/
void gortt_calc_epgap( gortt_parameters *p )
{
	int h;
	int t;
	int s;
	int maxs_p;
	double sum_Ps;
	double m_s;


	/*
	* n.b. t < nth, not nth-1. gortt_calc_kopen integrates this table
	* over t = 1..nth-1 inclusive, so leaving the last angle unfilled
	* built the final trapezium panel of k_openep out of a zero.
	*/
	for (h = 0; h < p->nh_probs; h++) {
		for (t = 0; t < p->nth; t++) {
			p->epgap[h][t] = 0.0;
			sum_Ps = 0.0;
			m_s = 0.0;

			/*
			* Integrate over path lengths.
			*/

			maxs_p = gortt_s_to_index( p, (p->z2_p - p->height_p[h]) / cos(p->theta_p[t]));

			for (s = 0; s <= maxs_p;  s++) {
				m_s += gortt_index_to_s( p, s ) * p->pd_s[h][t][s] ;
				sum_Ps += p->pd_s[h][t][s] ;
				p->epgap[h][t] += gortt_calc_pgap( p, gortt_index_to_s( p, s ), p->theta[t]) * p->pd_s[h][t][s] ;
			}

			/*
			if(h==0)   printf("\n %4.2f %8.4f %8.4f", RTOD(p->theta[t]),p->p_n0[h][t],p->epgap[h][t]);
			*/
		}

	}

	return ;

}



/*
* Get a gap probability given a pathlength and a zenith angle. At the
* moment, the zenith angle is not used because we are assuming a
* spherical LAD.
*/
double gortt_calc_pgap( gortt_parameters *p, double s, double th )
{
   return (exp(-s * p->tau_p));

}



/**
This is Lewis's functional approximation to the probability terms
as used in Quaife et al. (2008)
It is ONLY coded up for h=0 - use with caution.
**/
void gortt_gap_probabilities_Q08( gortt_parameters *p )
{


	int t; /*theta index*/
	double c; /*cover*/
	double l; /*tree lai*/
	double k1, k2, a ;

	double tmp1, tmp1_last;
	double tmp2, tmp2_last;


	c = M_PI * p->rr * p->lambda ;

	l = p->favd * p->b * 4./3. * c ;

	k2 = 0.348535 * pow( c , ( -1.08069 - 0.0874595 * c ) );
    k1 = 0.0014166;

	a = c * ( exp( k1*c*c ) - exp( -k2*l ) );

	p->k_open[0] = 0.0 ;
	p->k_openep[0] = 0.0;

	/*
	n.b. the trapezium is seeded AFTER the t=0 probabilities are
	computed. It used to be seeded before, so tmp1_last and tmp2_last
	were read out of the still-uninitialised arrays. That happened to
	be harmless -- theta[0] is exactly zero, so both were multiplied by
	sin(0) -- but it was harmless by accident, not by construction.
	*/

	p->p_n0[0][0]  = exp( -c / ( cos( p->theta_p[0] ) ) );
	p->epgap[0][0] = exp( -a / ( cos( p->theta_p[0] ) ) ) - p->p_n0[0][0] ;

	tmp1_last = p->p_n0[0][0] * sin(2.0*p->theta[0]);
	tmp2_last = p->epgap[0][0] * sin(2.0*p->theta[0]);


	for (t = 1; t < p->nth; t++) {

		p->p_n0[0][t]  = exp( -c / ( cos( p->theta_p[t] ) ) );
		p->epgap[0][t] = exp( -a / ( cos( p->theta_p[t] ) ) ) - p->p_n0[0][t] ;

		tmp1 = p->p_n0[0][t] * sin(2.0*p->theta[t]);
		p->k_open[0] += (tmp1 + tmp1_last) / 2.0 * p->dth;
		tmp1_last = tmp1;

		tmp2 = p->epgap[0][t] * sin(2.0*p->theta[t]);
		p->k_openep[0] += (tmp2 + tmp2_last) / 2.0 * p->dth;
		tmp2_last = tmp2;


	}

	return ;
}





