#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gortt.h>


void gortt_gap_probabilities(  p, g  )
/**
This is Lewis's functional approximation to the probability terms
it is ONLY coded up for h=0!
**/
gortt_parameters *p ;
gortt_geometry *g ;
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

	tmp1_last = p->p_n0[0][0] * sin(2.0*p->theta[0]);
	tmp2_last = p->epgap[0][0] * sin(2.0*p->theta[0]);

	p->p_n0[0][0]  = exp( -c / ( cos( p->theta_p[0] ) ) );
	p->epgap[0][0] = exp( -a / ( cos( p->theta_p[0] ) ) ) - p->p_n0[0][0] ;


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



/*
These are required by some of the other bits of the code
they are not used by Lewis's function
*/

int gortt_s_to_index( p, s )
gortt_parameters *p ;
double s;
{
	return ( (int) ( s / p->ds + 0.5 ) );
}

double gortt_index_to_s( p, index )
gortt_parameters *p ;
int index;
{
	return ( (double) index * p->ds );
}



