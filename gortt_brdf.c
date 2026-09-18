#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gortt.h>


/*
n.b. the relative azimuth is passed in rather than taken from g: gortt_kc
needs Kg on the principal plane (raa = 0 and 180 deg) as well as at the
actual geometry, and it used to get them by writing to g->raa and putting
it back afterwards.
*/
double gortt_kg( const gortt_parameters *p, const gortt_geometry *g, double raa )
{

	double Kg, overlap ;

	overlap = gortt_overlap( p, g, raa );

	Kg = exp( -( p->lambda*pow( p->r, 2 )*M_PI*( SEC(g->sza_prime)+SEC(g->vza_prime)-overlap ) ) );

	return( Kg );

}


/*
O-bar: the overlap between the illumination and view shadows of a single
crown, in units of pi r^2, so that Kg = exp(-lambda pi r^2 (sec i' + sec v' - O)).
All angles are primed.

For a crown centred at height h (Li & Strahler 1992 eqs 5-6):

	O(h)   = ( t - sin(t) cos(t) ) secs / pi
	cos(t) = (h/b) sep / secs
	secs   = sec i' + sec v'
	sep    = sqrt( D^2 + (tan i' tan v' sin raa)^2 )
	D^2    = tan^2 i' + tan^2 v' - 2 tan i' tan v' cos raa

D is the separation of the two shadow centres per unit height, as in the
linear kernels, and sep is the effective separation the overlap depends on.
Its sin(raa) term (the Ambrals/LiSparse form) matters only off the principal
plane; in the principal plane sep = D = |tan i' - tan v' cos raa|, which is
L&S eq. 6.  O(h) falls to zero when cos(t) reaches 1, where the shadows stop
touching, and at the hotspot (sep = 0) it is secs/2 whatever the height.

L&S eq. 3 wants the overlap averaged over the centre heights, which gortt
takes as uniform in [h1,h2].  The average is a Gauss-Legendre sum over the
NH_OVERLAP heights p->h_overlap and weights p->w_overlap, set up in
gortt_init_params; only h varies inside the loop.

Previously a single h = 2b + h2 - h1 (the canopy depth) was used, and two
alternative forms of sep and secs sat here behind if(TRUE)/if(FALSE).  All
three were scored against a Monte Carlo and are worse: KC-INVESTIGATION.md,
section 4, and figures/goms_kernels.py, where they survive as switches.
*/
double gortt_overlap( const gortt_parameters *p, const gortt_geometry *g, double raa )
{

	double overlap, t, D, d ;
	double sep, secs, cos_t, ratio ;
	int i ;


	/* As in the linear kernels*/
	d = pow( tan( g->sza_prime ),2 ) + pow( tan( g->vza_prime ),2 ) \
		- 2.0 * tan( g->sza_prime ) * tan( g->vza_prime ) * cos( raa ) ;

	D = sqrt( MAX( 0.0, d ) );

	sep  = sqrt( D*D + pow( ( tan(g->sza_prime)*tan(g->vza_prime)*sin(raa) ), 2) );
	secs = SEC(g->sza_prime) + SEC(g->vza_prime) ;

	/*cos(t) = h * ratio: everything but the height, hoisted out of the loop*/
	ratio = sep / ( secs * p->b ) ;

	overlap = 0.0 ;
	for( i=0; i<NH_OVERLAP; i++ ){

		cos_t = p->h_overlap[ i ] * ratio ;

		cos_t = MAX( -1.0, cos_t );
		cos_t = MIN(  1.0, cos_t );

		t = acos( cos_t );

		overlap += p->w_overlap[ i ] * MAX( 0.0, ( t - sin( t )*cos_t )*secs/M_PI ) ;
	}

	return( overlap );

}

double gortt_kc_ambrals( const gortt_parameters *p, const gortt_geometry *g )
{

	double Kc, cos_xi_p ;

	cos_xi_p = cos( g->sza_prime )*cos( g->vza_prime )+sin( g->sza_prime )*sin( g->vza_prime )*cos( g->raa );

	Kc = ( 1.0 - exp( -p->lambda*M_PI*p->rr*SEC( g->vza_prime ) ) ) * 0.5 * ( 1.0 + cos_xi_p );

	return( Kc );

}


/*
This version interpolates between values of
f calulated on the principal plane to produce
the actual Kc value. Based on the method used
by Wenge in her code, which is Li & Strahler
(1992) IEEE TGRS 30(2), sec. H, eq. (25).
*/
double gortt_kc( const gortt_parameters *p, const gortt_geometry *g, double Kg )
{

	double Kc ;
	double beta, junk ;
	double f, f0, f180 ;
	double F, F0, F180 ;
	double Kg0, Kg180 ;
	double frac ;

	/*at the actual geometry: only F and beta are used, f is discarded*/
	gortt_kc_fFbeta( p, g, g->raa, Kg, &f, &F, &beta );

	/*and on the two halves of the principal plane*/
	Kg0 = gortt_kg( p, g, DTOR( 0. ) ) ;
	gortt_kc_fFbeta( p, g, DTOR( 0. ), Kg0, &f0, &F0, &junk );

	Kg180 = gortt_kg( p, g, DTOR( 180. ) ) ;
	gortt_kc_fFbeta( p, g, DTOR( 180. ), Kg180, &f180, &F180, &junk );

	frac = g->raa / M_PI ;
	if( frac > 1.0 ) frac = 2.0 - frac ;

	/*slot in user defined beta*/
	if( p->use_user_beta == TRUE ) beta = p->beta ;

	/*
	Li & Strahler (1992) sec. H, eq. (25): f is computed at raa = 0, pi
	and interpolated to the actual azimuth. gortt_kc_fFbeta already
	returns f with its factor of F in it (eq. 17), so it is not applied
	again here. This used to read f0*F0 and f180*F180, i.e. F squared.
	*/
	f = (1.-frac)*f0 + frac*f180 ;

	f = beta*f+(1.0-beta)*F ;

	Kc = f * ( 1.0 - Kg );

	return( Kc ) ;
}

/*
Only calculates f, F and beta of the kc function
this is done so we can interpolate bewteen various
f's to get the true value of Kc for OPP cases.
*/
void gortt_kc_fFbeta( const gortt_parameters *p, const gortt_geometry *g, double raa, double Kg, double *f, double *F, double *beta )
{

	double phase_prime, Gamma, Gamma_v, Gamma_c, Gamma_i ;
	double overlap, Mv, Mi, M ;
	double PvMv, PiMi, Po ;
	double theta_Mi, D ;

	/*double f, F, beta ; */

	overlap = gortt_overlap( p, g, raa );
	phase_prime = cos( g->vza_prime )*cos( g->sza_prime )+sin( g->vza_prime )*sin( g->sza_prime )*cos( raa ) ;

	Mi = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->sza_prime) ) )  ;
	Mv = ( 1.0  - (1.0 - exp( -p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )/( p->lambda * M_PI * p->rr * SEC(g->vza_prime) ) )  ;

	Gamma   = M_PI * p->rr * ( SEC(g->sza_prime)+SEC(g->vza_prime)-overlap ) ;
	Gamma_c = M_PI * p->rr *   SEC(g->vza_prime) * 0.5 * ( 1.0 + phase_prime ) ;
	Gamma_v = M_PI * p->rr *   SEC(g->vza_prime) ;

	*F = Gamma_c / Gamma ;

	M = 1.0 - ( 1.0 - Kg )/( p->lambda * Gamma ) ;

	theta_Mi = acos( 1.0 - 2.0*Mi );

	/*IGARSS 92 paper: PvMv is its eq. 4, PiMi its eq. 5*/

	/*
	Gamma_i for beta: the shadow area a crown casts in the illumination
	direction. IGARSS '92 eq. 7 defines it as the shadowing area on a level,
	and Abuelgasim & Strahler (1994) call lambda*Gamma_i the "coverage and sun
	position". It used to be Gamma_v, which made beta depend on the view
	zenith. KC-INVESTIGATION.md section 8; OPEN-QUESTIONS #9.
	*/
	Gamma_i = M_PI * p->rr * SEC(g->sza_prime) ;

	PiMi = ( 1 - cos( theta_Mi * ( 1 - ( g->sza_prime - g->vza_prime * cos( raa ) )/M_PI ) ) ) / 2.0 ;
	PvMv = Mv -( 1.0 - cos( g->vza_prime * cos( raa ) - g->sza_prime ) ) / 2.0 ;


	/*
	Po: the part of the crown shaded from both directions, following Schaaf,
	Li & Strahler (1994) eq. 23A.  Beyond the hotspot in backscatter
	(vza' > sza') Po = PiMi, so the decrement is PvMv alone; the uniform-height
	Kc is then exactly 1 - Kg - Kz, all visible crown sunlit (L&S 1992 eq. 20).
	Everywhere else the decrement is the larger of PvMv and PiMi, i.e.
	Po = min(PvMv, PiMi): the smaller mutually shadowed set lies inside the
	larger, as Li & Strahler (1992, IGARSS) conclusion 1 and a Monte Carlo of
	the scene both show (KC-INVESTIGATION.md, sections 2 and 9).

	Two rules this replaces, both wrong somewhere:
	- PvMv throughout forward scatter (L&S 1992 p. 282) cancels the view
	  shadowing, and once PvMv overtakes PiMi at large forward zeniths Kc grows
	  like sec(vza') past 1 - Kg - Kz.
	- min(PvMv, PiMi) everywhere (0204222) subtracts PiMi in backscatter when
	  the sun is near overhead and PiMi > PvMv, and Kc goes negative.
	*/
	if( ( g->vza_prime > g->sza_prime ) && !( ( raa < DTOR( 270. ) ) && ( raa > DTOR( 90. ) ) ) )
		Po = PiMi ;
	else
		Po = MIN( PvMv, PiMi ) ;

	/*
	beta: Li & Strahler (1992, IGARSS) eq. 7.  D = r cot(sza'/2) is infinite
	for an overhead sun, so (h2-h1)/D -> 0 and beta -> 1; the formula itself
	would divide by zero there, hence the guard.  It used to return 0, which
	made beta jump from 0.99 at sza 0.1 deg to 0 at sza 0 (KC-INVESTIGATION.md
	section 9).
	*/
	if( g->sza_prime < 0.000000001 ){

		*beta = 1.0 ;

	}else{

		D = p->r * ( 1.0 / tan ( g->sza_prime/2.0 ) );
		*beta = (p->lambda*Gamma_i)/(p->lambda*Gamma_i+(p->h2-p->h1)/D)*(1.0-exp(-p->lambda*Gamma_i-(p->h2-p->h1)/D))/(1.0-exp(-p->lambda*Gamma_i)) ;

	}

	*f = *F * ( 1.0 - Gamma_v*( PvMv + PiMi - Po )/Gamma_c ) / ( 1.0 - M ) ;


	return ;
}


/*
Kuusk's bidirectional within-crown gap probability, as used in Ni et al.
(1999) eqs 3-7:

	P(n>0 | h=0, i, v) = P(n>0 | h=0, i) P(n>0 | h=0, v) H

	H         = exp( tau(i) * path_mean * hotspot )
	path_mean = sqrt( lsza * lvza )
	hotspot   = ( 1 - exp(-path_sep/r) ) / ( path_sep/r )
	path_sep  = sqrt( lsza^2 + lvza^2 - 2 lsza lvza cos(xi) )

tau(i) = k(i) FAVD.  lsza and lvza are the within-crown path lengths along the
sun and view directions, recovered by inverting Beer's law on the one-way
within-crown gap probabilities: l = -ln(P) / (k FAVD) (eqs 5-6).  path_sep is
the distance between the ends of the two paths (eq. 7), and the crown radius r
plays the role of Kuusk's correlation length (k_l = 1/r).

At the hotspot path_sep = 0 and hotspot = 1, so H = exp(tau path_mean), which
undoes the double counting of a shared path -- exactly, while k is the same for
the sun and view directions (true as long as the leaf angle distribution is
hardwired; OPEN-QUESTIONS #4-5).  As the directions separate, hotspot falls
towards 0 and H towards 1: the two gaps become independent.

n.b. cos(xi) uses the unprimed angles; a primed version is left commented out.
*/
double gortt_kuusk( const gortt_parameters *p, const gortt_geometry *g )
{

	double kuusk, lsza, lvza, path_sep, cos_xi, H ;
	double path_mean, hotspot ;

	cos_xi = cos( g->sza )*cos( g->vza )+sin( g->sza )*sin( g->vza )*cos( g->raa );
	/*05/06/06*/
	/*cos_xi = cos( g->sza_prime )*cos( g->vza_prime )+sin( g->sza_prime )*sin( g->vza_prime )*cos( g->raa );*/

	lsza = -log( g->p_ngt0_heq0_sza ) / ( p->k * p->favd ) ;
	lvza = -log( g->p_ngt0_heq0_vza ) / ( g->k_vza * p->favd ) ;


	/*sqrt trap for -ve; zero separation is the hotspot, where the limit is 1*/
	if( ( lsza*lsza + lvza*lvza - 2.*lsza*lvza*cos_xi ) > 0.0 ){
		path_sep = sqrt( lsza*lsza + lvza*lvza - 2.*lsza*lvza*cos_xi );
		hotspot = ( 1.0 - exp( -path_sep/p->r ) ) / ( path_sep/p->r ) ;
	}else{
		hotspot = 1.0 ;
	}

	if( ( lsza * lvza ) > 0.0 ){
		path_mean = sqrt( lsza * lvza ) ;
	}else{
		path_mean = 0.0 ;
	}

	H = exp( p->k * p->favd * path_mean * hotspot ) ;

	kuusk = g->p_ngt0_heq0_sza * g->p_ngt0_heq0_vza * H ;

	/*fprintf( stderr, "### %f\n", kuusk );*/

	return( kuusk );

}

/*
Fill in the radiative transfer terms for one geometry and one waveband.

These used to be a call tree of small functions, each recomputing its
dependencies from scratch every time it was asked for: one pass of the
waveband loop in gortt_rsurf evaluated the semi-infinite transmittance
T_inf_ff thirty-nine times for a value that does not change within it.
They are filled here instead, once, in dependency order -- which is the
order Ni et al. (1999) derive them in, so the block still reads against
the paper.
*/
/*
The two-stream path scattering parameters of a horizontally homogeneous
layer: Ni et al. (1999) eqs (8)-(11) and the four semi-infinite quantities
underneath them, for a layer of one-sided leaf area index L viewed at
incident zenith theta.

This used to be written out inline in gortt_rt_init at the single depth
k*elai. It is a function so that -rt_crown_path can ask for a second depth
-- the paper's eqs (15)-(16) set the crown path reflectances equal to the
homogeneous ones, and which homogeneous layer that means is exactly the
question those options exist to test. Called with L = p->elai and
theta = sza_prime it reproduces the published model bit for bit, which is
what the golden files pin.

The expressions are kept in the exact form they had inline -- SEC(theta)
rather than 1/mu, the pow() calls, the order of the products -- because at
-O0 there is no reassociation and the goldens are a byte comparison.
*/
static void gortt_two_stream( double w, double gamma, double k, double L, double theta,
                              int upscatter, double beta_up, gortt_slab *sl )
{

	double a, b, lambda, sbeam, det, P, Q ;

	sl->t_0 = exp( -( k*L*SEC( theta ) ) ) ;

	if( upscatter == TRUE ){

		/*
		The same four quantities with a general upscatter parameter.

		gamma is not passed through here because gamma IS the beta = 1/2
		assumption: the closure underneath eqs (8)-(11) is

			dF_down/dl = -a F_down + b F_up + (1-beta) w s exp(-s l)
			dF_up/dl   = -b F_down + a F_up -    beta  w s exp(-s l)

		with a = 2k[1 - w(1-beta)] and b = 2k w beta, whose eigenvalue
		is 2k sqrt( [1-w(1-beta)]^2 - (w beta)^2 ). At beta = 1/2 that
		collapses to 2k sqrt(1-w) = 2k gamma, and only then. See
		latex_notes/two_stream_derivation.tex section 3.2.

		The beam upscatter is taken equal to the diffuse one. That is a
		choice, not a result: the proper beam parameter depends on the
		incidence direction (Sellers 1985). This is the one-parameter
		version of the generalisation, and it is the same "half goes
		each way" assumption being relaxed in both places at once.
		*/

		a      = 2.0*k*( 1.0 - w*( 1.0 - beta_up ) ) ;
		b      = 2.0*k*w*beta_up ;
		lambda = 2.0*k*sqrt( ( 1.0 - w*( 1.0 - beta_up ) )*( 1.0 - w*( 1.0 - beta_up ) )
		                     - ( w*beta_up )*( w*beta_up ) ) ;
		sbeam  = k*SEC( theta ) ;
		det    = lambda*lambda - sbeam*sbeam ;

		/*
		b/(a+lambda), not the (a-lambda)/b the eigenvector gives
		directly: they are equal because lambda^2 = a^2 - b^2, and this
		form stays finite for a black leaf where b is zero and a-lambda
		is zero with it.
		*/
		sl->R_inf_ff = b / ( a + lambda ) ;
		sl->T_inf_ff = exp( -lambda*L ) ;

		if( det == 0.0 ){

			/*
			The beam decays at exactly the rate of the diffuse
			eigenmode, so the particular solution below is 0/0.
			This is the same removable singularity as the one in
			the published branch -- of which this branch is the
			beta_up != 1/2 generalisation -- and the limits are

			  R_inf_df -> w[ beta + (1-beta) R_inf_ff ] / 2
			  T_inf_df -> w L t_0 [ (1-beta)(a+s) + b beta ] / 2

			the first by l'Hopital in sbeam, the numerator of
			Q - P R_inf_ff vanishing on the same line because
			lambda^2 = a^2 - b^2; the second from
			T_inf_ff - t_0 -> -L (lambda - sbeam) t_0. At
			beta_up = 1/2 both reduce to the published limits.
			*/
			sl->R_inf_df = w*( beta_up + ( 1.0 - beta_up )*sl->R_inf_ff ) / 2.0 ;
			sl->T_inf_df = w*L*sl->t_0*( ( 1.0 - beta_up )*( a + sbeam ) + b*beta_up ) / 2.0 ;

		}else{

			/*the particular solution for the beam source*/
			P = w*sbeam*( ( 1.0 - beta_up )*( a + sbeam ) + b*beta_up ) / det ;
			Q = w*sbeam*( beta_up*( a - sbeam ) + b*( 1.0 - beta_up ) ) / det ;

			/*
			top-of-medium upward flux, and the internal downward
			diffuse flux at depth L, of the semi-infinite medium
			-- which is what these two are, notwithstanding their
			names
			*/
			sl->R_inf_df = Q - P*sl->R_inf_ff ;
			sl->T_inf_df = -P*( sl->T_inf_ff - sl->t_0 ) ;

		}

	}else{

	/*
	the closure's a, at the beta = 1/2 that gamma encodes:
	a = 2k[1 - w(1-beta)] = k(2-w). Not used by the expressions below --
	they are written in gamma -- but the conservative-scattering limit in
	the adding section is common to both branches and is written in a.
	*/
	a = k*( 2.0 - w ) ;

	/*hemispherical-hemispherical reflectance*/
	sl->R_inf_ff = ( 1.0 - gamma )/( 1.0 + gamma ) ;

	/*directional-hemispherical reflectance*/
	sl->R_inf_df = ( 1.0 - gamma )/( 1.0 + 2.0 *cos( theta )*gamma ) ;

	/*hemispherical-hemispherical transmittance*/
	sl->T_inf_ff = exp(-( 2.0*gamma*k*L )) ;

	/*directional-hemispherical transmittance*/
	if( ( 1. - pow( ( 2.*gamma*cos( theta ) ), 2 ) ) == 0.0 ){

		/*
		2 gamma mu == 1. The diffuse eigenvalue 2 gamma k and the beam
		extinction k sec(theta) are then the same number, so T_inf_ff
		and t_0 are equal and the quotient is 0/0. It is removable:

		  ( T_inf_ff - t_0 ) / ( 1 - (2 gamma mu)^2 ) -> k L sec(theta) t_0 / 2

		by expanding exp(-x u) - exp(-u) about x = 1 with u = k L sec(theta).
		Only the exact hit needs this -- a part in 1e11 off the line the
		quotient is already right to nine figures -- but omega = 0.75 with
		an overhead sun lands on it exactly.
		*/
		sl->T_inf_df  = ( w / 2.0 ) ;
		sl->T_inf_df *= ( 1. + 2. * cos( theta ) ) ;
		sl->T_inf_df *= k*L*SEC( theta )*sl->t_0 / 2.0 ;

	}else{

		sl->T_inf_df  = ( w / 2.0 ) ;
		sl->T_inf_df *= ( 1. + 2. * cos( theta ) ) / ( 1. - pow( ( 2.*gamma*cos( theta ) ), 2 ) )  ;
		sl->T_inf_df *= ( sl->T_inf_ff - sl->t_0 ) ;

	}

	}

	if( ( 1. - pow( sl->R_inf_ff*sl->T_inf_ff, 2 ) ) == 0.0 ){

		/*
		Conservative scattering, omega == 1. gamma is zero, and with it
		the eigenvalue lambda = 2 k gamma, so R_inf_ff and T_inf_ff are
		both exactly 1 and the adding denominator vanishes in both of
		the expressions below. The limits are the classical conservative
		two-stream pair, reached by expanding R_inf_ff -> 1 - lambda/a
		and T_inf_ff -> 1 - lambda L in the eigenvalue:

		  t_ff -> 1/(1 + a L),    p_ff -> a L/(1 + a L)

		i.e. 1/(1+kL) and kL/(1+kL) in the published branch, where
		a = k(2-w) is k at omega = 1. They satisfy t_ff + p_ff = 1, as
		a non-absorbing layer must, and with the T_inf_df and R_inf_df
		above they give p_df + t_df + t_0 = 1 to the last bit.
		*/
		sl->t_ff = 1.0 / ( 1.0 + a*L ) ;
		sl->p_ff = a*L / ( 1.0 + a*L ) ;

	}else{

		/*hemispherical hemispherical transmittance*/
		sl->t_ff  = sl->T_inf_ff ;
		sl->t_ff *= ( 1. - pow( sl->R_inf_ff, 2 ) ) ;
		sl->t_ff /= ( 1. - pow( sl->R_inf_ff*sl->T_inf_ff, 2 ) );

		/*hemispherical-hemispherical reflectance*/
		sl->p_ff  = sl->R_inf_ff ;
		sl->p_ff *= ( 1. - pow( sl->T_inf_ff, 2 ) ) ;
		sl->p_ff /= ( 1. - pow( sl->T_inf_ff*sl->R_inf_ff, 2 ) )  ;

	}

	/*directional hemispherical transmittance*/
	sl->t_df = sl->T_inf_df - sl->p_ff \
			* ( sl->t_0 * sl->R_inf_df + sl->T_inf_df * sl->R_inf_ff ) ;

	/*direction hemispherical reflectance*/
	sl->p_df = sl->R_inf_df - sl->t_ff * \
		 	(  sl->t_0     * sl->R_inf_df + \
				 sl->T_inf_df * sl->R_inf_ff ) ;

	return ;

}


/*
The recollision-probability rescaling, for -rt_recollision.

Spectral invariants (the "p-theory" of Knyazikhin and of Stenberg): of the
radiation a scattering body intercepts, a fraction omega is scattered, of
which a fraction p_c strikes the body again rather than escaping. Summing
that series, the fraction of intercepted radiation that eventually leaves is

	omega_c = omega ( 1 - p_c ) / ( 1 - omega p_c )

and the fraction of INCIDENT radiation that leaves is i0_c omega_c with i0_c
the interceptance. Both are angular-distribution-free, which is the point:
they say how much a crown scatters without saying where it goes.

p_c comes from GORT's own geometry rather than from a fitted formula. For a
convex body of mean chord l and extinction sigma, the first-flight escape
probability of a uniformly distributed isotropic source is 1/(1 + sigma l)
-- the Wigner rational approximation, exact in both limits -- so

	p_c = sigma l / ( 1 + sigma l ),	sigma l = k * p->crown_lai

since crown_lai is favd times the mean chord 4b/3. p_c -> 0 for an optically
thin crown and -> 1 for a thick one, as it must.

The two-stream is kept for the one thing p-theory cannot give, the split
between the upward and downward hemispheres: the slab's total escaping
scattered flux p_df + t_df is rescaled to i0_c omega_c and its up/down ratio
is left alone. For a closed canopy the factor is close to one, so the option
is near-inert where Ni et al. (1999) was validated -- worth checking whenever
this is touched.

n.b. the SCENE-level estimator p = 1 - (1 - K_open)/elai was tried and is not
used: it lowers the effective albedo and makes the sparse-canopy near
infrared worse, because a scene-level p describes absorption by the stand,
not the brightness of a crown.
*/
static double gortt_recollision_scale( const gortt_parameters *p, double w, double k,
                                       const gortt_slab *sl )
{

	double sigma_l, p_c, w_c, i0_c, escaping ;


	escaping = sl->p_df + sl->t_df ;

	/*black leaves scatter nothing and there is nothing to rescale*/
	if( w <= 0.0 || escaping <= 0.0 ) return( 1.0 );

	sigma_l = k * p->crown_lai ;
	p_c     = sigma_l / ( 1.0 + sigma_l ) ;
	w_c     = w * ( 1.0 - p_c ) / ( 1.0 - w * p_c ) ;
	i0_c    = 1.0 - sl->t_0 ;

	return( i0_c * w_c / escaping );

}


/*
The Henyey-Greenstein phase function, normalised so that its average over the
whole sphere is one: (1/4pi) int P dOmega = 1.

g is the asymmetry parameter, <cos Theta>, and it is rt->g_asym, i.e. Li &
Strahler's -(4/9)(rL - tL)/omega. That is negative when the leaf reflects more
than it transmits, which is the usual case and which makes the scattering
backward-peaked -- the right sign, and the reason no new user input is needed.
|g| <= 4/9 by construction, so the denominator cannot approach zero, and a
black leaf gives g = 0 and P = 1 rather than a NaN (see the guard on g_asym in
gortt_rt_init).
*/
static double gortt_hg_phase( double g, double cos_theta )
{

	double d ;

	d = 1.0 + g*g - 2.0*g*cos_theta ;

	return( ( 1.0 - g*g ) / ( d * sqrt( d ) ) );

}


/*
The scattering angle between the incident and exit directions, as a cosine.

Theta is measured between the two directions of *propagation*, so exact
backscatter -- the sensor looking down the sun's own beam -- is Theta = 180
degrees and cos Theta = -1, which is where a backward-peaked phase function
peaks. That is why this is the negative of the gortt_kuusk phase angle.

Everything here is in the primed frame, where the crown is a sphere and the
chord slab is horizontal, which is the frame gortt_two_stream is already
evaluated in. n.b. gortt_kuusk forms its own cos_xi from the *unprimed*
zeniths (the 05/06/06 comments); the two are inconsistent, and deliberately
left that way -- OPEN-QUESTIONS.md has the prime-transform question, and
CROWN-BRDF-NEXT.md section 6 says not to bundle it with this.
*/
static double gortt_cos_scatter( double mu_i, double mu_v, double raa )
{

	double si, sv ;

	si = sqrt( MAX( 0.0, 1.0 - mu_i*mu_i ) ) ;
	sv = sqrt( MAX( 0.0, 1.0 - mu_v*mu_v ) ) ;

	return( -( mu_i*mu_v + si*sv*cos( raa ) ) );

}


/*
Single-scattering reflectance factor of a horizontal slab of one-sided leaf
area index L, for one incident and one exit direction: the exact first-order
term, with no ground underneath.

	R_ss = omega P(Theta) / ( 4 ( mu_i + mu_v ) )
	       * [ 1 - exp( -k L ( 1/mu_i + 1/mu_v ) ) ]

Both limits are the point of using it. Optically thin, the bracket goes to
k L (1/mu_i + 1/mu_v) and R_ss becomes proportional to L, so a crown's
brightness tracks the chord and the profile across the projected disc is
sqrt(1 - u^2) -- the Lambert-sphere shape. Optically thick, the bracket goes
to one and R_ss becomes omega P / (4(mu_i + mu_v)), which is Lommel-Seeliger:
flat across the disc, and strongly view-dependent.

Both are what FLIGHT shows a single crown doing. Measured over a black scene
(flight_compare/run_crown_image.py), the ratio of the outermost to the
innermost ring of the projected disc at zero phase runs from 0.24 at a crown
optical radius of 0.23 to 0.94 at 15, against the 0.22 of a Lambert sphere and
the 1.00 of the uniform disc the published model assumes.
*/
static double gortt_r_ss( double w, double g_asym, double k, double L,
                          double mu_i, double mu_v, double raa )
{

	double P ;

	P = gortt_hg_phase( g_asym, gortt_cos_scatter( mu_i, mu_v, raa ) ) ;

	return( w * P / ( 4.0 * ( mu_i + mu_v ) )
	        * ( 1.0 - exp( -k*L*( 1.0/mu_i + 1.0/mu_v ) ) ) );

}


/*
The directional-hemispherical reflectance that gortt_r_ss implies, i.e.

	rho_ss = (1/pi) int R_ss mu_v dOmega_v

over the upward hemisphere. This is not a physical quantity anybody wants; it
is what has to be subtracted from rho_df so that the isotropic remainder
carries exactly the rest, and no energy is created by giving the crown a
direction.

The zenith half reuses the Gauss-Legendre rule in p->xq_gl on [-1,1] mapped to
mu_v on [0,1], which is the same rule and the same -nquad scaling as the chord
average above. The azimuth is a uniform midpoint rule, and is needed at all
only because P(Theta) depends on it.
*/
static double gortt_rho_ss( const gortt_parameters *p, double w, double g_asym,
                            double k, double L, double mu_i )
{

	int i, j, nphi ;

	double mu_v, wgt, raa, dphi, sum ;


	nphi = NQ_BRDF_PHI * p->nquad ;
	dphi = 2.0*M_PI / nphi ;
	sum  = 0.0 ;

	for( i=0; i<p->nq_gl; i++ ){

		mu_v = 0.5 * ( p->xq_gl[ i ] + 1.0 ) ;
		wgt  = 0.5 * p->wq_gl[ i ] * mu_v ;

		/*mu_v = 0 is a grazing exit direction: R_ss is finite there but
		the weight mu_v is zero, so nothing is lost by skipping it and a
		1/mu_v is kept away from zero*/
		if( mu_v <= 0.0 ) continue ;

		for( j=0; j<nphi; j++ ){

			raa = ( j + 0.5 ) * dphi ;
			sum += wgt * dphi
			     * gortt_r_ss( w, g_asym, k, L, mu_i, mu_v, raa ) ;

		}

	}

	/*(1/pi) int ... dOmega, and dOmega = dmu_v dphi*/
	return( sum / M_PI );

}


/*
The crown's bidirectional path reflectance, for -rt_crown_brdf: Option A of
CROWN-BRDF-NEXT.md, single scattering exactly and multiple scattering
isotropically, which is Hapke's IMSA construction.

	R = R_ss + R_ms,	R_ms chosen so that  (1/pi) int R mu_v dOmega = rho_df

rho_df is therefore the *normalisation*, not the answer: the two-stream keeps
doing the thing it is good at, which is how much of the light is scattered,
and the single-scattering term does the thing it is good at, which is where it
goes. Energy is conserved by construction, so switching this on can only move
radiation between directions and never add any.

The band structure flatters it. In the near infrared omega is 0.9 and multiple
scattering dominates the magnitude, but it is also genuinely the closer to
isotropic; in the red omega is 0.12, single scattering dominates, and that half
is treated exactly. The band that needs the anisotropy most is the band where
the approximation is best.

If single scattering alone already over-runs rho_df -- which it can, the
two-stream being an approximation and the two not being constrained to agree --
R_ss is scaled down to rho_df rather than R_ms being clamped at zero. That
keeps the angular shape instead of throwing it away, and keeps the integral
exact in both branches.
*/
static double gortt_crown_brdf( const gortt_parameters *p, const gortt_geometry *g,
                                const gortt_rt *rt, double k, double L,
                                double rho_df )
{

	double mu_i, mu_v, r_ss, rho_ss ;


	mu_i = cos( g->sza_prime ) ;
	mu_v = cos( g->vza_prime ) ;

	r_ss   = gortt_r_ss( rt->w, rt->g_asym, k, L, mu_i, mu_v, g->raa ) ;
	rho_ss = gortt_rho_ss( p, rt->w, rt->g_asym, k, L, mu_i ) ;

	if( rho_ss > rho_df ){

		/*nothing left for the isotropic part; rescale so the
		hemispherical integral is still exactly rho_df*/
		return( ( rho_ss > 0.0 ) ? r_ss * rho_df / rho_ss : rho_df );

	}

	/*an isotropic BRF of value x has directional-hemispherical
	reflectance x, so the remainder goes in unchanged*/
	return( r_ss + ( rho_df - rho_ss ) );

}


/*
The path reflectances of the discontinuous canopy, rho'_df and rho'_ff.

Ni et al. (1999) eqs (15)-(16) set them equal to those of the equivalent
homogeneous layer, i.e. to the layer holding all the foliage spread over the
whole ground area (optical depth k*elai). That is the "squeezing all the
plants together" of its section III-C, done vertically.

It is also the approximation that fails as the canopy opens up. rho_df falls
to zero with elai, so a stand of widely spaced crowns is predicted to have
black crowns -- although nothing about a crown has changed. Compared against
FLIGHT with a fixed crown and a thinning stand, the effective sunlit-crown
reflectance in the near infrared should be flat at about 0.32 and gortt runs
0.06 to 0.37 across covers of 5 to 65 per cent.

-rt_crown_path picks a different layer. -rt_path selects which:

	scene	 p->elai, the published choice; a control, not a change
	squeeze	 p->squeeze_lai = elai/(1-K_open), the same foliage spread
		 over the part of the hemisphere that is not open, which is
		 the reading of section III-C that keeps its own eq. (14)
		 self-consistent
	crown	 the crown's own depth, averaged over the chord distribution
		 of the spheroid (the default)

The chord average, not rho_df at the mean chord: the chord at impact
parameter s is 2b sqrt(1-(s/r)^2) with area weight 2s ds/r^2 over the
projected disc, and rho_df is concave in L, so averaging rho_df over chords
is smaller than rho_df at the mean chord -- 0.398 against 0.417 for the
sparse sphere. It reuses the Gauss-Legendre rule gortt_init_params already
built for the gap probabilities, so it adds no new quadrature and follows
-nquad like everything else.

-rt_crown_brdf then changes what is averaged rather than what it is averaged
over: rho_df(L) becomes a bidirectional reflectance with the same
hemispherical integral. It requires -rt_crown_path for the same reason
-rt_recollision does -- it describes a crown, and there is no chord to
evaluate it on without one.
*/
static void gortt_crown_reflectance( const gortt_parameters *p, const gortt_geometry *g,
                                     const gortt_rt *rt, double k,
                                     double *p_df, double *p_ff )
{

	int i ;

	double u, wgt, chord, f, L ;
	gortt_slab sl ;


	if( p->rt_crown_path == FALSE ){

		/*eqs (15)-(16) as published. -rt_recollision cannot get here:
		the parser refuses it without -rt_crown_path, because p_c is a
		property of a crown and rescaling a stand-sized layer by it
		would divide a crown quantity by a scene one.*/
		*p_df = rt->p_df ;
		*p_ff = rt->p_ff ;

		return ;

	}

	if( p->rt_crown_path_mode != RT_PATH_CROWN ){

		L = ( p->rt_crown_path_mode == RT_PATH_SQUEEZE ) ? p->squeeze_lai : p->elai ;

		gortt_two_stream( rt->w, rt->gamma, k, L, g->sza_prime,
		                  p->rt_upscatter, rt->beta_up, &sl ) ;

		f = ( p->rt_recollision == TRUE )
		    ? gortt_recollision_scale( p, rt->w, k, &sl ) : 1.0 ;

		*p_df = f * sl.p_df ;
		*p_ff = f * sl.p_ff ;

		return ;

	}

	/*
	chord average. u = s/r runs over [0,1] with weight 2u du, which
	integrates to one; the Gauss-Legendre rule is on [-1,1], so
	u = (x+1)/2 and the weight becomes w_i u_i.
	*/
	*p_df = 0.0 ;
	*p_ff = 0.0 ;

	for( i=0; i<p->nq_gl; i++ ){

		u     = 0.5 * ( p->xq_gl[ i ] + 1.0 ) ;
		wgt   = p->wq_gl[ i ] * u ;
		chord = 2.0 * p->b * sqrt( MAX( 0.0, 1.0 - u*u ) ) ;

		gortt_two_stream( rt->w, rt->gamma, k, p->favd * chord, g->sza_prime,
		                  p->rt_upscatter, rt->beta_up, &sl ) ;

		f = ( p->rt_recollision == TRUE )
		    ? gortt_recollision_scale( p, rt->w, k, &sl ) : 1.0 ;

		/*
		-rt_crown_brdf replaces rho_df(L) with a bidirectional
		reflectance of the same hemispherical integral, chord by chord.
		It goes inside this loop rather than outside it because the
		angular shape depends on the depth: a short chord at the limb
		is in the single-scattering limit while the chord through the
		middle need not be. Structurally it is one for one with the
		line it replaces -- same rule, same weights, same -nquad.

		rho'_ff is left alone. It is hemispherical-hemispherical, so
		there is no direction in it to get right, and it is what the
		crown-ground and neighbour-crown series are summed with.
		*/
		*p_df += wgt * f * ( ( p->rt_crown_brdf == TRUE )
		                     ? gortt_crown_brdf( p, g, rt, k,
		                                         p->favd * chord, sl.p_df )
		                     : sl.p_df ) ;
		*p_ff += wgt * f * sl.p_ff ;

	}

	return ;

}


void gortt_rt_init( const gortt_parameters *p, const gortt_geometry *g,
                    const gortt_spectra *s, int iw, gortt_rt *rt )
{

	double k, k_open ;

	gortt_slab slab ;


	/*
	leaf single scattering albedo, and the two quantities
	derived from it that are used throughout
	*/

	rt->w = *( s->rleaf + iw ) + *( s->tleaf + iw ) ;
	rt->gamma = sqrt( 1 - rt->w ) ;

	/*
	The phase function asymmetry, Li & Strahler's g = -(4/9)(rL - tL)/omega.

	A perfectly black leaf makes that 0/0, and the NaN propagates straight
	into CdC and out as the BRF: "gortt -refl_leaf 0 -trans_leaf 0" used to
	print -nan, where omega = 1e-6 prints 0.1521. Nothing is being asked for
	there -- with no scattered radiation the asymmetry of the scattering is
	undefined -- and the term it feeds is multiplied by omega anyway, so
	zero is both finite and right.
	*/
	rt->g_asym = ( rt->w > 0.0 )
	             ? -(4.0/9.0)*( *(s->rleaf+iw) - *(s->tleaf+iw) )/( rt->w )
	             : 0.0 ;

	/*
	Two-stream upscatter parameter, for -rt_upscatter.

	beta = 1/2 + (rL - tL)/(8 omega) is the standard canopy value for a
	spherical leaf angle distribution, which is the only one gortt has
	(LAD_05). It is bounded by [3/8, 5/8] because (rL-tL)/(rL+tL) cannot
	exceed one, so it needs no clamp; a black leaf is the 0/0 and takes
	the 1/2 that the published model assumes everywhere.

	It is identically 1/2 whenever rL = tL, so the option is inert in the
	near infrared of the FLIGHT comparison by construction. The red is
	what it is for: 0.5417 there, an 8 per cent error in the backscatter
	fraction that the published model does not carry.
	*/
	rt->beta_up = ( rt->w > 0.0 )
	              ? 0.5 + ( *(s->rleaf+iw) - *(s->tleaf+iw) ) / ( 8.0 * rt->w )
	              : 0.5 ;

	k = gortt_leaf_angle_distribution( p, g->sza );


	/*
	direct unscattered transmittance.

	n.b. the primed solar zenith, where the paper writes theta_i -- the
	05/06/06 change, which applies to every angle in this function. It is
	a no-op for spherical crowns and a deliberate deviation otherwise;
	OPEN-QUESTIONS.md records it.
	*/

	/*rt->t_0 = exp( -( k*p->favd*SEC(g->sza) ) ) ;*/
	/*05/06/06*/
	rt->t_0 = exp( -( k*p->elai*SEC( g->sza_prime ) ) ) ;


	/*
	direct-unscattered transmittance for a discontinuous
	plant canopy of height h at given geometry
	*/

	rt->t_prime_0 = g->p_neq0_heq0_sza + g->p_ngt0_heq0_sza ;


	/*
	-----------------------------------------------
	semi-infinite and finite homogeneous layer
	-----------------------------------------------

	Ni et al. (1999) eqs (8)-(11) and the semi-infinite quantities under
	them, at the published optical depth k*elai. The algebra now lives in
	gortt_two_stream so that -rt_crown_path can evaluate it at a second
	depth; with L = p->elai it is bit-for-bit what was written out here.
	*/

	gortt_two_stream( rt->w, rt->gamma, k, p->elai, g->sza_prime,
	                  p->rt_upscatter, rt->beta_up, &slab ) ;

	rt->R_inf_ff = slab.R_inf_ff ;
	rt->R_inf_df = slab.R_inf_df ;
	rt->T_inf_ff = slab.T_inf_ff ;
	rt->T_inf_df = slab.T_inf_df ;
	rt->t_ff     = slab.t_ff ;
	rt->p_ff     = slab.p_ff ;
	rt->t_df     = slab.t_df ;
	rt->p_df     = slab.p_df ;


	/*
	--------------------
	discontinuous canopy
	--------------------
	*/

	/*
	directional hemispherical path transmittance
	for a discontinuous plant canopy of height h
	at given geometry
	*/

	rt->t_prime_df = rt->t_df * ( 1 - rt->t_prime_0 ) ;

	/*
	hemispherical-hemispherical path transmittance
	for a discontinuous plant canopy of height h
	at given geometry
	*/

	k_open = p->k_open[ 0 ] + p->k_openep[ 0 ] ;
	rt->t_prime_ff = rt->t_ff * ( 1.0 - k_open ) + k_open ;

	/*
	the path reflectances of a discontinuous canopy are taken to be those
	of the equivalent homogeneous one (eqs 15-16), unless -rt_crown_path
	or -rt_recollision says otherwise -- see gortt_crown_reflectance
	*/

	gortt_crown_reflectance( p, g, rt, k, &rt->p_prime_df, &rt->p_prime_ff ) ;


	rt->kuusk = gortt_kuusk( p, g ) ;

	return ;

}


