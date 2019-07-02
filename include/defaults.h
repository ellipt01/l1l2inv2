#ifndef _DEFAULTS_H_
#define _DEFAULTS_H_

char	ifn[80] = "input.data";
char	tfn[80] = "terrain.data";
char	bfn[80];

int		type = TYPE_L0;
double	alpha = 1.;
double	log10_lambda_upper = 0.;
double	log10_lambda_lower = -4.;
double	log10_dlambda = 0.1;
double	tol = 1.e-5;
int		maxiter = 100000000;

// upper and lower bounds of solution
double	lower = -INFINITY;
double	upper = +INFINITY;

double	lambda = 0.;

/* scale factor of mgcal */
double	magscale = 100.;

/*** FLAGS ***/
bool	use_log10_lambda_upper = false;
// use initial beta
bool	use_initial_beta = false;
// use upper and lower limits of solutions
bool	constraint = false;
// output vector to file
bool	output_vector = false;
// output matrix to file
bool	output_matrix = false;
// use parallel CDA
bool	parallel = false;
// use stochastic CDA
bool	stochastic = false;
// verbose mode
bool	verbose = false;

#endif // _DEFAULTS_H_
