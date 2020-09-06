#ifndef _EXTERN_H_
#define _EXTERN_H_

extern char		*ifn;
extern char		*tfn;
extern char		*sfn;
extern char		*bfn;

extern int		type;
extern double	alpha;
extern double	log10_lambda_upper;
extern double	log10_lambda_lower;
extern double	log10_dlambda;
extern double	tol;
extern int		maxiter;

extern double	*w;

// upper and lower bounds of solution
extern double	lower;
extern double	upper;

extern double	lambda;

/* scale factor of mgcal */
extern double	magscale;

/*** FLAGS ***/
extern bool		use_log10_lambda_upper;
// use initial beta
extern bool		use_initial_beta;
// use upper and lower limits of solutions
extern bool		constraint;
// output vector to file
extern bool		output_vector;
// output weighted beta
extern bool		output_weighted;
// output matrix to file
extern bool		output_matrix;
// use parallel CDA
extern bool		parallel;
// use stochastic CDA
extern bool		stochastic;
// stretching the grid on the edge of the model space
extern bool		stretching_grid;
// verbose mode
extern bool		verbose;

extern bool		use_fixed_lambda1;
extern double	fixed_lambda1;
extern bool		use_fixed_lambda2;
extern double	fixed_lambda2;

#endif // _EXTERN_H_
