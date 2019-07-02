#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include "mgcal.h"
#include "cdescent.h"

#include "consts.h"
#include "simeq.h"
#include "utils.h"
#include "settings.h"
#include "defaults.h"

void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;

	fprintf (stderr, "\n");
	version_info (p);
	fprintf (stderr, "\n");

	fprintf (stderr, "DESCRIPTION:\n");
	fprintf (stderr, "This program derives the subsurface magnetic model\n");
	fprintf (stderr, "by performing L1-L2 norm regularized inversion.\n\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -a [alpha]\n");
	fprintf (stderr, "       -w [log10_lambda_lower:log10_dlambda]\n");

	fprintf (stderr, "[optional]\n");

	fprintf (stderr, "       -t [tolerance: default=1.e-5]\n");
	fprintf (stderr, "       -m [maximum iteration number: default=1000000]\n");
	fprintf (stderr, "       -n [lower:upper bounds of solutions]\n");
	fprintf (stderr, "       -s [parameter setting file: default=./settings]\n");
	fprintf (stderr, "       -b [initial beta: -b <filename>]\n");
	fprintf (stderr, "       -p (use parallel CDA: default is not use)\n");
	fprintf (stderr, "       -c (use stochastic CDA: default is not use)\n");
	fprintf (stderr, "       -v (verbose mode)\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	return;
}

static int
num_separator (char *str, const char c)
{
	char	*ptr = str;
	int		n = 0;
	while (1) {
		char	*p = strchr (ptr, c);
		if (!p) break;
		ptr = ++p;
		n++;
	}
	return n;
}

bool
read_input_params (int argc, char **argv)
{
	int		nsep;
	int		tmp;
	char	c;

	while ((c = getopt (argc, argv, ":a:w:t:m:l:n:x:s:k:b:z:dpcvh")) != EOF) {
		switch (c) {

			case 'a':
				alpha = (double) atof (optarg);
				break;

			case 'w':
				nsep = num_separator (optarg, ':');
				switch (nsep) {
					case 0:	// lower
						log10_lambda_lower = (double) atof (optarg);
						break;
					case 1:	// lower:dlambda
						sscanf (optarg, "%lf:%lf", &log10_lambda_lower, &log10_dlambda);
						break;
					case 2:	// lower:dlambda:upper
						sscanf (optarg, "%lf:%lf:%lf",
								&log10_lambda_lower, &log10_dlambda, &log10_lambda_upper);
						use_log10_lambda_upper = true;
						break;
					default:
						break;
				}
				break;

			case 't':
				tol = (double) atof (optarg);
				break;

			case 'm':
				maxiter = atoi (optarg);
				break;

			case 'n':
				constraint = true;
				nsep = num_separator (optarg, ':');
				if (nsep == 0) lower = (double) atof (optarg);
				else if (nsep == 1) sscanf (optarg, "%lf:%lf", &lower, &upper);
				else {
					fprintf (stderr, "ERROR: invalid parameter specification: -n %s\n", optarg);
					return false;
				}
				break;

			case 's':
				strcpy (sfn, optarg);
				break;

			case 'b':
				use_initial_beta = true;
				strcpy (bfn, optarg);
				break;

			case 'p':
				parallel = true;
				break;

			case 'c':
				stochastic = true;
				break;

			case 'v':
				verbose = true;
				break;

			case 'h':
				return false;

			case ':':
				fprintf (stdout, "%c needs value.\n", c);
     			return false;

			case '?':
				fprintf (stderr, "unknown option %c.\n", c);
     			return false;

			default:
     			return false;

		}
	}

	return true;
}

static bool
l1l2inv_constraint_func (cdescent *cd, const int j, const double etaj, double *val)
{
	double	scale = 1.;
	if (cd->lreg->xnormalized) scale = sqrt (cd->lreg->xtx[j]);

	*val = lower * scale;
	if (cd->beta->data[j] + etaj < *val) return false;
	*val = upper * scale;
	if (cd->beta->data[j] + etaj > *val) return false;
	return true;
}

bool
l1l2inv (void)
{
	simeq				*eq;
	linregmodel			*lreg;
	cdescent			*cd;

	if (verbose) fprintf (stderr, "preparing simeq object... ");
	if ((eq = read_input (type, ifn, tfn)) == NULL) return false;
	if (verbose) fprintf (stderr, "done\n");

	if (verbose) fprintf (stderr, "preparing linregmodel object... ");
	lreg = linregmodel_new (eq->y, eq->x, eq->d, DO_NORMALIZING_X);
	if (verbose) fprintf (stderr, "done\n");
#ifdef DEBUG
	{
		int			k;
		double		ret;
		mm_dense	*z;
		
		ret = mm_real_xj_trans_dot_yk (lreg->x, 2399, lreg->y, 0);
		fprintf (stdout, "ret = %f\n", ret);

		z = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, lreg->x->n, 1, lreg->x->n);
		mm_real_set_all (z, 0.);
		mm_real_x_dot_yk (true, 1., lreg->x, lreg->y, 0, 0., z);
		mm_real_fwrite (stdout, z, "%f");

		mm_real_axjpy (1., lreg->x, 810, lreg->y);
		mm_real_fwrite (stdout, lreg->y, "%f");

		mm_real_axjpy_atomic (1., lreg->x, 1230, lreg->y);
		mm_real_fwrite (stdout, lreg->y, "%f");

		for (k = 0; k < lreg->x->n; k++) {
			fprintf (stdout, "sx[%d] = %f\n", k, lreg->sx[k]);
			fprintf (stdout, "xtx[%d] = %f\n", k, lreg->xtx[k]);
		}	
		exit (1);
	}
#endif
	if (verbose) fprintf (stderr, "preparing cdescent object... ");
	cd = cdescent_new (alpha, lreg, tol, maxiter, parallel);
	if (use_initial_beta) {
		FILE	*fp = fopen (bfn, "r");
		if (fp) {
			mm_dense	*beta0 = mm_real_fread (fp);
			fclose (fp);
			cdescent_init_beta (cd, beta0);
			mm_real_free (beta0);
		}
	}
	if (verbose) fprintf (stderr, "done\n");
	fprintf (stderr, "[m, n] = [%d, %d]\n", *cd->m, *cd->n);

	if (stochastic) {
		time_t	t = time (NULL);
		cdescent_set_stochastic (cd, (unsigned int *) &t);
	}

	cdescent_not_use_intercept (cd);
	if (constraint) cdescent_set_constraint (cd, l1l2inv_constraint_func);
	if (verbose) cd->verbose = true;

	cdescent_set_outputs_fullpath (cd, NULL);
	cdescent_set_outputs_info (cd, NULL);
	cdescent_set_log10_lambda_lower (cd, log10_lambda_lower);
	cdescent_set_log10_dlambda (cd, log10_dlambda);
	if (use_log10_lambda_upper) cdescent_set_log10_lambda_upper (cd, log10_lambda_upper);

	fprintf (stderr, "regression start\n");
	if (!cdescent_do_pathwise_optimization (cd)) fprintf (stderr, "not converged.\n");
	fprintf (stderr, "total num of iter = %d\n", cd->total_iter);
	if (cd->use_intercept) fprintf (stderr, "intercept = %.4e\n", cd->b0);

	cdescent_free (cd);
	simeq_free (eq);

	return true;
}

int
main (int argc, char **argv)
{
	if (argc <= 1 || !read_input_params (argc, argv)) {
		usage (argv[0]);
		return EXIT_FAILURE;
	}
	if (!read_settings (sfn)) return EXIT_FAILURE;
	mgcal_set_scale_factor (magscale);

	fprintf_settings (stderr);
	l1l2inv ();

	return EXIT_SUCCESS;
}

