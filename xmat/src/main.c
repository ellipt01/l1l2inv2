#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include "mgcal.h"
#include "cdescent.h"

#include "consts.h"
#include "simeq_xmat.h"
#include "utils_xmat.h"
#include "l1l2inv_xmat.h"
#include "settings.h"
#include "defaults.h"

bool	create_xmat = true;
bool	penalty_for_actual_magnetization = false;

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
	fprintf (stderr, "by performing L1-L2 norm regularized inversion.\n");
	fprintf (stderr, "** This is a version that store matrix X in files **\n\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -a [alpha]\n");
	fprintf (stderr, "       -w [log10_lambda_lower:log10_dlambda]\n");

	fprintf (stderr, "[optional]\n");

	fprintf (stderr, "       -t [tolerance: default=1.e-5]\n");
	fprintf (stderr, "       -m [maximum iteration number: default=1000000]\n");
	fprintf (stderr, "       -n [lower:upper bounds of solutions]\n");
	fprintf (stderr, "       -s [parameter setting file: default=./settings]\n");
	fprintf (stderr, "       -b [initial beta: -b <filename>]\n");
	fprintf (stderr, "       -g [0(false) or 1(true): stretch the grid cells\n");
	fprintf (stderr, "           at the edge of the model space outward,\n");
	fprintf (stderr, "           default is 1]\n");
	fprintf (stderr, "       -p (use parallel CDA: default is not use)\n");
	fprintf (stderr, "       -c (use stochastic CDA: default is not use)\n");
	fprintf (stderr, "       -v (verbose mode)\n");
	fprintf (stderr, "       -x (use already exists xmat files, default=false)\n");
	fprintf (stderr, "       -h (show this message)\n\n");

	fprintf (stderr, "       -r [regression type: 0=L1,1=L1L2,2=L1TSV,3=L1L2TSV\n");
	fprintf (stderr, "           default is 1=L1L2]\n"); 
	fprintf (stderr, "       -d [wx:wy:wz(L1TSV) or w0:wx:wy:wz(L1L2TSV)\n");
	fprintf (stderr, "           when regression type is 2=L1TSV or 3=L1L2TSV,\n"); 
	fprintf (stderr, "           quadratic penalty is constructed as\n");
	fprintf (stderr, "           [wx*Dx;wy*Dy;wz*Dz] or [w0*E;wx*Dx;wy*Dy;wz*Dz]]\n");
	fprintf (stderr, "       -k (use quadratic penalty for actual magnetization)\n\n");

	return;
}

bool
read_input_params (int argc, char **argv)
{
	int		nsep;
	int		tmp;
	char	c;

	stretch_grid_at_edge = true;
	while ((c = getopt (argc, argv, ":r:d:a:w:t:m:n:s:b:g:kpcouxvh")) != EOF) {
		switch (c) {

			case 'r':
				type = atoi(optarg);
				break;

			case 'd':
				nsep = num_separator (optarg, ':');
				switch (nsep) {
					case 2:
						weight = (double *) malloc (3 * sizeof (double));
						sscanf (optarg, "%lf:%lf:%lf", &weight[0], &weight[1], &weight[2]);
						break;
					case 3:
						weight = (double *) malloc (4 * sizeof (double));
						sscanf (optarg, "%lf:%lf:%lf:%lf", &weight[0], &weight[1], &weight[2], &weight[3]);
						break;
					default:
						break;
				}
				break;

			case 'k':
				penalty_for_actual_magnetization = true;
				break;

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

			case 'o':
				output_vector = true;
				break;

			case 'u':
				output_weighted = true;
				break;

			case 'x':
				create_xmat = false;
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

	if (type != TYPE_L1 && type != TYPE_L1L2
		&& type != TYPE_L1TSV && type != TYPE_L1L2TSV) {
		fprintf (stderr, "ERROR: type must be 0, 1, 2, or 3\n");
		return false;	
	}

	return true;
}

static void
weight_d (mm_real *d, double *xtx)
{
	int		j;
	for (j = 0; j < d->n; j++) mm_real_xj_scale (d, j, 1. / sqrt (xtx[j]));
	return;
}

bool
run (void)
{
	simeq	*eq;

	if (verbose) fprintf (stderr, "preparing simeq object... ");
	if ((eq = read_input (type, ifn, tfn, weight, create_xmat)) == NULL) return false;
	if (verbose) fprintf (stderr, "done\n");

	if (penalty_for_actual_magnetization) {
		double	*xtx;
		int		n = eq->d->n;
		FILE	*fp = fopen ("xtx.vec", "rb");
		if (!fp) {
			fprintf (stderr, "ERROR: cannot open file xtx.vec.\nprogram abort!\n");
			exit (1);
		} else {
			int	ret;
			xtx = (double *) malloc (n * sizeof (double));
			ret = fread (xtx, sizeof (double), n, fp);
			fclose (fp);
		}
		weight_d (eq->d, xtx);
		free (xtx);
	}

	l1l2inv (eq, NULL, NULL);

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
	run ();

	return EXIT_SUCCESS;
}

