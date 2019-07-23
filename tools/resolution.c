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
#include "l1l2inv.h"
#include "settings.h"
#include "defaults.h"

int		lx, ly, lz;
int		incx, incy, incz;

extern void		dcopy_  (const int *n, const double *x, const int *incx, double *y, const int *incy);

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
	fprintf (stderr, "This program evaluates the model resolution\n");
	fprintf (stderr, "by performing L1-L2 norm inversion\n");
	fprintf (stderr, "using the input data of a synthetic anomaly\n");
	fprintf (stderr, "that produced by a subsurface block\n");
	fprintf (stderr, "(with dimension of [lx, ly, lz]).\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -l [lx:ly:lz: dimension of the synthetic block]\n");
	fprintf (stderr, "       -i [incx:incy:incz:\n");
	fprintf (stderr, "           interval to evaluate the model resolution]\n");
	fprintf (stderr, "       -a [alpha]\n");
	fprintf (stderr, "       -w [log10_lambda_lower:log10_dlambda]\n");

	fprintf (stderr, "[optional]\n");

	fprintf (stderr, "       -t [tolerance: default=1.e-5]\n");
	fprintf (stderr, "       -m [maximum iteration number: default=1000000]\n");
	fprintf (stderr, "       -n [lower:upper bounds of solutions]\n");
	fprintf (stderr, "       -s [parameter setting file: default=./settings]\n");
	fprintf (stderr, "       -c (use stochastic CDA: default is not use)\n");
	fprintf (stderr, "       -v (verbose mode)\n");
	fprintf (stderr, "       -x (use already exists xmat files, default=false)\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	return;
}

bool
read_input_params (int argc, char **argv)
{
	int		nsep;
	int		tmp;
	char	c;

	bool	set_l = false;
	bool	set_i = false;

	while ((c = getopt (argc, argv, ":l:i:a:w:t:m:n:s:pcvh")) != EOF) {

		switch (c) {

			case 'l':
				nsep = num_separator (optarg, ':');
				if (nsep != 2) {
					fprintf (stderr, "ERROR: parameter specification invalid: -l %s\n", optarg);
					return false;
				}
				sscanf (optarg, "%d:%d:%d", &lx, &ly, &lz);
				set_l = true;
				break;

			case 'i':
				nsep = num_separator (optarg, ':');
				if (nsep != 2) {
					fprintf (stderr, "ERROR: parameter specification invalid: -i %s\n", optarg);
					return false;
				}
				sscanf (optarg, "%d:%d:%d", &incx, &incy, &incz);
				set_i = true;
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

	if (!set_l || !set_i) return false;

	return true;
}

bool
resolution (void)
{
	int			i, j, k;
	int			m, n;
	simeq		*eq;
	mm_dense	*x;

	if (verbose) fprintf (stderr, "preparing simeq object... ");
	if ((eq = read_input (type, ifn, tfn)) == NULL) return false;
	if (verbose) fprintf (stderr, "done\n");

	m = eq->x->m;
	n = eq->x->n;

	x = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, 1, m);

	for (i = 0; i < ngrd[0]; i += incx) {

		for (j = 0; j < ngrd[1]; j += incy) {

			for (k = 0; k < ngrd[2]; k += incz) {
				int			p, q, r;

				char		path_fn[80];
				char		info_fn[80];

				// set output file name
				{
					int		ix0, iy0, iz0;

					ix0 = i + (int) ((double) lx / 2.);
					iy0 = j + (int) ((double) ly / 2.);
					iz0 = k + (int) ((double) lz / 2.);

					sprintf (path_fn, "beta_path%03d%03d%03d.data", ix0 , iy0, iz0);
					sprintf (info_fn, "regression_info%03d%03d%03d.data", ix0 , iy0, iz0);
				}

				// read X of columns in ranges of
				// x in [ix - lx / 2, ix + lx / 2]
				// y in [iy - ly / 2, iy + ly / 2]
				// z in [iz - lz / 2, iz + lz / 2].
				// sum of the above columns is store in eq->y
				mm_real_set_all (eq->y, 0.);
				for (p = 0; p < lx; p++) {
					int		ix = i + p;
					for (q = 0; q < ly; q++) {
						int		iy = j + q;
						for (r = 0; r < lz; r++) {
							int		one = 1;
							int		iz = k + r;
							int		l = ix + iy * ngrd[0] + iz * ngrd[0] * ngrd[1];
							dcopy_ (&m, eq->x->data + l, &one, x->data, &one);
							mm_real_axjpy (1., x, 0, eq->y);
						}
					}
				}

				l1l2inv (eq, path_fn, info_fn);

			}
		}
	}

	mm_real_free (x);
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
	resolution ();

	return EXIT_SUCCESS;
}

