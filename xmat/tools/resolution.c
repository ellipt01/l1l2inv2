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
#include "l1l2inv.h"
#include "settings.h"
#include "defaults.h"

bool	create_xmat = true;
int		num_xfiles;
int		xfile_len;
int		lx, ly, lz;
int		incx, incy, incz;
int		nx0, nx1, ny0, ny1, nz0, nz1;

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
	fprintf (stderr, "using the synthetic anomaly\n");
	fprintf (stderr, "produced by a magnetized block\n");
	fprintf (stderr, "(with dimension of [lx, ly, lz]).\n");
	fprintf (stderr, "The distance between true model (a magnetized block)\n");
	fprintf (stderr, "and estimated model indicates the \"resolution\"\n");
	fprintf (stderr, "at the point where the magnetized block located.\n");
	fprintf (stderr, "** This is a version that store matrix X in files **\n\n");

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
	fprintf (stderr, "       -r [range for calculation:\n");
	fprintf (stderr, "           specify start and end grid number\n");
	fprintf (stderr, "           to be estimated the resolution.\n");
	fprintf (stderr, "           format is, grid num of\n");
	fprintf (stderr, "           xstart:xend:ystart:yend:zstart:zend,\n");
	fprintf (stderr, "           default is 0:ngrd[0]:0:ngrd[1]:0:ngrd[2],\n");
	fprintf (stderr, "           if xend, yend, or zend is specified as -1,\n");
	fprintf (stderr, "           xend=ngrd[0], yend=ngrd[1], zend=ngrd[2] is used]\n");
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

	nx0 = 0;
	nx1 = ngrd[0];
	ny0 = 0;
	ny1 = ngrd[2];
	nz0 = 0;
	nz1 = ngrd[2];

	while ((c = getopt (argc, argv, ":l:i:a:w:t:m:n:r:s:pcxvh")) != EOF) {

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

			case 'r':
				nsep = num_separator (optarg, ':');
				if (nsep != 5) {
					fprintf (stderr, "ERROR: parameter specification invalid: -r %s\n", optarg);
					return false;
				}
				sscanf (optarg, "%d:%d:%d:%d:%d:%d", &nx0, &nx1, &ny0, &ny1, &nz0, &nz1);
				if (nx1 == -1) nx1 = ngrd[0];
				if (ny1 == -1) ny1 = ngrd[1];
				if (nz1 == -1) nz1 = ngrd[2];
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

	if (!set_l || !set_i) return false;

	return true;
}

// read xj from xmat file
static mm_dense *
read_xj_from_xmat (int j, int l, int m)
{
	mm_dense	*x;
	char		fn_xmat[80];
	FILE		*fp_xmat;
	sprintf (fn_xmat, "x%03d.mat", l);
	if ((fp_xmat = fopen (fn_xmat, "rb")) == NULL) {
		fprintf (stderr, "ERROR: cannot open file %s.\n", fn_xmat);
		exit (1);
	}
	x = mm_real_read_xj_xmatfile (fp_xmat, j, m);
	fclose (fp_xmat);
	return x;
}

/* read X of columns in ranges of
   x in [ix - lx / 2, ix + lx / 2]
   y in [iy - ly / 2, iy + ly / 2]
   z in [iz - lz / 2, iz + lz / 2].
   sum of the above columns is store in y */
static void
anomaly_by_grids (mm_dense *y, int i0, int j0, int k0, int lx, int ly, int lz, int m, double *w)
{
	int			i, j, k;
	mm_dense	*x;

	mm_real_set_all (y, 0.);
	for (i = 0; i < lx; i++) {
		int		ix = i0 + i;
		for (j = 0; j < ly; j++) {
			int		iy = j0 + j;
			for (k = 0; k < lz; k++) {
				int		iz = k0 + k;
				int		s = ix + iy * ngrd[0] + iz * ngrd[0] * ngrd[1];
				int		l = (int) (s / xfile_len);
				x = read_xj_from_xmat (s, l, m);
				// y = norm(x) * x + y
				mm_real_axjpy (sqrt (w[s]), x, 0, y);
				mm_real_free (x);
			}
		}
	}
	return;
}

bool
resolution (void)
{
	int			i, j, k;
	int			m, n;
	simeq		*eq;
	mm_dense	*xtx;

	if (verbose) fprintf (stderr, "preparing simeq object... ");
	if ((eq = read_input (type, ifn, tfn, create_xmat)) == NULL) return false;
	if (verbose) fprintf (stderr, "done\n");

	num_xfiles = ngrd[2];
	xfile_len = ngrd[0] * ngrd[1];
	m = eq->y->m;
	n = ngrd[0] * ngrd[1] * ngrd[2];

	{
		char	fn_xtx[] = "xtx.vec";
		FILE	*fp_xtx = fopen (fn_xtx, "rb");
		if (!fp_xtx) {
			fprintf (stderr, "ERROR: cannot open file %s.\n", fn_xtx);
			exit (1);
		}
		xtx = mm_real_read_xj_xmatfile (fp_xtx, 0, n);
	}

	for (k = nz0; k < nz1; k += incz) {
		for (j = ny0; j < ny1; j += incy) {
			for (i = nx0; i < nx1; i += incx) {

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

				anomaly_by_grids (eq->y, i, j, k, lx, ly, lz, m, xtx->data);

				l1l2inv (eq, path_fn, info_fn);

			}
		}
	}

	mm_real_free (xtx);
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

