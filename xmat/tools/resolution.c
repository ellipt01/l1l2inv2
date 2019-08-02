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
int		nx0 = 0, nx1 = -1;
int		ny0 = 0, ny1 = -1;
int		nz0 = 0, nz1 = -1;
double	bval = 1.;

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

	fprintf (stderr, "       -b [intensity of test sources: default is 1 A/m]\n");
	fprintf (stderr, "       -t [tolerance: default=1.e-5]\n");
	fprintf (stderr, "       -m [maximum iteration number: default=1000000]\n");
	fprintf (stderr, "       -n [lower:upper bounds of solutions]\n");
	fprintf (stderr, "       -r [range for evaluate resolution.\n");
	fprintf (stderr, "           format is, nx0:nx1:ny0:ny1:nz0:nz1,\n");
	fprintf (stderr, "           if nx1, ny1, or nz1 is set to -1,\n");
	fprintf (stderr, "           ngrd[0], ngrd[1], and ngrd[2] is used.\n");
	fprintf (stderr, "           the test sources are placed on grid numper of\n");
	fprintf (stderr, "           nx0 to nx1 with inclements of incx,\n");
	fprintf (stderr, "           ny0 to ny1 with inclements of incy,\n");
	fprintf (stderr, "           nz0 to nz1 with inclements of incz, respectively.\n");
	fprintf (stderr, "           default is 0:-1:0:-1:0:-1]\n");
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
	bool	set_a = false;
	bool	set_w = false;

	while ((c = getopt (argc, argv, ":l:i:a:w:b:t:m:n:r:s:pcxvh")) != EOF) {

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
				set_a = true;
				break;

			case 'w':
				nsep = num_separator (optarg, ':');
				switch (nsep) {
					case 0:	// lower
						log10_lambda_lower = (double) atof (optarg);
						set_w = true;
						break;
					case 1:	// lower:dlambda
						sscanf (optarg, "%lf:%lf", &log10_lambda_lower, &log10_dlambda);
						set_w = true;
						break;
					case 2:	// lower:dlambda:upper
						sscanf (optarg, "%lf:%lf:%lf",
								&log10_lambda_lower, &log10_dlambda, &log10_lambda_upper);
						use_log10_lambda_upper = true;
						set_w = true;
						break;
					default:
						break;
				}
				break;

			case 'b':
				bval = (double) atof (optarg);
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

	if (!set_l || !set_i || !set_a || !set_w) return false;

	return true;
}

// read xj from xmat file
static mm_dense *
read_xj_from_xmatfile (int j, int l, int m)
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
anomaly_by_grids_xmatfile (mm_dense *y, int i0, int j0, int k0, int lx, int ly, int lz, int m, double *w)
{
	int			i, j, k;
	mm_dense	*x;

	for (i = 0; i < lx; i++) {
		int		ix = i0 + i;
		for (j = 0; j < ly; j++) {
			int		iy = j0 + j;
			for (k = 0; k < lz; k++) {
				int		iz = k0 + k;
				int		s = ix + iy * ngrd[0] + iz * ngrd[0] * ngrd[1];
				int		l = (int) (s / xfile_len);
				x = read_xj_from_xmatfile (s, l, m);
				// y = norm(x) * x + y
				mm_real_axjpy (bval * sqrt (w[s]), x, 0, y);
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

	mm_real_set_all (eq->y, 0.);
	for (k = nz0; k < nz1; k += incz) {
		for (j = ny0; j < ny1; j += incy) {
			for (i = nx0; i < nx1; i += incx) {
				anomaly_by_grids_xmatfile (eq->y, i, j, k, lx, ly, lz, m, xtx->data);
			}
		}
	}

	{
		FILE	*fp = fopen ("y.vec", "w");
		if (fp) {
			mm_real_fwrite (fp, eq->y, "%f");
			fclose (fp);
		}
	}
	l1l2inv (eq, NULL, NULL);

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

	if (nx1 == -1) nx1 = ngrd[0];
	if (ny1 == -1) ny1 = ngrd[1];
	if (nz1 == -1) nz1 = ngrd[2];
	resolution ();

	return EXIT_SUCCESS;
}

