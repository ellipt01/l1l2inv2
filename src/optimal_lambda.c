#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_spline.h>

#include "mgcal.h"
#include "cdescent.h"

#include "consts.h"
#include "simeq.h"
#include "utils.h"
#include "settings.h"

char	ifn[80];
char	ofn[80] = "curvature_spl.data";
double	dt = 0.01;

typedef struct {
	int		n;
	double	*lambda;
	double	*curve;
} curvature;

curvature *
curvature_new (const int n)
{
	curvature	*c = (curvature *) malloc (sizeof (curvature));
	c->n = n;
	c->lambda = (double *) malloc (n * sizeof (double));
	c->curve = (double *) malloc (n * sizeof (double));
	return c;
}

void
curvature_free (curvature *c)
{
	free (c->lambda);
	free (c->curve);
	free (c);
	return;
}

static int
count_data (FILE *fp)
{
	int		count;
	char	buf[BUFSIZ];

	count = 0;
	while (fgets (buf, BUFSIZ, fp) != NULL) count++;
	fseek (fp, 0L, SEEK_SET);
	return count;
}

static curvature *
read_data (const char *fn)
{
	int			i, n;
	curvature	*cv;
	char		buf[BUFSIZ];
	FILE		*fp;

	if ((fp = fopen (fn, "r")) == NULL) {
		fprintf (stderr, "ERROR: file %s not exists.\n", fn);
		return NULL;
	}
	n = count_data (fp);

	cv = curvature_new (n);
	i = 0;
	while (fgets (buf, BUFSIZ, fp) != NULL) {
		double	l, c, v;
		sscanf (buf, "%lf %lf %lf", &l, &c, &v);
		cv->lambda[i] = l;
		cv->curve[i] = c;
		if (++i > n) break;
	}
	fclose (fp);
	return cv;
}

static double
interp (const int n, const double *x, const double *y)
{
	int					i;
	double				*t;
	double				ti;
	double				max[2] = {0., 0.};
	gsl_spline			*spl = gsl_spline_alloc (gsl_interp_cspline, n);
	gsl_interp_accel	*acc = gsl_interp_accel_alloc ();

	FILE				*fp = fopen (ofn, "w");

	t = (double *) malloc (n * sizeof (double));
	for (i = 0; i < n; i++) t[i] = log10 (x[i]);

	gsl_spline_init (spl, t, y, n);

	for (i = 0, ti = t[0]; ti <= t[n - 1]; i++, ti += dt) {
		double	xi = pow (10., ti);
		double	yi = gsl_spline_eval (spl, ti, acc);
		if (i == 0 || max[1] < yi) {
			max[0] = xi;
			max[1] = yi;
		}
		if (fp) fprintf (fp, "%.8e, %.8e\n", xi, yi);
	}
	if (fp) fclose (fp);
	gsl_spline_free (spl);
	gsl_interp_accel_free (acc);
	free (t);

	return max[0];
}

static void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;

	fprintf (stderr, "\n");
	version_info (p);
	fprintf (stderr, "\n");

	fprintf (stderr, "DESCRIPTION:\n");
	fprintf (stderr, "This program estimates the optimal regularization parameter lambda\n");
	fprintf (stderr, "which realizes the maximum curvature of the L-curve\n");
	fprintf (stderr, "with aid of the smooth B-spline.\n");
	fprintf (stderr, "NOTE: this program uses GSL(Gnu Scientific Library) version 2.0 or later.\n\n");

	fprintf (stderr, "OUTPUT:\n");
	fprintf (stderr, "[STDOUT]: most optimal lambda\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -f <name of file which stores curvature of L-curve,\n");
	fprintf (stderr, "          default is curvature.data,\n");
	fprintf (stderr, "          which is created by \"lcurve_interp\".\n");
	fprintf (stderr, "       -d <log10(incx) of spline: default=0.01>\n");
	fprintf (stderr, "       -o <output filename>\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	return;
}

static bool
read_input_params (int argc, char **argv)
{
	char	c;
	bool	infile_specified = false;

	while ((c = getopt (argc, argv, "f:d:o:h")) != EOF) {
		switch (c) {
			case 'f':
				strcpy (ifn, optarg);
				infile_specified = true;
				break;
			case 'd':
				dt = (double) atof (optarg);
				break;
			case 'o':
				strcpy (ofn, optarg);
				break;
			case 'h':
			case ':':
			case '?':
				return false;
				break;
		}
	}
	if (!infile_specified) {
		fprintf (stderr, "ERROR: curvature data file is not specified.\n");
		usage (argv[0]);
	}
	return true;
}

int
main (int argc, char **argv)
{
	double		optl;
	curvature	*cv;

	if (!read_input_params (argc, argv)) {
		usage (argv[0]);
		exit (EXIT_FAILURE);
	}

	cv = read_data (ifn);
	if (!cv) return EXIT_FAILURE;

	optl = interp (cv->n, cv->lambda, cv->curve);
	fprintf (stdout, "%.8e\n", optl);
	curvature_free (cv);

	return EXIT_SUCCESS;
}


