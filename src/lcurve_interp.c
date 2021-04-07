#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>

#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_version.h>

#include "mgcal.h"
#include "cdescent.h"

#include "simeq.h"
#include "utils.h"
#include "settings.h"

static bool		infile_specified = false;
static char		fn[80];
static double	alpha = 1.0;
static double	beta_val = 1.e-1;

void
fprintf_curvature (FILE *stream, const gsl_vector *l, const gsl_vector *cv, const gsl_vector *y)
{
	int		i;
	for (i = 0; i < cv->size; i++)
		fprintf (stream, "%.4e\t%.4e\t%.4e\n",
			pow (10., gsl_vector_get (l, i + 2)), gsl_vector_get (cv, i), pow (10, gsl_vector_get (y, i + 2)));
	return;
}

void
fprintf_lxy_vector (FILE *stream, const gsl_vector *l, const gsl_vector *x, const gsl_vector *y)
{
	int		i;
	for (i = 0; i < x->size; i++) {
		fprintf (stream, "%f\t%f\t%f\n",
			pow (10., gsl_vector_get (l, i)),
			pow (10., gsl_vector_get (x, i)),
			pow (10., gsl_vector_get (y, i)));
	}
	return;
}

void
fprintf_synth_vector (FILE *stream, const gsl_vector *l,
	const gsl_vector *x, const gsl_vector *y,
	const gsl_vector *dx, const gsl_vector *dy,
	const gsl_vector *ddx, const gsl_vector *ddy)
{
	int		i;
	for (i = 0; i < x->size; i++) {
		fprintf (stream, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
			pow (10., gsl_vector_get (l, i)),
			pow (10., gsl_vector_get (x, i)),
			pow (10., gsl_vector_get (y, i)),
			pow (10., gsl_vector_get (dx, i)),
			pow (10., gsl_vector_get (dy, i)),
			pow (10., gsl_vector_get (ddx, i)),
			pow (10., gsl_vector_get (ddy, i)));
	}
	return;
}

/* 基底関数の行列を作成 */
gsl_matrix *
bspline_fit_matrix (const gsl_vector *x, gsl_bspline_workspace *bw)
{
	int			i, j, k;
	size_t		ncoef = gsl_bspline_ncoeffs (bw);
	gsl_matrix	*a = gsl_matrix_alloc (x->size, ncoef);

	k = gsl_bspline_order (bw);
	gsl_matrix_set_zero (a);
	for (i = 0; i < x->size; i++) {
		size_t		istart, iend;
		double		xi = gsl_vector_get (x, i);
		gsl_vector	*bi = gsl_vector_alloc (k);
		gsl_bspline_eval_nonzero (xi, bi, &istart, &iend, bw);
		for (j = 0; j < bi->size; j++) gsl_matrix_set (a, i, j + istart, gsl_vector_get (bi, j));
		gsl_vector_free (bi);
	}
	return a;
}

/* 基底関数の微分行列を作成 */
gsl_matrix *
#ifndef GSL_MAJOR_VERSION
bspline_nderiv_matrix (size_t nderiv, const gsl_vector *x, gsl_bspline_workspace *bw, gsl_bspline_deriv_workspace *dw)
#else
bspline_nderiv_matrix (size_t nderiv, const gsl_vector *x, gsl_bspline_workspace *bw)
#endif
{
	int			i, j, k;
	size_t		ncoef;
	size_t		nderiv_max = 2;
	gsl_matrix	*a;

	if (nderiv < 0 || 2 < nderiv) {
		fprintf (stderr, "ERROR: nderiv must be 0, 1 or 2.\n");
		exit (EXIT_FAILURE);
	}

	ncoef = gsl_bspline_ncoeffs (bw);
	a = gsl_matrix_alloc (x->size, ncoef);

	k = gsl_bspline_order (bw);
	gsl_matrix_set_zero (a);
	for (i = 0; i < x->size; i++) {
		size_t		istart, iend;
		double		xi = gsl_vector_get (x, i);
		gsl_matrix	*db = gsl_matrix_alloc (k, nderiv_max + 1);
#ifndef GSL_MAJOR_VERSION
		gsl_bspline_deriv_eval_nonzero (xi, nderiv, db, &istart, &iend, bw, dw);
#else
		gsl_bspline_deriv_eval_nonzero (xi, nderiv, db, &istart, &iend, bw);
#endif
		for (j = 0; j < db->size1; j++) gsl_matrix_set (a, i, j + istart, gsl_matrix_get (db, j, nderiv));
		gsl_matrix_free (db);
	}
	return a;
}

/* QR分解で連立方程式 b = a * x を解く */
gsl_vector *
bspline_eval_coeffs_QR (gsl_matrix *a, gsl_vector *b)
{
	gsl_vector	*x;
	gsl_vector	*tau;
	gsl_vector	*res;

	tau = gsl_vector_alloc (GSL_MIN (a->size1, a->size2));
	gsl_linalg_QR_decomp (a, tau);
	x = gsl_vector_alloc (a->size2);
	res = gsl_vector_alloc (a->size1);
	gsl_linalg_QR_lssolve (a, tau, b, x, res);	   

	gsl_vector_free (tau);
	gsl_vector_free (res);

	return x;
}

gsl_matrix *
smooth_matrix (size_t n)
{
	int			i;
	gsl_matrix	*c;

	c = gsl_matrix_alloc (n - 2, n);
	gsl_matrix_set_all (c, 0.);
	for (i = 0; i < c->size1; i++) {
		gsl_matrix_set (c, i, i, 1.);
		gsl_matrix_set (c, i, i + 1, -2.);
		gsl_matrix_set (c, i, i + 2, 1.);
	}
	return c;
}

gsl_vector *
bspline_smooth_eval_coeffs (const double beta, gsl_matrix *a, gsl_vector *b)
{
	int			i, j;
	gsl_vector	*atb = gsl_vector_alloc (a->size2);
	gsl_matrix	*ata = gsl_matrix_alloc (a->size2, a->size2);
	gsl_matrix	*c;
	gsl_matrix	*ctc;
	gsl_vector	*e;

	gsl_blas_dgemv (CblasTrans, 1., a, b, 0., atb);

	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1., a, a, 0., ata);

	c = smooth_matrix (a->size2);
	ctc = gsl_matrix_alloc (c->size2, c->size2);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, beta, c, c, 0., ctc);
	gsl_matrix_free (c);

	for (i = 0; i < ata->size1; i++) {
		for (j = 0; j < ata->size2; j++) {
			gsl_matrix_set (ata, i, j, gsl_matrix_get (ata, i, j) + gsl_matrix_get (ctc, i, j));
		}
	}
	gsl_matrix_free (ctc);

	e = bspline_eval_coeffs_QR (ata, atb);

	gsl_vector_free (atb);
	gsl_matrix_free (ata);

	return e;
}

/* 得られたBスプライン係数からデータを再現 */
gsl_vector *
bspline_synthesize (const gsl_matrix *a, const gsl_vector *e)
{
	gsl_vector	*y = gsl_vector_alloc (a->size1);
	/* y = a * e */
	gsl_blas_dgemv (CblasNoTrans, 1., a, e, 0., y);
	return y;
}

int
count_data (char *fn)
{
	int		n;
	char	buf[BUFSIZ];
	FILE	*fp;

	if ((fp = fopen (fn, "r")) == NULL) {
		fprintf (stderr, "ERROR: cannot open file %s.\n", fn);
		exit (EXIT_FAILURE);
	}

	n = 0;
	while (fgets (buf, BUFSIZ, fp) != NULL) {
		if (buf[0] == '#' || buf[0] == ' '
			|| buf[0] == '\t' || buf[0] == '\n') continue;
		n++;
	}
	fclose (fp);
	return n;
}

void
read_data (char *fn, gsl_vector **_l, gsl_vector **_x, gsl_vector **_y)
{
	int			i;
	int			n;
	gsl_vector	*l;
	gsl_vector	*x;
	gsl_vector	*y;

	char		buf[BUFSIZ];
	FILE		*fp;

	n = count_data (fn);
	l = gsl_vector_alloc (n);
	x = gsl_vector_alloc (n);
	y = gsl_vector_alloc (n);

	if ((fp = fopen (fn, "r")) == NULL) {
		fprintf (stderr, "ERROR: cannot open file %s.\n", fn);
		exit (EXIT_FAILURE);
	}

	i = 0;
	while (fgets (buf, BUFSIZ, fp) != NULL) {
		double	nrm1, nrm2, rss2, lambda1, lambda2;
		if (buf[0] == '#' || buf[0] == ' '
			|| buf[0] == '\t' || buf[0] == '\n') continue;
		sscanf (buf, "%lf\t%lf\t%lf\t%lf\t%lf", &nrm1, &nrm2, &rss2, &lambda1, &lambda2);
		if (fabs (nrm1) <= 1.e-12) continue;
		gsl_vector_set (x, i, log10 (rss2));

		if (fabs (alpha) > 1.e-8) gsl_vector_set (l, i, log10 (lambda1 / alpha));
		else gsl_vector_set (l, i, log10 (lambda2));
		gsl_vector_set (y, i, log10 (alpha * nrm1 + (1. - alpha) * nrm2 / 2.));

		i++;
	}
	if (i != n) {
		gsl_vector	*l1 = gsl_vector_alloc (i);
		gsl_vector	*x1 = gsl_vector_alloc (i);
		gsl_vector	*y1 = gsl_vector_alloc (i);

		for (i = 0; i < l1->size; i++) {
			gsl_vector_set (l1, i, gsl_vector_get (l, i));
			gsl_vector_set (x1, i, gsl_vector_get (x, i));
			gsl_vector_set (y1, i, gsl_vector_get (y, i));
		}
		gsl_vector_free (l);
		gsl_vector_free (x);
		gsl_vector_free (y);
		if (_l) *_l = l1;
		else gsl_vector_free (l1);
		if (_x) *_x = x1;
		else gsl_vector_free (x1);
		if (_y) *_y = y1;
		else gsl_vector_free (y1);
	} else {
		if (_l) *_l = l;
		else gsl_vector_free (l);
		if (_x) *_x = x;
		else gsl_vector_free (x);
		if (_y) *_y = y;
		else gsl_vector_free (y);
	}
	return;
}

bool
bspline_eval_coeffs (gsl_vector *l, gsl_vector *x, gsl_vector *y, gsl_bspline_workspace	*bw,
	gsl_vector **_ex, gsl_vector **_ey)
{
	gsl_matrix	*a;
	gsl_vector	*ex;
	gsl_vector	*ey;

	a = bspline_fit_matrix (l, bw);

	ex = bspline_smooth_eval_coeffs (beta_val, a, x);	   
	ey = bspline_smooth_eval_coeffs (beta_val, a, y);	   
	gsl_matrix_free (a);

	if (_ex) *_ex = ex;
	else gsl_vector_free (ex);
	if (_ey) *_ey = ey;
	else gsl_vector_free (ey);

	return true;
}

#ifndef GSL_MAJOR_VERSION
bool
bspline_eval_derivs (gsl_vector *l, gsl_vector *ex, gsl_vector *ey,
	gsl_bspline_workspace *bw, gsl_bspline_deriv_workspace *dw,
	gsl_vector **_dx, gsl_vector **_dy, gsl_vector **_ddx, gsl_vector **_ddy)
#else
bool
bspline_eval_derivs (gsl_vector *l, gsl_vector *ex, gsl_vector *ey, gsl_bspline_workspace *bw,
	gsl_vector **_dx, gsl_vector **_dy, gsl_vector **_ddx, gsl_vector **_ddy)
#endif
{
	gsl_matrix	*a;
	gsl_vector	*x1;
	gsl_vector	*y1;
	gsl_vector	*dx1;
	gsl_vector	*dy1;
	gsl_vector	*ddx1;
	gsl_vector	*ddy1;

	FILE		*fp;

	a = bspline_fit_matrix (l, bw);
	x1 = bspline_synthesize (a, ex);
	y1 = bspline_synthesize (a, ey);
	gsl_matrix_free (a);

#ifndef GSL_MAJOR_VERSION
	a = bspline_nderiv_matrix (1, l, bw, dw);
#else
	a = bspline_nderiv_matrix (1, l, bw);
#endif
	dx1 = bspline_synthesize (a, ex);
	dy1 = bspline_synthesize (a, ey);
	gsl_matrix_free (a);

#ifndef GSL_MAJOR_VERSION
	a = bspline_nderiv_matrix (2, l, bw, dw);
#else
	a = bspline_nderiv_matrix (2, l, bw);
#endif
	ddx1 = bspline_synthesize (a, ex);
	ddy1 = bspline_synthesize (a, ey);
	gsl_matrix_free (a);

	fp = fopen ("lcurve_splined.data", "w");
	if (fp) {
		fprintf_synth_vector (fp, l, x1, y1, dx1, dy1, ddx1, ddy1);
		fclose (fp);
	}
	gsl_vector_free (x1);
	gsl_vector_free (y1);

	if (_dx) *_dx = dx1;
	else gsl_vector_free (dx1);
	if (_dy) *_dy = dy1;
	else gsl_vector_free (dy1);
	if (_ddx) *_ddx = ddx1;
	else gsl_vector_free (ddx1);
	if (_ddy) *_ddy = ddy1;
	else gsl_vector_free (ddy1);

	return true;
}

gsl_vector *
eval_curvature (const gsl_vector *dx, const gsl_vector *dy, const gsl_vector *ddx, const gsl_vector *ddy)
{
	int			i;
	gsl_vector	*cv = gsl_vector_alloc (dx->size - 4);
	for (i = 0; i < cv->size; i++) {
		double	dxi = gsl_vector_get (dx, i + 2);
		double	dyi = gsl_vector_get (dy, i + 2);
		double	ddxi = gsl_vector_get (ddx, i + 2);
		double	ddyi = gsl_vector_get (ddy, i + 2);
		double	c1 = dxi * ddyi - dyi * ddxi;
		double	c2 = pow (dxi, 2.) + pow (dyi, 2.);
		gsl_vector_set (cv, i, c1 / pow (c2, 1.5));
	}
	return cv;
}

typedef struct {
	double	l;
	double	cv;
	double	y;
} cv_data;

static int
compare_cv (const void *_a, const void *_b)
{
	cv_data	*a = (cv_data *) _a;
	cv_data	*b = (cv_data *) _b;

	if (a->l < b->l) return -1;
	if (a->l > b->l) return 1;
	return 0;
}

bool
bspline_curvature (gsl_vector *l, gsl_vector *x, gsl_vector *y)
{
	int				k;
	size_t			nbreaks;

	double			lmax;
	double			lmin;

	gsl_vector		*ex;
	gsl_vector		*ey;

	gsl_vector		*dx;
	gsl_vector		*dy;
	gsl_vector		*ddx;
	gsl_vector		*ddy;

	gsl_vector		*cv;
	FILE			*fp;

	gsl_bspline_workspace	*bw;
#ifndef GSL_MAJOR_VERSION
	gsl_bspline_deriv_workspace	*dw;
#endif
	/* 3次Bスプライン関数(ブレイクポイント数 = 8) */
	k = 4;
	nbreaks = 16;
	bw = gsl_bspline_alloc (k, nbreaks);
#ifndef GSL_MAJOR_VERSION
	dw = gsl_bspline_deriv_alloc (k);
#endif
	lmin = gsl_vector_min (l);
	lmax = gsl_vector_max (l);
	gsl_bspline_knots_uniform (lmin, lmax, bw);

	/* Bスプラインの係数決定 */
	bspline_eval_coeffs (l, x, y, bw, &ex, &ey);

	/* 微係数の推定 */
#ifndef GSL_MAJOR_VERSION
	bspline_eval_derivs (l, ex, ey, bw, dw, &dx, &dy, &ddx, &ddy);
#else
	bspline_eval_derivs (l, ex, ey, bw, &dx, &dy, &ddx, &ddy);
#endif
	gsl_vector_free (ex);
	gsl_vector_free (ey);
	gsl_bspline_free (bw);
#ifndef GSL_MAJOR_VERSION
	gsl_bspline_deriv_free (dw);
#endif

	cv = eval_curvature (dx, dy, ddx, ddy);
	gsl_vector_free (dx);
	gsl_vector_free (dy);
	gsl_vector_free (ddx);
	gsl_vector_free (ddy);

	// sort results
	fp = fopen ("cv.data", "w");
	if (fp) {
		int		i;
		int		n = cv->size;
		cv_data	cdata[n];

		for (i = 0; i < n; i++) {
			cdata[i].l = pow (10., gsl_vector_get (l, i + 2));
			cdata[i].cv = gsl_vector_get (cv, i);
			cdata[i].y = pow (10, gsl_vector_get (y, i + 2));
		}
		qsort ((void *) &cdata, n, sizeof (cv_data), compare_cv);

		// do not output first and last elements of cv
		// because on the edge of the range,
		// estimation of curvature is likely to unstable.
		for (i = 1; i < n - 1; i++) fprintf (fp, "%.4e\t%.4e\t%.4e\n", cdata[i].l, cdata[i].cv, cdata[i].y);

	}

	/*
	if (fp) {
		fprintf_curvature (fp, l, cv, y);
		fclose (fp);
	}
	*/
	gsl_vector_free (cv);

	return true;
}

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
	fprintf (stderr, "This program performs the interpolation of the discrete L-curve\n");
	fprintf (stderr, "obtained by \"l1l2inv\" program,\n");
	fprintf (stderr, "and estimates the curvature of the interpolated L-curve,\n");
	fprintf (stderr, "by reading a file which stores the regression infomations.\n");
	fprintf (stderr, "(e.g., L1 norm, L2 norm, RSS, etc. for the models\n");
	fprintf (stderr, "derived by \"l1l2inv\" with sequential reqularization parameters.)\n");
	fprintf (stderr, "NOTE: this program uses GSL(Gnu Scientific Library) version 2.0 or later.\n\n");

	fprintf (stderr, "OUTPUT:\n");
	fprintf (stderr, "[FILE]: lcurve_splined.data\n");
	fprintf (stderr, "        L-curve interpolated by smooth B-spline.\n");
	fprintf (stderr, "        format is <lambda> <RSS> <magnitude of penalty>\n");
	fprintf (stderr, "        and first and second derivatives of <RSS> and <mag. of penalty>.\n"); 
	fprintf (stderr, "[FILE]: cv.data\n");
	fprintf (stderr, "        curvature of the interpolated L-curve.\n\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -f <input file name>\n");
	fprintf (stderr, "       -a <alpha>\n");
	fprintf (stderr, "       -b <beta: dumping for bspline>\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	return;
}

bool
read_input_params (int argc, char **argv)
{
	char	c;

	while ((c = getopt (argc, argv, "f:a:b:h")) != EOF) {
		switch (c) {
			case 'f':
				infile_specified = true;
				strcpy (fn, optarg);
				break;
			case 'a':
				alpha = (double) atof (optarg);
				break;
			case 'b':
				beta_val = (double) atof (optarg);
				break;
			case 'h':
				return false;
		}
	}
	if (!infile_specified) {
		fprintf (stderr, "ERROR: input file not specified.\n");
		return false;
	}
	return true;
}

int
main (int argc, char **argv)
{
	gsl_vector		*l;
	gsl_vector		*x;
	gsl_vector		*y;

	if (!read_input_params (argc, argv)) {
		usage (argv[0]);
		return EXIT_FAILURE;
	}

	/* データ読み込み */
	read_data (fn, &l, &x, &y);

	bspline_curvature (l, x, y);

	gsl_vector_free (l);
	gsl_vector_free (x);
	gsl_vector_free (y);

	return 0;
}
