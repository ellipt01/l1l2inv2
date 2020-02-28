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
#include "settings.h"
#include "extern.h"

extern bool	create_xmat;
int			num_xfiles;
int			xfile_len;

int
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

static void
fprintf_vectors (linregmodel *lreg)
{
	FILE	*fp;

	fp = fopen ("y.data", "w");
	if (fp) {
		int		i;
		for (i = 0; i < lreg->m; i++) fprintf (fp, "%f\n", lreg->y->data[i]);
		fclose (fp);
	}

	fp = fopen ("xtx.data", "w");
	if (fp) {
		int		j;
		for (j = 0; j < lreg->n; j++) fprintf (fp, "%f\n", lreg->xtx[j]);
		fclose (fp);
	}

	return;
}

bool
l1l2inv (simeq *eq, char *path_fn, char *info_fn)
{
	int					m, n;
	linregmodel			*lreg;
	cdescent			*cd;

	m = eq->y->m;
	n = ngrd[0] * ngrd[1] * ngrd[2];
	num_xfiles = ngrd[2];
	xfile_len = ngrd[0] * ngrd[1];

	if (verbose) fprintf (stderr, "preparing linregmodel object... ");
	lreg = linregmodel_new (m, n, num_xfiles, eq->y, eq->d);

	if (type == TYPE_L1D1) {
		int		j;
		for (j = 0; j < lreg->d->n; j++)
			mm_real_xj_scale (lreg->d, j, 1. / sqrt (lreg->xtx[j]));
	}

	if (verbose) fprintf (stderr, "done\n");
	if (output_vector) fprintf_vectors (lreg);

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
	fprintf (stderr, "[m, n] = [%d, %d]\n", cd->lreg->m, cd->lreg->n);

	if (stochastic) {
		time_t	t = time (NULL);
		cdescent_set_stochastic (cd, (unsigned int *) &t);
	}

	cdescent_not_use_intercept (cd);
	if (constraint) cdescent_set_constraint (cd, l1l2inv_constraint_func);
	if (output_weighted) cd->output_rescaled = false;
	if (verbose) cd->verbose = true;

	cdescent_set_outputs_fullpath (cd, path_fn);
	cdescent_set_outputs_info (cd, info_fn);
	cdescent_set_log10_lambda_lower (cd, log10_lambda_lower);
	cdescent_set_log10_dlambda (cd, log10_dlambda);
	if (use_log10_lambda_upper) cdescent_set_log10_lambda_upper (cd, log10_lambda_upper);

	fprintf (stderr, "regression start\n");
	if (!cdescent_do_pathwise_optimization (cd)) fprintf (stderr, "not converged.\n");
	fprintf (stderr, "total num of iter = %d\n", cd->total_iter);
	if (cd->use_intercept) fprintf (stderr, "intercept = %.4e\n", cd->b0);

	cdescent_free (cd);
	linregmodel_free (lreg);

	return true;
}

