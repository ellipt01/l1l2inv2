#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "mgcal.h"
#include "cdescent.h"

#include "simeq_xmat.h"
#include "smooth.h"
#include "consts.h"

simeq *
simeq_new (void)
{
	simeq	*eq = (simeq *) malloc (sizeof (simeq));
	eq->y = NULL;
	eq->d = NULL;
	return eq;
}

void
simeq_free (simeq *eq)
{
	if (eq) {
		if (eq->y) mm_real_free (eq->y);
		if (eq->d) mm_real_free (eq->d);
	}
	return;
}

static mm_dense *
create_observation (const data_array *array)
{
	int			i;
	int			m = array->n;
	mm_dense	*y = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, 1, m);
	for (i = 0; i < m; i++) y->data[i] = array->data[i];
	return y;
}

static void
create_kernel_matrix_xmatfile (const double inc, const double dec,
	const data_array *array, const grid *gsrc, const mgcal_func *func)
{
	int			i, l;
	int			m = array->n;
	int			n = gsrc->n;
	int			nnz = m * n;

	int			num_xfiles = ngrd[2];
	int			xfile_len = ngrd[0] * ngrd[1];

	vector3d	*e;

	mm_real		*w = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n);
	mm_real		*s = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, n, 1, n);

	bool		nextfile;
	FILE		*fp = NULL;

	e = vector3d_new_with_geodesic_poler (1., inc, dec);

#pragma omp parallel for private (i)
	for (l = 0; l < num_xfiles; l++) {
		int		k;
		char	fn[80];
		FILE	*fp_xmat;

		mm_real	*xj = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, 1, m);

		vector3d	*obs = vector3d_new (0., 0., 0.);
		source		*src = source_new (0., 0.);
		src->exf = vector3d_copy (e);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (0., 0., 0.);
		src->begin->mgz = vector3d_copy (e);

		sprintf (fn, "x%03d.mat", l);
		fp_xmat = fopen (fn, "wab");
		if (!fp_xmat) {
			fprintf (stderr, "ERROR: cannot open file %s.\nprogram abort!\n", fn);
			exit (1);
		}
		for (k = 0; k < xfile_len; k++) {
			int		j = k + l * xfile_len;
			grid_get_nth (gsrc, j, src->begin->pos, src->begin->dim);
			for (i = 0; i < m; i++) {
				vector3d_set (obs, array->x[i], array->y[i], array->z[i]);
				xj->data[i] = func->function (obs, src, func->parameter);
			}
			s->data[j] = mm_real_xj_sum (xj, 0);
			w->data[j] = mm_real_xj_ssq (xj, 0);
			// normalize
			mm_real_xj_scale (xj, 0, 1. / sqrt (w->data[j]));
			fwrite (xj->data, sizeof (double), m, fp_xmat);
		}
		fclose (fp_xmat);
		mm_real_free (xj);
		vector3d_free (obs);
		source_free (src);
	}

	fp = fopen("xs.vec", "wb");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file xs.vec.\nprogram abort!\n");
		exit (1);
	}
	fwrite (s->data, sizeof (double), n, fp);
	fclose (fp);
	mm_real_free (s);

	fp = fopen("xtx.vec", "wb");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file xtx.vec\nprogram abort!\n");
		exit (1);
	}
	fwrite (w->data, sizeof (double), n, fp);
	fclose (fp);
	mm_real_free (w);

	vector3d_free (e);

	return;
}

simeq *
create_simeq (const int type, const double inc, const double dec,
	const data_array *array, const grid *gsrc, const mgcal_func *func, bool create_xmat)
{
	simeq	*eq;
	FILE	*fp;

	eq = simeq_new ();
	if (array) eq->y = create_observation (array);
	if (create_xmat) create_kernel_matrix_xmatfile (inc, dec, array, gsrc, func);

	switch (type) {
		case TYPE_L1L2:
			eq->d = mm_real_eye (MM_REAL_SPARSE, gsrc->n);
			fprintf (stderr, "d = eye\n");
			break;
		case TYPE_L1D1:
			eq->d = mm_real_smooth_1 (MM_REAL_SPARSE, gsrc->nx, gsrc->ny, gsrc->nz);
			fprintf (stderr, "d = smooth_L1\n");
			break;
		case TYPE_L1D01:
			eq->d = mm_real_smooth_l01 (MM_REAL_SPARSE, gsrc->nx, gsrc->ny, gsrc->nz);
			fprintf (stderr, "d = smooth_L0_L1\n");
			break;
		default:
			break;
	}

	return eq;
}

