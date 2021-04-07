#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <mgcal.h>
#include <cdescent.h>

#include "smooth.h"
#include "simeq.h"

simeq *
simeq_new (void)
{
	simeq	*eq = (simeq *) malloc (sizeof (simeq));
	eq->x = NULL;
	eq->y = NULL;
	eq->d = NULL;
	return eq;
}

void
simeq_free (simeq *eq)
{
	if (eq) {
		if (eq->x) mm_real_free (eq->x);
		if (eq->y) mm_real_free (eq->y);
		if (eq->d) mm_real_free (eq->d);
	}
	return;
}

static void
centering (mm_dense *x)
{
	int		j;
#pragma omp parallel for
	for (j = 0; j < x->n; j++) {
		double	xm;
		xm = mm_real_xj_sum (x, j) / (double) x->m;
		mm_real_xj_add_const (x, j, - xm);
	}	
	return;
}

static void
normalizing (mm_dense *x)
{
	int		j;
#pragma omp parallel for
	for (j = 0; j < x->n; j++) {
		double	xnrm;
		xnrm = mm_real_xj_nrm2 (x, j);
		mm_real_xj_scale (x, j, 1. / xnrm);
	}	
	return;
}

void
simeq_centering_y (simeq *eq)
{
	centering (eq->y);
	return;
}

void
simeq_centering_x (simeq *eq)
{
	centering (eq->x);
	return;
}

void
simeq_normalizing_x (simeq *eq)
{
	normalizing (eq->x);
	return;
}

void
simeq_standardizing_x (simeq *eq)
{
	centering (eq->x);
	normalizing (eq->x);
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

static mm_dense *
create_kernel_matrix_dense (const double exf_inc, const double exf_dec,
	const double mag_inc, const double mag_dec,
	const data_array *array, const grid *gsrc, const mgcal_func *func)
{
	int			i;
	int			m = array->n;
	int			n = gsrc->n;
	int			nnz = m * n;
	vector3d	*exf;
	vector3d	*mag;
	mm_dense	*a = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, nnz);

	exf = vector3d_new_with_geodesic_poler (1., exf_inc, exf_dec);
	mag = vector3d_new_with_geodesic_poler (1., mag_inc, mag_dec);
	kernel_matrix_set (a->data, array, gsrc, mag, exf, func);
	vector3d_free (exf);
	vector3d_free (mag);

	return a;
}

static double
hdist (const vector3d *pos0, const vector3d *pos1)
{
	return sqrt (pow (pos0->x - pos1->x, 2.) + pow (pos0->y - pos1->y, 2.));
}

// obsolete
static mm_sparse *
create_kernel_matrix_sparse (const double exf_inc, const double exf_dec,
	const double mag_inc, const double mag_dec,
	const data_array *array, const grid *gsrc, const mgcal_func *func, const double rthres)
{
	int			j;
	int			m = array->n;
	int			n = gsrc->n;
	int			nnz = m * n;
	vector3d	*exf;
	vector3d	*mag;
	mm_dense	*a = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);

	exf = vector3d_new_with_geodesic_poler (1., exf_inc, exf_dec);
	mag = vector3d_new_with_geodesic_poler (1., mag_inc, mag_dec);

#pragma omp parallel
	{
		int			i, k;
		vector3d	*obs = vector3d_new (0., 0., 0.);
		source		*src = source_new (0., 0.);
		src->exf = vector3d_copy (exf);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (0., 0., 0.);
		src->begin->mgz = vector3d_copy (mag);

		k = 0;
#pragma omp for
		for (j = 0; j < n; j++) {
			grid_get_nth (gsrc, j, src->begin->pos, src->begin->dim);

			for (i = 0; i < m; i++) {
				vector3d_set (obs, array->x[i], array->y[i], array->z[i]);

				if (hdist (obs, src->begin->pos) <= rthres) {
					a->data[k] = func->function (obs, src, func->parameter);
					a->i[k] = i;
					k++;
				}
			}
			a->p[j + 1] = k;
		}

		mm_real_realloc (a, k);

		vector3d_free (obs);
		source_free (src);
	}			

	vector3d_free (exf);
	vector3d_free (mag);

	return a;
}

simeq *
create_simeq (const int type, const double exf_inc, const double exf_dec,
	const double mag_inc, const double mag_dec,
	const data_array *array, const grid *gsrc, const mgcal_func *func, const double *w)
{
	simeq	*eq;
	FILE	*fp;

	eq = simeq_new ();
	if (array) eq->y = create_observation (array);

	eq->x = create_kernel_matrix_dense (exf_inc, exf_dec, mag_inc, mag_dec, array, gsrc, func);

	switch (type) {
		case TYPE_L1L2:
			eq->d = mm_real_eye (MM_REAL_SPARSE, eq->x->n);
			break;
		case TYPE_L1TSV:
			eq->d = mm_real_smooth_1 (MM_REAL_SPARSE, gsrc->nx, gsrc->ny, gsrc->nz, w);
			break;
		case TYPE_L1L2TSV:
			eq->d = mm_real_smooth_l01_1 (MM_REAL_SPARSE, gsrc->nx, gsrc->ny, gsrc->nz, w);
			break;
		default:
			break;
	}

	return eq;
}

