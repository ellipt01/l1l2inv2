/*
 * kernel.c
 *
 *  Created on: 2015/03/15
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../include/vector3d.h"
#include "source.h"
#include "data_array.h"
#include "grid.h"
#include "scattered.h"
#include "kernel.h"
#include "private/util.h"

static mgcal_func *
mgcal_func_alloc (void)
{
	mgcal_func	*f = (mgcal_func *) malloc (sizeof (mgcal_func));
	f->function = NULL;
	f->parameter = NULL;
	return f;
}

mgcal_func *
mgcal_func_new (const mgcal_theoretical func, void *data)
{
	mgcal_func	*f = mgcal_func_alloc ();
	f->function = func;
	f->parameter = data;
	return f;
}

void
mgcal_func_free (mgcal_func *f)
{
	if (f) free (f);
	return;
}

void
kernel_matrix_set (double *a, const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		m;
	int		nx;
	int		ny;
	int		nz;
	int		nh;

	if (!a) error_and_exit_mgcal ("kernel_matrix_set", "double *a is empty.", __FILE__, __LINE__);

	m = array->n;
	nx = g->nx;
	ny = g->ny;
	nz = g->nz;
	nh = g->nh;

#pragma omp parallel
	{
		int			i, j, k, l;
		double		*z1 = NULL;
		vector3d	*obs = vector3d_new (0., 0., 0.);
		source		*src = source_new (0., 0.);
		if (exf) src->exf = vector3d_copy (exf);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (0., 0., 0.);
		if (mgz) src->begin->mgz = vector3d_copy (mgz);

#pragma omp for
		for (k = 0; k < nz; k++) {
			double	*zk = g->z + k;	// for parallel calculation
			double	*dzk = g->dz + k;	// for parallel calculation
			double	*yj = g->y;
			double	*dyj = g->dy;
			for (j = 0; j < ny; j++) {
				double	*xi = g->x;
				double	*dxi = g->dx;
				if (g->z1) z1 = g->z1 + j * nx;
				for (i = 0; i < nx; i++) {
					double	*al = a + (k * nh + j * nx + i) * m;
					double	*xl = array->x;
					double	*yl = array->y;
					double	*zl = array->z;
					double	z1k = *zk;
					if (z1) z1k += z1[i];
					vector3d_set (src->begin->pos, *xi, *yj, z1k);
					vector3d_set (src->begin->dim, *dxi, *dyj, *dzk);
					for (l = 0; l < m; l++) {
						vector3d_set (obs, *xl, *yl, *zl);
						*al = f->function (obs, src, f->parameter);
						al++;
						xl++;
						yl++;
						zl++;
					}
					xi++;
					dxi++;
				}
				yj++;
				dyj++;
			}
		}
		vector3d_free (obs);
		source_free (src);
	}
	return;
}

double *
kernel_matrix (const data_array *array, const grid *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		m, n;
	double	*a;

	m = array->n;
	n = g->n;
	a = (double *) malloc (m * n * sizeof (double));
	if (!a) error_and_exit_mgcal ("kernel_matrix", "failed to allocate memory of *a.", __FILE__, __LINE__);
	kernel_matrix_set (a, array, g, mgz, exf, f);
	return a;
}

void
kernel_matrix_scattered_set (double *a, const data_array *array, const scattered *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		m;
	int		n;

	if (!a) error_and_exit_mgcal ("kernel_matrix_set", "double *a is empty.", __FILE__, __LINE__);

	m = array->n;
	n = g->n;

#pragma omp parallel
	{
		int			i, j;
		vector3d	*obs = vector3d_new (0., 0., 0.);
		source		*src = source_new (0., 0.);
		if (exf) src->exf = vector3d_copy (exf);
		source_append_item (src);
		src->begin->pos = vector3d_new (0., 0., 0.);
		src->begin->dim = vector3d_new (1., 1., 1.);
		if (mgz) src->begin->mgz = vector3d_copy (mgz);

#pragma omp for
		for (j = 0; j < n; j++) {
			vector3d_set (src->begin->pos, g->x[j], g->y[j], g->z[j]);
			if (g->dx && g->dy && g->dz) vector3d_set (src->begin->dim, g->dx[j], g->dy[j], g->dz[j]);
			for (i = 0; i < m; i++) {
				vector3d_set (obs, array->x[i], array->y[i], array->z[i]);
				a[i + j * m] = f->function (obs, src, f->parameter);
			}
		}
		vector3d_free (obs);
		source_free (src);
	}
	return;
}

double *
kernel_matrix_scattered (const data_array *array, const scattered *g, const vector3d *mgz, const vector3d *exf, const mgcal_func *f)
{
	int		m, n;
	double	*a;

	m = array->n;
	n = g->n;
	a = (double *) malloc (m * n * sizeof (double));
	if (!a) error_and_exit_mgcal ("kernel_matrix_scattered", "failed to allocate memory of *a.", __FILE__, __LINE__);
	kernel_matrix_scattered_set (a, array, g, mgz, exf, f);
	return a;
}
