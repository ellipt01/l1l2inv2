/*
 * scattered_scattered.c
 *
 *  Created on: 2016/08/31
 *      Author: utsugi
 */

#include <stdlib.h>
#include <stdbool.h>

#include "../include/vector3d.h"
#include "private/util.h"
#include "scattered.h"

static scattered *
scattered_alloc (void)
{
	scattered	*g = (scattered *) malloc (sizeof (scattered));
	g->n = 0;

	g->x = NULL;
	g->y = NULL;
	g->z = NULL;

	g->dx = NULL;
	g->dy = NULL;
	g->dz = NULL;

	g->data = NULL;

	return g;
}

static scattered *
scattered_new_0 (const int n, const double x[], const double y[], const double z[], const double *dx, const double *dy, const double *dz)
{
	scattered	*g = scattered_alloc ();

	g->n = n;

	if (dx) {
		g->dx = (double *) malloc (n * sizeof (double));
		array_copy (n, g->dx, dx);
	}
	g->x = (double *) malloc (n * sizeof (double));
	array_copy (n, g->x, x);

	if (dy) {
		g->dy = (double *) malloc (n * sizeof (double));
		array_copy (n, g->dy, dy);
	}
	g->y = (double *) malloc (n * sizeof (double));
	array_copy (n, g->y, y);

	if (dz) {
		g->dz = (double *) malloc (n * sizeof (double));
		array_copy (n, g->dz, dz);
	}
	g->z = (double *) malloc (n * sizeof (double));
	array_copy (n, g->z, z);

	return g;
}

scattered *
scattered_new (const int n, const double x[], const double y[], const double z[])
{
	scattered	*g;
	if (n <= 0) error_and_exit_mgcal ("scattered_new", "n must be >= 1.", __FILE__, __LINE__);
	if (!x || !y || !z) error_and_exit_mgcal ("scattered_new", "x, y and z must be not empty.", __FILE__, __LINE__);
	g = scattered_new_0 (n, x, y, z, NULL, NULL, NULL);
	if (!g) error_and_exit_mgcal ("scattered_new", "failed to create scattered object.", __FILE__, __LINE__);
	return g;
}

scattered *
scattered_new_full (const int n, const double x[], const double y[], const double z[], const double *dx, const double *dy, const double *dz)
{
	scattered	*g;
	if (n <= 0) error_and_exit_mgcal ("scattered_new_full", "n must be >= 1.", __FILE__, __LINE__);
	if (!x || !y || !z) error_and_exit_mgcal ("scattered_new_full", "x, y and z must be not empty.", __FILE__, __LINE__);
	g = scattered_new_0 (n, x, y, z, dx, dy, dz);
	if (!g) error_and_exit_mgcal ("scattered_new_full", "failed to create scattered object.", __FILE__, __LINE__);
	return g;
}

void
scattered_free (scattered *g)
{
	if (g) {
		if (g->x) free (g->x);
		if (g->y) free (g->y);
		if (g->z) free (g->z);
		if (g->dx) free (g->dx);
		if (g->dy) free (g->dy);
		if (g->dz) free (g->dz);
		free (g);
	}
	return;
}

void
scattered_get_nth (const scattered *g, const int n, vector3d *center, vector3d *dim)
{
	if (!g) error_and_exit_mgcal ("scattered_get_nth", "scattered *g is empty.", __FILE__, __LINE__);
	if (n < 0 || n >= g->n) error_and_exit_mgcal ("scattered_get_nth", "index invalid.", __FILE__, __LINE__);
	if (center) vector3d_set (center, g->x[n], g->y[n], g->z[n]);
	if (dim && g->dx && g->dy && g->dz) vector3d_set (dim, g->dx[n], g->dy[n], g->dz[n]);
	return;
}
