/*
 * calc.c
 *
 *  Created on: 2017/08/21
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "../include/vector3d.h"
#include "private/util.h"
#include "source.h"
#include "calc.h"

#define SIGN(a) ((a) < 0. ? -1. : +1.)

extern double scale_factor;

static double
dipole_tf_kernel (
	vector3d *exf, vector3d *mgz,
	double xobs, double yobs, double zobs, double xsrc, double ysrc, double zsrc
)
{
	/*
		Calculate magnetic total force due to a dipole
		according to eq.(5) after Bhattacharya (1964).

		In his paper, the coordinate system is (Fig. 1)

		x-axis(EW) directs to the East,
		y-axis(NS) directs to the South,
		z-axis(Up-Down) directs downward.

		However coordinates of this program is

		x-axis(EW) directs to the East(same),
		y-axis(NS) directs to the North,
		z-axis(Up-Down) directs upward.

		So, the sign of y and z coordinates, components have to be changed.
	*/

	// coordinates
	double	x = xsrc - xobs;
	double	y = - (ysrc - yobs);
	double	z = - (zsrc - zobs);

	// polarization vector 
	double	al = mgz->x;
	double	am = - mgz->y;
	double	an = - mgz->z;
	
	// geomagnetic field unit vector 
	double	bL = exf->x;
	double	bM = - exf->y;
	double	bN = - exf->z;

	double	x2 = pow (x, 2.);
	double	y2 = pow (y, 2.);
	double	z2 = pow (z, 2.);
	
	double	r = sqrt (x2 + y2 + z2);
	double	r3 = pow (r, 3.);
	double	r5 = pow (r, 5.);

	double	a1 = al * bL;
	double	a2 = am * bM;
	double	a3 = an * bN;
	double	a12 = bL * am + bM * al;
	double	a13 = bL * an + bN * al;
	double	a23 = bM * an + bN * am;

	double	ct = a1 + a2 + a3;

	double	ret = (
		- ct / r3
		+ 3. * (
			a1 * x2 + a2 * y2 + a3 * z2
			+ a12 * x * y + a13 * x * z + a23 * y * z
		) / r5
	);

	return ret;
}

static double
prism_tf_kernel (
	vector3d *exf, vector3d *mgz,
	double xobs, double yobs, double zobs, double xsrc, double ysrc, double zsrc
)
{
	/*
		Calculate magnetic total force due to a magnetized prism
		according to eq.(10) after Bhattacharya (1964).

		In his paper, the coordinate system is (Fig. 1)

		x-axis(EW) directs to the East,
		y-axis(NS) directs to the South,
		z-axis(Up-Down) directs downward.

		However coordinates of this program is

		x-axis(EW) directs to the East(same),
		y-axis(NS) directs to the North,
		z-axis(Up-Down) directs upward.

		So, the sign of y and z coordinates, components have to be changed.
	*/

	// coordinates
	double	x = xsrc - xobs;
	double	y = - (ysrc - yobs);
	double	z = - (zsrc - zobs);

	// polarization vector 
	double	al = mgz->x;
	double	am = - mgz->y;
	double	an = - mgz->z;
	
	// geomagnetic field unit vector 
	double	bL = exf->x;
	double	bM = - exf->y;
	double	bN = - exf->z;

	double	xy = x * y;
	double	x2 = pow (x, 2.);
	double	y2 = pow (y, 2.);
	double	z2 = pow (z, 2.);
	
	double	r2 = x2 + y2 + z2;
	double	r = sqrt (r2);

	double	a1 = al * bL;
	double	a2 = am * bM;
	double	a3 = an * bN;
	double	a12 = bL * am + bM * al;
	double	a13 = bL * an + bN * al;
	double	a23 = bM * an + bN * am;

	double	ct = a1 + a2 + a3;

	double	ret = (
		  0.5 * a23 * (log (r - x) - log (r + x))
		+ 0.5 * a13 * (log (r - y) - log (r + y))
		- a12 * log (r + z)
		- a1 * atan (xy / (x2 + r * z + z2))
		- a2 * atan (xy / (r2 + r * z - x2))
		+ a3 * atan (xy / (r * z))
	);

	return ret;
}

double
dipole_tf (const vector3d *obs, const source *s)
{
	double		xobs, yobs, zobs;
	double		xsrc, ysrc, zsrc;
	double		f;
	vector3d	*exf;
	source_item	*cur;

	if (!obs) error_and_exit_mgcal ("dipole", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit_mgcal ("dipole", "source *s is empty.", __FILE__, __LINE__);

	xobs = obs->x;
	yobs = obs->y;
	zobs = obs->z;

	exf = s->exf;

	f = 0.;
	cur = s->begin;
	while (cur) {
		double	dx, dy, dz;
		double	flag;

		if (!cur->pos) error_and_exit_mgcal ("dipolem", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit_mgcal ("dipole", "magnetization of source item is empty.", __FILE__, __LINE__);

		xsrc = cur->pos->x;
		ysrc = cur->pos->y;
		zsrc = cur->pos->z;

		f += dipole_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, xsrc, ysrc, zsrc);

		cur = cur->next;
	}
	f *= scale_factor;
	return f;
}

double
prism_tf (const vector3d *obs, const source *s)
{
	double		a[2], b[2], c[2];
	double		xobs, yobs, zobs;
	double		xsrc, ysrc, zsrc;
	double		f;
	vector3d	*exf;
	double		tmp[8];
	source_item	*cur;

	if (!obs) error_and_exit_mgcal ("prism", "vector3d *obs is empty.", __FILE__, __LINE__);
	if (!s) error_and_exit_mgcal ("prism", "source *s is empty.", __FILE__, __LINE__);

	xobs = obs->x;
	yobs = obs->y;
	zobs = obs->z;

	exf = s->exf;

	f = 0.;
	cur = s->begin;
	while (cur) {
		double	dx, dy, dz;
		double	flag;

		if (!cur->pos) error_and_exit_mgcal ("prism", "position of source item is empty.", __FILE__, __LINE__);
		if (!cur->dim) error_and_exit_mgcal ("prism", "dimension of source item is empty.", __FILE__, __LINE__);
		if (!cur->mgz) error_and_exit_mgcal ("prism", "magnetization of source item is empty.", __FILE__, __LINE__);

		dx = cur->dim->x;
		dy = cur->dim->y;
		dz = cur->dim->z;
		flag = SIGN (dx) * SIGN (dy) * SIGN (dz);

		xsrc = cur->pos->x;
		ysrc = cur->pos->y;
		zsrc = cur->pos->z;

		a[0] = xsrc - 0.5 * dx;
		b[0] = ysrc - 0.5 * dy;
		c[0] = zsrc - 0.5 * dz;

		a[1] = a[0] + dx;
		b[1] = b[0] + dy;
		c[1] = c[0] + dz;

		tmp[0] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[1], b[1], c[1]);
		tmp[2] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[1], b[0], c[1]);
		tmp[4] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[0], b[1], c[1]);
		tmp[6] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[0], b[0], c[1]);

		if (fabs (dz) < DBL_EPSILON) {
			tmp[1] = 0.;
			tmp[3] = 0.;
			tmp[5] = 0.;
			tmp[7] = 0.;
		} else {
			tmp[1] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[1], b[1], c[0]);
			tmp[3] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[1], b[0], c[0]);
			tmp[5] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[0], b[1], c[0]);
			tmp[7] = prism_tf_kernel (exf, cur->mgz, xobs, yobs, zobs, a[0], b[0], c[0]);
		}

		f += - flag * (
			+ tmp[0] - tmp[1] - tmp[2] + tmp[3]
			- tmp[4] + tmp[5] + tmp[6] - tmp[7]
		);

		cur = cur->next;
	}
	f *= scale_factor;
	return f;
}

double
total_force_dipole (const vector3d *obs, const source *src, void *data)
{
	return dipole_tf (obs, src);
}

double
total_force_prism (const vector3d *obs, const source *src, void *data)
{
	return prism_tf (obs, src);
}


