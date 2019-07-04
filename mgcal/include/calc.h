/*
 * calc.h
 *
 *  Created on: 2015/03/14
 *      Author: utsugi
 */

#ifndef CALC_H_
#define CALC_H_

#ifdef __cplusplus
extern "C" {
#endif

vector3d	*dipole (const vector3d *obs, const source *s);
vector3d	*prism (const vector3d *obs, const source *s);
double		total_force_dipole (const vector3d *obs, const source *src, void *data);
double		total_force_prism (const vector3d *obs, const source *src, void *data);

double		dipole_tf (const vector3d *obs, const source *s);
double		prism_tf (const vector3d *obs, const source *s);
double		total_force_dipole_bh (const vector3d *obs, const source *src, void *data);
double		total_force_prism_bh (const vector3d *obs, const source *src, void *data);

#ifdef __cplusplus
}
#endif

#endif /* CALC_H_ */
