#ifndef _const_H_
#define _const_H_

/* shift the grid coordinates by half of a grid */
bool	shift_grid;

/* declination and inclination of external field */
double	exf_dec;
double	exf_inc;

/* grid setting */
int		ngrd[3];
double	xgrd[2];
double	ygrd[2];
double	zgrd[2];

bool	stretch_grid_at_edge;
bool	use_dz_array;
double	*dz;

#endif /* _const_H_ */

