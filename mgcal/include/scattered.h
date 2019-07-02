#ifndef SCATTERED_H
#define SCATTERED_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_scattered	scattered;

struct s_scattered {
	int			n;

	double		*x;
	double		*y;
	double		*z;

	double		*dx;
	double		*dy;
	double		*dz;

	void		*data;
};

scattered	*scattered_new (const int n, const double x[], const double y[], const double z[]);
scattered	*scattered_new_full (const int n, const double x[], const double y[], const double z[], const double *dx, const double *dy, const double *dz);
void		scattered_free (scattered *s);
void		scattered_get_nth (const scattered *g, const int n, vector3d *center, vector3d *dim);

#ifdef __cplusplus
}
#endif

#endif /* SCATTERED_H */
