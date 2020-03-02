#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cdescent.h>
#include <mgcal.h>

#include "simeq_xmat.h"
#include "utils_xmat.h"
#include "_version_.h"

extern const double	exf_dec;
extern const double	exf_inc;

extern const int	ngrd[];
extern const double	xgrd[];
extern const double	ygrd[];
extern const double	zgrd[];
extern bool			stretch_grid_at_edge;
extern const bool	use_dz_array;
extern double		dz[];

double *
fread_z (FILE *fp, const int n)
{
	int		k;
	double	xval, yval, zval;
	double	*z = (double *) malloc (n * sizeof (double));

	k = 0;
	while (fscanf (fp, "%lf\t%lf\t%lf", &xval, &yval, &zval) != EOF) {
		z[k] = zval;
		if (++k > n) break;
	}
	return z;
}

/* skip one line:
   go ahead file descripter until fgetc reacheas '\n' */
static bool
fskip_one_line (FILE *fp)
{
	bool	read = false;
	char	c;
	while ((c = fgetc (fp)) != EOF) {
		if (c == '\n') break;
		read = true;
	}
	return read;
}


/* read data contained in n-th line and return it */
double *
get_nth_data (FILE *fp, const int n, int *len)
{
	int		i;
	int		l;
	double	*data;

	for (i = 0; i < n; i++) {
		if (!fskip_one_line (fp)) break;
	}
	if (i != n) return NULL;
	l = fread_ndata_of_one_line (fp);
	if (l <= 0) return NULL;
	data = fread_one_line (fp, l);
	if (len) *len = l;

	return data;
}

mm_dense *
extract_beta (const char *fn, const int j)
{
	int			i, m;
	double		*data;
	mm_dense	*beta;
	FILE		*fp;

	fp = fopen (fn, "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s.\n", fn);
		return NULL;
	}
	data = get_nth_data (fp, j, &m);
	fclose (fp);
	if (!data) return NULL;
	beta = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m - 2, 1, m - 2);
	for (i = 0; i < m - 2; i++) beta->data[i] = data[i + 2];
	free (data);

	return beta;
}

simeq *
read_input (const int type, const char *ifn, const char *tfn, const double *w, bool create_xmat)
{
	simeq		*eq;
	grid		*g;
	double		*z;
	data_array	*array;
	mgcal_func	*func;
	FILE		*fp;

	mgcal_theoretical	f;

	array = NULL;
	if (ifn) {
		/* read input data */
		if ((fp = fopen (ifn, "r")) == NULL) {
			fprintf (stderr, "ERROR: cannot open file %s.\n", ifn);
			return NULL;
		}
		array = fread_data_array (fp);
		fclose (fp);
		{
			FILE	*fp = fopen ("array.data", "w");
			fwrite_data_array (fp, array, NULL);
			fclose (fp);
		}
	}
	/* create simeq */
	if (!use_dz_array) {
		g = grid_new (ngrd[0], ngrd[1], ngrd[2], xgrd, ygrd, zgrd);
	} else {
		int		i;
		double	zgrd0[2];
		zgrd0[0] = zgrd[0];
		zgrd0[1] = zgrd[0];
		for (i = 0; i < ngrd[2]; i++) zgrd0[1] += dz[i];
		g = grid_new_full (ngrd[0], ngrd[1], ngrd[2], xgrd, ygrd, zgrd0, NULL, NULL, dz, NULL);
	}
	if (tfn && (fp = fopen (tfn, "r")) != NULL) {
		z = fread_z (fp, g->nh);
		grid_set_surface (g, z);
		free (z);
		fclose (fp);
	}

	if (stretch_grid_at_edge) {
		double	l = 1000.;
		g->x[0] -= l / 2.;
		g->x[g->nx - 1] += l / 2.;
		g->y[0] -= l / 2.;
		g->y[g->ny - 1] += l / 2.;
		g->z[g->nz - 1] -= l / 2.;

		g->dx[0] += l;
		g->dx[g->nx - 1] += l;
		g->dy[0] += l;
		g->dy[g->ny - 1] += l;
		g->dz[g->nz - 1] -= l;
	}

	{
		FILE	*fp;
		fp = fopen ("grid.data", "w");
		fwrite_grid (fp, g);
		fclose (fp);
		fp = fopen ("grid_xyz.data", "w");
		fwrite_grid_to_xyz (fp, g, "%f\t%f\t%f\t%f");
		fclose (fp);
	}

	f = total_force_prism;
	func = mgcal_func_new (f, NULL);
	eq = create_simeq (type, exf_inc, exf_dec, array, g, func, w, create_xmat);

	grid_free (g);
	data_array_free (array);
	mgcal_func_free (func);

	return eq;
}

/*********************************************
 *  入力データの作成
 *  ファイル model.par からパラメータを読み取り
 *********************************************/
source *
read_model_par (FILE *fp, const double exf_inc, const double exf_dec)
{
	int		i;
	source	*s;
	char	buf[BUFSIZ];

	double	inc, dec;

	s = source_new (exf_inc, exf_dec);
	i = 0;
	while (fgets (buf, BUFSIZ, fp) != NULL) {
		double	x, y, z, mgz, l, w, h;
		char	*p = buf;
		while (p[0] == ' ' || p[0] == '\t' || p[0] == '\r' || p[0] == '\n') p++;
		if (p[0] == '#') continue;
		if (strlen (p) <= 1) continue;
		if (i == 0) {
			sscanf (p, "%lf %lf", &inc, &dec);
			i++;
		} else {
			sscanf (p, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &l, &w, &h, &mgz);
			source_append_item (s);
			source_set_position (s, x, y, z);
			source_set_dimension (s, l, w, h);
			source_set_magnetization (s, mgz, inc, dec);
		}
	}

	return s;
}

void
printf_beta (const mm_dense *beta)
{
	int			i;
	grid		*g;
	FILE		*fp;

	fp = fopen ("beta_opt.data", "w");
	if (!fp) return;

	g = grid_new (ngrd[0], ngrd[1], ngrd[2], xgrd, ygrd, zgrd);
	fwrite_grid_with_data (fp, g, beta->data, NULL);
	grid_free (g);
	fclose (fp);
	return;
}

void
printf_mm_real (const char *fn, const mm_real *x, const char *format)
{
	FILE	*fp;
	/* output matrix to file */
	fp = fopen (fn, "w");
	if (fp) {
		mm_real_fwrite (fp, x, format);
		fclose (fp);
	}
	return;
}

void
printf_array (const char *fn, const int m, const int n, const double *a, const char *format)
{
	FILE	*fp;
	/* output matrix to file */
	fp = fopen (fn, "w");
	if (fp) {
		int		i;
		mm_real	*x = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, m * n);
		for (i = 0; i < x->m; i++) x->data[i] = a[i];
		mm_real_fwrite (fp, x, format);
		fclose (fp);
		mm_real_free (x);
	}
	return;
}

void
printf_estimated (const cdescent *cd, const mm_dense *beta)
{
	int			m;
	double		zobs;
	double		*z;
	grid		*g;
	mm_dense	*c;
	FILE		*fp;

	m = cd->lreg->m;
	c = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, 1, m);
	mm_real_set_all (c, 0.);
	if ((fp = fopen ("estimated.data", "w")) == NULL) return;
	mm_real_x_dot_yk_xmatfile (cd->lreg->fp_xmat, cd->lreg->m, cd->lreg->n, false, 1., beta, 0, 0., c);
	mm_real_fwrite (fp, c, NULL);
	fclose (fp);
	mm_real_free (c);
	return;
}

/* count num data which contained in one line
   count separatir (' ', '\t' or ',' and also '\n') of a line */
int
fread_ndata_of_one_line (FILE *fp)
{
	int		len;	
	char	c;
	fpos_t	pos;

	if (fgetpos (fp, &pos) != 0) return -1;

	len = 0;
	while ((c = fgetc (fp)) != EOF) {
		if (c == ' ' || c == '\t' || c == ',' || c == '\n') len++;
		if (c == '\n') break;
	}

	fsetpos (fp, &pos);

	return len;
}

/* read data contained in a line and return it */
double *
fread_one_line (FILE *fp, int len)
{
	char	c;
	int		i, k;
	char	tmp[BUFSIZ];
	double	*data;

	if (len <= 0) return NULL;
	data = (double *) malloc (len * sizeof (double));

	i = 0;
	k = 0;
	while ((c = fgetc (fp)) != EOF) {
		if (c == ' ' || c == '\t' || c == ',' || c == '\n') {
			if (k > 0) {
				tmp[k] = '\0';
				k = 0;
				data[i] = (double) atof (tmp);
				if (++i >= len) break;
			}
			if (c == '\n') break;
		} else tmp[k++] = c;
	}
	if (i <= 0) {
		free (data);
		data = NULL;
	}
	return data;
}

/* printout Z-file of the cross section
   int dir indicates the direction of the cross section (whether x or y). */
void
fprintf_zfile_2d (FILE *stream, grid *g, const int dir, const double *v,
	const double *t, const double *offset)
{
	int			k, n;
	int			nn = (dir == CROSS_SECTION_X) ? g->nx : g->ny;
	vector3d	*pos0 = vector3d_new (0., 0., 0.);
	vector3d	*pos1 = vector3d_new (0., 0., 0.);
	double		*p = (dir == CROSS_SECTION_X) ? &pos0->x : &pos0->y;
	vector3d	*dim0 = vector3d_new (0., 0., 0.);
	vector3d	*dim1 = vector3d_new (0., 0., 0.);
	double		*d = (dir == CROSS_SECTION_X) ? &dim0->x : &dim0->y;

	k = 0;
	for (n = 0; n < g->n; n++) {
		double	val = v[n];
		grid_get_nth (g, n, pos0, dim0);
		if (++k == nn) {
			vector3d_set (pos1, pos0->x, pos0->y, pos0->z);
			vector3d_set (dim1, dim0->x, dim0->y, dim0->z);
			k = 0;
		} else {
			grid_get_nth (g, n + 1, pos1, dim1);
		}
		if (t && fabs (val) < *t) {
			if (offset) val += *offset;
			else continue;
		}
		fprintf (stream, "> -Z %.8f\n", val);
		fprintf (stream, "%.4f\t%.4f\n", (*p) - 0.5 * (*d), pos0->z + 0.5 * dim0->z);
		fprintf (stream, "%.4f\t%.4f\n", (*p) - 0.5 * (*d), pos0->z - 0.5 * dim0->z);
		fprintf (stream, "%.4f\t%.4f\n", (*p) + 0.5 * (*d), pos1->z - 0.5 * dim1->z);
		fprintf (stream, "%.4f\t%.4f\n", (*p) + 0.5 * (*d), pos1->z + 0.5 * dim1->z);
	}
	vector3d_free (pos0);
	vector3d_free (pos1);
	vector3d_free (dim0);
	vector3d_free (dim1);
	return;
}

void
fprintf_zfile_3d (FILE *stream, grid *g, const double *v,
	const double *t, const double *offset)
{
	int			n;
	vector3d	*pos = vector3d_new (0., 0., 0.);
	vector3d	*dim = vector3d_new (0., 0., 0.);

	for (n = 0; n < g->n; n++) {
		double	val = v[n];
		if (t && fabs (val) < *t) {
			if (offset) val += *offset;
			else continue;
		}
		grid_get_nth (g, n, pos, dim);
		fprintf (stream, "> -Z %.8f\n", val);
		fprintf (stream, "%.4f\t%.4f\t%.4f\n", pos->x - 0.5 * dim->x, pos->y + 0.5 * dim->y, pos->z);
		fprintf (stream, "%.4f\t%.4f\t%.4f\n", pos->x - 0.5 * dim->x, pos->y - 0.5 * dim->y, pos->z);
		fprintf (stream, "%.4f\t%.4f\t%.4f\n", pos->x + 0.5 * dim->x, pos->y - 0.5 * dim->y, pos->z);
		fprintf (stream, "%.4f\t%.4f\t%.4f\n", pos->x + 0.5 * dim->x, pos->y + 0.5 * dim->y, pos->z);
	}
	vector3d_free (pos);
	vector3d_free (dim);
	return;
}

void
fprintf_sources (FILE *stream, const source *s)
{
	source_item	*i;
	for (i = s->item->next; i; i = i->next) {
		fprintf (stderr, "(%.4f, %.4f, %.4f), (%.4f, %.4f, %.4f), (%.4f, %.4f, %.4f)\n",
			i->pos->x, i->pos->y, i->pos->z,
			i->dim->x, i->dim->y, i->dim->z,
			i->mgz->x, i->mgz->y, i->mgz->z);
	}
	return;
}

void
version_info (const char *toolname)
{
	fprintf (stderr, "%s -- VERSION:%s, LAST_MODIFIED:%s\n",
		toolname, _VERSION_, _LAST_MODIFIED_);
	return;
}

