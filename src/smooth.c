#include <stdlib.h>
#include <cdescent.h>

#include <omp.h>

static void
col2ijk (const int nx, const int nh, const int col, int *i, int *j, int *k)
{
	int		h;
	*k = col / nh;
	h = col % nh;
	*j = h / nx;
	*i = h % nx;
	return;
}

static mm_sparse *
mm_real_smooth_sparse_x (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			n;
	int			nh;
	int			nnz;
	mm_sparse	*s;

	double		val = 1.;

	nh = nx * ny;
	m = (nx - 1) * ny * nz;
	n = nh * nz;
	nnz = 2 * m;
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);

#pragma omp parallel for
	for (col = 0; col < s->n; col++) {
		int		i, j, k, p;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = col - j - k * ny;

		p = 2 * index - 1;
		if (i == 0) p++;

		if (i > 0) {
			s->i[p] = index - 1;
			s->data[p++] = val;
		}
		if (i < nx - 1) {
			s->i[p] = index;
			s->data[p++] = - val;
		}
		s->p[col + 1] = p;
	}
	return s;
}

static mm_dense *
mm_real_smooth_dense_x (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			n;
	int			nh;
	int			nnz;
	mm_dense	*d;

	double		val = 1.;

	nh = nx * ny;
	m = (nx - 1) * ny * nz;
	n = nh * nz;
	nnz = m * n;
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, nnz);
	mm_real_set_all (d, 0.);

#pragma omp parallel for
	for (col = 0; col < d->n; col++) {
		int		i, j, k;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = col - j - k * ny;
		if (i > 0) d->data[index - 1 + col * d->m] = val;
		if (i < nx - 1) d->data[index + col * d->m] = - val;
	}
	return d;
}

static mm_real *
mm_real_smooth_x (MMRealFormat format, const int nx, const int ny, const int nz)
{
	return (format == MM_REAL_SPARSE) ? mm_real_smooth_sparse_x (nx, ny, nz) : mm_real_smooth_dense_x (nx, ny, nz);
}

static mm_sparse *
mm_real_smooth_sparse_y (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			n;
	int			nh;
	int			nnz;
	mm_sparse	*s;

	double		val = 1.;

	nh = nx * ny;
	m = nx * (ny - 1) * nz;
	n = nh * nz;
	nnz = 2 * m;
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);

#pragma omp parallel for
	for (col = 0; col < s->n; col++) {
		int		i, j, k, p;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = col - k * nx;

		p = 2 * index - nx;
		if (j == 0) p -= i - nx;
		else if (j == ny - 1) p -= i;

		if (j > 0) {
			s->i[p] = index - nx;
			s->data[p++] = val;
		}
		if (j < ny - 1) {
			s->i[p] = index;
			s->data[p++] = - val;
		}
		s->p[col + 1] = p;
	}
	return s;
}

static mm_dense *
mm_real_smooth_dense_y (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			n;
	int			nh;
	int			nnz;
	mm_dense	*d;

	double		val = 1.;

	nh = nx * ny;
	m = nx * (ny - 1) * nz;
	n = nh * nz;
	nnz = m * n;
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, nnz);
	mm_real_set_all (d, 0.);

#pragma omp parallel for
	for (col = 0; col < d->n; col++) {
		int		i, j, k;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = col - k * nx;
		if (j > 0) d->data[index - nx + col * d->m] = val;
		if (j < ny - 1) d->data[index + col * d->m] = - val;
	}
	return d;
}

static mm_real *
mm_real_smooth_y (MMRealFormat format, const int nx, const int ny, const int nz)
{
	return (format == MM_REAL_SPARSE) ? mm_real_smooth_sparse_y (nx, ny, nz) : mm_real_smooth_dense_y (nx, ny, nz);
}

static mm_sparse *
mm_real_smooth_sparse_z (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			n;
	int			nnz;
	int			nh;
	mm_sparse	*s;

	double		val = 1.;

	nh = nx * ny;
	m = nh * (nz - 1);
	n = nh * nz;
	nnz = 2 * m;
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, n, nnz);

#pragma omp parallel for
	for (col = 0; col < s->n; col++) {
		int		i, j, k, p;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = col;

		p = 2 * index - nh;
		if (k == 0) p -= i + j * nx - nh;
		else if (k == nz - 1) p -= i + j * nx;

		if (k > 0) {
			s->i[p] = index - nh;
			s->data[p++] = val;
		}
		if (k + 1 < nz) {
			s->i[p] = index;
			s->data[p++] = - val;
		}
		s->p[col + 1] = p;
	}
	return s;
}

static mm_dense *
mm_real_smooth_dense_z (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			n;
	int			nh;
	int			nnz;
	mm_dense	*d;

	double		val = 1.;

	nh = nx * ny;
	m = nh * (nz - 1);
	n = nh * nz;
	nnz = m * n;
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, n, nnz);
	mm_real_set_all (d, 0.);

#pragma omp parallel for
	for (col = 0; col < d->n; col++) {
		int		i, j, k;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = col;
		if (k > 0) d->data[index - nh + col * d->m] = val;
		if (k + 1 < nz) d->data[index + col * d->m] = - val;
	}
	return d;
}

static mm_real *
mm_real_smooth_z (MMRealFormat format, const int nx, const int ny, const int nz)
{
	return (format == MM_REAL_SPARSE) ? mm_real_smooth_sparse_z (nx, ny, nz) : mm_real_smooth_dense_z (nx, ny, nz);
}

mm_real *
mm_real_smooth (MMRealFormat format, const int nx, const int ny, const int nz)
{
	mm_real	*dx;
	mm_real	*dy;
	mm_real	*dz;
	mm_real	*dxy;
	mm_real	*d;

	/* [dx; dy] */
	dx = mm_real_smooth_x (format, nx, ny, nz);
	dy = mm_real_smooth_y (format, nx, ny, nz);
	dxy = mm_real_vertcat (dx, dy);
	mm_real_free (dx);
	mm_real_free (dy);

	/* [dx; dy; dz] */
	dz = mm_real_smooth_z (format, nx, ny, nz);
	d = mm_real_vertcat (dxy, dz);
	mm_real_free (dxy);
	mm_real_free (dz);

	return d;
}

/*
	D = [ 
			E
			DX
			DY
			DZ
		]
*/
mm_real *
mm_real_smooth_l01 (MMRealFormat format, const int nx, const int ny, const int nz, const double *w)
{
	mm_real	*d0;
	mm_real	*dx;
	mm_real	*dy;
	mm_real	*dz;
	mm_real	*d0x;
	mm_real	*d0xy;
	mm_real	*d;

	d0 = mm_real_eye (MM_REAL_SPARSE, nx * ny * nz);
	dx = mm_real_smooth_x (format, nx, ny, nz);
	dy = mm_real_smooth_y (format, nx, ny, nz);
	dz = mm_real_smooth_z (format, nx, ny, nz);

	if (w) {
		int		j;
		for (j = 0; j < d0->n; j++) mm_real_xj_scale (d0, j, w[0]);
		for (j = 0; j < dx->n; j++) mm_real_xj_scale (dx, j, w[1]);
		for (j = 0; j < dy->n; j++) mm_real_xj_scale (dy, j, w[2]);
		for (j = 0; j < dz->n; j++) mm_real_xj_scale (dz, j, w[3]);

	}

	/* [E; dx] */
	d0x = mm_real_vertcat (d0, dx);
	mm_real_free (d0);
	mm_real_free (dx);

	/* [E; dx; dy] */
	d0xy = mm_real_vertcat (d0x, dy);
	mm_real_free (d0x);
	mm_real_free (dy);

	/* [E; dx; dy; dz] */
	d = mm_real_vertcat (d0xy, dz);
	mm_real_free (d0xy);
	mm_real_free (dz);

	return d;
}

static mm_sparse *
mm_real_smooth_sparse_x_1 (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			nh;
	int			nnz;
	mm_sparse	*s;

	double		val = 1.;

	nh = nx * ny;
	m = nh * nz;
	nnz = 2 * m - ny * nz;
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, m, nnz);

#pragma omp parallel for
	for (col = 0; col < s->n; col++) {
		int		i, j, k, p;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = i + j * nx + k * nh;
		p = 2 * index - j - k * ny;
		s->i[p] = index;
		s->data[p++] = val;
		if (i < nx - 1) {
			s->i[p] = index + 1;
			s->data[p++] = - val;
		}
		s->p[col + 1] = p;
	}
	return s;
}

static mm_dense *
mm_real_smooth_dense_x_1 (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			nh;
	mm_dense	*d;

	double		val = 1.;

	nh = nx * ny;
	m = nh * nz;
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, m, m * m);

#pragma omp parallel for
	for (col = 0; col < d->n; col++) {
		int		i, j, k;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = i + j * nx + k * nh;
		d->data[index + col * d->m] = val;
		if (i < nx - 1) d->data[index + 1 + col * d->m] = - val;
	}
	return d;
}

static mm_real *
mm_real_smooth_x_1 (MMRealFormat format, const int nx, const int ny, const int nz)
{
	return (format == MM_REAL_SPARSE) ? mm_real_smooth_sparse_x_1 (nx, ny, nz) : mm_real_smooth_dense_x_1 (nx, ny, nz);
}

static mm_sparse *
mm_real_smooth_sparse_y_1 (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			nh;
	int			nnz;
	mm_sparse	*s;

	double		val = 1.;

	nh = nx * ny;
	m = nh * nz;
	nnz = 2 * m - nx * nz;
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, m, nnz);

#pragma omp parallel for
	for (col = 0; col < s->n; col++) {
		int		i, j, k, p;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = i + j * nx + k * nh;
		p = 2 * index - k * nx;
		if (j == ny - 1) p -= i;
		s->i[p] = index;
		s->data[p++] = val;
		if (j < ny - 1) {
			s->i[p] = index + nx;
			s->data[p++] = - val;
		}
		s->p[col + 1] = p;
	}
	return s;
}

static mm_dense *
mm_real_smooth_dense_y_1 (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			nh;
	mm_dense	*d;

	double		val = 1.;

	nh = nx * ny;
	m = nh * nz;
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, m, m * m);
	mm_real_set_all (d, 0.);

#pragma omp parallel for
	for (col = 0; col < d->n; col++) {
		int		i, j, k;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = i + j * nx + k * nh;
		d->data[index + col * d->m] = val;
		if (j < ny - 1) d->data[index + nx + col * d->m] = - val;
	}
	return d;
}

static mm_real *
mm_real_smooth_y_1 (MMRealFormat format, const int nx, const int ny, const int nz)
{
	return (format == MM_REAL_SPARSE) ? mm_real_smooth_sparse_y_1 (nx, ny, nz) : mm_real_smooth_dense_y_1 (nx, ny, nz);
}

static mm_sparse *
mm_real_smooth_sparse_z_1 (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			nh;
	int			nnz;
	mm_sparse	*s;

	double		val = 1.;

	nh = nx * ny;
	m = nh * nz;
	nnz = 2 * m - nh;
	s = mm_real_new (MM_REAL_SPARSE, MM_REAL_GENERAL, m, m, nnz);

#pragma omp parallel for
	for (col = 0; col < s->n; col++) {
		int		i, j, k, p;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = i + j * nx + k * nh;
		p = 2 * index;
		if (k == nz - 1) p -= i + j * nx;
		s->i[p] = index;
		s->data[p++] = val;
		if (k < nz - 1) {
			s->i[p] = index + nh;
			s->data[p++] = - val;
		}
		s->p[col + 1] = p;
	}
	return s;
}

static mm_dense *
mm_real_smooth_dense_z_1 (const int nx, const int ny, const int nz)
{
	int			col;
	int			m;
	int			nh;
	mm_dense	*d;

	double		val = 1.;

	nh = nx * ny;
	m = nh * nz;
	d = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, m, m * m);
	mm_real_set_all (d, 0.);

#pragma omp parallel for
	for (col = 0; col < d->n; col++) {
		int		i, j, k;
		int		index;
		col2ijk (nx, nh, col, &i, &j, &k);
		index = i + j * nx + k * nh;
		d->data[index + col * d->m] = val;
		if (k < nz - 1) d->data[index + nh + col * d->m] = - val;
	}
	return d;
}

static mm_real *
mm_real_smooth_z_1 (MMRealFormat format, const int nx, const int ny, const int nz)
{
	return (format == MM_REAL_SPARSE) ? mm_real_smooth_sparse_z_1 (nx, ny, nz) : mm_real_smooth_dense_z_1 (nx, ny, nz);
}

/*
	D = [ 
			DX
			DY
			DZ
		]
*/
mm_real *
mm_real_smooth_1 (MMRealFormat format, const int nx, const int ny, const int nz, const double *w)
{
	mm_real	*dx;
	mm_real	*dy;
	mm_real	*dz;
	mm_real	*dxy;
	mm_real	*d;

	/* [dx; dy] */
	dx = mm_real_smooth_x_1 (format, nx, ny, nz);
	dy = mm_real_smooth_y_1 (format, nx, ny, nz);
	dz = mm_real_smooth_z_1 (format, nx, ny, nz);
	if (w) {
		int		j;
		fprintf (stderr, "w[0] = %f, w[1] = %f, w[2] = %f\n", w[0], w[1], w[2]);
		for (j = 0; j < dx->n; j++) mm_real_xj_scale (dx, j, w[0]);
		for (j = 0; j < dy->n; j++) mm_real_xj_scale (dy, j, w[1]);
		for (j = 0; j < dz->n; j++) mm_real_xj_scale (dz, j, w[2]);
	}
	dxy = mm_real_vertcat (dx, dy);
	mm_real_free (dx);
	mm_real_free (dy);

	/* [dx; dy; dz] */
	d = mm_real_vertcat (dxy, dz);
	mm_real_free (dxy);
	mm_real_free (dz);

	return d;
}

