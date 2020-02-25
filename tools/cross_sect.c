#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#include "cdescent.h"
#include "mgcal.h"

#include "consts.h"
#include "simeq.h"
#include "utils.h"
#include "settings.h"
#include "defaults.h"

char		fn[80];
bool		terrain_file_specified = false;
char		fn_ter[80];
int			iter = -1;
int			direction = -1;
double		coordinate;
bool		zfile = true;

static double *
read_z (FILE *fp, const int n)
{
	int		i;
	double	xval, yval, zval;
	double	*z = (double *) malloc (n * sizeof (double));

	i = 0;
	while (fscanf (fp, "%lf\t%lf\t%lf", &xval, &yval, &zval) != EOF) {
		z[i] = zval;
		if (++i >= n) break;
	}
	return z;
}

static void
fprintf_cross_section (FILE *stream, bool zfile, const int d, const double c,
	const double *data, const double *z1)
{
	int			i, j, k, l, n;
	grid		*g;
	double		*v;

	bool		terrain_specified = false;
	double		*z = NULL;

	int			start = -1;

	g = grid_new (ngrd[0], ngrd[1], ngrd[2], xgrd, ygrd, zgrd);

	if (z1) terrain_specified = true;

	switch (d) {
		case CROSS_SECTION_X:
			n = g->nx * g->nz;
			if (terrain_specified) z = (double *) malloc (g->nx * sizeof (double));
			break;
		case CROSS_SECTION_Y:
			n = g->ny * g->nz;
			if (terrain_specified) z = (double *) malloc (g->ny * sizeof (double));
			break;
	}

	v = (double *) malloc (n * sizeof (double));

	switch (d) {
		case CROSS_SECTION_X:
			for (k = 0; k < g->ny; k++) {
				double	yk = g->y[k];
				if (c <= yk) {
					start = k;
					break;
				}
			}
			if (start < 0) {
				fprintf (stderr, "ERROR: no valid grid found.\n");
				exit (1);
			}
			l = 0;
			for (j = 0; j < g->nz; j++) {
				for (i = 0; i < g->nx; i++) {
					k = i + start * g->nx + j * g->nh;
					v[l++] = data[k];
				}
			}
			break;

		case CROSS_SECTION_Y:
			for (k = 0; k < g->nx; k++) {
				double	xk = g->x[k];
				if (c <= xk) {
					start = k;
					break;
				}
			}
			if (start < 0) {
				fprintf (stderr, "ERROR: no valid grid found.\n");
				exit (1);
			}
			l = 0;
			for (j = 0; j < g->nz; j++) {
				for (i = 0; i < g->ny; i++) {
					k = i * g->nx + start + j * g->nh;
					v[l++] = data[k];
				}
			}
			break;
	}
	if (l == 0) {
		fprintf (stderr, "ERROR: no valid grid found.\n");
		exit (1);
	}

	if (terrain_specified) {
		switch (d) {
			case CROSS_SECTION_X:
				k = start * g->nx;
				for (i = 0; i < g->nx; i++) {
					z[i] = z1[k];
					k++;
					if (k >= g->n) break;
				}
				break;
			case CROSS_SECTION_Y:
				k = start;
				for (i = 0; i < g->ny; i++) {
					z[i] = z1[k];
					k += g->nx;
					if (k >= g->n) break;
				}
				break;
		}
	}

	switch (d) {
		case CROSS_SECTION_X:
			fprintf (stderr, "OUTPUT: cross section through y = %.4f\n", g->y[start]);
			break;

		case CROSS_SECTION_Y:
			fprintf (stderr, "OUTPUT: cross section through x = %.4f\n", g->x[start]);
			break;
	}
	grid_free (g);

	switch (d) {
		case CROSS_SECTION_X:
			g = grid_new (ngrd[0], 1, ngrd[2], xgrd, &c, zgrd);
			break;
		case CROSS_SECTION_Y:
			g = grid_new (1, ngrd[1], ngrd[2], &c, ygrd, zgrd);
			break;
	}

	if (terrain_specified) grid_set_surface (g, z);
	if (zfile) fprintf_zfile_2d (stream, g, d, v);
	else fwrite_grid_with_data (stream, g, v, NULL);

	grid_free (g);
	free (v);
	return;
}

static bool
check_range (const int d, const double c)
{
	if (d == CROSS_SECTION_X) {
		return (ygrd[0] <= c && c <= ygrd[1]);
	} else if (d == CROSS_SECTION_Y) {
		return (xgrd[0] <= c && c <= xgrd[1]);
	}
	return false;	
}

static void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;

	fprintf (stderr, "\n");
	version_info (p);
	fprintf (stderr, "\n");

	fprintf (stderr, "DESCRIPTION:\n");
	fprintf (stderr, "This program derives the cross-section of a model\n");
	fprintf (stderr, "obtained by \"l1l2inv\" program,\n\n");
	fprintf (stderr, "by reading solution path file (e.g. beta_path.data).\n");
	fprintf (stderr, "Users have to specify\n");
	fprintf (stderr, "   number of model (i.e. iteration number)\n");
	fprintf (stderr, "   direction of the profile\n");
	fprintf (stderr, "   coordinate of the profile\n\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -f <solution path file name (e.g. beta_path.data)>\n");
	fprintf (stderr, "       -i <iteration number>\n");
	fprintf (stderr, "       -d <direction: x = 0, y = 1>\n");
	fprintf (stderr, "       -c <x or y val of perpendicular coordinate>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -o (output to grid (x,y,z, and magnetization) format,\n");
	fprintf (stderr, "           defaukt is GMT closed polygon format)\n");
	fprintf (stderr, "       -t <terrain file>\n");
	fprintf (stderr, "          if this option is not specified,\n");
	fprintf (stderr, "          terrain is assumed to be a flat plane\n");
	fprintf (stderr, "       -s <parameter setting file: default=./settings>\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	exit (1);
}

static bool
read_input_params (int argc, char **argv)
{
	char	c;
	bool	file_specified = false;
	bool	iter_specified = false;
	bool	direction_specified = false;
	bool	coordinate_specified = false;

	while ((c = getopt (argc, argv, ":f:i:d:c:z:t:s:oh")) != -1) {
		switch (c) {
			case 'f':
				strcpy (fn, optarg);
				file_specified = true;
				break;
			case 'i':
				iter = atoi (optarg);
				iter_specified = true;
				break;
			case 'd':
				direction = atoi (optarg);
				direction_specified = true;
				break;
			case 'c':
				coordinate = (double) atof (optarg);
				coordinate_specified = true;
				break;
			case 'o':
				zfile = false;
				break;
			case 't':
				strcpy (fn_ter, optarg);
				terrain_file_specified = true;
				break;
			case 's':
				strcpy (sfn, optarg);
				break;
			case 'h':
				usage (argv[0]);
				break;
			case ':':
				fprintf (stdout, "%c needs value.\n", c);
     			break;
			case '?':
				fprintf (stderr, "unknown option.\n");
				break;
			default:
				break;
		}
	}

	if (!file_specified) {
		fprintf (stderr, "ERROR: input file name is not specified (-f).\n");
		return false;
	}
	if (!iter_specified) {
		fprintf (stderr, "ERROR: iteration number is not specified (-i).\n");
		return false;
	}
	if (!direction_specified) {
		fprintf (stderr, "ERROR: direction is not specified (-d)\n");
		return false;
	}
	if (!coordinate_specified) {
		fprintf (stderr, "ERROR: coordinate is not specified (-c)\n");
		return false;
	}
	if (strlen (fn) <= 0) {
		fprintf (stderr, "ERROR: specified input file name invalid.\n");
		return false;
	}
	if (iter < 0) {
		fprintf (stderr, "ERROR: specified iteration number invalid.\n");
		return false;
	}
	if (direction != 0 && direction != 1) {
		fprintf (stderr, "ERROR: specified direction invalid.\n");
		return false;
	}
	return true;
}

int
main (int argc, char **argv)
{
	double		*z1 = NULL;
	mm_dense	*beta;

	if (!read_input_params (argc, argv)) usage (argv[0]);
	if (!read_settings (sfn)) return EXIT_FAILURE;

	if (!check_range (direction, coordinate)) {
		fprintf (stderr, "ERROR: specified coordinate out of range.\n");
		exit (EXIT_FAILURE);
	}

	beta = extract_beta (fn, iter);
	if (!beta) {
		fprintf (stderr, "ERROR: failed to extracu %dth beta.\n", iter);
		exit (EXIT_FAILURE);
	}

	if (terrain_file_specified) {
		FILE	*fp;
		if ((fp = fopen (fn_ter, "r")) == NULL) {
			fprintf (stderr, "ERROR: cannot open file %s.\n", fn_ter);
			exit (1);	
		}
		z1 = read_z (fp, beta->nnz);
	}

	fprintf_cross_section (stdout, zfile, direction, coordinate, beta->data, z1);
	mm_real_free (beta);

	return EXIT_SUCCESS;
}

