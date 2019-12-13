#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cdescent.h>
#include <mgcal.h>

#include "consts.h"
#include "simeq.h"
#include "utils.h"
#include "settings.h"
#include "defaults.h"

#define	TYPE_GRID	0
#define	TYPE_ARRAY	1

int		obs_points_type;

static void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;

	fprintf (stderr, "\n");
	version_info (p);
	fprintf (stderr, "\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -f <model parameter file>\n");
	fprintf (stderr, "       -z <altitude of observation points>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -d <std of noise: default=0>\n");
	fprintf (stderr, "       -r <seed of rand-generator>\n");
	fprintf (stderr, "       -s <parameter setting file: default=./settings>\n");
	fprintf (stderr, "       -o <output file name>\n");
	fprintf (stderr, "       -h (show this message)\n");
	exit (1);
}

char		ifn[80];
bool		pfn_specified = false;
char		pfn[80];
char		ofn[80];
double		std = 0.;
int			seed = 100;
double		*zobs = NULL;
double		zobs_val;

static bool
read_input_params (int argc, char **argv)
{
	char	c;

	obs_points_type = TYPE_GRID;
	strcpy (ofn, "input.data");
	while ((c = getopt (argc, argv, "f:z:i:d:r:s:o:h")) != EOF) {
		switch (c) {
			case 'f':
				strcpy (pfn, optarg);
				pfn_specified = true;
				break;
			case 'i':
				obs_points_type = TYPE_ARRAY;
				strcpy (ifn, optarg);
				break;
			case 'z':
				zobs_val = (double) atof (optarg);
				zobs = &zobs_val;
				break;
			case 'd':
				std = (double) atof (optarg);
				break;
			case 'r':
				seed = atoi (optarg);
				break;
			case 's':
				strcpy (sfn, optarg);
				break;
			case 'o':
				strcpy (ofn, optarg);
				break;
			case 'h':
			case ':':
			case '?':
				return false;
			default:
				break;		
		}
	}
	if (!pfn_specified) {
		fprintf (stderr, "ERROR: please specify input parameter file.\n");
		return false;
	}
	if (obs_points_type == TYPE_GRID && !zobs) {
		fprintf (stderr, "ERROR: please specify altitude of observation.\n");
		return false;
	}
	return true;
}

static void
create_input_data_by_grid (FILE *stream, const grid *g, const source *s, mgcal_theoretical func)
{
	int			n;
	gsl_rng		*rng = gsl_rng_alloc (gsl_rng_default);
	double		*a = (double *) malloc (g->n * sizeof (double));
	vector3d	obs;

	FILE		*fp = fopen ("noise.data", "w");

	gsl_rng_set (rng, seed);

	for (n = 0; n < g->n; n++) {
		// variance sigma^2 = (mgcal_get_scale_factor () / 5.)^2
		double	r = (std > 0.) ? gsl_ran_gaussian (rng, std) : 0.;
		grid_get_nth (g, n, &obs, NULL);
		a[n] = func (&obs, s, NULL) + r;
		if (fp) fprintf (fp, "%.4e\n", r);
	}
	if (fp) fclose (fp);
	fprintf (stderr, "standard deviation of noise: sigma = %f\n", std);
	fwrite_grid_with_data (stream, g, a, NULL);
	gsl_rng_free (rng);
	free (a);

	return;
}

static void
create_input_data_by_array (FILE *stream, const data_array *array, const source *s, mgcal_theoretical func)
{
	int			n;
	gsl_rng		*rng = gsl_rng_alloc (gsl_rng_default);
	double		*a = (double *) malloc (array->n * sizeof (double));
	vector3d	obs;

	FILE		*fp = fopen ("noise.data", "w");

	gsl_rng_set (rng, seed);

	for (n = 0; n < array->n; n++) {
		// variance sigma^2 = (mgcal_get_scale_factor () / 5.)^2
		double	r = (std > 0.) ? gsl_ran_gaussian (rng, std) : 0.;
		obs.x = array->x[n];
		obs.y = array->y[n];
		obs.z = array->z[n];
		a[n] = func (&obs, s, NULL) + r;
		if (fp) fprintf (fp, "%.4e\n", r);
	}
	if (fp) fclose (fp);
	fprintf (stderr, "standard deviation of noise: sigma = %f\n", std);
	fwrite_data_array_with_data (stream, array, a, "%.4f\t%.4f\t%.4f\t%.4f");
	gsl_rng_free (rng);
	free (a);

	return;
}

int
main (int argc, char **argv)
{
	int		n = 1;
	source	*s;
	FILE	*fpi;
	FILE	*fpo;

	mgcal_theoretical	f;

	if (!read_input_params (argc, argv)) {
		usage (argv[0]);
		return EXIT_FAILURE;
	}
	if (!read_settings (sfn)) return EXIT_FAILURE;

	if ((fpi = fopen (pfn, "r")) == NULL) {
		fprintf (stderr, "ERROR: cannot open file %s.\nAbort.\n", pfn);
		return EXIT_FAILURE;
	}
	if ((fpo = fopen (ofn, "w")) == NULL) {
		fprintf (stderr, "ERROR: cannot open file %s.\nAbort.\n", ofn);
		return EXIT_FAILURE;
	}

	mgcal_set_scale_factor (magscale);
	s = read_model_par (fpi, exf_inc, exf_dec);
	fclose (fpi);
	fprintf_sources (stderr, s);

	f = total_force_prism;

	if (obs_points_type == TYPE_GRID) {
		grid	*g = grid_new (ngrd[0], ngrd[1], n, xgrd, ygrd, zobs);
		create_input_data_by_grid (fpo, g, s, f);
		grid_free (g);
	} else if (obs_points_type == TYPE_ARRAY) {
		FILE		*fp;
		data_array *array;

		fp = fopen (ifn, "r");
		if (!fp) {
			fprintf (stderr, "ERROR: cannot open file %s.\nAbort.\n", ifn);
			return EXIT_FAILURE;
		}
		array = fread_data_array (fp);
		fclose (fp);
		create_input_data_by_array (fpo, array, s, f);
	} else {
		fprintf (stderr, "ERROR: observation points type is invalid.\nAbort.\n");
		return EXIT_FAILURE;
	}
	fclose (fpo);

	source_free (s);

	return EXIT_SUCCESS;
}

