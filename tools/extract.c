#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "cdescent.h"
#include "mgcal.h"

#include "consts.h"
#include "simeq.h"
#include "utils.h"
#include "settings.h"
#include "defaults.h"

extern char	*optarg;

int			idx = -1;
char		*fn = NULL;
char		*ter_fn = NULL;

int			len_fn = 256;

typedef enum {
	OUTPUT_TO_GRID,
	OUTPUT_TO_MMREAL
} OutputFormat;

OutputFormat	oformat = OUTPUT_TO_GRID;

static void
fprintf_grid (FILE *stream, const double *data, char *fn_ter)
{
	grid	*g = grid_new (ngrd[0], ngrd[1], ngrd[2], xgrd, ygrd, zgrd);
	if (fn_ter) {
		// if terrain file name "fn_ter" is specified, try to read it 
		FILE	*fp = fopen (fn_ter, "r");
		if (fp) {
			double	*z = fread_z (fp, g->nh);
			grid_set_surface (g, z);
			free (z);
			fclose (fp);
		} else {
			fprintf (stderr, "WARNING: terrain file %s not exists.\n", fn_ter);
		}
	} else {
		// if fn_ter is not specified,
		// try to read default terrain file tfn = "./terrain.data" 
		FILE	*fp = fopen (tfn, "r");
		if (fp) {
			double	*z = fread_z (fp, g->nh);
			grid_set_surface (g, z);
			free (z);
			fclose (fp);
		}
		// else, terrain is ignored and assume to be a flat plane
	}
	fwrite_grid_with_data (stdout, g, data, NULL);
	grid_free (g);
	return;
}

static bool
read_input_params (int argc, char **argv)
{
	char	c;

	while ((c = getopt (argc, argv, "f:i:o:t:s:h")) != EOF) {
		switch (c) {
			case 'f':
				fn = (char *) malloc (len_fn * sizeof (char));
				strcpy (fn, optarg);
				break;
			case 'i':
				idx = atoi (optarg);
				break;
			case 'o':
				oformat = (atoi (optarg) == 0) ? OUTPUT_TO_GRID : OUTPUT_TO_MMREAL;
				break;
			case 't':
				ter_fn = (char *) malloc (len_fn * sizeof (char));
				strcpy (ter_fn, optarg);
				break;
			case 's':
				strcpy (sfn, optarg);
				break;
			case 'h':
			case ':':
			case '?':
				return false;
			default:
				break;
		}
	}
	if (!fn || idx < 0) return false;

	return true;
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
	fprintf (stderr, "This program extract specified number of model\n");
	fprintf (stderr, "derived by \"l1l2inv\",\n");
	fprintf (stderr, "by reading solution path file (e.g. beta_path.data).\n");
	fprintf (stderr, "Users have to specify number of model (i.e. iteration number).\n");
	fprintf (stderr, "Users can choose output file format either\n");
	fprintf (stderr, "grid (x, y, z, and magnetization) or matrix-market format.\n\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -f <solution path file name (e.g. beta_path.data)>\n");
	fprintf (stderr, "       -i <iteration number>\n");

	fprintf (stderr, "[optional]\n");

	fprintf (stderr, "       -o <output format: grid=0,mm_real=1, default is 0>\n");
	fprintf (stderr, "       -t <terrain file>\n");
	fprintf (stderr, "          if this option is not specified,\n");
	fprintf (stderr, "          terrain is assumed to be a flat plane.\n");
	fprintf (stderr, "          When \"-o 1\" is specified, this option is ignored.\n");
	fprintf (stderr, "       -s <parameter setting file: default=./settings>\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	exit (1);
}

int
main (int argc, char **argv)
{
	mm_dense	*beta;
	FILE		*fp;

	if (!read_input_params (argc, argv)) usage (argv[0]);
	if (!read_settings (sfn)) return EXIT_FAILURE;

	if ((beta = extract_beta (fn, idx)) == NULL) {
		fprintf (stderr, "ERROR: failed to read %d'th beta.\n", idx);
		return EXIT_FAILURE;
	}
	if (beta->m != ngrd[0] * ngrd[1] * ngrd[2]) {
		fprintf (stderr, "ERROR: data length mismatch %d / %d.\n", beta->n, ngrd[0] * ngrd[1] * ngrd[2]);
		return EXIT_FAILURE;
	}

	if (oformat == OUTPUT_TO_GRID) fprintf_grid (stdout, beta->data, ter_fn);
	else mm_real_fwrite (stdout, beta, "%.8e");

	//fprintf (stderr, "|beta[%d]| = %.4e\n", idx, mm_real_xj_nrm2 (beta, 0));
	mm_real_free (beta);

	return EXIT_SUCCESS;
}

