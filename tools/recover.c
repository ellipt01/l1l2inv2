#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

#include "cdescent.h"
#include "mgcal.h"

#include "simeq.h"
#include "utils.h"
#include "settings.h"
#include "defaults.h"

static int	iter = -1;
static bool	beta_file_specified = false;
static char	fn[80];
static bool	input_file_specified = false;

void
fprintf_estimated (FILE *stream, const char *fn, const double *val, const char *format)
{
	data_array	*data;
	FILE		*fp;

	if ((fp = fopen (fn, "r")) == NULL) return;
	data = fread_data_array (fp);
	fclose (fp);

	fwrite_data_array_with_data (stream, data, val, format);
	data_array_free (data);

	return;
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
	fprintf (stderr, "This program estimates recoverd anomaly,\n");
	fprintf (stderr, "by reading solution path file (e.g. beta_path.data).\n");
	fprintf (stderr, "Users have to specify number of model (i.e. iteration number).\n\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -f <solution path file name (e.g. beta_path.data)>\n");
	fprintf (stderr, "       -i <iteration number>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -d <input data filename, default is input.data>\n");
	fprintf (stderr, "       -s <parameter setting file: default=./settings>\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	return;
}

static bool
read_input_params (int argc, char **argv)
{
	char	c;

	while ((c = getopt (argc, argv, "f:i:d:s:h")) != EOF) {
		switch (c) {
			case 'f':
				strcpy (fn, optarg);
				beta_file_specified = true;
				break;
			case 'i':
				iter = atoi (optarg);
				break;
			case 'd':
				strcpy (ifn, optarg);
				input_file_specified = true;
				break;
			case 's':
				strcpy (sfn, optarg);
				break;
			case 'h':
			case ':':
			case '?':
				return false;
		}	
	}
	if (!beta_file_specified || iter < 0) return false;
	return true;
}

int
main (int argc, char **argv)
{
	simeq		*eq;
	mm_dense	*beta;
	mm_dense	*f;

	if (!read_input_params (argc, argv)) {
		usage (argv[0]);
		return EXIT_FAILURE;
	}
	if (!read_settings (sfn)) return EXIT_FAILURE;
	mgcal_set_scale_factor (magscale);

	// extract beta
	beta = extract_beta (fn, iter);
	if (!beta) {
		fprintf (stderr, "ERROR: failed to read %d-th beta.\n", iter);
		exit (EXIT_FAILURE);
	}

	// create simeq
	if ((eq = read_input (-1, ifn, tfn, NULL)) == NULL) return EXIT_FAILURE;

	// calc and output z = X * beta
	f = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, eq->x->m, 1, eq->x->m);
	mm_real_set_all (f, 0.);
	mm_real_x_dot_yk (false, 1., eq->x, beta, 0, 0., f);
	mm_real_free (beta);
	simeq_free (eq);

	if (!input_file_specified) strcpy (ifn, "input.data");
	fprintf_estimated (stdout, ifn, f->data, NULL);
	mm_real_free (f);

	return EXIT_SUCCESS;
}

