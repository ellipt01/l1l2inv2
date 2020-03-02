#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

#include "cdescent.h"
#include "mgcal.h"

#include "consts.h"
#include "simeq_xmat.h"
#include "utils_xmat.h"
#include "settings.h"
#include "defaults.h"

static int	iter = -1;
static bool	beta_file_specified = false;
static char	fn[80];
static bool	input_file_specified = false;

bool		create_xmat = true;
int			num_xfiles;
int			xfile_len;


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
	fprintf (stderr, "Users have to specify number of model (i.e. iteration number).\n");
	fprintf (stderr, "** This is a version that store matrix X in files **\n\n");

	fprintf (stderr, "USAGE: %s\n", p);
	fprintf (stderr, "       -f <solution path file name (e.g. beta_path.data)>\n");
	fprintf (stderr, "       -i <iteration number>\n");
	fprintf (stderr, "[optional]\n");
	fprintf (stderr, "       -d <input data filename, default is input.data>\n");
	fprintf (stderr, "       -s <parameter setting file: default=./settings>\n");
	fprintf (stderr, "       -x (xmat files are already exist, and use them, default=false)\n");
	fprintf (stderr, "       -h (show this message)\n\n");
	return;
}

static bool
read_input_params (int argc, char **argv)
{
	char	c;

	while ((c = getopt (argc, argv, "f:i:d:s:xh")) != EOF) {
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
			case 'x':
				create_xmat = false;
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
	int			i;
	int			m, n, ret;

	double		*xtx;

	simeq		*eq;
	mm_dense	*beta;
	mm_dense	*f;

	FILE		**fp_xmat;
	FILE		*fp_xtx;

	if (!read_input_params (argc, argv)) {
		usage (argv[0]);
		return EXIT_FAILURE;
	}
	if (!read_settings (sfn)) return EXIT_FAILURE;
	mgcal_set_scale_factor (magscale);

	// create simeq
	if ((eq = read_input (-1, ifn, tfn, NULL, create_xmat)) == NULL) return EXIT_FAILURE;

	m = eq->y->m;
	n = ngrd[0] * ngrd[1] * ngrd[2];
	num_xfiles = ngrd[2];
	xfile_len = ngrd[0] * ngrd[1];

	fp_xmat = (FILE **) malloc (num_xfiles * sizeof (FILE *));
	for (i = 0; i < num_xfiles; i++) {
		char	fn[80];
		sprintf (fn, "x%03d.mat", i);
		fp_xmat[i] = fopen (fn, "rb");
		if (!fp_xmat[i]) {
			fprintf (stderr, "ERROR: cannot open file %s.program abort!\n", fn);
			exit (1);
		}
	}

	fp_xtx = fopen ("xtx.vec", "rb");
	if (!fp_xtx) {
		fprintf (stderr, "ERROR: cannot open file xtx.vec.\nprogram abort!\n");
		exit (1);
	}
	xtx = (double *) malloc (n * sizeof (double));
	ret = fread (xtx, sizeof (double), n, fp_xtx);
	fclose (fp_xtx);

	// extract beta
	beta = extract_beta (fn, iter);
	if (!beta) {
		fprintf (stderr, "ERROR: failed to read %d-th beta.\n", iter);
		exit (EXIT_FAILURE);
	}
	for (i = 0; i < n; i++) beta->data[i] *= sqrt (xtx[i]);

	// calc and output z = X * beta
	f = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, m, 1, m);
	mm_real_set_all (f, 0.);
	mm_real_x_dot_yk_xmatfile (fp_xmat, m, n, false, 1., beta, 0, 0., f);
	mm_real_free (beta);
	simeq_free (eq);
	for (i = 0; i < num_xfiles; i++) fclose (fp_xmat[i]);
	free (fp_xmat);

	if (!input_file_specified) strcpy (ifn, "input.data");
	fprintf_estimated (stdout, ifn, f->data, NULL);
	mm_real_free (f);

	return EXIT_SUCCESS;
}

