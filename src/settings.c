#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <float.h>
#include <math.h>

#include "consts.h"

static bool
check_num (int n[3])
{
	return (n[0] > 0 && n[1] > 0 && n[2] > 0);
}

static bool
check_range (double r[2])
{
	return (fabs (r[0] - r[1]) > DBL_EPSILON);
}

static char *
skip_headers (char *buf)
{
	char	*p = strrchr (buf, ':');
	if (!p) return NULL;
	return ++p;
}

static void
init (void)
{
	ngrd[0] = 0;
	ngrd[1] = 0;
	ngrd[2] = 0;

	xgrd[0] = 0.;
	xgrd[1] = 0.;

	ygrd[0] = 0.;
	ygrd[1] = 0.;

	zgrd[0] = 0.;
	zgrd[1] = 0.;

	use_dz_array = false;
	dz = NULL;

	return;
}

bool
read_settings (char *fn)
{
	char	buf[BUFSIZ];
	FILE	*fp;

	if (!fn) {
		fprintf (stderr, "ERROR: parameter setting file ./settings is not exists.");
		return false;
	}

	fp = fopen (fn, "r");
	if (!fp) {
		fprintf (stderr, "ERROR: cannot open file %s.\n", fn);
		return false;
	}

	init ();
	while (fgets (buf, BUFSIZ, fp) != NULL) {
		char	*p = skip_headers (buf);
		if (!p) continue;
		switch (buf[0]) {
			case '#':
			case ' ':
			case '\t':
				continue;
			case '1':
				sscanf (p, "%lf,%lf", &exf_inc, &exf_dec);
				break;
			case '2':
				sscanf (p, "%d,%d,%d", &ngrd[0], &ngrd[1], &ngrd[2]);
				break;
			case '3':
				sscanf (p, "%lf,%lf", &xgrd[0], &xgrd[1]);
				break;
			case '4':
				sscanf (p, "%lf,%lf", &ygrd[0], &ygrd[1]);
				break;
			case '5':
				sscanf (p, "%lf,%lf", &zgrd[0], &zgrd[1]);
				break;

			default:
				break;
		}
	}
	fclose (fp);

	if (!check_num (ngrd)) {
		fprintf (stderr, "ERROR: specified ngrd is invalid.\n");
		return false;
	}
	if (!check_range (xgrd)) {
		fprintf (stderr, "ERROR: specified xgrd is invalid.\n");
		return false;
	}
	if (!check_range (ygrd)) {
		fprintf (stderr, "ERROR: specified ygrd is invalid.\n");
		return false;
	}
	if (!check_range (zgrd)) {
		fprintf (stderr, "ERROR: specified zgrd is invalid.\n");
		return false;
	}

	return true;
}

void
fprintf_settings (FILE *stream)
{
	fprintf (stream, "1. inclination, declination: %f, %f\n", exf_inc, exf_dec);
	fprintf (stream, "2. num of grids: %d, %d, %d\n", ngrd[0], ngrd[1], ngrd[2]);
	fprintf (stream, "3. x range: %f, %f\n", xgrd[0], xgrd[1]);
	fprintf (stream, "4. y range: %f, %f\n", ygrd[0], ygrd[1]);
	fprintf (stream, "5. z range: %f, %f\n", zgrd[0], zgrd[1]);

	return;
}

