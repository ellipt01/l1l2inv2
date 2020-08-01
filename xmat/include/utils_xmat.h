#ifndef _IO_H_
#define _IO_H_

#ifndef CROSS_SECTION_X
#define CROSS_SECTION_X		0
#endif

#ifndef CROSS_SECTION_Y
#define CROSS_SECTION_Y		1
#endif

double		*fread_z (FILE *fp, const int n);
double		*get_nth_data (FILE *fp, const int n, int *len);
mm_dense	*extract_beta (const char *fn, const int j);

// read
simeq		*read_input (const int type, const char *ifn, const char *tfn, const double *w, bool create_xmat);
source		*read_model_par (FILE *fp, const double exf_inc, const double exf_dec);

// output
void		printf_beta (const mm_dense *beta);
void		printf_mm_real (const char *fn, const mm_real *x, const char *format);
void		printf_array (const char *fn, const int m, const int n, const double *a, const char *format);
void		printf_estimated (const cdescent *cd, const mm_dense *beta);
int			fread_ndata_of_one_line (FILE *fp);
double		*fread_one_line (FILE *fp, int len);
void		fprintf_zfile_2d (FILE *stream, grid *g, const int dir, const double *v,
							  const double *t, const double *offset);
void		fprintf_zfile_3d (FILE *stream, grid *g, const double *v,
							 const double *t, const double *offset);
void		fprintf_regularization_type (FILE *stream, const int type);
void		fprintf_sources (FILE *stream, const source *s);
void		version_info (const char *toolname);

#endif // _IO_H_
