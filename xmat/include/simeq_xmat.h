#ifndef _SIMEQ_H_
#define _SIMEQ_H_

typedef struct {
	mm_dense	*y;
	mm_sparse	*d;
} simeq;

enum {
	TYPE_NONE = -1,
	TYPE_L1 = 0,
	TYPE_L1L2 = 1,
	TYPE_L1D1 = 2,
	TYPE_L1D01 = 3
};

simeq	*simeq_new (void);
void	simeq_free (simeq *eq);
void	simeq_centering_y (simeq *eq);
void	simeq_centering_x (simeq *eq);
void	simeq_normalizing_x (simeq *eq);
void	simeq_standardizing_x (simeq *eq);
simeq	*create_simeq (
			const int type, const double inc, const double dec,
			const data_array *array, const grid *gsrc, const mgcal_func *func,
			bool create_xmat
		);

#endif // _SIMEQ_H_
