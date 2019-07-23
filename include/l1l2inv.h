#ifndef _L1L2INV_H_
#define _L1L2INV_H_

int		num_separator (char *str, const char c);
bool	l1l2inv_constraint_func (cdescent *cd, const int j, const double etaj, double *val);
bool	l1l2inv (simeq *eq, char *path_fn, char *info_fn);

#endif // _L1L2INV_H_
