/*
 * linregmodel.h
 *
 *  linear regression model
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef LINREGMODEL_H_
#define LINREGMODEL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <objects.h>

/* linregmodel.c */
linregmodel	*linregmodel_new (const int m, const int n, int nxfiles, mm_dense *y, const mm_real *d);
void		linregmodel_free (linregmodel *l);

#ifdef __cplusplus
}
#endif

#endif /* LINREGMODEL_H_ */
