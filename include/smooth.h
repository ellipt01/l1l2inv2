#ifndef _SMOOTH_H_
#define _SMOOTH_H_

mm_real	*mm_real_smooth (MMRealFormat format, const int nx, const int ny, const int nz);
mm_real	*mm_real_smooth_l01 (MMRealFormat format, const int nx, const int ny, const int nz);
mm_real	*mm_real_smooth_1 (MMRealFormat format, const int nx, const int ny, const int nz);

#endif // _SMOOTH_H_
