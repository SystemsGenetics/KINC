#ifndef _STATS_
#define _STATS_

#include <stdlib.h>
#include <float.h>
#include <math.h>

/**
 * The code contained in this file and the stats.c file was borrowed from the
 * R package and is reproduced here in accordance with the GNU General Public
 * License.
 */

#define ML_POSINF INFINITY
#define ML_NEGINF -INFINITY
#define ML_NAN NAN

#define FALSE 0
#define TRUE 1

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)    \
  if (log_p) {                                    \
    if(p > 0)                                     \
      return ML_NAN;                              \
    if(p == 0) /* upper bound*/                   \
      return lower_tail ? _RIGHT_ : _LEFT_;       \
    if(p == ML_NEGINF)                            \
      return lower_tail ? _LEFT_ : _RIGHT_;       \
  }                                               \
  else { /* !log_p */                             \
    if(p < 0 || p > 1)                            \
      return ML_NAN;                              \
    if(p == 0)                                    \
      return lower_tail ? _LEFT_ : _RIGHT_;       \
    if(p == 1)                                    \
      return lower_tail ? _RIGHT_ : _LEFT_;       \
  }



#define R_D_Lval(p)     (lower_tail ? (p) : (0.5 - (p) + 0.5))  /*  p  */
#define R_D_Cval(p)     (lower_tail ? (0.5 - (p) + 0.5) : (p))  /*  1 - p */

#define R_DT_qIv(p) (log_p ? (lower_tail ? exp(p) : - expm1(p)) \
                           : R_D_Lval(p))

#define R_DT_CIv(p) (log_p ? (lower_tail ? -expm1(p) : exp(p)) \
                           : R_D_Cval(p))

#define R_D__0  (log_p ? ML_NEGINF : 0.)    /* 0 */
#define R_D__1  (log_p ? 0. : 1.)     /* 1 */
#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)    /* 0 */
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0)    /* 1 */

//#ifndef min
//# define min(a, b) ((a) > (b) ? (b) : (a))
//#endif

#define M_SQRT_32 5.656854249492380195206754896838  /* sqrt(32) */
#define M_1_SQRT_2PI  0.398942280401432677939946059934  /* 1/sqrt(2pi) */

#define R_FINITE(x)    R_finite(x)
#define SIXTEEN  16 /* Cutoff allowing exact "*" and "/" */

int R_finite(double);

double qnorm(double p, double mu, double sigma, int lower_tail, int log_p);
double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
double * ppoints(int n, float a);
double sign(double x);

#endif
