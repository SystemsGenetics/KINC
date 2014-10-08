#include "stats.h"


/**
 *
 */
int R_finite(double x) {
#ifdef HAVE_WORKING_ISFINITE
    return isfinite(x);
#else
    return (!isnan(x) & (x != ML_POSINF) & (x != ML_NEGINF));
#endif
}

/**
 * This function computes the  'signum(.)' function:
 */
double sign(double x) {
  if (isnan(x)) {
    return x;
  }
  return ((x > 0) ? 1 : ((x == 0)? 0 : -1));
}

/**
 * Compute the quantile function for the normal distribution.
 */
double qnorm(double p, double mu, double sigma, int lower_tail, int log_p) {
  double p_, q, r, val;

#ifdef IEEE_754
    if (isnan(p) || isnan(mu) || isnan(sigma))
  return p + mu + sigma;
#endif

  R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

  p_ = R_DT_qIv(p); // real lower_tail prob. p
  q = p_ - 0.5;

  /*-- use AS 241 --- *
   * double ppnd16_(double *p, long *ifault)*
   * ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
   *
   * Produces the normal deviate Z corresponding to a given lower
   * tail area of P; Z is accurate to about 1 part in 10**16.
   *
   * (original fortran code used PARAMETER(..) for the coefficients
   * and provided hash codes for checking them...)
   */
  /* 0.075 <= p <= 0.925 */
  if (fabs(q) <= .425) {

    r = .180625 - q * q;
    val = q * (((((((r * 2509.0809287301226727 +
                     33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                     1971.5909503065514427) * r + 133.14166789178437745) * r +
                     3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                     21213.794301586595867) * r + 5394.1960214247511077) * r +
                     687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
  }
  /* closer than 0.075 from {0,1} boundary */
  else {
    /* r = min(p, 1-p) < 0.075 */
    if (q > 0) {
      r = R_DT_CIv(p); /* 1-p */
    }
    else {
      r = p_; /* = R_DT_Iv(p) ^=  p */
    }
    r = sqrt(- ((log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ? p : /* else */ log(r)));

    if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
      r += -1.6;
      val = (((((((r * 7.7454501427834140764e-4 +
                   .0227238449892691845833) * r + .24178072517745061177) *
                    r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                    1.42343711074968357734)
              / (((((((r *
                       1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                       r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                       r + 1.6763848301838038494) * r +
                       2.05319162663775882187) * r + 1.);
    }
    else { /* very close to  0 or 1 */
      r += -5.;
      val = (((((((r * 2.01033439929228813265e-7 +
                   2.71155556874348757815e-5) * r +
                   .0012426609473880784386) * r + .026532189526576123093) *
                   r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772)
             / (((((((r *
                   2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                   r + 1.8463183175100546818e-5) * r +
                   7.868691311456132591e-4) * r + .0148753612908506148525)
                   * r + .13692988092273580531) * r +
                  .59983220655588793769) * r + 1.);
    }

    if (q < 0.0) {
      val = -val;
    }
    /* return (q >= 0.)? r : -r ;*/
  }
  return mu + sigma * val;
}

/**
 *
 */
double pnorm(double x, double mu, double sigma, int lower_tail, int log_p) {
    double p, cp;

  /* Note: The structure of these checks has been carefully thought through.
   * For example, if x == mu and sigma == 0, we get the correct answer 1.
   */
#ifdef IEEE_754
  if(isnan(x) || isnan(mu) || isnan(sigma)) {
    return x + mu + sigma;
  }
#endif
  if (!R_FINITE(x) && mu == x) {
    return ML_NAN; /* x-mu is NaN */
  }
  if (sigma <= 0) {
    if(sigma < 0) {
      return ML_NAN;
    }
    /* sigma = 0 : */
    return (x < mu) ? R_DT_0 : R_DT_1;
  }
  p = (x - mu) / sigma;
  if (!R_FINITE(p)) {
    return (x < mu) ? R_DT_0 : R_DT_1;
  }
  x = p;

  pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

  return (lower_tail ? p : cp);
}

/**
 *
 */
void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p) {
  /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
     if (lower) return  *cum := P[X <= x]
     if (upper) return *ccum := P[X >  x] = 1 - P[X <= x]
  */
  const static double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const static double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };
  const static double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const static double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };
  const static double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const static double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  double xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
  double min = DBL_MIN;
#endif
  int i, lower, upper;

#ifdef IEEE_754
  if (isnan(x)) {
    *cum = *ccum = x;
    return;
  }
#endif

  /* Consider changing these : */
  eps = DBL_EPSILON * 0.5;

  /* i_tail in {0,1,2} =^= {lower, upper, both} */
  lower = i_tail != 1;
  upper = i_tail != 0;

  y = fabs(x);
  if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
    if (y > eps) {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (i = 0; i < 3; ++i) {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
      }
    }
    else {
      xnum = xden = 0.0;
    }

    temp = x * (xnum + a[3]) / (xden + b[3]);
    if (lower) {
      *cum = 0.5 + temp;
    }
    if (upper) {
      *ccum = 0.5 - temp;
    }
    if (log_p) {
      if (lower) {
        *cum = log(*cum);
      }
      if (upper) {
        *ccum = log(*ccum);
      }
    }
  }
  else if (y <= M_SQRT_32) {
    /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)             \
    xsq = trunc(X * SIXTEEN) / SIXTEEN;       \
    del = (X - xsq) * (X + xsq);          \
    if(log_p) {             \
      *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp); \
      if((lower && x > 0.) || (upper && x <= 0.))     \
      *ccum = log1p(-exp(-xsq * xsq * 0.5) *    \
        exp(-del * 0.5) * temp);    \
    }               \
    else {                \
      *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;  \
      *ccum = 1.0 - *cum;           \
    }

#define swap_tail           \
    if (x > 0.) {/* swap  ccum <--> cum */      \
        temp = *cum; if(lower) *cum = *ccum; *ccum = temp;  \
    }

    do_del(y);
    swap_tail;
  }

  /* else   |x| > sqrt(32) = 5.657 :
   * the next two case differentiations were really for lower=T, log=F
   * Particularly  *not*  for  log_p !

   * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
   *
   * Note that we do want symmetry(0), lower/upper -> hence use y
   */
  else if((log_p && y < 1e170) /* avoid underflow below */
    /*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
     * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

     xsq = x*x;

     if(xsq * DBL_EPSILON < 1.)
        del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
     else
        del = 0.;
     *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
     *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

     swap_tail;

     [Yes, but xsq might be infinite.]

    */
      || (lower && -37.5193 < x  &&  x < 8.2924)
      || (upper && -8.2924  < x  &&  x < 37.5193)
  ) {

    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
    xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i) {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (M_1_SQRT_2PI - temp) / y;

    do_del(x);
    swap_tail;
  }
  else { /* large x such that probs are 0 or 1 */
    if (x > 0) {
      *cum = R_D__1;
      *ccum = R_D__0;
    }
    else {
      *cum = R_D__0;
      *ccum = R_D__1;
    }
  }

#ifdef NO_DENORMS
  /* do not return "denormalized" -- we do in R */
  if(log_p) {
    if(*cum > -min) {
      *cum = -0.;
    }
    if(*ccum > -min) {
      *ccum = -0.;
    }
  }
  else {
    if(*cum < min) {
      *cum = 0.;
    }
    if(*ccum < min) {
      *ccum = 0.;
    }
  }
#endif
  return;
}

/**
 * Generates a sequence of probability points.
 *
 * This function was adapted from the ppoints R function
 * http://stat.ethz.ch/R-manual/R-devel/library/stats/html/ppoints.html
 *
 * @param int n
 *   The number of points generated
 * @param int a
 *   The offset fraction to be used; typically in (0,1). Recommend values
 *   are 3/8 if n < 10 and 1/2 otherwise.
 */
double * ppoints(int n, float a) {
  double * points = (double *) malloc(sizeof(double) * n);

  if (n > 0) {
    int i;
    for (i = 0; i < n; i++) {
      points[i] = (i + 1 - a) / (n + 1 - 2 * a);
    }
  }
  return points;
}
