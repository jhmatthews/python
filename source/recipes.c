
/***********************************************************/
/** @file  recipes.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  These are collection of Numberical Recipes routines
 * used in Python
 *
 *
 * The routines have been modified to make them ansi compatable
 *
 ***********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include "atomic.h"
#include <math.h>
#include "recipes.h"
#include "log.h"



#define EPS 1.0e-6
#define JMAX 100
#define JMAXP JMAX+1
#define K 5


/**********************************************************/
/**
 * @brief      Integrate a function between a and b
 *
 *
 * @param [in] func  The function which one wishes a root for
 * @param [in] a  one end of the interval
 * @param [in] b  the other end of the interval
 * @param [in] eps   A tolerence
 * @return   The root of the function f                           .
 *
 * @details
 * A Numerical Recipes routine to integrates the double precision function func from a to b.
 *
 * ### Notes ###
 * JMAX limits the total number number of times the
 * trapezoidal rule can be called, and K determines the order of polint
 *
 *
 * This routine has been modified from that that appears
 * in Numberical recipes so that the tolerance could be specified.
 * Additioanlly, the convergence criterion for exiting the
 * exiting trapzd loop
 * to include = sign.  This change was in the ansi version
 * of qromb routine, and was needed to integrate functions
 * which were in some cases zero throughout the range. Not
 * only did this seem to eliminate a number of error
 * returns, it sped up some portions of the program.
 *
 *
 **********************************************************/


double
qromb (func, a, b, eps)
     double a, b;
     double (*func) (double);
     double eps;
{
  double ss, dss, trapzd ();
  double s[JMAXP + 1], h[JMAXP + 1];
  int j;
  void polint ();

  ss = 0.0;
  if (a >= b)
  {
    Error ("Error qromb: a %e>=b %e\n", a, b);
  }

  h[1] = 1.0;
  for (j = 1; j <= JMAX; j++)
  {
    s[j] = trapzd (func, a, b, j);
    if (j >= K)
    {
      polint (&h[j - K], &s[j - K], K, 0.0, &ss, &dss);
      // ksl --56g --05Oct -- Changed to include = sign
      if (fabs (dss) <= eps * fabs (ss))
      {
        recipes_error = 0;
        return ss;
      }
    }
    s[j + 1] = s[j];
    h[j + 1] = 0.25 * h[j];
  }
  Error ("Too many steps in routine QROMB\n");
  recipes_error = -1;
  return ss;                    /* I set this to the best value but the user should beware */
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define FUNC(x) ((*func)(x))

double
trapzd (func, a, b, n)
     double a, b;
     double (*func) (double);   /* ANSI: double (*func)(double); */
     int n;
{
  double x, tnm, sum, del;
  static double s;
  static int it;
  double ssss;
  int j;

  if (n == 1)
  {
    it = 1;
    return (s = 0.5 * (b - a) * (FUNC (a) + FUNC (b)));
  }
  else
  {
    tnm = it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    for (sum = 0.0, j = 1; j <= it; j++, x += del)
      sum += ssss = FUNC (x);
    it *= 2;
    s = 0.5 * (s + (b - a) * sum / tnm);
    return s;
  }
}


/*******************************************
 * Given arrays xa[] and y[a] defined from elements 1 to n and
   given x, polint returns a value y and an error estimate dy.  If
   P(x) is an polynomial of degree n-1, the results will be exact. Note
   that polint will return a value outside the specified range without
   bothering to warn you
**************************************************/

void
polint (xa, ya, n, x, y, dy)
     double xa[], ya[], x, *y, *dy;
     int n;
{
  int i, m, ns = 1;
  double den, dif, dift, ho, hp, w;
  double *c, *d, *vector ();
  void free_vector ();


  dif = fabs (x - xa[1]);
  c = vector (1, n);
  d = vector (1, n);
  for (i = 1; i <= n; i++)
  {
    if ((dift = fabs (x - xa[i])) < dif)
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  *y = ya[ns--];
  for (m = 1; m < n; m++)
  {
    for (i = 1; i <= n - m; i++)
    {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      if ((den = ho - hp) == 0.0)
        Error ("Error in routine POLINT\n");
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    *y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
  }
  free_vector (d, 1, n);
  free_vector (c, 1, n);
}



#define ITMAX 100
#define EPS 3.0e-8

/**********************************************************/
/**
 * @brief      A routine that finds the rooot of  a function f
 *
 *
 * @param [in] *funcd  The function which one wishes a root for
 * @param [in] x1  one end of the interval
 * @param [in] x2  the other end of the interval
 * @param [in] tol   A tolerence for the root
 * @return   The root of the function f                           .
 *
 * @details
 * A Numerical Recipes routine to find the root of a function known to
 * be between x1 and x2
 *
 * ### Notes ###
 *
 **********************************************************/


double
zbrent (func, x1, x2, tol)
     double x1, x2, tol;
     double (*func) (double);   /* ANSI: double (*func)(double); */
{
  int iter;
  double a = x1, b = x2, c, d, e, min1, min2;
  double fa = (*func) (a), fb = (*func) (b), fc, p, q, r, s, tol1, xm;

  c = d = e = 0;                // to avoid -03 warning


  if (fb * fa > 0.0)
  {
    Log ("ZBRENT: Min %e & Max %e must bracket zero, but got %e & %e\n", x1, x2, fa, fb);
  }
  fc = fb;
  for (iter = 1; iter <= ITMAX; iter++)
  {
    if (fb * fc > 0.0)
    {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs (fc) < fabs (fb))
    {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * EPS * fabs (b) + 0.5 * tol;
    xm = 0.5 * (c - b);
    if (fabs (xm) <= tol1 || fb == 0.0)
      return b;
    if (fabs (e) >= tol1 && fabs (fa) > fabs (fb))
    {
      s = fb / fa;
      if (a == c)
      {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0)
        q = -q;
      p = fabs (p);
      min1 = 3.0 * xm * q - fabs (tol1 * q);
      min2 = fabs (e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2))
      {
        e = d;
        d = p / q;
      }
      else
      {
        d = xm;
        e = d;
      }
    }
    else
    {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs (d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs (tol1) : -fabs (tol1));
    fb = (*func) (b);
  }
  Error ("Maximum number of iterations exceeded in ZBRENT\n");
  return b;
}

#undef ITMAX
#undef EPS




#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER

/******************************
 * The next two routines were written by ksl.  They were not part of
   the recipes programs which I had but I think they are what was intended
********************************/

double *
vector (i, j)
     int i, j;
{
  double dummy, *d;
  d = calloc (sizeof (dummy), (j - i + 1) + 1);
  return (d);
}

void
free_vector (a, i, j)
     double *a;
     int i, j;
{
  free (a);
}




#include <math.h>
#define MAXIT 100

/**********************************************************/
/**
 * @brief      A routine that mimimizes a function f
 *
 *
 * @param [in] *funcd
 * @param [in] x1
 * @param [in] x2
 * @param [in] xacc
 * @return   The minimum of the function f                           .
 *
 * @details
 * A Numerical Recipes routine to find the zero of an equation.  A function
 * which returns the value of the function and and the derivative at x is required.
 * The root must be bracketed by x1 and x2.
 *
 * ### Notes ###
 *
 **********************************************************/

double
rtsafe (void (*funcd) (double, double *, double *), double x1, double x2, double xacc)
{
  int j;
  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;

  (*funcd) (x1, &fl, &df);
  (*funcd) (x2, &fh, &df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
  {
    Error ("rtsafe: x1 %8.3e x2 %8.3e fl %8.3e fh %8.3e \n", x1, x2, fl, fh);
    Error ("rtsafe: Root must be bracketed in RTSAFE\n");
  }
  if (fl == 0.0)
    return x1;
  if (fh == 0.0)
    return x2;
  if (fl < 0.0)
  {
    xl = x1;
    xh = x2;
  }
  else
  {
    xh = x1;
    xl = x2;
  }
  rts = 0.5 * (x1 + x2);
  dxold = fabs (x2 - x1);
  dx = dxold;
  (*funcd) (rts, &f, &df);
  for (j = 1; j <= MAXIT; j++)
  {
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) || (fabs (2.0 * f) > fabs (dxold * df)))
    {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts)
        return rts;
    }
    else
    {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }
    if (fabs (dx) < xacc)
      return rts;
    (*funcd) (rts, &f, &df);
    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  Error (" Maximum number of iterations exceeded in RTSAFE\n");
  return 0.0;
}



#undef MAXIT

#define R 0.61803399
#define CC (1.0-R)
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

/**********************************************************/
/**
 * @brief      A routine that mimimizes a function f
 *
 *
 * @param [in] ax
 * @param [in] bx
 * @param [in] cx
 * @param [in] f  the function to minimize
 * @param [in] tol A tolerance
 * @param [out] * xmin - The place where the minimum occors
 * @return   The minimum of the function f                           .
 *
 * @details
 * This is a Numerical Recipes routine.  Fiven a funYtion f,
 * and bracketing triplet of abscissas a,b,c, such that it is known
 * that f(b) < f(a) and f(b)<f(c), golden returns the minimum value
 * of the function, and the value xmin where the minimum occurs
 *
 * ### Notes ###
 *
 **********************************************************/

double
golden (ax, bx, cx, f, tol, xmin)
     double ax, bx, cx, tol, *xmin;
     double (*f) (double);      /* ANSI: double (*f)(double); */
{
  double f0, f1, f2, f3, x0, x1, x2, x3;

  x0 = ax;
  x3 = cx;
  if (fabs (cx - bx) > fabs (bx - ax))
  {
    x1 = bx;
    x2 = bx + CC * (cx - bx);
  }
  else
  {
    x2 = bx;
    x1 = bx - CC * (bx - ax);
  }
  f1 = (*f) (x1);
  f2 = (*f) (x2);
  while (fabs (x3 - x0) > tol * (fabs (x1) + fabs (x2)))
  {
    if (f2 < f1)
    {
    SHFT (x0, x1, x2, R * x1 + CC * x3) SHFT (f0, f1, f2, (*f) (x2))}
    else
    {
    SHFT (x3, x2, x1, R * x2 + CC * x0) SHFT (f3, f2, f1, (*f) (x1))}
  }
  if (f1 < f2)
  {
    *xmin = x1;
    return f1;
  }
  else
  {
    *xmin = x2;
    return f2;
  }
}

#undef CC
#undef R
