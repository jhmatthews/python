

/* These are numerical recipes routines used in the Monte Carlo programs 
04mar 	ksl	modified to make all of the calls ansi compatible*/

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

/* qromb integrates the double precision function func from a to b.  EPS limits the
   fractional accuracy; JMAX limits the total number number of times the
   trapezoidal rule can be called, and K determines the order of polint 
	01oct	ksl	Modified qromb so that the accuracy EPS could be 
			specified.
	05oct	ksl	Modified convergence condition for exiting trapzd loop
			to include = sign.  This change was in the ansi version
			of qromb routine, and was needed to integrate functions
			which were in some cases zero throughout the range. Not
			only did this seem to eliminate a number of error 
			returns, it sped up some portions of the program.	
*/


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
  return ss;			/* I set this to the best value but the user should beware */
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define FUNC(x) ((*func)(x))

double
trapzd (func, a, b, n)
     double a, b;
     double (*func) (double);	/* ANSI: double (*func)(double); */
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


/* Given arrays xa[] and y[a] defined from elements 1 to n and
   given x, polint returns a value y and an error estimate dy.  If
   P(x) is an polynomial of degree n-1, the results will be exact. Note
   that polint will return a value outside the specified range without 
   bothering to warn you */

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

double
zbrent (func, x1, x2, tol)
     double x1, x2, tol;
     double (*func) (double);	/* ANSI: double (*func)(double); */
{
  int iter;
  double a = x1, b = x2, c, d, e, min1, min2;
  double fa = (*func) (a), fb = (*func) (b), fc, p, q, r, s, tol1, xm;

  c = d = e = 0;		// to avoid -03 warning


  if (fb * fa > 0.0)
    {
      Log ("ZBRENT: Min %e & Max %e must bracket zero, but got %e & %e\n",
	      x1, x2, fa, fb);
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


void
spline (x, y, n, yp1, ypn, y2)
     double x[], y[], yp1, ypn, y2[];
     int n;
{
  int i, k;
  double p, qn, sig, un, *u, *vector ();
  void free_vector ();

  u = vector (1, n - 1);
  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0.0;
  else
    {
      y2[1] = -0.5;
      u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
    }
  for (i = 2; i <= n - 1; i++)
    {
      sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
      p = sig * y2[i - 1] + 2.0;
      y2[i] = (sig - 1.0) / p;
      u[i] =
	(y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] -
								     x[i -
								       1]);
      u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else
    {
      qn = 0.5;
      un =
	(3.0 / (x[n] - x[n - 1])) * (ypn -
				     (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
    }
  y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);
  for (k = n - 1; k >= 1; k--)
    y2[k] = y2[k] * y2[k + 1] + u[k];
  free_vector (u, 1, n - 1);
}

void
splint (xa, ya, y2a, n, x, y)
     double xa[], ya[], y2a[], x, *y;
     int n;
{
  int klo, khi, k;
  double h, b, a;

  klo = 1;
  khi = n;
  while (khi - klo > 1)
    {
      k = (khi + klo) >> 1;
      if (xa[k] > x)
	khi = k;
      else
	klo = k;
    }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
    {
      Error ("Bad XA input to routine SPLINT");
      exit (0);
    }
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  *y =
    a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] +
				 (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}



#include <math.h>
#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPS 1.0e-7

double
expint (int n, double x)
{
  int i, ii, nm1;
  double a, b, c, d, del, fact, h, psi, ans;

  nm1 = n - 1;
  ans = 0;			// To avoid -O3 warning -- ksl 04dec
  if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
    {
      Error ("bad arguments in expint");
    }
  else
    {
      if (n == 0)
	ans = exp (-x) / x;
      else
	{
	  if (x == 0.0)
	    ans = 1.0 / nm1;

	  else
	    {
	      if (x > 1.0)
		{
		  b = x + n;
		  c = 1.0 / FPMIN;
		  d = 1.0 / b;
		  h = d;
		  for (i = 1; i <= MAXIT; i++)
		    {
		      a = -i * (nm1 + i);
		      b += 2.0;
		      d = 1.0 / (a * d + b);
		      c = b + a / c;
		      del = c * d;
		      h *= del;
		      if (fabs (del - 1.0) < EPS)
			{
			  ans = h * exp (-x);
			  return ans;
			}
		    }
		  Error ("continued fraction failed in expint");
		}
	      else
		{
		  ans = (nm1 != 0 ? 1.0 / nm1 : -log (x) - EULER);
		  fact = 1.0;
		  for (i = 1; i <= MAXIT; i++)
		    {
		      fact *= -x / i;
		      if (i != nm1)
			del = -fact / (i - nm1);
		      else
			{
			  psi = -EULER;
			  for (ii = 1; ii <= nm1; ii++)
			    psi += 1.0 / ii;
			  del = fact * (-log (x) + psi);
			}
		      ans += del;
		      if (fabs (del) < fabs (ans) * EPS)
			return ans;
		    }
		  Error ("series failed in expint");
		}
	    }
	}
    }
  return ans;
}

#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER

/* The next two routines were written by ksl.  They were not part of
   the recipes programs which I had but I think they are what was intended */

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



/* A Numerical Recipes routine to find the zero of an equation.  A function
   which returns the value of the function and and the derivative at x is required.
  The root must be bracketed by x1 and x2.
 

History
	05jul	ksl	Replaced old version of rtsafe with newer ansi version
*/

#include <math.h>
#define MAXIT 100

double
rtsafe (void (*funcd) (double, double *, double *), double x1, double x2,
	double xacc)
{
  int j;
  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;

  (*funcd) (x1, &fl, &df);
  (*funcd) (x2, &fh, &df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    {
      Error ("rtsafe: x1 %8.3e x2 %8.3e fl %8.3e fh %8.3e \n", x1, x2, fl,
	     fh);
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
      if ((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0)
	  || (fabs (2.0 * f) > fabs (dxold * df)))
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
/*
#define MAXIT 100

double
rtsafe (funcd, x1, x2, xacc)
     double x1, x2, xacc;
     // void (*funcd) ();               // original 
     void (*funcd) (double, double *, double *);	// ANSI
{
  int j;
  double df, dx, dxold, f, fh, fl;
  double swap, temp, xh, xl, rts;

  (*funcd) (x1, &fl, &df);
  (*funcd) (x2, &fh, &df);
  if (fl * fh >= 0.0)
    {
      Error ("rtsafe: x1 %8.3e x2 %8.3e fl %8.3e fh %8.3e \n", x1, x2, fl,
	     fh);
      Error ("rtsafe: Root must be bracketed in RTSAFE\n");
      return (-INFINITY);
    }

  if (fl < 0.0)
    {
      xl = x1;
      xh = x2;
    }
  else
    {
      xh = x1;
      xl = x2;
      swap = fl;
      fl = fh;
      fh = swap;
    }
  rts = 0.5 * (x1 + x2);
  dxold = fabs (x2 - x1);
  dx = dxold;
  (*funcd) (rts, &f, &df);
  for (j = 1; j <= MAXIT; j++)
    {
      if ((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0)
	  || (fabs (2.0 * f) > fabs (dxold * df)))
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
	{
	  xl = rts;
	  fl = f;
	}
      else
	{
	  xh = rts;
	  fh = f;
	}
    }
  Error (" Maximum number of iterations exceeded in RTSAFE\n");
  mytrap ();
  exit (0);
}

*/

/* This is a Numerical Recipes routine.  Fiven a function f,
and bracketing triplet of abscissas a,b,c, such that it is known
that f(b) < f(a) and f(b)<f(c), golden returns the minimum value 
of the function, and the value xmin where the minimum occurs */

#undef MAXIT

#define R 0.61803399
#define CC (1.0-R)
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double
golden (ax, bx, cx, f, tol, xmin)
     double ax, bx, cx, tol, *xmin;
     double (*f) (double);	/* ANSI: double (*f)(double); */
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
