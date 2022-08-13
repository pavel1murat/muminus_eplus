//-----------------------------------------------------------------------------
// closure approximation:
//
// P[0] : normalization constant
// P[1] : kMax
//-----------------------------------------------------------------------------
#ifndef __fun_closure__
#define __fun_closure__

double fun_closure(double* X, double* P) {
  double f{0};

  double kmax = P[1];

  double x = X[0]/kmax;
  
  if (x <= 1) {
    f = P[0]*(1-2*x+2*x*x)*x*(1-x)*(1-x);
  }

  return f;
}

#endif
