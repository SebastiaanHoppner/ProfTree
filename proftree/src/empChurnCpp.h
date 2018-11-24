#include <Rcpp.h> // R/C++ interface class library
using namespace Rcpp;

double B(double x, double a, double b);
void FilterRch(NumericVector F1, NumericVector F0, int LIP, int RIP);
int FurthestPoint(NumericVector F1, NumericVector F0, int LIP, int RIP);
IntegerVector order_asc(NumericVector v);
List ROCconvexhull(NumericVector scores, IntegerVector classes);
double empChurnCpp(NumericVector scores, IntegerVector classes, double alpha, double beta, double clv, double d, double f);