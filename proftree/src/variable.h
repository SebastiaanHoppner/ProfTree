#include <Rcpp.h>
#include <iostream>
#include <R.h>
#include <Rinternals.h>
using namespace std;

class variable{
public:
  bool isCat;
  double* sortedValues;
  int nCats;
  variable(int independendVar, int dependendVar, int nInst, double** data, int nCats);
  ~variable();
  void sortValues();
};
