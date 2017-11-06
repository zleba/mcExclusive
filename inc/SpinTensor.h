#ifndef _SpinTensor_
#define _SpinTensor_


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <complex>

#include "Basics.h"

using namespace std;


class SpinTensor {


public:
  SpinTensor() { Reset(); }

  //inline complex<Double> Element(int index, int k, int l) const;
  inline complex<Double> Element(int index, int k, int l) const
  { return P[2*k+l][index]; }

  inline complex<Double> Element(int in, int jn, int i1, int j1) const
  { return conj(P[2*i1+j1][2*in+jn]); }


  void Reset();
  void PrintPnow() const;

  void AddEmission(int idIn, int idOut, Double z, Double ph);

  static Double ClassicalSplitting(int idIn, int idOut, Double z);

public:
  complex<Double> P[4][4]; //Matrix

  complex<Double> Pnow[4][4];

  inline bool In(int x, int a, int b) {return (x >= a && x <= b); }

  void CalculateMatrixGG(Double z, Double ph);
  void CalculateMatrixGQ(Double z, Double ph);
  void CalculateMatrixQQ(Double z, Double ph);
  void CalculateMatrixQG(Double z, Double ph);


};


#endif
