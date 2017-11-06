#ifndef _KMRlumi_
#define _KMRlumi_

#include "Basics.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <complex>
#include "Pythia8/Pythia.h"
#include "KMRlum.h"

using namespace std;


class KMRlumi {

public:

  KMRlumi(Double sqrtS, Double qMin, Double alphaS, Double mc, Double mb, Pythia8::BeamParticle *_beamPtr );

  static Double Uniform(Double a, Double b);
  Double gluon(int id, Double x, Double q2 ) const;

  pair<Double,Double> gluonPairChic(int id, Double x, Double q) const;

  Double LeifGluonChicUnintegrated(int id, Double x, Double q, Double y, Double Int) const;

  vector< complex<Double> > LumSqrtRandom(int id1, Double M, Double y, Double Mju, Double MjuLoop, const Vec2 &p1, const Vec2 &p2);
  void LuminosityRand(int id1, Double M, Double y, Double Mju, Double MjuLoop, Vec2 &p1, Vec2 &p2, complex<Double> lum[], complex<Double> lumC[] );

  Double LuminosityInc(int id1, Double M, Double y, Double Mju);


  Double SplittingAlphaS(int id, Double  LnScale2) const;

  Double alphaStrong(Double LnScale2) const;
  Double Splitting(int id, Double Kt) const;

  Double InsideSudakov(int id, Double LnQ2, Double LnM2) const;
  Double NoEmProbability(Double _x1, Double _x2, Double LnQ2, Double LnM2, bool backward=true);

  Double integral( Double (KMRlumi::*func)(int id, Double x) const, int id, Double xmin, Double xmax, Double &err ) const;
  void SaveIntegral(STATUS &st, Double xmin, Double xmax, Double Int) const;

  static  void GetYintGluon(Double x, Double &y, Double &Int);
  static  void GetYintQuark(Double x, Double &y, Double &Int);

private:

  Double Lambda5QCD, LnFreeze2;
  Double MinQt2;
  Double bSlope;
  Double LnCharmMass2,  LnBottomMass2;
  Double CharmMass,  BottomMass;

  Double x1,x2,s;

  Pythia8::BeamParticle *beamPtr;
  int idLeft, idRight;


};


#endif
