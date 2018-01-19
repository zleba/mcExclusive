#ifndef _KMRlum_
#define _KMRlum_

#include "Basics.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <complex>
#include "Pythia8/Pythia.h"

using namespace std;

class Vec2;
struct STATUS;


class KMRlum {

public:

  KMRlum(Double sqrtS, Double qMin, Double alphaS, Double mc, Double mb, Pythia8::BeamParticle *_beamPtr );

  static Double Uniform(Double a, Double b);
  Double gluon(Double x, Double q2 ) const;

  pair<Double,Double> gluonPairChic(Double x, Double q) const;

  Double LeifGluonChicUnintegrated(Double x, Double q, Double y, Double Int) const;

  vector< complex<Double> > LumSqrtRandom(Double M, Double y, Double Mju, Double MjuLoop, const Vec2 &p1, const Vec2 &p2);
  void LuminosityRand(Double M, Double y, Double Mju, Double MjuLoop, Vec2 &p1, Vec2 &p2, complex<Double> lum[], complex<Double> lumC[] );

  Double LuminosityInc(Double M, Double y, Double Mju);


  Double SplittingAlphaS(Double  LnScale2) const;

  Double alphaStrong(Double LnScale2) const;
  Double Splitting(Double Kt) const;

  Double InsideSudakov(Double LnQ2, Double LnM2) const;
  Double NoEmProbability(Double _x1, Double _x2, Double LnQ2, Double LnM2, bool backward=true);

  Double integral( Double (KMRlum::*func)(Double x) const, Double xmin, Double xmax, Double &err ) const;
  void SaveIntegral(STATUS &st, Double xmin, Double xmax, Double Int) const;

  static  void GetYint(Double x, Double &y, Double &Int);

private:

  Double Lambda5QCD, LnFreeze2;
  Double MinQt2;
  Double bSlope;
  Double LnCharmMass2,  LnBottomMass2;
  Double CharmMass,  BottomMass;

  Double x1,x2,s;

  Pythia8::BeamParticle *beamPtr;
  static Pythia8::Rndm *rndmPtr;


};


inline bool close(double a, double b)
{
  return abs(a - b) < 1e-12;
}


struct DATA {
  Double xmin, xmax, Int;

  DATA(Double _xmin, Double _xmax, Double _Int)
             : xmin(_xmin), xmax(_xmax), Int(_Int) {}

  inline bool isSame(Double _xmin, Double _xmax)
  { return ( close(xmin,_xmin) && close(xmax,_xmax) ); }

};

struct STATUS {
  Double Mred;
  vector<DATA> points;
  void print() {
    cout << "Start" << endl;
    for(unsigned int i=0; i < points.size(); ++i)
      cout << points[i].xmin<<" "<<points[i].xmax<<" : "<< points[i].Int << endl;
    cout << "End" << endl;
  }
};





class Vec2 {

public:

  Vec2(Double xIn = 0., Double yIn = 0.)
      : xx(xIn), yy(yIn) { }

  Vec2(const Vec2& v) : xx(v.xx), yy(v.yy) { }

  void setRPhi(Double  r, Double phi )
    { xx=r*cos(phi); yy=r*sin(phi); }
  void setXY(Double  x, Double y )
    { xx=x; yy=y; }



  void x(Double xIn) {xx = xIn;}
  void y(Double yIn) {yy = yIn;}

  Double x() const {return xx;}
  Double y() const {return yy;}
  Double px() const {return xx;}
  Double py() const {return yy;}

  Double norm()  const {return sqrt(xx*xx+yy*yy);}
  Double norm2() const {return (xx*xx+yy*yy);}

  Vec2 operator-()
    { return Vec2(-xx,-yy);}
  Vec2& operator+=(const Vec2& v)
    {xx += v.xx; yy += v.yy; return *this;}
  Vec2& operator-=(const Vec2& v)
    {xx -= v.xx; yy -= v.yy; return *this;}
  Vec2& operator*=(Double f)
    {xx *= f; yy *= f; return *this;}
  Vec2& operator/=(Double f)
    {xx /= f; yy /= f; return *this;}

  // Operator overloading with friends
  friend Vec2 operator+(const Vec2& v1, const Vec2& v2);
  friend Vec2 operator-(const Vec2& v1, const Vec2& v2);
  friend Vec2 operator*(Double f, const Vec2& v1);
  friend Vec2 operator*(const Vec2& v1, Double f);
  friend Vec2 operator/(const Vec2& v1, Double f);
  friend Double operator*(const Vec2& v1, const Vec2& v2);
  friend Double cross(const Vec2& v1, const Vec2& v2);


private:
  Double xx, yy;


};

inline Vec2 operator+(const Vec2& v1, const Vec2& v2)
  {Vec2 v = v1 ; return v += v2;}

inline Vec2 operator-(const Vec2& v1, const Vec2& v2)
  {Vec2 v = v1 ; return v -= v2;}

inline Vec2 operator*(Double f, const Vec2& v1)
  {Vec2 v = v1; return v *= f;}

inline Vec2 operator*(const Vec2& v1, Double f)
  {Vec2 v = v1; return v *= f;}

inline Vec2 operator/(const Vec2& v1, Double f)
  {Vec2 v = v1; return v /= f;}

inline Double operator*(const Vec2& v1, const Vec2& v2)
  {return  v1.xx*v2.xx + v1.yy*v2.yy;}


inline Double cross(const Vec2& v1, const Vec2& v2)
  {return  v1.xx*v2.yy - v1.yy*v2.xx;}




#endif
