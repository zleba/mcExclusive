#ifndef _Sigma2Process_
#define _Sigma2Process_

#include "ColorTensor.h"
#include "SpinTensor.h"
#include "Pythia8/Pythia.h"
#include "TopoSelector.h"

#include <string>
#include <complex>
#include <cmath>
#include <vector>
#include <map>


using namespace std;


class ColArr {
  public:
    void setColAcol( int cc1, int cc2, int cc3, int cc4, int cc5, int cc6, int cc7, int cc8) 
    { cols[0] = cc1; cols[1] = cc2; cols[2] = cc3; cols[3] = cc4;
      cols[4] = cc5; cols[5] = cc6; cols[6] = cc7; cols[7] = cc8; }

    void print()const{ cout << "("<<cols[0]<<","<<cols[1]<<") + ("<<cols[2]<<","<<cols[3]<<") -> " 
                            << "("<<cols[4]<<","<<cols[5]<<") + ("<<cols[6]<<","<<cols[7]<<")" << endl; }

    int c1() const { return cols[0]; }
    int c2() const { return cols[1]; }
    int c3() const { return cols[2]; }
    int c4() const { return cols[3]; }
    int c5() const { return cols[4]; }
    int c6() const { return cols[5]; }
    int c7() const { return cols[6]; }
    int c8() const { return cols[7]; }

    int calcHash() const 
      { int hash = 0; for(int i=0; i<8; ++i) hash += cols[i] << (3*i);
      return hash; }

    void substitute(map<int,int> rule) {
      for(int i = 0; i < 8; ++i) { 
        if(rule.count(cols[i]) > 0)
          cols[i] = rule[cols[i]];
      }
    }
    int maximum() const
    { int m = 0; for(int i =0; i < 8; ++i) m = max(m,cols[i]); return m; }

    void swapColAcol() 
    { swap(cols[0], cols[1]); swap(cols[2], cols[3]);
      swap(cols[4], cols[5]); swap(cols[6], cols[7]); }

    void swapCol1234()
    { swap(cols[0],cols[2]); swap(cols[1],cols[3]);
      swap(cols[4],cols[6]); swap(cols[5],cols[7]);}
    void addoffset(int ofs) {
      for(int i=0; i <8 ; ++i) if(cols[i] > 0) cols[i] += ofs;  }

    friend bool operator==(const ColArr &a, const ColArr &b);

  private:
    int cols[8];

};


class Sigma2Process {
public:

  Sigma2Process(); 

  Double CrossSection(bool isColorSinglet, bool isSpinSinglet);
  Double CrossSectionFull(bool isColorSinglet, complex<Double> lum[3],  complex<Double>  lumC[3]  );
  Double CrossSectionFullFast(bool isColorSinglet, complex<Double> lum[3],  complex<Double>  lumC[3]  );

  Double GetCrossSectionBase();

  virtual void Amplitudes() = 0;
  virtual void ColorMatrix() = 0;

#ifdef MADGRAPH
  virtual void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg)  = 0;
#endif

  virtual Double PyCrossSection() = 0;
  virtual ColArr PyColorFlow() = 0;
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr) = 0;

#ifdef MADGRAPH
  void MadGraphResult();
#endif

  string getName() const { return name; }
  unsigned int getNcolFlows() const { return nColFlow; }

	void SetEThetaPhi(Double E, Double Theta, Double Phi);

  void SetSTUPhi(Double _s, Double _t, Double _u,  Double Phi);
	void SetMasses(Double m1, Double m2, Double m3, Double m4) {mP1=m1;mP2=m2;mP3=m3;mP4=m4;}
  void SetScaleAlphas(Double _scale, Double _alphaS, Double _alphaEm);

  void SetIds   (int id1, int id2);
  void SetIdsOut(int id1, int id2) {idOut1 = id1; idOut2 = id2; }

  void SetColors(ColorTensor *cL, ColorTensor *cR) { colL = cL; colR = cR; }
  void SetSpins(SpinTensor *sL, SpinTensor *sR) { spinL = sL; spinR = sR; }

  void AddEmissionL(int idIn, int idOut, Double z, Double phi);
  void AddEmissionR(int idIn, int idOut, Double z, Double phi);

  void PrintInfo() const;
  void PrintAmplitudes() const;
  void PrintColMat() const;
  void Init();
	void RemoveShowers() { SetIds(idLeft,idRight); }

  enum dirs {pp = 0, pm = 1, mp = 2, mm = 3};

	complex<Double> GetAmplitude(int i, int j) const { return JAMP[i][j][0]; }

  TopoSelector topo;

protected:
  int idPythia; //Pythia process ID
  
  int nJamp;    //size of color matrix

  int nPol;     //number of final state polarisations

  int idLeft, idRight; //id of left and right incoming particle

  int idOut1, idOut2;  //Out1 <-> Left, Out2 <-> Right

  unsigned int nColFlow; //number of possible color flows

  string name; //process name


  vector<vector<Double> > colMatrix;  //color matrix;

  vector<vector<vector< complex<Double> > > > JAMP;
  vector<vector<vector< complex<Double> > > > JAMPmg; //MadGraph amplitudes
  vector<vector< bool > > NonZero;


  Double mP1, mP2, mP3, mP4;
  Double s,t,u, phi;
  int Nsym;

  const complex<Double> i_;

  Double alphaS;  //default 0.118
  Double alphaEm; 
  Double scale;   //scale in GeV
  Double g2S, e2;


  void DataQQ(Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1);
  void DataGQ(Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c);
  void DataGG( Double &A12, Double &A13, Double &A14,
                            Double &B234, Double &B243, Double &B342, Double &B432);




  void CopyConj(dirs d1, dirs d2, int par = 1)
  {
    dirs d1Inv, d2Inv;
    d1Inv = dirs(3 - d1);
    d2Inv = dirs(3 - d2);
    
    for(int i=0; i < nJamp; ++i)
      JAMP[d1Inv][d2Inv][i] = Double(par)*conj( JAMP[d1][d2][i] );

    NonZero[d1Inv][d2Inv] = NonZero[d1][d2];
  }



  void SetZero(dirs d1, dirs d2)
  {

    dirs d1Inv, d2Inv;
    d1Inv = dirs(3 - d1);
    d2Inv = dirs(3 - d2);

    for(int i=0; i < nJamp; ++i) {
      JAMP[d1][d2][i] = 0;
      JAMP[d1Inv][d2Inv][i] = 0;
    }

    NonZero[d1Inv][d2Inv] = false;
    NonZero[d1][d2] = false;

  }

private:

  complex<Double> PartialSigma(int finPol, int colID1, int colID2, bool isSinglet ) const;
  void GeneralAmplitudes(complex<Double> arr[][4], int finPol, int colID1, int colID2) const;
  complex<Double> MM(complex<Double> arr[][4], int a, int b) const;

  void GeneralSpinAmplitudes(complex<Double> SpinAmp[][4]) const;
  void GeneralSplitting(int a, int b, complex<Double> SplittingAmp[][4]) const;

  void GeneralSplittingAlt(int a, int b, complex<Double> SplittingAmp[][4]) const;
  void GeneralSpinAmplitudesAlt(complex<Double> SpinAmp[][4]) const;


  void AddEmissionGeneral(SpinTensor *spin, ColorTensor *col, 
                          int idIn, int idOut, Double z, Double phi);


  bool isColorSinglet;
  ColorTensor *colL,  *colR;
  SpinTensor  *spinL, *spinR;


};



#endif
