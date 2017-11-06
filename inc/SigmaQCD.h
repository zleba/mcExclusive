#ifndef _SigmaQCD_
#define _SigmaQCD_

#include "Sigma2Process.h"

#ifdef MADGRAPH
extern "C" void amplitudesgg2gg_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesuubar2ddbar_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesuubar2uubar_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesud2ud_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesuu2uu_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesgg2uu_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesgg2tt_(double *P, int *NHEL, complex<double> *JAMP, double *MQ);
extern "C" void amplitudesuu2gg_(double *P, int *NHEL, complex<double> *JAMP);

extern "C" void amplitudesgu2gu_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesuxg2uxg_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesug2ug_(double *P, int *NHEL, complex<double> *JAMP);

extern "C" void amplitudesuxux2uxux_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesudx2udx_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesuxdx2uxdx_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesuux2ttx_(double *P, int *NHEL, complex<double> *JAMP, double *MQ);
#endif





class Sigma2gg2gg : public Sigma2Process
{
  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2gg2gg() {
    name = "g g -> g g";
    idLeft = idRight = 21;
    idPythia = 111;
    nJamp = 6;
    nPol = 4;
    Nsym = 2;
    nColFlow = 6;
	  topo.SetPs(1, 74, 3.4, 5.7 );
    Init();
  }
    #ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    amplitudesgg2gg_(P, NHEL, _JAMPmg);
  }
    #endif
  private:
    Double sigTS, sigUS, sigTU;
    Double sigSum;

};


class Sigma2gg2qqbar: public Sigma2Process
{
  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  void printInv() {cout << "sigTS + sigUS = sigSum  : " << sigTS<<" + "<<sigUS <<" = "<<sigSum<< endl; }

  Sigma2gg2qqbar() {
    name = "g g -> q qbar (uds)";
    idLeft = idRight = 21;
    idPythia = 112;
    nJamp = 2;
    nPol = 4;
    nColFlow = 2;
	  topo.SetPs(1, 0, 0, 3 );
    Init();
  }
    #ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    amplitudesgg2uu_(P, NHEL, _JAMPmg);
  }
    #endif

  private:
    Double sigTS, sigUS;
    Double sigSum;

};

class  Sigma2qg2qg : public Sigma2Process
{
  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2qg2qg() {
    name = "q g -> q g";
    //idLeft = idRight = 21;
    idPythia = 113;
    nJamp = 2;
    nPol = 4;
    nColFlow = 2;
	  topo.SetPs(1, 0, 8, 3 );
    Init();
  }

		#ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    //amplitudesgu2gu_(P, NHEL, _JAMPmg);
    //return;
    if(idLeft > 0)
      amplitudesug2ug_(P,NHEL, _JAMPmg);
    else
      amplitudesuxg2uxg_(P,NHEL, _JAMPmg);
  }
		#endif

  private:
    Double sigTS, sigTU;
    Double sigSum;

};


class  Sigma2qq2qq : public Sigma2Process
{
  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2qq2qq() {
    name = "q q(bar) -> q q(bar)";
    //idLeft = idRight = 21;
    idPythia = 114;
    nJamp = 2;
    nPol = 4;
    nColFlow = 2;
	  topo.SetPs(1, 0, 0, 0 );
    Init();
  }
		#ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {

    if(idLeft == idRight && idLeft > 0)
      amplitudesuu2uu_(P, NHEL, _JAMPmg);
    else if(idLeft == idRight && idLeft < 0)
      amplitudesuxux2uxux_(P,NHEL,_JAMPmg);
    else if(idLeft ==-idRight)
      amplitudesuubar2uubar_(P, NHEL, _JAMPmg);
    else if(idLeft != idRight && idLeft > 0 && idRight > 0)
      amplitudesud2ud_(P, NHEL, _JAMPmg);
    else if(idLeft != idRight && idLeft > 0 && idRight < 0)
      amplitudesudx2udx_(P, NHEL, _JAMPmg);
    else if(idLeft != idRight && idLeft < 0 && idRight < 0)
      amplitudesuxdx2uxdx_(P, NHEL, _JAMPmg);
  }
		#endif

  private:
    Double sigT, sigU, sigS, sigTU, sigST;
    Double sigSum;

};


class  Sigma2qqbar2gg : public Sigma2Process
{
  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2qqbar2gg() {
    name = "q qbar -> g g";
    //idLeft = idRight = 21;
    idPythia = 115;
    nJamp = 2;
    nPol = 4;
    Nsym = 2;
    nColFlow = 2;
	  topo.SetPs(1, 0, 0, 0 );
    Init();
  }
		#ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    amplitudesuu2gg_(P, NHEL, _JAMPmg);
  }
		#endif

  private:
    Double sigTS, sigUS;
    Double sigSum;

};

class Sigma2qqbar2qqbarNew : public Sigma2Process
{
  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2qqbar2qqbarNew() {
    name = "q qbar -> q' qbar' (uds)";
    //idLeft = idRight = 21;
    idPythia = 116;
    nJamp = 2;
    nPol = 4;
    nColFlow = 1;
	  topo.SetPs(1, 0, 0, 0 );
    Init();
  }
		#ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    amplitudesuubar2ddbar_(P, NHEL, _JAMPmg);
  }
		#endif

  private:
    Double sigS, sigT, sigST;
};


//Pythia needs to be add
class Sigma2gg2QQbar : public Sigma2Process
{
  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2gg2QQbar() {
    name = "g g -> Q Qbar";
    //idLeft = idRight = 21;
    idPythia = -1; //will be know later
    nJamp = 2;
    nPol = 4;
    mP3=mP4 = 170;
    nColFlow = 2;
	  topo.SetPs(1, 5, 1, 2 );
    Init();
  }
		#ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    double mass = mP3;
    amplitudesgg2tt_(P, NHEL, _JAMPmg, &mass);
  }
		#endif

  private:
    Double sigTS, sigUS;
    Double sigSum;

};


class Sigma2qqbar2QQbar : public Sigma2Process
{

  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2qqbar2QQbar() {
    name = "q qbar -> Q Qbar";
    //idLeft = idRight = 21;
    idPythia = -1; //will be know later
    nJamp = 2;
    nPol = 4;
    mP3=mP4 = 170;
    nColFlow = 1;
	  topo.SetPs(1, 0, 0, 0 );
    Init();
  }
		#ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    double mass = mP3;
    amplitudesuux2ttx_(P, NHEL, _JAMPmg, &mass);
  }
		#endif

  private:
    Double sigS;

};



#endif
