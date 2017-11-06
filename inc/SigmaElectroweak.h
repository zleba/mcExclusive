#include "Sigma2Process.h"

#define MADGRAPH

#ifdef MADGRAPH
extern "C" void amplitudesuubar2aa_(double *P, int *NHEL, complex<double> *JAMP);

extern "C" void amplitudesuu2mumu_(double *P, int *NHEL, complex<double> *JAMP);
extern "C" void amplitudesdd2mumu_(double *P, int *NHEL, complex<double> *JAMP);

extern "C" void amplitudesbb2h2bb_(double *P, int *NHEL, complex<double> *JAMP);
#endif



class Sigma2gg2gammagamma : public Sigma2Process
{

  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2gg2gammagamma() {
    name = "g g -> gamma gamma";
    idLeft = idRight = 21;
    idPythia = 205; 
    nJamp = 1;
    nPol = 4;
    Nsym = 2;
    nColFlow = 1;

	  topo.SetPs(1, 1, 1, 1 );
    Init();
  }
  void MgAmplitudes(double *, int *, complex<double> *) { }

};


class Sigma2qqbar2gammagamma: public Sigma2Process
{

  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);

  Sigma2qqbar2gammagamma() {
    name = "q qbar -> gamma gamma";
    idLeft  = 1;
    idRight =-1;
    idPythia = 204; 
    nJamp = 1;
    nPol = 4;
    mP3=mP4 = 0;
    Nsym = 2;
    nColFlow = 1;

    Init();
  }

    #ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    amplitudesuubar2aa_(P,NHEL, _JAMPmg);
  }
    #endif
  private:
    Double sigTU;

};





class Sigma2ffbar2ffbarsgmZ: public Sigma2Process
{

  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);


  Sigma2ffbar2ffbarsgmZ() {
    name = "q qbar -> f fbar";
    idLeft  = 2;
    idRight =-2;
    idPythia = 224; 
    nJamp = 1;
    nPol = 4;
    mP3=mP4 = 0;
    Nsym = 1;
    nColFlow = 1;

	  topo.SetPs(1, 1, 1, 1 );
    Init();
  }




		#ifdef MADGRAPH
  void MgAmplitudes(double *P, int *NHEL, complex<double> *_JAMPmg) {
    if(abs(idLeft) == 1 || abs(idLeft) == 3 || abs(idLeft) == 5)
      amplitudesdd2mumu_(P,NHEL, _JAMPmg);
    if(abs(idLeft) == 2 || abs(idLeft) == 4 || abs(idLeft) == 6)
      amplitudesuu2mumu_(P,NHEL, _JAMPmg);
  }
		#endif


};


class Sigma1gg2H: public Sigma2Process
{

  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);


  Sigma1gg2H() {
    name = "g g -> b bbar (higgs)";
    idLeft  = 21;
    idRight = 21;
    idPythia = 902; 
    nJamp = 1;
    nPol = 4;
    mP3=mP4 = 4.2;
    Nsym = 1;
    nColFlow = 1;
	  topo.SetPs(1, 2000, 1, 3 );

    Init();
  }




		#ifdef MADGRAPH
  void MgAmplitudes(double *, int *, complex<double> *) {

  }
		#endif


};


class  Sigma1ffbar2H: public Sigma2Process
{

  public:

  virtual void Amplitudes();
  virtual void ColorMatrix();
  virtual Double PyCrossSection();
  virtual ColArr PyColorFlow();
	virtual void PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr);


  Sigma1ffbar2H() {
    name = "q qbar -> b bbar (higgs)";
    idLeft  =  5;
    idRight = -5;
    idPythia = 901; 
    nJamp = 1;
    nPol = 4;
    mP3=mP4 = 4.7;
    Nsym = 1;
    nColFlow = 1;
	  topo.SetPs(1, 1, 1, 1 );

    Init();
  }


		#ifdef MADGRAPH
  void MgAmplitudes( double *P, int *NHEL, complex<double> *_JAMPmg ) {
    amplitudesbb2h2bb_(P, NHEL, _JAMPmg);
  }
		#endif


};
