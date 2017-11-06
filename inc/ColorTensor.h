#ifndef _ColorTensor_
#define _ColorTensor_


#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>

#include "Basics.h"
#include "ColorData.h"


using namespace std;
//using std::cout;
//using std::endl;

class ColorTensor {
  
public:

  void AddG();
  void AddQ();
  void AddQx();

  void Reset();
  
  void Print() const;
  string GetOrg() const {return Org;}
  string GetNow() const {return Now;}

  ColorTensor() { Reset(); }



  ColorTensor(const char *cParton) {
    SetStartingParton(cParton);
  }

  void SetStartingParton(const char *cParton)
  {
    Reset();
    string parton = cParton;
    Now = Org = parton;
    if(parton=="G") 
      A12 = 4;
    else if(parton=="Q" || parton=="Qx" ) 
      Ki1i2 = 1;
    else {
      cout << "Wrong parton type! : " << parton <<  endl;
      exit(1);
    }
  }



  vector<Double> GetGGVec() const;
  vector<Double> GetGQVec() const;
  vector<Double> GetQQVec() const;
  

private:
  void Flip(const ColorTensor &t);

  Double A12, A13, A14;

  Double B234, B243, B324;
  Double B342, B423, B432;

  Double Ki1j1, Ki1j2, Ki1i2;

  Double Dij;
  Double T1,  T2;
  Double T1c, T2c;

  string Now;
  string Org;


};


#endif
