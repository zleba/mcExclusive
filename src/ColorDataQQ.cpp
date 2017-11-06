#include "ColorData.h"
#include "ColorTensor.h"
#include <vector>
using std::vector;

void SingletQQ::DataQQ(const ColorTensor *colL, const ColorTensor *colR, 
                       Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1)
{

  const int BaseDl1l1c[3][3] = {
    { 9, 3, 3 },
    { 3, 1, 1 },
    { 3, 1, 1 }
  };
  const Double factBaseDl1l1c = 3.;


  const int BaseDl1r1c[3][3] = {
    { 3, 1, 1 },
    { 1, 1, 3 },
    { 1, 3, 1 }
  };
  const Double factBaseDl1r1c = 3.;

 
  const int BaseDl1r1[3][3] = {
    { 3, 1, 1 },
    { 1, 3, 1 },
    { 1, 1, 3 }
  };
  const Double factBaseDl1r1  = 3.;


  vector<Double> colLar, colRar;

  colLar = colL->GetQQVec(); //quark -> quark
  colRar = colR->GetQQVec(); //quark -> quark

  Dl1l1c = Dl1r1c = Dl1r1 = 0;


  for(int i=0; i < 3; ++i)
    for(int j=0; j < 3; ++j) {
      Double colIJ = colLar[i]*colRar[j];
      Dl1l1c += BaseDl1l1c[i][j] * colIJ;
      Dl1r1  += BaseDl1r1[i][j]  * colIJ;
      Dl1r1c += BaseDl1r1c[i][j] * colIJ;
    }


  Dl1l1c *= factBaseDl1l1c;
  Dl1r1  *= factBaseDl1r1;
  Dl1r1c *= factBaseDl1r1c;

}



void InclusiveQQ::DataQQ(const ColorTensor *colL, const ColorTensor *colR, 
                         Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1)
{

  const int BaseDl1l1c[3][3] = {
    { 9, 3, 3 },
    { 3, 1, 1 },
    { 3, 1, 1 }
  };

  const Double factBaseDl1l1c = 9.;
  const Double factBaseDl1r1c = 3.;
  const Double factBaseDl1r1  = 3.;



  vector<Double> colLar, colRar;

  colLar = colL->GetQQVec(); //quark -> quark
  colRar = colR->GetQQVec(); //quark -> quark


  Dl1l1c = 0;


  for(int i=0; i < 3; ++i)
    for(int j=0; j < 3; ++j) {
      Double colIJ = colLar[i]*colRar[j];
      Dl1l1c += BaseDl1l1c[i][j] * colIJ;
    }

  
  Dl1r1c = Dl1r1 = Dl1l1c;

  Dl1l1c *= factBaseDl1l1c;
  Dl1r1  *= factBaseDl1r1;
  Dl1r1c *= factBaseDl1r1c;

}




void SingletQQ::DataGG(const ColorTensor *colL, const ColorTensor *colR, 
                       Double &A12,  Double &A13,  Double &A14,
                       Double &B234, Double &B243, Double &B342, Double &B432)
{

  const int BaseA12[5][5] = {
    { 36, 12, 12, 12, 12 },
    { 12,  2, 11, -4, 32 },
    { 12, 11,  2, 32, -4 },
    { 12, -4, 32,  2, 11 },
    { 12, 32, -4, 11,  2 }
  };
  const Double factBaseA12 = 1./24;



  const int BaseA13[5][5] = {
    { 9, 3, 3, 3, 3 },
    { 3, 1, 1, 1, 1 },
    { 3, 1, 1, 1, 1 },
    { 3, 1, 1, 1, 1 },
    { 3, 1, 1, 1, 1 }
  };
  const Double factBaseA13 = 4./3;



  const int BaseA14[5][5] = {
    { 36, 12, 12, 12, 12 },
    { 12, 11,  2, 32, -4 },
    { 12,  2, 11, -4, 32 },
    { 12, 32, -4, 11,  2 },
    { 12, -4, 32,  2, 11 }
  };
  const Double factBaseA14 = 1./24;


  const int BaseB234[5][5] = {
    { -12, -4, -4, -4, -4 },
    {  -4, -3, -3, -8, -8 },
    {  -4, -3, -3, -8, -8 },
    {  -4, -8, -8, -3, -3 },
    {  -4, -8, -8, -3, -3 }
  };
  const Double factBaseB234 = 1./24;


  const int BaseB243[5][5] = {
    { 96, 32, 32, 32, 32 },
    { 32,  6, 27, -8, 76 },
    { 32, 12,  6, 16, -8 },
    { 32, -8, 76,  6, 27 },
    { 32, 16, -8, 12,  6 }
  };
  const Double factBaseB243 = 1./24;


  const int BaseB342[5][5] = {
    { 96, 32, 32, 32, 32 },
    { 32,  6, 12, -8, 16 },
    { 32, 27,  6, 76, -8 },
    { 32, -8, 16,  6, 12 },
    { 32, 76, -8, 27,  6 }
  };
  const Double factBaseB342 = 1./24;


  const int BaseB432[5][5] = {
    { -12, -4, -4, -4, -4 },
    {  -4, -3, -3, -8, -8 },
    {  -4, -3, -3, -8, -8 },
    {  -4, -8, -8, -3, -3 },
    {  -4, -8, -8, -3, -3 }
  };
  const Double factBaseB432 = 1./24;





  vector<Double> colLar, colRar;

  colLar = colL->GetGQVec(); 
  colRar = colR->GetGQVec();

  A12=0; A13=0; A14=0;
  B234=0; B243=0;
  B342=0; B432=0;


  for(int i=0; i < 5; ++i)
    for(int j=0; j < 5; ++j) {
      Double colIJ = colLar[i]*colRar[j];
      A12  += BaseA12[i][j]  * colIJ;
      A13  += BaseA13[i][j]  * colIJ;
      A14  += BaseA14[i][j]  * colIJ;
      B234 += BaseB234[i][j] * colIJ;
      B243 += BaseB243[i][j] * colIJ;
      B342 += BaseB342[i][j] * colIJ;
      B432 += BaseB432[i][j] * colIJ;
    }

  A12 *= factBaseA12;
  A13 *= factBaseA13;
  A14 *= factBaseA14;
  B234*= factBaseB234;
  B243*= factBaseB243;
  B342*= factBaseB342;
  B432*= factBaseB432;

}


void InclusiveQQ::DataGG(const ColorTensor *colL, const ColorTensor *colR, 
                Double &A12,  Double &A13,  Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432)
{

  const int BaseA12[5][5] = {
    { 9, 3, 3, 3, 3 },
    { 3, 1, 1, 1, 1 },
    { 3, 1, 1, 1, 1 },
    { 3, 1, 1, 1, 1 },
    { 3, 1, 1, 1, 1 }
  };

  const Double factBaseA12 = 1./2;
  const Double factBaseA13 = 4.;
  const Double factBaseA14 = 1./2;


  const Double factBaseB234 = -1./6;
  const Double factBaseB243 =  4./3;
  const Double factBaseB342 =  4./3;
  const Double factBaseB432 = -1./6;



  vector<Double> colLar, colRar;

  colLar = colL->GetGQVec(); 
  colRar = colR->GetGQVec();

  A12=0;

  for(int i=0; i < 5; ++i)
    for(int j=0; j < 5; ++j) {
      A12  += BaseA12[i][j]  * colLar[i]*colRar[j];
    }

  A13=A14 = A12;
  B234=B243=B342=B432 = A12;


  A12 *= factBaseA12;
  A13 *= factBaseA13;
  A14 *= factBaseA14;
  B234*= factBaseB234;
  B243*= factBaseB243;
  B342*= factBaseB342;
  B432*= factBaseB432;

}



void InclusiveQQ::DataGQ(const ColorTensor *colL, const ColorTensor *colR, 
                         Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c)
{

  const int BaseDelta[5][3] = {
    { 9, 3, 3 },
    { 3, 1, 1 },
    { 3, 1, 1 },
    { 3, 1, 1 },
    { 3, 1, 1 }
  };
  const Double factBaseDelta   = 6;
  const Double factBaseTl1cl1C = 2;
  const Double factBaseTl1l1c  = 2;

  vector<Double> colGar, colQar;

  if(colL->GetOrg() == "G" && colR->GetOrg() != "G") {
    colGar = colL->GetGQVec();
    colQar = colR->GetQQVec();
  }
  else if(colL->GetOrg() != "G" && colR->GetOrg() == "G") {
    colQar = colL->GetQQVec();
    colGar = colR->GetGQVec();
  }

  Delta=0;

  for(int i=0; i < 5; ++i)
    for(int j=0; j < 3; ++j) {
      Delta  += BaseDelta[i][j]  * colGar[i]*colQar[j];
    }

  Tl1cl1C=Tl1l1c = Delta;

  Delta   *= factBaseDelta;
  Tl1cl1C *= factBaseTl1cl1C;
  Tl1l1c  *= factBaseTl1l1c;

}



void SingletQQ::DataGQ(const ColorTensor *colL, const ColorTensor *colR, 
                       Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c)
{

  const int BaseDelta[5][3] = {
    { 9, 3, 3 },
    { 3, 1, 1 },
    { 3, 1, 1 },
    { 3, 1, 1 },
    { 3, 1, 1 }
  };
  const Double factBaseDelta   = 2;


  const int BaseTl1cl1C [5][3] = {
    { 36, 12, 12 },
    { 12, 32, 11 },
    { 12, -4,  2 },
    { 12, 11, 32 },
    { 12,  2, -4 }
  };
  const Double factBaseTl1cl1C = 1./6;


  const int BaseTl1l1c[5][3] = {
    { 36, 12, 12 },
    { 12,  2, -4 },
    { 12, 11, 32 },
    { 12, -4,  2 },
    { 12, 32, 11 }
  };
  const Double factBaseTl1l1c  = 1./6;

  vector<Double> colGar, colQar;

  if(colL->GetOrg() == "G" && colR->GetOrg() != "G") {
    colGar = colL->GetGQVec();
    colQar = colR->GetQQVec();
  }
  else if(colL->GetOrg() != "G" && colR->GetOrg() == "G") {
    colQar = colL->GetQQVec();
    colGar = colR->GetGQVec();
  }

  Delta=0;
  Tl1cl1C=0;
  Tl1l1c =0;

  for(int i=0; i < 5; ++i)
    for(int j=0; j < 3; ++j) {
      Double colIJ = colGar[i]*colQar[j];
      Delta   += BaseDelta[i][j]   * colIJ;
      Tl1cl1C += BaseTl1cl1C[i][j] * colIJ;
      Tl1l1c  += BaseTl1l1c[i][j]  * colIJ;
    }

  Delta   *= factBaseDelta;
  Tl1cl1C *= factBaseTl1cl1C;
  Tl1l1c  *= factBaseTl1l1c;

}
