#include "ColorData.h"
#include "ColorTensor.h"
#include <vector>
using std::vector;

void SingletGG::DataQQ( const ColorTensor *colL, const ColorTensor *colR, 
                        Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1)
//CORRECTED
{

  const int BaseDl1l1c[5][5] = {
    { 9, 3, 3, 3, 3, },
    { 3, 1, 1, 1, 1, },
    { 3, 1, 1, 1, 1, },
    { 3, 1, 1, 1, 1, },
    { 3, 1, 1, 1, 1, }

  };
  const Double factBaseDl1l1c = 2;


  const int BaseDl1r1c[5][5] = {
    { 36, 12, 12, 12, 12, },
    { 12, - 4, 32, 2, 11, },
    { 12, 32, - 4, 11, 2, },
    { 12, 2, 11, - 4, 32, },
    { 12, 11, 2, 32, - 4, }

  };
  const Double factBaseDl1r1c = 1./6;


  const int BaseDl1r1[5][5] = {
    { 36, 12, 12, 12, 12, },
    { 12, 2, 11, - 4, 32, },
    { 12, 11, 2, 32, - 4, },
    { 12, - 4, 32, 2, 11, },
    { 12, 32, - 4, 11, 2, }


  };
  const Double factBaseDl1r1  = 1./6;

  vector<Double> colLar, colRar;

  colLar = colL->GetGQVec(); 
  colRar = colR->GetGQVec(); 


  Dl1l1c = Dl1r1c = Dl1r1 = 0;


  for(int i=0; i < 5; ++i)
    for(int j=0; j < 5; ++j) {
      Double colIJ = colLar[i]*colRar[j];
      Dl1l1c += BaseDl1l1c[i][j] * colIJ;
      Dl1r1  += BaseDl1r1[i][j]  * colIJ;
      Dl1r1c += BaseDl1r1c[i][j] * colIJ;
    }


  Dl1l1c *= factBaseDl1l1c;
  Dl1r1  *= factBaseDl1r1;
  Dl1r1c *= factBaseDl1r1c;

}


void InclusiveGG::DataQQ( const ColorTensor *colL, const ColorTensor *colR, 
                    Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1)
{
//CORRECTED

  const int BaseDl1l1c[5][5] = {
  { 9, 3, 3, 3, 3, },
  { 3, 1, 1, 1, 1, },
  { 3, 1, 1, 1, 1, },
  { 3, 1, 1, 1, 1, },
  { 3, 1, 1, 1, 1, }

  };
  const Double factBaseDl1l1c = 16;

  const Double factBaseDl1r1c = 16/3.;
  const Double factBaseDl1r1  = 16/3.;

  vector<Double> colLar, colRar;

  colLar = colL->GetGQVec();
  colRar = colR->GetGQVec();


  Dl1l1c = 0;


  for(int i=0; i < 5; ++i)
    for(int j=0; j < 5; ++j) {
      Double colIJ = colLar[i]*colRar[j];
      Dl1l1c += BaseDl1l1c[i][j] * colIJ;
    }
  Dl1r1c = Dl1r1 = Dl1l1c;

  Dl1l1c *= factBaseDl1l1c;
  Dl1r1  *= factBaseDl1r1;
  Dl1r1c *= factBaseDl1r1c;

}


void SingletGG::DataGQ(const ColorTensor *colL, const ColorTensor *colR,
                       Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c)
{
//CORRECTED

  const int BaseDelta[9][5] = {
    { 9, 3, 3, 3, 3, },
    { 72, 24, 24, 24, 24, },
    { 9, 3, 3, 3, 3, },
    { - 3, - 1, - 1, - 1, - 1, },
    { 24, 8, 8, 8, 8, },
    { 24, 8, 8, 8, 8, },
    { 24, 8, 8, 8, 8, },
    { 24, 8, 8, 8, 8, },
    { - 3, - 1, - 1, - 1, - 1, }

  };
  const Double factBaseDelta = 1./6;


  const int BaseTl1cl1C[9][5] = {
    { 12, 32, - 4, 11, 2, },
    { 96, 32, 32, 32, 32, },
    { 12, - 4, 32, 2, 11, },
    { - 4, - 8, - 8, - 3, - 3, },
    { 32, 16, - 8, 12, 6, },
    { 32, - 8, 76, 6, 27, },
    { 32, 76, - 8, 27, 6, },
    { 32, - 8, 16, 6, 12, },
    { - 4, - 8, - 8, - 3, - 3, }

  };
  const Double factBaseTl1cl1C = 1./24;


  const int BaseTl1l1c[9][5] = {
    { 12, 2, 11, - 4, 32, },
    { 96, 32, 32, 32, 32, },
    { 12, 11, 2, 32, - 4, },
    { - 4, - 3, - 3, - 8, - 8, },
    { 32, 6, 27, - 8, 76, },
    { 32, 12, 6, 16, - 8, },
    { 32, 6, 12, - 8, 16, },
    { 32, 27, 6, 76, - 8, },
    { - 4, - 3, - 3, - 8, - 8, }

  };
  const Double factBaseTl1l1c  =  1./24;


  vector<Double> colGar, colQar;

  if(colL->GetOrg() == "G" && colR->GetOrg() != "G") {
    colGar = colL->GetGGVec();
    colQar = colR->GetGQVec();
  }
  else if(colL->GetOrg() != "G" && colR->GetOrg() == "G") {
    colQar = colL->GetGQVec();
    colGar = colR->GetGGVec();
  }



  Delta   = 0;
  Tl1cl1C = 0;
  Tl1l1c  = 0;

  for(int i=0; i < 9; ++i)
    for(int j=0; j < 5; ++j) {
      Double colIJ = colGar[i]*colQar[j];
      Delta   += BaseDelta[i][j]   * colIJ;
      Tl1cl1C += BaseTl1cl1C[i][j] * colIJ;
      Tl1l1c  += BaseTl1l1c[i][j]  * colIJ;
    }

  Delta   *= factBaseDelta;
  Tl1cl1C *= factBaseTl1cl1C;
  Tl1l1c  *= factBaseTl1l1c;

}




void InclusiveGG::DataGQ(const ColorTensor *colL, const ColorTensor *colR,
                         Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c)
{
//CORRECTED

  const int Base[9][5] = {
  { 9, 3, 3, 3, 3, },
  { 72, 24, 24, 24, 24, },
  { 9, 3, 3, 3, 3, },
  { - 3, - 1, - 1, - 1, - 1, },
  { 24, 8, 8, 8, 8, },
  { 24, 8, 8, 8, 8, },
  { 24, 8, 8, 8, 8, },
  { 24, 8, 8, 8, 8, },
  { - 3, - 1, - 1, - 1, - 1, }


  };
  const Double factBaseDelta = 4./3;

  const Double factBaseTl1cl1C = 4./9;

  const Double factBaseTl1l1c  = 4./9;


  vector<Double> colGar, colQar;

  if(colL->GetOrg() == "G" && colR->GetOrg() != "G") {
    colGar = colL->GetGGVec();
    colQar = colR->GetGQVec();
  }
  else if(colL->GetOrg() != "G" && colR->GetOrg() == "G") {
    colQar = colL->GetGQVec();
    colGar = colR->GetGGVec();
  }



  Delta   = 0;

  for(int i=0; i < 9; ++i)
    for(int j=0; j < 5; ++j) {
      Delta   += Base[i][j] * colGar[i]*colQar[j] ;
    }

  Tl1cl1C = Tl1l1c = Delta;

  Delta   *= factBaseDelta;
  Tl1cl1C *= factBaseTl1cl1C;
  Tl1l1c  *= factBaseTl1l1c;
}



void SingletGG::DataGG( const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432)
{

const int BaseA12[9][9] = {
{ 24, 3, 3, 8, 8, - 1, 8, - 1, 8, },
{ 3, 24, 3, - 1, 8, 8, 8, 8, - 1, },
{ 3, 3, 24, 8, - 1, 8, - 1, 8, 8, },
{ 8, - 1, 8, 4, - 2, - 2, - 2, - 2, 19, },
{ 8, 8, - 1, - 2, 4, - 2, 19, - 2, - 2, },
{ - 1, 8, 8, - 2, - 2, 4, - 2, 19, - 2, },
{ 8, 8, - 1, - 2, 19, - 2, 4, - 2, - 2, },
{ - 1, 8, 8, - 2, - 2, 19, - 2, 4, - 2, },
{ 8, - 1, 8, 19, - 2, - 2, - 2, - 2, 4, } };
const Double factBaseA12 = 1./24;


const int BaseA13[9][9] = {
{ 9, 72, 9, - 3, 24, 24, 24, 24, - 3, },
{ 72, 576, 72, - 24, 192, 192, 192, 192, - 24, },
{ 9, 72, 9, - 3, 24, 24, 24, 24, - 3, },
{ - 3, - 24, - 3, 1, - 8, - 8, - 8, - 8, 1, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ - 3, - 24, - 3, 1, - 8, - 8, - 8, - 8, 1, } };
const Double factBaseA13 = 1./72;

const int BaseA14[9][9] = {
{ 3, 3, 24, 8, - 1, 8, - 1, 8, 8, },
{ 3, 24, 3, - 1, 8, 8, 8, 8, - 1, },
{ 24, 3, 3, 8, 8, - 1, 8, - 1, 8, },
{ 8, - 1, 8, 19, - 2, - 2, - 2, - 2, 4, },
{ - 1, 8, 8, - 2, - 2, 4, - 2, 19, - 2, },
{ 8, 8, - 1, - 2, 4, - 2, 19, - 2, - 2, },
{ - 1, 8, 8, - 2, - 2, 19, - 2, 4, - 2, },
{ 8, 8, - 1, - 2, 19, - 2, 4, - 2, - 2, },
{ 8, - 1, 8, 4, - 2, - 2, - 2, - 2, 19, } };
const Double factBaseA14 = 1./24;



const int BaseB234[9][9] = {
{ 72, - 9, 72, 171, - 18, - 18, - 18, - 18, 36, },
{ - 9, - 72, - 9, 3, - 24, - 24, - 24, - 24, 3, },
{ 72, - 9, 72, 36, - 18, - 18, - 18, - 18, 171, },
{ 36, 3, 171, 83, 11, 11, 11, 11, 83, },
{ - 18, - 24, - 18, 11, 2, 2, - 43, - 43, 11, },
{ - 18, - 24, - 18, 11, 2, 2, - 43, - 43, 11, },
{ - 18, - 24, - 18, 11, - 43, - 43, 2, 2, 11, },
{ - 18, - 24, - 18, 11, - 43, - 43, 2, 2, 11, },
{ 171, 3, 36, 83, 11, 11, 11, 11, 83, } };
const Double factBaseB234 = 1./216;

const int BaseB243[9][9] = {
{ 72, 72, - 9, - 18, 171, - 18, 36, - 18, - 18, },
{ 72, 576, 72, - 24, 192, 192, 192, 192, - 24, },
{ - 9, 72, 72, - 18, - 18, 36, - 18, 171, - 18, },
{ - 18, - 24, - 18, 11, - 43, 2, 2, - 43, 11, },
{ 36, 192, - 18, 2, 92, 56, 92, - 34, 2, },
{ - 18, 192, 171, - 43, - 34, 92, - 34, 407, - 43, },
{ 171, 192, - 18, - 43, 407, - 34, 92, - 34, - 43, },
{ - 18, 192, 36, 2, - 34, 92, 56, 92, 2, },
{ - 18, - 24, - 18, 11, - 43, 2, 2, - 43, 11, } };
const Double factBaseB243 = 1./216;

const int BaseB342[9][9] = {
{ 72, 72, - 9, - 18, 36, - 18, 171, - 18, - 18, },
{ 72, 576, 72, - 24, 192, 192, 192, 192, - 24, },
{ - 9, 72, 72, - 18, - 18, 171, - 18, 36, - 18, },
{ - 18, - 24, - 18, 11, 2, - 43, - 43, 2, 11, },
{ 171, 192, - 18, - 43, 92, - 34, 407, - 34, - 43, },
{ - 18, 192, 36, 2, 56, 92, - 34, 92, 2, },
{ 36, 192, - 18, 2, 92, - 34, 92, 56, 2, },
{ - 18, 192, 171, - 43, - 34, 407, - 34, 92, - 43, },
{ - 18, - 24, - 18, 11, 2, - 43, - 43, 2, 11, } };
const Double factBaseB342 = 1./216;

const int BaseB432[9][9] = {
{ 72, - 9, 72, 36, - 18, - 18, - 18, - 18, 171, },
{ - 9, - 72, - 9, 3, - 24, - 24, - 24, - 24, 3, },
{ 72, - 9, 72, 171, - 18, - 18, - 18, - 18, 36, },
{ 171, 3, 36, 83, 11, 11, 11, 11, 83, },
{ - 18, - 24, - 18, 11, 2, 2, - 43, - 43, 11, },
{ - 18, - 24, - 18, 11, 2, 2, - 43, - 43, 11, },
{ - 18, - 24, - 18, 11, - 43, - 43, 2, 2, 11, },
{ - 18, - 24, - 18, 11, - 43, - 43, 2, 2, 11, },
{ 36, 3, 171, 83, 11, 11, 11, 11, 83, } };
const Double factBaseB432 = 1./216;

vector<Double> colLar, colRar;

colLar = colL->GetGGVec(); //gluon -> gluon
colRar = colR->GetGGVec(); //gluon -> gluon


A12=0; A13=0; A14=0;
B234=0; B243=0;
B342=0; B432=0;

for(int i=0; i < 9; ++i)
  for(int j=0; j < 9; ++j) {
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

void InclusiveGG::DataGG( const ColorTensor *colL, const ColorTensor *colR,
                Double &A12, Double &A13, Double &A14,
                Double &B234, Double &B243, Double &B342, Double &B432)
{

const int BaseA12[9][9] = {
{ 9, 72, 9, - 3, 24, 24, 24, 24, - 3, },
{ 72, 576, 72, - 24, 192, 192, 192, 192, - 24, },
{ 9, 72, 9, - 3, 24, 24, 24, 24, - 3, },
{ - 3, - 24, - 3, 1, - 8, - 8, - 8, - 8, 1, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ 24, 192, 24, - 8, 64, 64, 64, 64, - 8, },
{ - 3, - 24, - 3, 1, - 8, - 8, - 8, - 8, 1, }
};
const Double factBaseA12 = 1./72;


const Double factBaseA13  = 1./9;
const Double factBaseA14  = 1./72;

const Double factBaseB234 =-1./216;
const Double factBaseB243 = 1./27;
const Double factBaseB342 = 1./27;
const Double factBaseB432 =-1./216;


vector<Double> colLar, colRar;

colLar = colL->GetGGVec(); //gluon -> gluon
colRar = colR->GetGGVec(); //gluon -> gluon



A12=0;

for(int i=0; i < 9; ++i)
  for(int j=0; j < 9; ++j) {
    A12  += BaseA12[i][j]  * colLar[i]*colRar[j];
  }

A13=A14=B234=B243=B342=B432 = A12;


A12 *= factBaseA12;
A13 *= factBaseA13;
A14 *= factBaseA14;
B234*= factBaseB234;
B243*= factBaseB243;
B342*= factBaseB342;
B432*= factBaseB432;


}
