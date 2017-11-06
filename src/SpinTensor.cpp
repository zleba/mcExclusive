#include "SpinTensor.h"

void SpinTensor::Reset()
{
  for(int i=0; i < 4; ++i)
    for(int j=0; j < 4; ++j) {
      P[i][j] = (i==j) ? 1 : 0;
    }
}

void SpinTensor::PrintPnow() const
{
  for(int i=0; i < 4; ++i) {
    for(int j=0; j < 4; ++j) 
      cout << Pnow[i][j] << " ";
    cout << endl;
  }

}


Double SpinTensor::ClassicalSplitting(int idIn, int idOut, Double z)
{

  if(idIn == 21 && idOut == 21) 
    return 6*(1/z + 1/(1-z) -2 +z*(1-z) );

  if(idIn == 21 && idOut != 21)
    return 1./2*(z*z + (1-z)*(1-z) );

  if(idIn != 21 && idOut == 21)
    return 4./3 * ( 1 + (1-z)*(1-z) )/z ;

  if(idIn != 21 && idOut != 21)
    return 4./3 * ( 1 + z*z)/(1-z) ;

  return 0.0;
}




void SpinTensor::CalculateMatrixGG(Double z, Double ph)
{
 //0  X(1,1) 
 //1  X(1,2) 
 //2  X(2,1) 
 //3  X(2,2) 

static const complex<Double> i_(0,1);

complex<Double> ExpP = exp( + Double(2)*i_*ph);
complex<Double> ExpM = conj(ExpP);


//Correct solution
//Step 11
Pnow[0][0] = - 1 + 1/z - z - z*z + 2/(1 - z);

Pnow[0][1] = -z*(1-z) *ExpP;

Pnow[0][2] = -z*(1-z)* ExpM;

Pnow[0][3] = - 3 + 1/z + 3*z - z*z;

//Step 12
Pnow[1][0] =  -(1-z)/z * ExpM;

Pnow[1][1] =  2.*z/(1-z);

Pnow[1][2] = 0;

Pnow[1][3] =  -(1-z)/z * ExpM;

//Step 21
Pnow[2][0] = -(1-z)/z * ExpP;

Pnow[2][1] = 0;

Pnow[2][2] =  2.*z/(1-z);

Pnow[2][3] =  -(1-z)/z * ExpP;


//Step 22
Pnow[3][0] = - 3 + 1/z + 3*z - z*z;

Pnow[3][1] = -z*(1-z) * ExpP;

Pnow[3][2] = -z*(1-z) * ExpM;

Pnow[3][3] = - 1 + 1/z - z - z*z + 2/(1 - z);

}

//From gluon to quark
void SpinTensor::CalculateMatrixGQ(Double z, Double ph)
{

   //0  X(1,1) 
   //1  X(1,2) 
   //2  X(2,1) 
   //3  X(2,2) 


  static const complex<Double> i_(0,1);

  complex<Double> ExpP = exp( + Double(2)*i_*ph);
  complex<Double> ExpM = conj(ExpP);



  //Aplitudes matrix

  Pnow[0][0] = z*z;
  Pnow[0][1] = -z*(1-z) * ExpP;
  Pnow[0][2] = -z*(1-z) * ExpM;
  Pnow[0][3] = 1 - 2*z + z*z;


  Pnow[1][0] = 0;
  Pnow[1][1] = 0;
  Pnow[1][2] = 0;
  Pnow[1][3] = 0;


  Pnow[2][0] = 0;
  Pnow[2][1] = 0;
  Pnow[2][2] = 0;
  Pnow[2][3] = 0;


  Pnow[3][0] = 1 - 2*z + z*z;
  Pnow[3][1] = -z*(1-z) * ExpP;
  Pnow[3][2] = -z*(1-z) * ExpM;
  Pnow[3][3] = z*z;

}


void SpinTensor::CalculateMatrixQQ(Double z, Double)
{

  Pnow[0][0] = - 1 - z + 2/(1 - z);
  Pnow[0][1] = 0;
  Pnow[0][2] = 0;
  Pnow[0][3] = 0;

  Pnow[1][0] = 0;
  Pnow[1][1] = - 2 + 2/(1 - z);
  Pnow[1][2] = 0;
  Pnow[1][3] = 0;

  Pnow[2][0] = 0;
  Pnow[2][1] = 0;
  Pnow[2][2] = - 2 + 2/(1 - z);
  Pnow[2][3] = 0;

  Pnow[3][0] = 0;
  Pnow[3][1] = 0;
  Pnow[3][2] = 0;
  Pnow[3][3] = - 1 - z + 2/(1 - z);


}


void SpinTensor::CalculateMatrixQG(Double z, Double ph)
{
  static const complex<Double> i_(0,1);
  complex<Double> ExpP = exp( + Double(2)*i_*ph);
  complex<Double> ExpM = conj(ExpP);

  Pnow[0][0] = 1/z;
  Pnow[0][1] = 0;
  Pnow[0][2] = 0;
  Pnow[0][3] = - 2 + 1/z + z;

  Pnow[1][0] = -(1-z)/z * ExpM;
  Pnow[1][1] = 0;
  Pnow[1][2] = 0;
  Pnow[1][3] = -(1-z)/z * ExpM;

  Pnow[2][0] = -(1-z)/z * ExpP;
  Pnow[2][1] = 0;
  Pnow[2][2] = 0;
  Pnow[2][3] = -(1-z)/z * ExpP;

  Pnow[3][0] = - 2 + 1/z + z;
  Pnow[3][1] = 0;
  Pnow[3][2] = 0;
  Pnow[3][3] = 1/z;


}



void SpinTensor::AddEmission(int idIn, int idOut, Double z, Double ph)
{

  complex<Double> Pold[4][4];


  if     (idIn == 21 && idOut == 21)
      CalculateMatrixGG(z, ph);
  else if(idIn == 21 && In(idOut,-6,6) )
      CalculateMatrixGQ(z, ph);
  else if(In(idIn,-6,6) && idOut == 21 )
      CalculateMatrixQG(z, ph);
  else if( In(idIn,-6,6) && idIn == idOut)
      CalculateMatrixQQ(z, ph);
  else {
    cout << "Something wrong with flavours :"<<endl;
    cout << "In->Out "<< idIn <<" -> "<<idOut << endl;
    exit(1);
  }



  for(int i=0; i < 4; ++i)
    for(int j=0; j < 4; ++j) {
      Pold[i][j] = P[i][j];
    }
      

  for(int i=0; i < 4; ++i)
    for(int j=0; j < 4; ++j) {
      P[i][j] = 0;
      for(int k=0; k < 4; ++k)
        P[i][j] += Pold[i][k] * Pnow[k][j];
    }

}


