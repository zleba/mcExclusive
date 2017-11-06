#include "ColorTensor.h"

void ColorTensor::Print() const
{

  cout << Now<<" -> .. -> "<<Org<<" -> hadr process" << endl;

  if(Now == "G" && Org == "G" ) {
    cout << "A12="<<A12<<", A13="<<A13<<", A14="<<A14<<endl;
    cout << "B234="<<B234<<", B243="<<B243<<", B324="<<B324<<endl;
    cout << "B342="<<B342<<", B423="<<B423<<", B432="<<B432<<endl;
  }
  else if( (Now == "Q" || Now == "Qx") && (Org == "Q" || Org == "Qx") ) {
    cout << "Ki1j1="<<Ki1j1<<" ,Ki1j2="<<Ki1j2<<" , Ki1i2="<<Ki1i2<< endl;
  }
  else if( ( Now == "G" && (Org == "Q" || Org == "Qx") ) ||
      ( Org == "G" && (Now == "Q" || Now == "Qx") ) ) {
    cout << "Dij="<<Dij<< endl;
    cout << "T1 ="<<T1<<", T2 ="<< T2<<endl;
    cout << "T1c="<<T1<<", T2c="<< T2<<endl;
  }

}


vector<Double> ColorTensor::GetGGVec() const
{
  vector<Double> colArr(9);

  colArr[0] = A12;
  colArr[1] = A13;
  colArr[2] = A14; 

  colArr[3] = B234;
  colArr[4] = B243;
  colArr[5] = B324;

  colArr[6] = B342;
  colArr[7] = B423;
  colArr[8] = B432;

  return colArr;
}

vector<Double> ColorTensor::GetGQVec() const
{

  vector<Double> colArr(5);

  colArr[0] = Dij;

  colArr[1] = T1;
  colArr[2] = T2;
  colArr[3] = T1c;
  colArr[4] = T2c;

  return colArr;

}

vector<Double> ColorTensor::GetQQVec() const
{
  vector<Double> colArr(3);

  colArr[0] = Ki1j1;
  colArr[1] = Ki1i2;
  colArr[2] = Ki1j2;

  return colArr;

}



void   ColorTensor::Reset()
{
  A12=0; A13=0; A14=0;
  B234=0; B243=0; B324=0;
  B342=0; B423=0; B432=0;

  Ki1j1=0; Ki1j2=0; Ki1i2=0;

  Dij=0;
  T1=0;  T2=0;
  T1c=0; T2c=0;

}

void  ColorTensor::Flip(const ColorTensor &t)
{
  A12=t.A12; A13=t.A13; A14=t.A14;
  B234=t.B234; B243=t.B243; B324=t.B324;
  B342=t.B342; B423=t.B423; B432=t.B432;

  Ki1j1=t.Ki1j1; Ki1j2=t.Ki1j2; Ki1i2=t.Ki1i2;

  Dij=t.Dij;
  T1 =t.T1;  T2 =t.T2;
  T1c=t.T1c; T2c=t.T2c;
}

void ColorTensor::AddG()
{
  ColorTensor t;
  if( Now == "Q"  && (Org == "Q" || Org == "Qx") ) { //CORRECTED
    //Dij * T(i3,j3)*d_(i1,j1)
    t.Dij = Ki1j1;
    t.T1 = Ki1i2;
    t.T1c= Ki1j2;

  }
  else if( Now == "Qx" && (Org == "Q" || Org == "Qx") ) { //CORRECTED
    t.Dij = Ki1j1;
    t.T2  = Ki1j2;
    t.T2c = Ki1i2;

  }
  else if( Now == "G" && (Org == "Q" || Org == "Qx") ) { //CORRECTED

    t.Dij= 1./2*T2c + 1./2*T1c + 1./2*T2 + 1./2*T1 + 3.*Dij;
    t.T1 = 3./2*T1;
    t.T1c= 3./2*T1c;
    t.T2 = 3./2*T2;
    t.T2c= 3./2*T2c;

  }
  else if( Now=="Q"  && Org == "G") { //CORRECTED

    t.A13  = Dij;

    t.B243 = T2c;
    t.B324 = T1c;

  }
  else if( Now=="Qx" && Org == "G") { //CORRECTED

    t.A13 = Dij;

    t.B423 = T2;
    t.B342 = T1;

  }
  else if( Now=="G" && Org == "G") { //OK

    t.A12  = - 1./2*B432 - 1./2*B234;

    t.A13  = 1./2*B423 + 1./2*B342 + 1./2*B324 + 1./2*B243 + 3*A13;

    t.A14  = -1./2*B432 - 1./2*B234;


    t.B234 = - 1./2* A14 - 1./2*A12;

    t.B243 = 3./2*B243 + 1./2*A12;

    t.B324 = 3./2*B324 + 1./2*A14;


    t.B342 = 3./2*B342 + 1./2*A12;

    t.B423 = + 3./2* B423 + 1./2*A14;

    t.B432 = -1./2*A14 - 1./2*A12;

  }
  else {
    cout << "Something Wrong in AddG " << endl;
    cout << "Now : "<< Now << ", Org : " << Org << endl;
    exit(1);
  }

  Now = "G";
  Flip(t);

}



void ColorTensor::AddQ()
{
  ColorTensor t;

  if(Now == "Q" && (Org == "Q" || Org == "Qx") ) { //CORRECTED
    //d_(i1,j1)*d_(i3,j3)
    t.Ki1j1 = 1./2*Ki1i2 + 1./2*Ki1j2 + 4./3*Ki1j1;

    //d_(i1,i2)*d_(j1,j2)
    t.Ki1i2 = - 1./6*Ki1i2;

    //d_(i1,j2)*d_(i2,j1)
    t.Ki1j2 = - 1./6*Ki1j2;
  }
  else if(Now == "Q" && Org == "G") { //CORRECTED

    t.Dij = 1./2*T2c + 1./2*T1c + 4./3*Dij;

    t.T1c = - 1./6*T1c;
    t.T2c = - 1./6*T2c;

  }
  else if(Now == "G" && (Org == "Q" || Org == "Qx") ) { //CORRECTED

    t.Ki1j1 =  5./18*T2c + 1./36*T1c + 5./18*T2 + 1./36*T1 + 2./3*Dij;;
    t.Ki1i2 = - 1./6*T2 + 7./12*T1;
    t.Ki1j2 = - 1./6*T2c + 7./12*T1c;

  }
  else if(Now == "G" && Org == "G") { //CORRECTED

    //T(i1,j1)*d_(i3,j3)
    t.Dij=  1./36*B432 + 5./18*B423 + 5./18*B342 + 1./36*B324 + 1./36*B243 + 1./36*B234 + 
          2./3*A13;
    //L(i1,j3,x2)*L(j1,x2,i3)
    t.T1c = - 1./12*B432 - 1./6*B342 + 7./12*B324 - 1./12*B234 + 1./4*A14;

    //L(j1,j3,x1)*L(i1,x1,i3)
    t.T2c = - 1./12*B432 - 1./6*B423 + 7./12*B243 - 1./12*B234 + 1./4*A12;

  }
  else {
    cout << "Something Wrong in AddQ " << endl;
    cout << "Now : "<< Now << ", Org : " << Org << endl;
    exit(1);
  }

  Now = "Q";
  Flip(t);

}

void ColorTensor::AddQx()
{
  ColorTensor t;

  if(Now == "Qx" && (Org == "Q" || Org == "Qx") ) { //CORRECTED
    //d_(i1,j1)*d_(i3,j3)
    t.Ki1j1 = 1./2*Ki1i2 + 1./2*Ki1j2 + 4./3*Ki1j1;

    //d_(i1,i2)*d_(j1,j2)
    t.Ki1i2 = - 1./6*Ki1i2;

    //d_(i1,j2)*d_(i2,j1)
    t.Ki1j2 = - 1./6*Ki1j2;
  }
  else if(Now == "Qx" && Org == "G") { //CORRECTED

    t.Dij = 1./2*T2 + 1./2*T1 + 4./3*Dij;

    t.T1 = - 1./6*T1;
    t.T2 = - 1./6*T2;
  } 
  else if(Now == "G" && (Org == "Q" || Org == "Qx") ) { //CORRECTED

    t.Ki1j1= 1./36*T2c + 5./18*T1c + 1./36*T2 + 5./18*T1 + 2./3*Dij;
    t.Ki1i2= 7./12*T2c - 1./6*T1c;
    t.Ki1j2= 7./12*T2 - 1./6*T1;
  }
  else if(Now == "G" && Org == "G") { //CORRECTED

    //T(i1,j1)*d_(i3,j3)
    t.Dij= 1./36*B432 + 1./36*B423 + 1./36*B342 + 5./18*B324 + 5./18*B243 + 1./36*B234 + 
          2./3*A13;

    //L(i1,i3,x3)*L(j1,x3,j3)
    t.T1 = - 1./12*B432 + 7./12*B342 - 1./6*B324 - 1./12*B234 + 1./4*A12;

    //L(j1,i3,x3)*L(i1,x3,j3)
    t.T2 = - 1./12*B432 + 7./12*B423 - 1./6*B243 - 1./12*B234 + 1./4*A14;

  }
  else {
    cout << "Something Wrong in AddQx " << endl;
    cout << "Now : "<< Now << ", Org : " << Org << endl;
    exit(1);
  }

  Now = "Qx";
  Flip(t);

}

