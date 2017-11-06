#include "Sigma2Process.h"


static double GetMomentum(double s, double m1, double m2)
{
	double temp = s*s + m1*m1*m1*m1 + m2*m2*m2*m2 - 2*s*(m1*m1+m2*m2) - 2*m1*m1*m2*m2;

	return sqrt(temp) / ( 2.*sqrt(s) );

}




Sigma2Process::Sigma2Process()
: Nsym(1), i_(complex<Double>(0.,1.)),
alphaS(0.118), alphaEm(1/132.507), g2S(4*M_PI*alphaS), e2(4*M_PI*alphaEm)  
{
	//set pointers to zero
	colL = colR = 0;
	spinL = spinR = 0;

}


void Sigma2Process::SetScaleAlphas(Double _scale, Double _alphaS, Double _alphaEm)
{
  scale = _scale;

  alphaS  = _alphaS;
  alphaEm = _alphaEm;

  g2S = 4*M_PI * alphaS;
  e2  = 4*M_PI * alphaEm;
}

void Sigma2Process::SetIds(int id1, int id2)
{
    idLeft = id1; idRight = id2;       

		if(colL) {
			if(1 <= idLeft && idLeft <= 6)
				colL->SetStartingParton("Q");
			else if(1 <=-idLeft &&-idLeft <= 6)
				colL->SetStartingParton("Qx");
			else if(idLeft == 21)
				colL->SetStartingParton("G");
			else {
				cout << "Something wrong with input parton type" << endl;
				exit(1);
			}
		}

		if(colR) {
			if(1 <= idRight && idRight <= 6)
				colR->SetStartingParton("Q");
			else if(1 <=-idRight &&-idRight <= 6)
				colR->SetStartingParton("Qx");
			else if(idRight == 21)
				colR->SetStartingParton("G");
			else {
				cout << "Something wrong with input parton type" << endl;
				exit(1);
			}
		}

    if(spinL) spinL->Reset();
    if(spinR) spinR->Reset();

}







Double Sigma2Process::GetCrossSectionBase()
{
  Amplitudes();

  //Double mat[2][2] = { {16./3, -2./3},
                       //{-2./3, 16./3} };

  Double mat[1][1] = { 3 };
	
  complex<Double> res = 0;
  for(int i=0; i < 4; ++i)
    for(int j=0; j < 4; ++j) {
      for(int k=0; k < nJamp; ++k)
      for(int l=0; l < nJamp; ++l)
        res += mat[k][l] * JAMP[i][j][k]* conj(JAMP[i][j][l]);
    }
  
  return res.real();
}

void Sigma2Process::AddEmissionGeneral(SpinTensor *spin, ColorTensor *col, 
                                       int idIn, int idOut, Double z, Double phiEm)
{

  spin->AddEmission(idIn, idOut, z, phiEm);

  if( idIn == 21 )
    col->AddG();
  else if(  1 <= idIn && idIn <=  5 )
    col->AddQ();
  else if( -5 <= idIn && idIn <= -1 )
    col->AddQx();
  else {
    cout << "Unknown PDG id of parton shower particle! " <<idIn << endl;
    exit(1);
  }

}


void Sigma2Process::AddEmissionL(int idIn, int idOut, Double z, Double phiEm)
{
  AddEmissionGeneral(spinL, colL, idIn, idOut, z, phiEm);
}

void Sigma2Process::AddEmissionR(int idIn, int idOut, Double z, Double phiEm)
{
  AddEmissionGeneral(spinR, colR, idIn, idOut, z, phiEm);
}



void Sigma2Process::SetEThetaPhi(Double E, Double Theta, Double Phi)
{
  Theta *= M_PI/180; //to rad
  Phi *= M_PI/180;//to rad


  s = E*E;

  double pIn  = GetMomentum(s, mP1, mP2);
  double pOut = GetMomentum(s, mP3, mP4);

	double E1 = sqrt(pIn*pIn+mP1*mP1);
	double E3 = sqrt(pOut*pOut+mP3*mP3);


  t = mP1*mP1 + mP3*mP3 - 2*( E1*E3 - pIn*pOut*cos(Theta) );

	//Double m;
	//mP3 = mP4 = m = _mQ;
  //Double p = sqrt(E*E - m*m);
  //t = m*m - 2*(E*E - E*p*cos(Theta));

  //t = - s/2 * (1-cos(Theta) );

  u = mP1*mP1+mP2*mP2+mP3*mP3+mP4*mP4 - s - t;

  phi = Phi;

}


void Sigma2Process::SetSTUPhi(Double _s, Double _t, Double _u,  Double Phi)
{
  s   = _s;
  t   = _t;
  u   = _u;
  phi = Phi;

  //mP3 = mP4 = sqrt( (s + t + u)/2 );

}



void Sigma2Process::PrintInfo() const
{

    cout << "name : "<<  name   << endl;
    cout << "idLeft :"<< idLeft <<", idRight : "<<  idRight << endl;
    cout << "idPythia: "<< idPythia << endl;
    cout << "nJamp : " <<  nJamp << endl;
    cout << "nPol  : " <<  nPol  << endl;

}

void Sigma2Process::PrintAmplitudes() const
{
  string names[4];
  names[mm] = "- -";
  names[pm] = "+ -";
  names[mp] = "- +";
  names[pp] = "+ +";


  cout << "In Program" << endl;
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      cout << names[i] <<" -> "<<names[j]<<" : ";
      for(int k=0; k < nJamp; ++k)
        cout << JAMP[i][j][k] <<" ";
      cout << endl;
    }
    cout << endl;
  }

  cout << endl;
  cout << "In MadGraph" << endl;
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      cout << names[i] <<" -> "<<names[j]<<" : ";
      for(int k=0; k < nJamp; ++k)
        cout << JAMPmg[i][j][k] <<" ";
      cout << endl;
    }
    cout << endl;
  }

  cout << "Compare both" << endl;
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      cout << names[i] <<" -> "<<names[j]<<" : ";
      for(int k=0; k < nJamp; ++k) {
        if( (abs(JAMPmg[i][j][k]-JAMP[i][j][k]) < 1e-11 && abs(JAMP[i][j][k]) < 1e-11) ||
            abs(JAMPmg[i][j][k]/JAMP[i][j][k] - static_cast<Double>(1) ) < 1e-11 )
          cout << "OK ";
        else
          cout << "!! ";
      }
      cout << endl;
    }
    cout << endl;
  }

  cout << "Ratio" << endl;
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      cout << names[i] <<" -> "<<names[j]<<" : ";
      for(int k=0; k < nJamp; ++k) {
        cout << JAMPmg[i][j][k]/JAMP[i][j][k] <<" : r-1 " <<  JAMPmg[i][j][k]/JAMP[i][j][k]-1.<< endl;
        //if( (abs(JAMPmg[i][j][k]-JAMP[i][j][k]) < 1e-11 && abs(JAMP[i][j][k]) < 1e-11) ||
            //abs(JAMPmg[i][j][k]/JAMP[i][j][k] - static_cast<Double>(1) ) < 1e-11 )
          //cout << "OK ";
        //else
          //cout << "!! ";
      }
      cout << endl;
    }
    cout << endl;
  }


}




#ifdef MADGRAPH
void Sigma2Process::MadGraphResult()
{
  const int NEXTERNAL = 4;

  double P[4*NEXTERNAL];
  int NHEL[NEXTERNAL];
  complex<double> *JAMPtemp = new complex<double>[nJamp];

  int Hel[4][2] ={ { 1, 1 },
                   { 1,-1 },
                   {-1,+1 },
                   {-1,-1 } };

  {

    int i;

	double pIn  = GetMomentum(s, mP1, mP2);
	double pOut = GetMomentum(s, mP3, mP4);

	double E1 = sqrt(pIn*pIn+mP1*mP1);
	double E2 = sqrt(pIn*pIn+mP2*mP2);
	double E3 = sqrt(pOut*pOut+mP3*mP3);
	double E4 = sqrt(pOut*pOut+mP4*mP4);

    //double E = sqrt(s)/2;
    //double th = acos( 1+ 2*t/s);

    //double m2 = mP3*mP3;

	double th = acos( ( E1*E3 - (-t+mP1*mP1+mP3*mP3)/2. )/ pIn/pOut );

    //double pOut = sqrt(E*E-m2);

    //double th = acos( (2*E*E - m2 + t)/(2*E*p) );

    //first IN
    i = 0;
    P[i*4+0] = E1;
    P[i*4+1] = 0;
    P[i*4+2] = 0;
    P[i*4+3] = pIn;

    //second IN
    i = 1;
    P[i*4+0] = E2;
    P[i*4+1] = 0;
    P[i*4+2] = 0;
    P[i*4+3] =-pIn;

    //first OUT
    i = 2;
    P[i*4+0] = E3;
    P[i*4+1] = pOut*sin(th)*cos(phi);
    P[i*4+2] = pOut*sin(th)*sin(phi);
    P[i*4+3] = pOut*cos(th);

    //second OUT
    i = 3;
    P[i*4+0] = E4;
    P[i*4+1] =-pOut*sin(th)*cos(phi);
    P[i*4+2] =-pOut*sin(th)*sin(phi);
    P[i*4+3] =-pOut*cos(th);

  }



  for(int i=0; i < 4; ++i)
    for(int j=0; j < 4; ++j) {
      NHEL[0] = Hel[i][0];
      NHEL[1] = Hel[i][1];
      NHEL[2] = Hel[j][0];
      NHEL[3] = Hel[j][1];

      MgAmplitudes(P, NHEL, JAMPtemp);

      for(int k=0; k < nJamp; ++k)
        JAMPmg[i][j][k] = JAMPtemp[k] / sqrt(Nsym);

    }
}
#endif

void Sigma2Process::PrintColMat() const
{
  for(int i = 0; i < nJamp; ++i) {
    for(int j = 0; j < nJamp; ++j)
      cout << colMatrix[i][j] << " ";
    cout << endl;
  }


}

void Sigma2Process::Init()
{

  colMatrix.resize(nJamp);
  for(int i = 0; i < nJamp; ++i)
    colMatrix[i].assign(nJamp, 0);

  NonZero.resize(4);
  for(int i = 0; i < 4; ++i)
    NonZero[i].assign(4,false);

  JAMP.resize(4);
  JAMP[0].resize(4);
  JAMP[1].resize(4);
  JAMP[2].resize(4);
  JAMP[3].resize(4);

  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
      JAMP[i][j].assign(nJamp,0);

  JAMPmg = JAMP;

  isColorSinglet=false;

}


void  Sigma2Process::DataQQ(Double &Dl1l1c, Double &Dl1r1c, Double &Dl1r1)
{
  GeneralColorData::DataQQ(isColorSinglet, colL, colR, Dl1l1c, Dl1r1c, Dl1r1);
}

void  Sigma2Process::DataGQ(Double &Delta, Double  &Tl1cl1C, Double &Tl1l1c)
{
  GeneralColorData::DataGQ(isColorSinglet, colL, colR, Delta, Tl1cl1C, Tl1l1c);
}

void  Sigma2Process::DataGG( Double &A12, Double &A13, Double &A14,
                             Double &B234, Double &B243, Double &B342, Double &B432)
{
  GeneralColorData::DataGG(isColorSinglet, colL, colR,
                          A12, A13, A14, B234, B243, B342, B432);
}



Double Sigma2Process::CrossSection(bool _isColorSinglet, bool isSpinSinglet)
{
  isColorSinglet = _isColorSinglet;

  //Initiale Color Matrix
  ColorMatrix();

	//PrintColMat();

  //Initiale JAMP
  Amplitudes();


  complex<Double> res=0;


  for(int pol = 0;  pol < nPol; ++pol) {
    for(int colID1 = 0; colID1 < nJamp; ++colID1)
    for(int colID2 = 0; colID2 < nJamp; ++colID2) {
      res += colMatrix[colID1][colID2] * PartialSigma(pol, colID1, colID2, isSpinSinglet );
			//cout << "partialSigma "<< PartialSigma(pol, colID1, colID2, isSpinSinglet ) <<" "<< res<<" "<<isSpinSinglet << endl;
		}
  }

  //cout <<"ResOrg " <<res<< endl;

  /*
  for(int pol = 0;  pol < nPol; ++pol) 
    for(int polOut = 0;  polOut < nPol; ++polOut) {
      for(int colID1 = 0; colID1 < nJamp; ++colID1)
      for(int colID2 = 0; colID2 < nJamp; ++colID2)
        res += colMatrix[colID1][colID2] * JAMP[pol][polOut][colID1] * conj(JAMP[pol][polOut][colID2]);
    }
  */

  return res.real();// * 1./64 * 1./4 * 1./(2*s);
}


void GetIndexes(int a, int &l, int &r)
{
  switch(a) {
    case 0:
      l=0;r=0; break;
    case 1:
      l=1;r=1; break;
    case 2:
      l=0;r=1; break;
    case 3:
      l=1;r=0; break;
  }
	
}



complex<Double> Sigma2Process::MM(complex<Double> arr[][4], int a, int b) const
{
  int l,r, lC,rC;
	GetIndexes(a, l, r);
	GetIndexes(b, lC,rC);


  //return arr[2*lC+l][2*rC+r];
  return arr[2*lC+l][2*r+rC];
}


Double Sigma2Process::CrossSectionFull(bool _isColorSinglet, complex<Double> lum[3],  complex<Double>  lumC[3]  )
{
  isColorSinglet = _isColorSinglet;

  //Initiale Color Matrix
  ColorMatrix();

  //Initiale JAMP
  Amplitudes();


  complex<Double> arr[4][4];

  complex<Double> res=0;

  complex<Double> lumArr[4] = { lum[0] +lum[1] , lum[0] -lum[1] , lum[2] , conj(lum[2])  };
  complex<Double> lumArrC[4]= { lumC[0]+lumC[1], lumC[0]-lumC[1], lumC[2], conj(lumC[2]) };

  //complex<Double> lumArr[4] = { 1 , 2 , complex<Double>(3,3) ,complex<Double>(3,-3)  };
  //complex<Double> lumArrC[4]= { 3 , 2, complex<Double>(1,-1), complex<Double>(1, 1) };


  for(int pol = 0;  pol < nPol; ++pol) 
    for(int colID1 = 0; colID1 < nJamp; ++colID1)
    for(int colID2 = 0; colID2 < nJamp; ++colID2) {

      GeneralAmplitudes(arr, pol, colID1, colID2);
       

      for(int k=0; k < 4; ++k)
      for(int l=0; l < 4; ++l) {
        res += lumArr[k]*lumArrC[l] *colMatrix[colID1][colID2]* MM(arr,k,l);
      }


    }




  //cout <<"ResNew "<< res  << endl;


  return res.real();
}

Double Sigma2Process::CrossSectionFullFast(bool _isColorSinglet, complex<Double> lum[3],  complex<Double>  lumC[3]  )
{
  isColorSinglet = _isColorSinglet;

  //Initiale Color Matrix
  ColorMatrix();

  //Initiale JAMP
  Amplitudes();


  complex<Double> lumArr[4] = { lum[0] +lum[1] , lum[0] -lum[1] , lum[2] , conj(lum[2])  };
  complex<Double> lumArrC[4]= { lumC[0]+lumC[1], lumC[0]-lumC[1], lumC[2], conj(lumC[2]) };


  complex<Double> SpinAmp[4][4];
  complex<Double> SplittingAmp[4][4];

  GeneralSpinAmplitudes(SpinAmp);

  complex<Double> res = 0;
  complex<Double> sing = 0;

	//cout << "CrossSec ";
  for(int k=0; k < 4; ++k)
  for(int l=0; l < 4; ++l) {
    GeneralSplitting(k, l, SplittingAmp);


    complex<Double> temp=0;
    for(int i=0; i<4; ++i)
    for(int j=0; j<4; ++j)
      temp += SpinAmp[i][j] * SplittingAmp[i][j];

    if(k <= 1 && l <= 1)
      sing += temp;


		//cout << temp <<" ";
    res += lumArr[k]*lumArrC[l] * temp;

  }
	//cout << endl;
  //cout << "Singlet " << sing << endl;

  return res.real();
}







complex<Double> Sigma2Process::PartialSigma(int finPol, int colID1, int colID2, bool isSinglet ) const
{
  complex<Double> res = 0;

  complex<Double> Amp1[2][2], Amp2[2][2];

  Amp1[0][0] = JAMP[pp][finPol][colID1];
  Amp1[0][1] = JAMP[pm][finPol][colID1];
  Amp1[1][0] = JAMP[mp][finPol][colID1];
  Amp1[1][1] = JAMP[mm][finPol][colID1];
  
  Amp2[0][0] = JAMP[pp][finPol][colID2];
  Amp2[0][1] = JAMP[pm][finPol][colID2];
  Amp2[1][0] = JAMP[mp][finPol][colID2];
  Amp2[1][1] = JAMP[mm][finPol][colID2];

  const int (*Indx)[2];
  const static int IndxSinglet[2][2]   = { {1,2}, {2,1} };
  const static int IndxInclusive[2][2] = { {3,0}, {0,3} };

  if(isSinglet)
    Indx = IndxSinglet;
  else
    Indx = IndxInclusive;


  for(int i1=0; i1 < 2; ++i1)
  for(int i2=0; i2 < 2; ++i2)
  for(int i3=0; i3 < 2; ++i3)
  for(int i4=0; i4 < 2; ++i4) {

    res +=       spinL->Element(0,i1,i2) * Amp1[i2][i3] *
           conj( spinR->Element(0,i3,i4) * Amp2[i1][i4] );


    res +=       spinL->Element(3,i1,i2) * Amp1[i2][i3] *
           conj( spinR->Element(3,i3,i4) * Amp2[i1][i4] );


    //Interference terms


    res +=       spinL->Element(Indx[0][0],i1,i2) * Amp1[i2][i3] *
           conj( spinR->Element(Indx[0][1],i3,i4) * Amp2[i1][i4] );

    res +=       spinL->Element(Indx[1][0],i1,i2) * Amp1[i2][i3] *
           conj( spinR->Element(Indx[1][1],i3,i4) * Amp2[i1][i4] );


  }
  //cout << "r00Old " << finPol << " "<< colID1 <<" "<< colID2 <<" "<< r00 << endl;

  return res;

}


void Sigma2Process::GeneralAmplitudes(complex<Double> arr[][4], int finPol, int colID1, int colID2) const
{
  complex<Double> tempL[4], tempR[4];

  complex<Double> Amp1[2][2], Amp2[2][2];

  Amp1[0][0] = JAMP[pp][finPol][colID1];
  Amp1[0][1] = JAMP[pm][finPol][colID1];
  Amp1[1][0] = JAMP[mp][finPol][colID1];
  Amp1[1][1] = JAMP[mm][finPol][colID1];
  
  Amp2[0][0] = JAMP[pp][finPol][colID2];
  Amp2[0][1] = JAMP[pm][finPol][colID2];
  Amp2[1][0] = JAMP[mp][finPol][colID2];
  Amp2[1][1] = JAMP[mm][finPol][colID2];

  for(int k=0; k < 4; ++k)
  for(int l=0; l < 4; ++l)
    arr[k][l] = 0;


  for(int i1=0; i1 < 2; ++i1)
  for(int i3=0; i3 < 2; ++i3) {

    for(int i = 0; i < 4; ++i) {
      tempL[i] =       spinL->Element(i,i1,0) * Amp1[ 0][i3] +  spinL->Element(i,i1,1) * Amp1[ 1][i3] ;
      tempR[i] = conj( spinR->Element(i,i3,0) * Amp2[i1][ 0] +  spinR->Element(i,i3,1) * Amp2[i1][ 1] );
    }
    for(int k=0; k < 4; ++k)
    for(int l=0; l < 4; ++l)
      arr[k][l] += tempL[k]*tempR[l];
  }

  //cout << "r00New " << finPol << " "<< colID1 <<" "<< colID2 <<" "<< arr[0][0] << endl;


}



void Sigma2Process::GeneralSplitting(int a, int b, complex<Double> SplittingAmp[][4]) const
{
  int l,r, lC,rC;

	GetIndexes(a, l,  r);
	GetIndexes(b, lC, rC);


	for(int il1=0; il1<2; ++il1)
	for(int jl1=0; jl1<2; ++jl1)
	for(int ir1=0; ir1<2; ++ir1)
	for(int jr1=0; jr1<2; ++jr1) {
    SplittingAmp[2*il1+jl1][2*ir1+jr1] = spinL->Element(l,lC, il1,jl1) *
																		     spinR->Element(r,rC, ir1,jr1) ;
  }

}


void Sigma2Process::GeneralSpinAmplitudes(complex<Double> SpinAmp[][4]) const
{

	for(int il1=0; il1<2; ++il1)
	for(int jl1=0; jl1<2; ++jl1)
	for(int ir1=0; ir1<2; ++ir1)
	for(int jr1=0; jr1<2; ++jr1) {

		SpinAmp[2*il1+jl1][2*ir1+jr1] = 0;

    for(int pol = 0;  pol < nPol; ++pol) 
      if(NonZero[2*il1+ir1][pol] && NonZero[2*jl1+jr1][pol] ) 
        for(int colID2 = 0; colID2 < nJamp; ++colID2) {

          complex<Double> tempSum = 0;
          for(int colID1 = 0; colID1 < nJamp; ++colID1)
            tempSum += colMatrix[colID1][colID2] * JAMP[2*il1+ir1][pol][colID1]; 

          SpinAmp[2*il1+jl1][2*ir1+jr1]+=tempSum*conj(JAMP[2*jl1+jr1][pol][colID2]);

        }
  }

}


/*
void Sigma2Process::GeneralSpinAmplitudesOld(complex<Double> SpinAmp[][4]) const
{
  //SpinAmp[2*i1+i2][2*i3+i4]

  for(int i1=0; i1<2; ++i1)
  for(int i2=0; i2<2; ++i2)
  for(int i3=0; i3<2; ++i3)
  for(int i4=0; i4<2; ++i4) {

    SpinAmp[2*i1+i2][2*i3+i4] = 0;
    for(int pol = 0;  pol < nPol; ++pol) 
      if(NonZero[2*i2+i3][pol] && NonZero[2*i1+i4][pol] ) 
        for(int colID2 = 0; colID2 < nJamp; ++colID2) {

          complex<Double> tempSum = 0;
          for(int colID1 = 0; colID1 < nJamp; ++colID1)
            tempSum += colMatrix[colID1][colID2] * JAMP[2*i2+i3][pol][colID1]; 

          SpinAmp[2*i1+i2][2*i3+i4]+=tempSum*conj(JAMP[2*i1+i4][pol][colID2]);

        }
  }

//  res +=       spinL->Element(0,i1,i2) * Amp1[i2][i3] *
//         conj( spinR->Element(0,i3,i4) * Amp2[i1][i4] );

}



void Sigma2Process::GeneralSplittingOld(int a, int b, complex<Double> SplittingAmp[][4]) const
{
  int l,r, lC,rC;

	GetIndexes(a, l,  r);
	GetIndexes(b, lC, rC);


  //index of left and right splitting matrix
  int ind1 = 2*lC+l;
  int ind2 = 2*r+rC;


  for(int i1=0; i1<2; ++i1)
  for(int i2=0; i2<2; ++i2)
  for(int i3=0; i3<2; ++i3)
  for(int i4=0; i4<2; ++i4) {
    SplittingAmp[2*i1+i2][2*i3+i4] = spinL->Element(ind1,i1,i2) *
                               conj( spinR->Element(ind2,i3,i4) );
  }

}
*/
