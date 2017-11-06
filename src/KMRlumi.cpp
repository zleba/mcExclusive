#include "KMRlumi.h"


//KMR initialisation
KMRlumi::KMRlumi(Double sqrtS, Double qMin, Double alphaS, Double mc, Double mb, Pythia8::BeamParticle *_beamPtr )
{
  s= sqrtS*sqrtS;

  bSlope = 4;

  MinQt2 = qMin*qMin;   //Minimum for the Integration


  //pdfName = PDFname;
  //LHAPDF::initPDFSetByName(PDFname);
  //LHAPDF::initPDF(0);

  //Get masses from PDF information
  CharmMass  = mc;
  BottomMass = mb;
  //CharmMass  = 1.4;
  //BottomMass = 4.75;

  const Double MZ = 91.18800;
  const Double ASConst5     = 6. / (33. - 2*5);

  Lambda5QCD   =  MZ * exp( -M_PI*ASConst5/ alphaS ) ;

  //cout <<"Lambda "<<setprecision(10)<< Lambda5QCD << endl;
  //exit(0);

  //log(scale^2) where alphaS freeze
  LnFreeze2 = log(MinQt2);
  LnCharmMass2   = 2*log(CharmMass);
  LnBottomMass2  = 2*log(BottomMass);


  //MinQ2PDF = 1; //good to know
 //LnMinQ2PDF = log(MinQ2PDF);
  beamPtr = _beamPtr;
}



inline Double KMRlumi::Uniform(Double a, Double b)
{
  return (  a +  (b-a) * rand() / (RAND_MAX + 0.0)  );
  //return (  a +  (b-a) * MyRand() );

}


inline Double KMRlumi::gluon(int id, Double x, Double q2 ) const
{
  // Standard parton distributions.
    return beamPtr->xf(id, x, q2);
}

//PDG independent - PDG + diggPDG
pair<Double,Double> KMRlumi::gluonPairChic(int id, Double x, Double q) const
{
  
  if(x > 0.9999)
     return make_pair(0.0,0.0);


  const Double q0=1.5; //was 1.5

  const Double eps=1e-3;
  const Double eps1=1e-4;
  Double qsq = q*q;

  if(qsq  > q0) {
    Double fPlus  = gluon(id, x, qsq+eps1);
    Double fMinus = gluon(id, x, qsq-eps1);

    Double Df = qsq*(fPlus-fMinus)/(2*eps1);
    Double  f = (fPlus+fMinus)/2;

    return make_pair(f, Df);
  }
  else {
     Double Lq0 = log(q0);


     const Double qp=q0+eps1;
     const Double qm=q0-eps1;

     Double glu1p=  gluon(id, x, qp); 
     Double glu1pp= gluon(id, x, qp+eps);
     Double glu1m=  gluon(id, x, qm);
     Double glu1mp= gluon(id, x, qm+eps);

		 if( glu1p < 1e-6 || glu1m < 1e-6 )
		 	 return make_pair(0, 0);

     Double diffp=qp*(glu1pp-glu1p)/(eps*glu1p);
     Double diffm=qm*(glu1mp-glu1m)/(eps*glu1m);

     Double ddiff=q0*(diffp-diffm)/(2*eps1);


     //Double fPlus   = LHAPDF::xfx(x, sqrt(qsq+eps1), 0);
     //Double fMinus  = LHAPDF::xfx(x, sqrt(qsq-eps1), 0);
     //Double fCenter = LHAPDF::xfx(x, sqrt(qsq), 0);
     //Double diffp   = q0*(fPlus-fMinus)/(2*eps1*fCenter);
     //Double ddiff   = diffp - diffp*diffp + q0*q0*( fPlus  -2*fCenter + fMinus )/(eps1*eps1*fCenter);
     //Double glu1p = fCenter;


     //Double lamn=((diffp+2*Lq0)/(1+2*Lq0)-(ddiff+2* Lq0+4)/4.0/(1+Lq0));
     //lamn=lamn/((1+Lq0)/(1+2*Lq0)- (2+Lq0)/4.0/(1+Lq0)) ;


     Double lamn=(4*(diffp+2*Lq0)*(1+Lq0) - (ddiff+2*Lq0+4)*(1+2*Lq0)  )   /  (4*(1+Lq0)*(1+Lq0)- (2+Lq0)*(1+2*Lq0));

     
     //Double lamt1=(diffp+2*Lq0-lamn*(1+Lq0))  /  (1+2*Lq0);
     Double lamt2=(ddiff+2*Lq0+4-lamn*(2+Lq0)) /4.0/(1+Lq0);

     Double rat = qsq / q0;
     Double power=2+(lamn-2)*rat + lamt2* rat*rat;

     if(power < 0)
         return make_pair(glu1p, 0);
     else {
        Double a= glu1p/pow(q0,lamn+lamt2);
        //return qsq*( ((lamn-2)/q0+2*lamt2*qsq/q0/q0)*log(qsq) + power/qsq ) * a*pow(qsq,power);
         //dxg=dxg*glu

        Double f  = a*pow(qsq, power) ;
        Double Df = qsq*( ((lamn-2)/q0+2*lamt2*qsq/q0/q0)*log(qsq) + power/qsq ) * f;
				//cout << a << glu1p<<" "<<q0<<" "<<lamn<<" "<<lamt2<<endl;
				//cout << glu1p<<" "<<glu1pp<<" "<<glu1m<<" "<<glu1mp<<endl;
        return make_pair(f,Df);

     }

  }
    
}



//PDG dependent
Double KMRlumi::LeifGluonChicUnintegrated(int id, Double x, Double q, Double y, Double Int) const
{
  Double weight = 1;

	if(id == 21) { //gluon
		weight *= Int/ (y*y*sqrt(1-y) );
		weight *= 16*y/M_PI *sqrt(y*(1-y));
	}
	else {         //quark or anti-quark
		weight *= Int/ ( y*y / sqrt(1-y) );
		weight *= 8*y/M_PI/x *sqrt(y/(1-y));
	}

  Double split = SplittingAlphaS(id, 2*log(q) );

  Double g,Dg;
  pair<Double,Double> res;
  res = gluonPairChic(id, x/4/y, q);
  g  = res.first;
  Dg = res.second;
  
  return  weight*(Dg + 0.5*g * split );

}


void KMRlumi::GetYintGluon(Double x, Double &y, Double &Int)
{
  Int = (4.-x)*sqrt(4.-x) * ( 3.*x*(5.*x+16.) + 128.) / 6720.;

  static const Double Max = 16/25./sqrt(5);

  Double wg, temp;
  do {
    y = Uniform(x/4,1);
    temp = Uniform(0,Max);
    wg = y*y * sqrt(1-y);

  } while(temp > wg); 

}

void KMRlumi::GetYintQuark(Double x, Double &y, Double &Int)
{
  Int = sqrt(4.-x) * ( x*(3.*x+16.) + 128.) / 240.;

  Double wg, temp;
  do {
    y = Uniform(0,1);
		y = 1- (1-x/4)*y*y;
    temp = Uniform(0, 1./sqrt(1-y)  );
    wg = y*y / sqrt(1-y);

  } while(temp > wg); 

}








vector< complex<Double> > KMRlumi::LumSqrtRandom(int id1, Double M, Double y, Double Mju, Double MjuLoop, const Vec2 &p1, const Vec2 &p2) 
{
	int id2 = (id1 == 21) ? id1 : -id1;
  vector< complex<Double> > lum(3,0);
  const Double pt2max = 2;
  if(p1.norm2() > pt2max|| p2.norm2() > pt2max) return lum;


  Double weight = M_PI/4;     //magic factor (from t and/or q integration)
  weight *= 1/(bSlope); // ??



  //double M =  1/(MInv); //exp(0.5*LnM2);
  Double sSqrt = sqrt(s);

  x1 = exp(y)   * M/sSqrt; 
  x2 = exp(-y)  * M/sSqrt; 


  if(x1 > 1 || x2 > 1 || x1 < 0 || x2 < 0) {
    cout << "x1,x2 out of range " << x1<<" "<<x2 << endl;
    cout << "M, y, sSqrt " << M<<" "<<y<<" "<<sSqrt << endl;
    exit(1);
  }

  
  //Setting of hadr scale

  Double Qt2InvA;


  //Importance sampling
  const Double q0i = 1/2.;
  //const Double Norm = q0i + 2 * q0i*sqrt(q0i)*(1/sqrt(q0i) - sqrt(MinQt2) );
  const Double Norm = q0i +  q0i*q0i*(1/q0i - MinQt2 );
  double r1 = Norm * Uniform(0,1);

  if(r1< q0i) {
    Qt2InvA = r1;
    weight *=1;
  }
  else {
    //Qt2InvA = 1/sqrt(q0i) - (r1-q0i)/(2*q0i*sqrt(q0i));
    //Qt2InvA = 1/(Qt2InvA*Qt2InvA);
    //weight *= pow(q0i/Qt2InvA,-1.5);
    Qt2InvA = 1/q0i  - (r1-q0i)/(q0i*q0i);
    Qt2InvA = 1/Qt2InvA;
    weight *= (Qt2InvA*Qt2InvA)/(q0i*q0i);
  } 

  weight *= Norm;

  Double qn = 1./sqrt(Qt2InvA);

  weight *= qn*qn; //Due to importance sampling in 1/q2



  //Integration over loop momenta
  //Qt2InvA =  rand()/(Double(RAND_MAX) *  MinQt2) ;
  //Qt2InvB =  rand()/(Double(RAND_MAX) *  MinQt2);

  //weight *= 1./MinQt2 * 1./MinQt2;


  if( qn > 10 || qn > MjuLoop ) //add due to SuperChic
    return lum;

  //Mju = M/2;
  Double LnMju2 = 2*log(Mju);


  //Possibility to generate events without collinear approximation
  Double q1n, q2n;
  Vec2  q, q1, q2;
  
  Double phi = Uniform(0,2*M_PI);


  //Importance sampling in y^2*(1-y)^0.5
  Double y1, y2, Int1, Int2;

	if( id1 == 21 ) {
		GetYintGluon(x1, y1, Int1);
		GetYintGluon(x2, y2, Int2);
	}
	else {
		GetYintQuark(x1, y1, Int1);
		GetYintQuark(x2, y2, Int2);
	}

  complex<Double> lumTemp[3] = {0,0,0};
  //Loop over symetric Luminosity configurations
  for(int i=0; i < 4; ++i) {
    Double weightNow = 1;
    
    q.setRPhi( qn, phi + i*M_PI/2 );

    q1 =  q - p1;
    q2 = -q - p2;

    q1n = min(qn, q1.norm() );
    q2n = min(qn, q2.norm() );

    //cut applied to remove non-perturbative events
    if( min(q1n,q2n) < sqrt(MinQt2) ) 
      weightNow = 0;

    weightNow /=  q1.norm2() * q2.norm2() ;


    //cout << "RADEK "<< y1 << " " << Int1 << endl;

    weightNow *=  max(LeifGluonChicUnintegrated(id1, x1, q1n, y1, Int1), 0.0) *
                  max(LeifGluonChicUnintegrated(id2, x2, q2n, y2, Int2), 0.0);

		//cout << "leif2 params "<<x2<<" "<< q2n<<" "<< y2<<" "<< Int2 << endl;
		//cout << "weightNow "<<setprecision(10)<< LeifGluonChicUnintegrated(x1, q1n, y1, Int1) << " " << LeifGluonChicUnintegrated(x2, q2n, y2, Int2) << endl;
		//cout << "FunNow "<<LeifGluonChicUnintegrated(0.5021384231,1.1548492578,0.1281635546,0.1517533879) << endl; 
		//gluonPairChic(0.9794875475,1.1548492578);
		//exit(0);

    Double LnQ21 =  2*log(q1n);
    Double LnQ22 =  2*log(q2n);

    weightNow *= exp(-0.5*( InsideSudakov(id1, LnQ21, LnMju2)+InsideSudakov(id2, LnQ22, LnMju2) ));

    const complex<Double> i_(0,1);

    lumTemp[0] +=weightNow* -1./2   *  q1*q2 ;
    lumTemp[1] +=weightNow* -i_/2.*  cross(q1,q2) ;

    lumTemp[2] +=weightNow* 1./2* ( ( q1.x()*q2.x() - q1.y()*q2.y() )   - i_* ( q1.x()*q2.y() + q1.y()*q2.x() ) );

		//cout << "LumTemp "<<i<<" "<<lumTemp[0]<<" "<<lumTemp[1]<<" "<<lumTemp[2]<<endl;
    /*
    qpA2 = 1./2* ( ( q1.x()*q2.x() - q1.y()*q2.y() )   + Icmp* ( q1.x()*q2.y() + q1.y()*q2.x() ) );

    */


  }

  for(int j=0; j <3; ++j)
    lum[j] = 1./4* lumTemp[j] * weight;

  return lum;
}

void KMRlumi::LuminosityRand(int id1, Double M, Double y, Double Mju, Double MjuLoop, Vec2 &p1, Vec2 &p2, complex<Double> lum[], complex<Double> lumC[] ) 
{
  const int Niter = 12;

  
  vector< complex<Double> > lum1(3,0.0);
  vector< complex<Double> > lum2(3,0.0);

  Double rnd1R, rnd1Phi;
  Double rnd2R, rnd2Phi;

  do {
  rnd1R = sqrt( -1./bSlope * log(Uniform(0,1)) );
  rnd1Phi = Uniform(0,2*M_PI);

  rnd2R = sqrt( -1./bSlope * log(Uniform(0,1)) );
  rnd2Phi = Uniform(0,2*M_PI);
  } while( rnd1R > 2 || rnd2R > 2);

  p1.setRPhi( rnd1R, rnd1Phi);
  p2.setRPhi( rnd2R, rnd2Phi);


  for(int k=0; k < Niter; ++k) {
    vector< complex<Double> > lum1Temp;
    vector< complex<Double> > lum2Temp;

    lum1Temp = LumSqrtRandom(id1, M, y, Mju, MjuLoop, p1, p2);
    lum2Temp = LumSqrtRandom(id1, M, y, Mju, MjuLoop, p1, p2);
    for(int i=0; i< 3; ++i) {
      lum1[i] += lum1Temp[i]/Double(Niter);
      lum2[i] += lum2Temp[i]/Double(Niter);
    }

		if(lum1[0].real() !=  lum1[0].real() )
			exit(1);

  }


  for(int i=0; i< 3; ++i) {
    lum[i]  = lum1[i];
    lumC[i] = conj(lum2[i]);
  }

}

Double KMRlumi::LuminosityInc(int id1, Double M, Double y, Double Mju)
{
	int id2 = (id1 == 21) ? id1 : -id1;

  Double sSqrt = sqrt(s);

  x1 = exp(y)   * M/sSqrt; 
  x2 = exp(-y)  * M/sSqrt; 

  return gluon(id1,x1,Mju*Mju) * gluon(id2,x2,Mju*Mju); 

}

Double KMRlumi::SplittingAlphaS(int id, Double LnScale2) const
{
  //cout <<"alphaS "<<  2*M_PI*alphaStrong(2*log(91) ) << endl;
  return alphaStrong(LnScale2) * Splitting(id, exp(0.5*LnScale2));

}

//Alpha Strong over 2pi
Double KMRlumi::alphaStrong(Double LnScale2) const
{
  //const Double ASConst     = 2* 12.0 * M_PI/22.0/2.0; //need to be check

  const Double ASConst3     = 6. / (33. - 2*3);
  const Double ASConst4     = 6. / (33. - 2*4);
  const Double ASConst5     = 6. / (33. - 2*5);

  //const Double Lambda5QCD   = 0.001 * 146;


  const Double LnLambda5QCD2= 2*log(Lambda5QCD);

  const Double LnLambda4QCD2= LnBottomMass2 - ASConst4/ASConst5 * (LnBottomMass2 - LnLambda5QCD2);
  const Double LnLambda3QCD2= LnCharmMass2  - ASConst3/ASConst4 * (LnCharmMass2  - LnLambda4QCD2);

  
  if(LnScale2 > LnBottomMass2)
    return ASConst5 / ( LnScale2  - LnLambda5QCD2 );
  else if(LnScale2 > LnCharmMass2)
    return ASConst4 / ( LnScale2  - LnLambda4QCD2 );
  else if(LnScale2 > LnFreeze2)
    return ASConst3 / ( LnScale2  - LnLambda3QCD2 );
  else
    return ASConst3 / ( LnFreeze2 - LnLambda3QCD2 );

}


Double KMRlumi::Splitting(int id, Double Kt) const
{

	Double M = sqrt(x1*x2*s);
	Double D = 1 - Kt/M * ( sqrt(1+ 0.25*Kt*Kt/M/M) - 0.5*Kt/M   ); //Pythia choice
	//Double D = 1 - Kt/M; //KMR choice


	if( id == 21) {
		int Nf;

				 if(Kt < CharmMass ) Nf = 3;
		else if(Kt < BottomMass) Nf = 4;
		else                     Nf = 5;
			
			
		Double D2 = D * D;
		Double D3 = D2 * D;

		return(
					6.0*(-D2 +  1./3 * D3 - 0.25 *D3 * D - log(1 - D)) +
					Nf *(0.5 * (D - D2) + 1./3 * D3 )
				);
	}
	else {
		Double zmax = D;
		Double zmin = 1-D;
	  Double result = 2* log(zmax/zmin) - 2*(zmax-zmin) +1./2*(zmax*zmax-zmin*zmin);
	  return 4./3* result;
	}

}



Double KMRlumi::InsideSudakov(int id, Double LnQ2, Double LnM2) const
{
  Double norm;
  if(LnQ2 < LnM2)
    norm = 1;
  else {
    norm = -1;
    swap(LnQ2, LnM2);
  }

  #define Includes(i) ((LnQ2 < breakPoints[i] && breakPoints[i] < LnM2 ) ? true : false )

  Double breakPoints[] = {LnFreeze2,  LnCharmMass2,  LnBottomMass2 };
  int nBP  = sizeof(breakPoints)/sizeof(breakPoints[0]);

  int start=-1, end=-1;

  for(int i=0; i < nBP; ++i) {
    if(Includes(i) && start==-1) start = i;
    if(Includes(i) && (i==nBP-1 || !Includes(i+1)) ) end = i;
  }
  
  Double result = 0;
  Double err    = 0;

  if(start < end) {
    for(int i= start; i<nBP-1 && i<end; ++i)
      result += integral(&KMRlumi::SplittingAlphaS,id, breakPoints[i], breakPoints[i+1], err );
  }

  if(start != -1) {
      result += integral(&KMRlumi::SplittingAlphaS,id, LnQ2, breakPoints[start], err );
      result += integral(&KMRlumi::SplittingAlphaS,id, breakPoints[end],   LnM2, err );
  }
  else {
      result += integral(&KMRlumi::SplittingAlphaS,id, LnQ2, LnM2, err );
  }

  //cout <<"Error "<<  err/result << endl;
  if(err/result> 1e-5)
    cout << "Relative error larger than 1e-05" << endl;
  return norm*result;
}







//Error is added to the original value
Double KMRlumi::integral( Double (KMRlumi::*func)(int id, Double x) const, int id, Double xmin, Double xmax, Double &err ) const
{

  static STATUS st;

  Double Mred = x1*x2;
  //st.print();

  //cout << "Size " << st.points.size() << " "<< 1e6*Mred <<endl;
  //If the integral exist
  if( abs(Mred - st.Mred) < 1e-12 * st.Mred ) {
    int size = st.points.size();
    for(int i = 0; i < size; ++i)
      if( st.points[i].isSame(xmin, xmax) )
        return st.points[i].Int;
  }


  static const Double nodsG7[] = {0,             0.405845151377397, 0.741531185599394, 0.949107912342759 };
  static const Double wG7[]= {0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870 };
  static const Double wK7[]= {0.209482141084728, 0.190350578064785, 0.140653259715525, 0.063092092629979 };

  static const Double nodsK15[] = {0.207784955007898, 0.586087235467691, 0.864864423359769, 0.991455371120813 };
  static const Double wK15[]= {    0.204432940075298, 0.169004726639267, 0.104790010322250, 0.022935322010529 };

  //cout <<"Masses "<<LnFreeze2<<" "<< LnCharmMass2  <<" "<< LnBottomMass2  << endl;
  //cout <<"Limits "<< 1e6*x1*x2 <<" "<<xmin <<" "<< xmax << endl;

  Double valuesG[7];
  Double valuesK[7];

  Double xmean = 0.5*(xmin+xmax);
  Double xdelta = xmean - xmin;

  valuesG[0] = (this->*func)(id, xmean);
  for(int i=1; i < 4; ++i) {
    valuesG[i] = (this->*func)(id, xmean + nodsG7[i]*xdelta) +  (this->*func)(id, xmean - nodsG7[i]*xdelta);
  }
  for(int i=0; i < 4; ++i) {
    valuesK[i] = (this->*func)(id, xmean + nodsK15[i]*xdelta) + (this->*func)(id, xmean - nodsK15[i]*xdelta);
  }


  Double intG = 0;
  Double intK = 0;

  for(int i=0; i < 4; ++i) {
    intG += wG7[i] * valuesG[i] ;

    intK += wK7[i]  * valuesG[i];
    intK += wK15[i] * valuesK[i];
  }
  //cout << intK << " "<< intG << endl;

  Double errTemp;
  errTemp = 200*abs(intK - intG);
  errTemp *= sqrt(err);
  errTemp *= xdelta;

  err += errTemp;

  Double Int = intK*xdelta;
  
  SaveIntegral(st, xmin, xmax, Int);



  return Int;
}

void KMRlumi::SaveIntegral(STATUS &st, Double xmin, Double xmax, Double Int) const
{
  Double Mred = x1*x2;

  if( abs(Mred - st.Mred) < 1e-12 * st.Mred ) {

    int size = st.points.size();
    int i;
    for(i = 0; i < size; ++i) {
      if( close(st.points[i].xmin,xmin) && !close(st.points[i].xmax,LnCharmMass2)
     &&  !close(st.points[i].xmax,LnBottomMass2) && !close(st.points[i].xmax,LnFreeze2) ) {
        st.points[i] = DATA(xmin, xmax, Int);
        break;
      }
      if( close(st.points[i].xmax,xmax) &&  !close(st.points[i].xmin,LnCharmMass2)
     &&  !close(st.points[i].xmin,LnBottomMass2) && !close(st.points[i].xmin,LnFreeze2) ) {
        st.points[i] = DATA(xmin, xmax, Int);
        break;
      }
    }
    if(i == size) 
      st.points.push_back( DATA(xmin, xmax, Int) );
  }
  else {
    st.Mred = Mred;
    st.points.clear();
    st.points.push_back( DATA(xmin, xmax, Int) );
  }

}



Double KMRlumi::NoEmProbability(Double _x1, Double _x2, Double LnQ2, Double LnM2, bool backward) 
{
  x1 = _x1;
  x2 = _x2;

  Double res = exp( - 2.*InsideSudakov(21, LnQ2, LnM2) );

  if(backward) {
		res *=
		gluon(21, x1, exp(LnQ2) ) / gluon(21, x1, exp(LnM2) ) *
		gluon(21, x2, exp(LnQ2) ) / gluon(21, x2, exp(LnM2) );
	}

  if(res > 1 || res < 0) {
    cout << "Strange value of prob " << res <<", Scales: "<< exp(0.5*LnQ2) <<" "<< exp(0.5*LnM2) << endl;
		//res = min(0.98, res);
    exit(1);
  }

  return res;

}
