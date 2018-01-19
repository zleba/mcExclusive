#include "SigmaQCD.h"
#include "KMRlumi.h"




inline Double Uniform(Double a, Double b)
{
  //return (  a +  (b-a) * rand() / (RAND_MAX + 0.0)  );
  return KMRlumi::Uniform(a,b);

}






//=========================================================================

//Derived class g g -> g g

void Sigma2gg2gg::ColorMatrix()
{

  Double A12, A13, A14;
  Double B234, B243, B342, B432;

  DataGG( A12, A13, A14,
          B234, B243, B342, B432);

  colMatrix[0][0] = 1./36*A12 + 7./12*B243;
  colMatrix[0][1] = 5./18*A12 - 1./6*B243;
  colMatrix[0][2] = 1./36*A12 - 1./12*B243 - 1./12*B432;
  colMatrix[0][3] = 1./36*A12 + 7./12*B432;
  colMatrix[0][4] = 1./36*A12 - 1./12*B243 - 1./12*B432;
  colMatrix[0][5] = 5./18*A12 - 1./6*B432;

  colMatrix[1][0] = 5./18*A12 - 1./6*B243;
  colMatrix[1][1] = 1./36*A12 + 7./12*B243;
  colMatrix[1][2] = 1./36*A12 - 1./12*B243 - 1./12*B432;
  colMatrix[1][3] = 5./18*A12 - 1./6*B432;
  colMatrix[1][4] = 1./36*A12 - 1./12*B243 - 1./12*B432;
  colMatrix[1][5] = 1./36*A12 + 7./12*B432;

  colMatrix[2][0] = 1./36*A12 - 1./12*B234 - 1./12*B243;
  colMatrix[2][1] = 1./36*A12 - 1./12*B234 - 1./12*B243;
  colMatrix[2][2] = 1./36*A12 - 1./12*B243 + 1./4*A13 - 1./12*B342;
  colMatrix[2][3] = 1./36*A12 - 1./12*B342 - 1./12*B432;
  colMatrix[2][4] = 1./36*A12 - 1./12*B234 + 1./4*A14 - 1./12*B432;
  colMatrix[2][5] = 1./36*A12 - 1./12*B342 - 1./12*B432;

  colMatrix[3][0] = 1./36*A12 + 7./12*B234;
  colMatrix[3][1] = 5./18*A12 - 1./6*B234;
  colMatrix[3][2] = 1./36*A12 - 1./12*B234 - 1./12*B342;
  colMatrix[3][3] = 1./36*A12 + 7./12*B342;
  colMatrix[3][4] = 1./36*A12 - 1./12*B234 - 1./12*B342;
  colMatrix[3][5] = 5./18*A12 - 1./6*B342;

  colMatrix[4][0] = 1./36*A12 - 1./12*B234 - 1./12*B243;
  colMatrix[4][1] = 1./36*A12 - 1./12*B234 - 1./12*B243;
  colMatrix[4][2] = 1./36*A12 - 1./12*B234 + 1./4*A14 - 1./12*B432;
  colMatrix[4][3] = 1./36*A12 - 1./12*B342 - 1./12*B432;
  colMatrix[4][4] = 1./36*A12 - 1./12*B243 + 1./4*A13 - 1./12*B342;
  colMatrix[4][5] = 1./36*A12 - 1./12*B342 - 1./12*B432;

  colMatrix[5][0] = 5./18*A12 - 1./6*B234;
  colMatrix[5][1] = 1./36*A12 + 7./12*B234;
  colMatrix[5][2] = 1./36*A12 - 1./12*B234 - 1./12*B342;
  colMatrix[5][3] = 5./18*A12 - 1./6*B342;
  colMatrix[5][4] = 1./36*A12 - 1./12*B234 - 1./12*B342;
  colMatrix[5][5] = 1./36*A12 + 7./12*B342;

}


void Sigma2gg2gg::Amplitudes()
{
  complex<Double> Exp = exp( Double(2)*i_*phi );

  Double norm = g2S/sqrt(Nsym) *(-4.);

  ///////////////////////
  //++ to ++
  JAMP[pp][pp][0] = norm*(s/u);
  JAMP[pp][pp][1] = norm*(s/t);
  JAMP[pp][pp][2] = norm*(s*s/u/t);
  JAMP[pp][pp][3] = norm*(s/t);
  JAMP[pp][pp][4] = norm*(s*s/u/t);
  JAMP[pp][pp][5] = norm*(s/u);

  NonZero[pp][pp] = true;
  CopyConj(pp,pp);

  ///////////////////////
  //+- to +-
  JAMP[pm][pm][0] = norm*Exp*( u/s );
  JAMP[pm][pm][1] = norm*Exp*(  u*u/(t*s) );
  JAMP[pm][pm][2] = norm*Exp*( u/t   );
  JAMP[pm][pm][3] = norm*Exp*(  u*u/(t*s) );
  JAMP[pm][pm][4] = norm*Exp*( u/t   );
  JAMP[pm][pm][5] = norm*Exp*( u/s );

  NonZero[pm][pm] = true;
  CopyConj(pm,pm);

  ///////////////////////
  //+- to -+ 
  JAMP[pm][mp][0] = norm*Exp*(  t*t/(u*s) );
  JAMP[pm][mp][1] = norm*Exp*(  s*t/(s*s) );
  JAMP[pm][mp][2] = norm*Exp*(  t/u       );
  JAMP[pm][mp][3] = norm*Exp*(  s*t/(s*s) );
  JAMP[pm][mp][4] = norm*Exp*(  t/u       );
  JAMP[pm][mp][5] = norm*Exp*(  t*t/(u*s) );

  NonZero[pm][mp] = true;
  CopyConj(pm,mp);

}


Double Sigma2gg2gg::PyCrossSection()
{
  Double t2 = t*t;
  Double s2 = s*s;
  Double u2 = u*u;
  // Calculate kinematics dependence.
  sigTS  = (9./4.) * (t2 / s2 + 2. * t / s + 3. + 2. * s / t
      + s2 / t2);
  sigUS  = (9./4.) * (u2 / s2 + 2. * u / s + 3. + 2. * s / u
      + s2 / u2);
  sigTU  = (9./4.) * (t2 / u2 + 2. * t / u + 3. + 2. * u / t
      + u2 / t2);
  sigSum = sigTS + sigUS + sigTU;

  // Answer with identical factor 0.5
  return  0.5* g2S*g2S *  sigSum;
}

ColArr Sigma2gg2gg::PyColorFlow()
{
  ColArr cl;
  // Three colour flow topologies, each with two orientations.
  Double sigRand = sigSum * Uniform(0,1);
  if (sigRand < sigTS) cl.setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
  else if (sigRand < sigTS + sigUS)
    cl.setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
  else                 cl.setColAcol( 1, 2, 3, 4, 1, 4, 3, 2);

  if (Uniform(0,1) > 0.5) cl.swapColAcol();

  return cl;
}

void Sigma2gg2gg::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
	pArr[0] = 0.5* sigTS / sigSum;

	clArr[1].setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
	pArr[1] = 0.5* sigUS / sigSum;

	clArr[2].setColAcol( 1, 2, 3, 4, 1, 4, 3, 2);
	pArr[2] = 0.5* sigTU / sigSum;

	//swap colors and anti-colors for all color topologies
	clArr[3] = clArr[0]; clArr[3].swapColAcol();
	clArr[4] = clArr[1]; clArr[4].swapColAcol();
	clArr[5] = clArr[2]; clArr[5].swapColAcol();

	pArr[3] = pArr[0];
	pArr[4] = pArr[1];
	pArr[5] = pArr[2];

}




//=========================================================================

//Derived class g g -> q qbar (uds)
//For one particular flavour


void Sigma2gg2qqbar::ColorMatrix()
{


  Double A12, A13, A14;
  Double B234, B243, B342, B432;

  DataGG( A12, A13, A14,
          B234, B243, B342, B432);



  colMatrix[0][0] = B243;
  colMatrix[0][1] = B234;

  colMatrix[1][0] = B432;
  colMatrix[1][1] = B342;


}


void Sigma2gg2qqbar::Amplitudes()
{

  complex<Double> Exp = exp( Double(2)*i_*phi );

  ///////////////////////
  //+- to +-
  JAMP[pm][pm][0] =-g2S*2.* u/s * sqrt(u/t)*Exp;
  JAMP[pm][pm][1] =+g2S*2.*1./s * sqrt(t*u)*Exp;

  NonZero[pm][pm] = true;
  CopyConj(pm,pm);


  ///////////////////////
  //+- to -+
  JAMP[pm][mp][0] =-g2S*2.*1./s * sqrt(t*u) *Exp;
  JAMP[pm][mp][1] =+g2S*2.* t/s * sqrt(t/u) *Exp;

  NonZero[pm][mp] = true;
  CopyConj(pm,mp);


}

Double Sigma2gg2qqbar::PyCrossSection()
{

  Double t2 = t*t;
  Double s2 = s*s;
  Double u2 = u*u;



  sigTS = (1./6.) * u / t - (3./8.) * u2 / s2;
  sigUS = (1./6.) * t / u - (3./8.) * t2 / s2;

  sigSum = sigTS + sigUS;

  return g2S*g2S* sigSum;


}


ColArr Sigma2gg2qqbar::PyColorFlow()
{
  ColArr cl;


  Double sigRand = sigSum * Uniform(0,1);
  if (sigRand < sigTS) cl.setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 cl.setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);

  return cl;
}


void Sigma2gg2qqbar::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
	pArr[0] = sigTS / sigSum;

	clArr[1].setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);
	pArr[1] = sigUS / sigSum;

}




//=========================================================================

//Derived class q g -> q g (uds c b)

void Sigma2qg2qg::ColorMatrix()
{

  Double Delta, Tl1cl1C, Tl1l1c;


  DataGQ( Delta, Tl1cl1C, Tl1l1c);


  colMatrix[0][0] =  1./2 * Delta - 1./6 * Tl1cl1C;
  colMatrix[0][1] = -1./6 * Tl1cl1C;

  colMatrix[1][0] = -1./6 * Tl1cl1C;
  colMatrix[1][1] =  4./3 * Tl1cl1C;

  //colMatrix[0][0] =  4./3 * Tl1cl1C;


}

void Sigma2qg2qg::Amplitudes()
{
  //reverse order ?
  complex<Double> ExpP = exp(+i_*phi );
  complex<Double> ExpM = conj(ExpP);

  /*
  //g u -> g u
  //+- to +-
  JAMP[pm][pm][0] =+g2S*2.*Exp*Exp* sqrt(-u/s)*  s / t ;
  JAMP[pm][pm][1] =+g2S*2.*Exp*Exp* sqrt(-u/s)*  u / t ;
  CopyConj(pm,pm);

  //++ to ++
  JAMP[pp][pp][0] =+g2S*2.* sqrt(-s/u) * s / t ;
  JAMP[pp][pp][1] =+g2S*2.* sqrt(-s/u) * u / t ;
  CopyConj(pp,pp);

  return;
  */

  //u g -> u g
  if(idLeft > 0) {
    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] =+g2S*2.* sqrt(-u/s)* s / t *ExpP;
    JAMP[pm][pm][1] =+g2S*2.* sqrt(-u/s)* u / t *ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);

    ///////////////////////
    //++ to ++
    JAMP[pp][pp][0] =+g2S*2.* sqrt(-s/u) * s / t *ExpM;
    JAMP[pp][pp][1] =+g2S*2.* sqrt(-s/u) * u / t *ExpM;

    NonZero[pp][pp] = true;
    CopyConj(pp,pp);

  }


  //ux g -> ux g
  else{

    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] =+g2S*2.* sqrt(-u/s)* u / t *ExpP;
    JAMP[pm][pm][1] =+g2S*2.* sqrt(-u/s)* s / t *ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);

    ///////////////////////
    //++ to ++
    JAMP[pp][pp][0] =+g2S*2.* sqrt(-s/u) * u / t *ExpM;
    JAMP[pp][pp][1] =+g2S*2.* sqrt(-s/u) * s / t *ExpM;

    NonZero[pp][pp] = true;
    CopyConj(pp,pp);

  }

}

Double  Sigma2qg2qg::PyCrossSection()
{

  Double t2 = t*t;
  Double s2 = s*s;
  Double u2 = u*u;

  // Calculate kinematics dependence.
  sigTS  = u2 / t2 - (4./9.) * u / s;
  sigTU  = s2 / t2 - (4./9.) * s / u;
  sigSum = sigTS + sigTU;

  // Answer.
  return  g2S*g2S * sigSum;

}


ColArr Sigma2qg2qg::PyColorFlow()
{
  ColArr cl;

  // Two colour flow topologies. Swap if first is gluon, or when antiquark.
  Double sigRand = sigSum * Uniform(0,1);
  if (sigRand < sigTS) cl.setColAcol( 1, 0, 2, 1, 3, 0, 2, 3);
  else                 cl.setColAcol( 1, 0, 2, 3, 2, 0, 1, 3);

  if (idLeft == 21) cl.swapCol1234();
  if (idLeft < 0 || idRight < 0) cl.swapColAcol();

  return cl;
}

void Sigma2qg2qg::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 0, 2, 1, 3, 0, 2, 3);
	pArr[0] = sigTS / sigSum;

	clArr[1].setColAcol( 1, 0, 2, 3, 2, 0, 1, 3);
	pArr[1] = sigTU / sigSum;

  if (idLeft == 21) {
		clArr[0].swapCol1234(); clArr[1].swapCol1234();
	}
  if (idLeft < 0 || idRight < 0) {
		clArr[0].swapColAcol(); clArr[1].swapColAcol();
	}

}




//=========================================================================

//Derived class q q(bar) -> q q(bar) (uds c b)

void Sigma2qq2qq::ColorMatrix()
{

  Double Dl1l1c, Dl1r1c;

  Double Dl1r1; 

  DataQQ( Dl1l1c, Dl1r1c, Dl1r1);


  colMatrix[0][0] = Dl1l1c;
  colMatrix[0][1] = Dl1r1c;

  colMatrix[1][0] = Dl1r1c;
  colMatrix[1][1] = Dl1l1c;

}


void Sigma2qq2qq::Amplitudes()
{
  
  complex<Double> ExpP = exp(+i_*phi);
  complex<Double> ExpM = conj(ExpP);

  //uu -> uu || ux ux -> ux ux
  if(idLeft == idRight ) { 
    Nsym = 2;
    Double norm = 1/sqrt(Nsym);

    ///////////////////////
    //++ to ++
    JAMP[pp][pp][0] =-2*g2S*norm*(-1./6*  s/t + 1./2 * s/u )*  ExpM;
    JAMP[pp][pp][1] =-2*g2S*norm*(+1./2*  s/t - 1./6 * s/u )*  ExpM;

    NonZero[pp][pp] = true;
    CopyConj(pp,pp);

    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] = 2*g2S*norm*(-1./6)* u/t *  ExpP;
    JAMP[pm][pm][1] = 2*g2S*norm* 1./2  * u/t *  ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);

    ///////////////////////
    //+- to -+
    JAMP[pm][mp][0] = 2*g2S*norm* 1./2  * t/u *  ExpP;
    JAMP[pm][mp][1] = 2*g2S*norm*(-1./6)* t/u *  ExpP;

    NonZero[pm][mp] = true;
    CopyConj(pm,mp);

  }
  //u ux -> u ux
  else if(idLeft == - idRight ) {
    Nsym = 1;
    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] = 2*g2S*(+1./6* u/s + -1./2 * u/t  ) * ExpP;
    JAMP[pm][pm][1] = 2*g2S*(-1./2* u/s + +1./6 * u/t  ) * ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);

    ///////////////////////
    //+- to -+
    JAMP[pm][mp][0] =-1./6* 2*g2S* t/s * ExpP;
    JAMP[pm][mp][1] =+1./2* 2*g2S* t/s * ExpP;

    NonZero[pm][mp] = true;
    CopyConj(pm,mp);

    ///////////////////////
    //++ to ++
    JAMP[pp][pp][0] =+1./2* 2*g2S* s/t * ExpM;
    JAMP[pp][pp][1] =-1./6* 2*g2S* s/t * ExpM;

    NonZero[pp][pp] = true;
    CopyConj(pp,pp);
  }

  //ud -> ud || ux dx -> ux dx
  else if(idLeft != idRight  && idLeft*idRight > 0 ) {
    Nsym = 1;

    ///////////////////////
    //++ to ++
    JAMP[pp][pp][0] =+1./6* 2*g2S* s/t * ExpM;
    JAMP[pp][pp][1] =-1./2* 2*g2S* s/t * ExpM;

    NonZero[pp][pp] = true;
    CopyConj(pp,pp);

    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] =-1./6* 2*g2S* u/t * ExpP;
    JAMP[pm][pm][1] =+1./2* 2*g2S* u/t * ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);


    ///////////////////////
    //+- to -+
    NonZero[pm][mp] = false;
    SetZero(pm,mp);

  }

  //udx -> udx
  else if(idLeft != idRight  && idLeft*idRight < 0 ) {
    Nsym = 1;
    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] = 2*g2S*( -1./2 * u/t  ) * ExpP;
    JAMP[pm][pm][1] = 2*g2S*( +1./6 * u/t  ) * ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);


    ///////////////////////
    //++ to ++
    JAMP[pp][pp][0] =+1./2* 2*g2S* s/t * ExpM;
    JAMP[pp][pp][1] =-1./6* 2*g2S* s/t * ExpM;

    NonZero[pp][pp] = true;
    CopyConj(pp,pp);

    ///////////////////////
    //+- to -+
    NonZero[pm][mp] = false;
    SetZero(pm,mp);
  }
  else {
    cout << "Problem with process identification "<<
            "in Sigma2qq2qq::Amplitudes()" << endl;
    exit(1);

  }
  
}


Double Sigma2qq2qq::PyCrossSection()
{
  Double t2 = t*t;
  Double s2 = s*s;
  Double u2 = u*u;

  sigT   = (4./9.) * (s2 + u2) / t2;
  sigU   = (4./9.) * (s2 + t2) / u2;
  sigTU  = - (8./27.) * s2 / (t * u);
  sigST  = - (8./27.) * u2 / (s * t);

  sigS   = (4./9.) * (t2 + u2) / s2;

  if (idLeft == idRight)
    nColFlow = 2;
  else
    nColFlow = 1;


  if      (idLeft ==  idRight) sigSum = 0.5 * (sigT + sigU + sigTU);
  else if (idLeft == -idRight) sigSum = sigT + sigST + sigS;
  else                         sigSum = sigT;


  return g2S*g2S * sigSum;
}


ColArr Sigma2qq2qq::PyColorFlow()
{
  ColArr cl;

  // Colour flow topologies. Swap when antiquarks.
  if (idLeft * idRight > 0)  cl.setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else										 	 cl.setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  if (idLeft == idRight && (sigT + sigU) * Uniform(0,1) > sigT)
                      cl.setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  if (idLeft < 0) cl.swapColAcol();


  return cl;
}


void Sigma2qq2qq::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);


	if (idLeft * idRight > 0)  clArr[0].setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
	else											 clArr[0].setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);

	if (idLeft < 0) clArr[0].swapColAcol();


	if(idLeft != idRight) {
		pArr[0] = 1;
	}
	else {
		pArr[0] = sigT/(sigT+sigU);

		clArr[1].setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
		pArr[1] = sigU/(sigT+sigU);

		if (idLeft < 0) clArr[1].swapColAcol();
		
	}

}





//=========================================================================

//Derived class q qbar -> g g


void Sigma2qqbar2gg::ColorMatrix()
{
  Double Dl1l1c, Dl1r1c;

  Double Dl1r1; 

  DataQQ( Dl1l1c, Dl1r1c, Dl1r1);


  colMatrix[0][0] =  7/12.*Dl1l1c + 1/36.*Dl1r1;

  colMatrix[0][1] = - 1/6.*Dl1l1c + 5/18.*Dl1r1;

  colMatrix[1][1] =  7/12.*Dl1l1c + 1/36.*Dl1r1;

  colMatrix[1][0] = - 1/6.*Dl1l1c + 5/18.*Dl1r1;

}



void Sigma2qqbar2gg::Amplitudes()
{
  
  complex<Double> Exp = exp(+i_*phi);
  Double norm = 1/sqrt(Nsym);


  ///////////////////////
  //+- to +-
  JAMP[pm][pm][1] =-g2S*norm*2.*  u/s *sqrt(u/t) *  Exp;
  JAMP[pm][pm][0] =+g2S*norm*2.* 1./s *sqrt(t*u) *  Exp;

  NonZero[pm][pm] = true;
  CopyConj(pm,pm);

  ///////////////////////
  //+- to -+
  JAMP[pm][mp][0] =+g2S*norm*2.* t/s* sqrt(t/u)  *  Exp;
  JAMP[pm][mp][1] =-g2S*norm*2.*1./s* sqrt(t*u)  *  Exp;

  NonZero[pm][mp] = true;
  CopyConj(pm,mp);


}


Double Sigma2qqbar2gg::PyCrossSection()
{

  Double t2 = t*t;
  Double s2 = s*s;
  Double u2 = u*u;


  // Calculate kinematics dependence.
  sigTS  = (32./27.) * u / t - (8./3.) * u2 / s2;
  sigUS  = (32./27.) * t / u - (8./3.) * t2 / s2;
  sigSum = sigTS + sigUS;

  // Answer contains factor 1/2 from identical gluons.
  return  g2S*g2S * 0.5 * sigSum;

}


ColArr Sigma2qqbar2gg::PyColorFlow()
{
  ColArr cl;

  Double sigRand = sigSum * Uniform(0,1);
  if (sigRand < sigTS) cl.setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                 cl.setColAcol( 1, 0, 0, 2, 3, 2, 1, 3);
  if (idLeft < 0) cl.swapColAcol();


  return cl;
}

void Sigma2qqbar2gg::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
	pArr[0] = sigTS / sigSum;

	clArr[1].setColAcol( 1, 0, 0, 2, 3, 2, 1, 3);
	pArr[1] = sigUS / sigSum;

  if (idLeft < 0) {
		clArr[0].swapColAcol();
		clArr[1].swapColAcol();
	}

}







//=========================================================================

//Derived class q qbar -> q' qbar' (q'=uds, q != q')



void Sigma2qqbar2qqbarNew::ColorMatrix()
{
  Double Dl1l1c, Dl1r1c;
  Double Dl1r1; 

  DataQQ( Dl1l1c, Dl1r1c, Dl1r1);

  if(abs(idOut1) != abs(idLeft)) {
    colMatrix[0][0] =  3*Dl1r1;
    colMatrix[0][1] =  Dl1r1;

    colMatrix[1][1] =  Dl1l1c;
    colMatrix[1][0] =  Dl1r1;
  }
  else {
    colMatrix[0][0] = Dl1l1c;
    colMatrix[0][1] = Dl1r1c;

    colMatrix[1][0] = Dl1r1c;
    colMatrix[1][1] = Dl1l1c;
  }


}


void Sigma2qqbar2qqbarNew::Amplitudes()
{

  complex<Double> ExpP = exp(+i_*phi);
  complex<Double> ExpM = conj(ExpP);

  //q qbar -> q' qbar'
  if(abs(idOut1) != abs(idLeft)) {
    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] =+1./6* 2*g2S* u/s * ExpP;
    JAMP[pm][pm][1] =-1./2* 2*g2S* u/s * ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);

    ///////////////////////
    //+- to -+
    JAMP[pm][mp][0] =-1./6* 2*g2S* t/s * ExpP;
    JAMP[pm][mp][1] =+1./2* 2*g2S* t/s * ExpP;

    NonZero[pm][mp] = true;
    CopyConj(pm,mp);

    ///////////////////////
    //++ to ++
    NonZero[pp][pp] = false;
    SetZero(pp,pp);

  }
  //q qbar -> q qbar
  else {
    ///////////////////////
    //+- to +-
    JAMP[pm][pm][0] = 2*g2S*(+1./6* u/s + -1./2 * u/t  ) * ExpP;
    JAMP[pm][pm][1] = 2*g2S*(-1./2* u/s + +1./6 * u/t  ) * ExpP;

    NonZero[pm][pm] = true;
    CopyConj(pm,pm);

    ///////////////////////
    //+- to -+
    JAMP[pm][mp][0] =-1./6* 2*g2S* t/s * ExpP;
    JAMP[pm][mp][1] =+1./2* 2*g2S* t/s * ExpP;

    NonZero[pm][mp] = true;
    CopyConj(pm,mp);

    ///////////////////////
    //++ to ++
    JAMP[pp][pp][0] =+1./2* 2*g2S* s/t * ExpM;
    JAMP[pp][pp][1] =-1./6* 2*g2S* s/t * ExpM;

    NonZero[pp][pp] = true;
    CopyConj(pp,pp);

  }

}



Double Sigma2qqbar2qqbarNew::PyCrossSection()
{

  Double t2 = t*t;
  Double s2 = s*s;
  Double u2 = u*u;

  sigS  = (4./9.) * (t2 + u2) / s2;

  sigT  = (4./9.) * (s2 + u2) / t2;

  sigST = - (8./27.) * u2 / (s * t);

  Double sum;
  if(abs(idOut1) != abs(idLeft)) 
    sum = sigS;
  else
    sum = sigS + sigT + sigST;


  // Answer is proportional to number of outgoing flavours.
  return  g2S*g2S *  sum;

}


ColArr Sigma2qqbar2qqbarNew::PyColorFlow()
{
  ColArr cl;

  // Colour flow topologies. Swap when antiquarks.
  cl.setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
	if (idLeft < 0) cl.swapColAcol();

  return cl;
}


void Sigma2qqbar2qqbarNew::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
	pArr[0] = 1;

	if (idLeft < 0) clArr[0].swapColAcol();

}




//=========================================================================

//Derived class g g -> Q Qbar (c, b, t)



void Sigma2gg2QQbar::ColorMatrix()
{


  Double A12, A13, A14;
  Double B234, B243, B342, B432;

  DataGG( A12, A13, A14,
          B234, B243, B342, B432);



  colMatrix[0][0] = B243;
  colMatrix[0][1] = B234;

  colMatrix[1][0] = B432;
  colMatrix[1][1] = B342;

}


void Sigma2gg2QQbar::Amplitudes()
{

  Double M = mP3;

  Double beta = sqrt(1.0-4.0*M*M/s);

  Double E = sqrt(s)/2;


  Double  pz = E - (M*M-t)/2./E;
  Double  pT = sqrt( beta*E*beta*E - pz*pz);


  Double Sin = pT/(E*beta);
  


  complex<Double> Exp3P = exp( Double(3)*i_*phi );
  complex<Double> Exp2P = exp( Double(2)*i_*phi );
  complex<Double> Exp1P = exp( Double(1)*i_*phi );
  complex<Double> Exp1M = conj(Exp1P);




  //++ to --
  JAMP[pp][mm][0] = g2S * 1./(t-M*M) *M*sqrt(s)*(1-beta) *Exp1P;
  JAMP[pp][mm][1] = g2S * 1./(u-M*M) *M*sqrt(s)*(1-beta) *Exp1P;

  NonZero[pp][mm] = true;
  CopyConj(pp,mm, -1);




  //++ to ++
  JAMP[pp][pp][0] = g2S * 1./(t-M*M) * M*sqrt(s)*(1+beta) *Exp1M;
  JAMP[pp][pp][1] = g2S * 1./(u-M*M) * M*sqrt(s)*(1+beta) *Exp1M;

  NonZero[pp][pp] = true;
  CopyConj(pp,pp, -1);


  //+- to ++
  JAMP[pm][pp][0] = -g2S * 1./(t-M*M) * M*sqrt(s) * beta *Sin*Sin *Exp1P;
  JAMP[pm][pp][1] = -g2S * 1./(u-M*M) * M*sqrt(s) * beta *Sin*Sin *Exp1P;

  NonZero[pm][pp] = true;
  CopyConj(pm,pp, -1);



  //+- to +-
  JAMP[pm][pm][0] = -g2S * ( 1 + s/2*(1+beta)/(t-M*M) ) * Sin *Exp2P;
  JAMP[pm][pm][1] =  g2S * ( 1 + s/2*(1-beta)/(u-M*M) ) * Sin *Exp2P;

  NonZero[pm][pm] = true;
  CopyConj(pm,pm);




  //+- to -+
  JAMP[pm][mp][0] = -g2S * ( 1 + s/2*(1-beta)/(t-M*M) ) * Sin *Exp2P;
  JAMP[pm][mp][1] =  g2S * ( 1+  s/2*(1+beta)/(u-M*M) ) * Sin *Exp2P;

  NonZero[pm][mp] = true;
  CopyConj(pm,mp);




  //+- to --
  JAMP[pm][mm][0] = g2S * 1./(t-M*M) * M*sqrt(s)*beta *Sin*Sin *Exp3P;
  JAMP[pm][mm][1] = g2S * 1./(u-M*M) * M*sqrt(s)*beta *Sin*Sin *Exp3P;

  NonZero[pm][mm] = true;
  CopyConj(pm,mm, -1);



}


Double Sigma2gg2QQbar::PyCrossSection()
{
  Double s2 = s*s;

  Double m32 = mP3*mP3;
  Double m42 = mP3*mP3;

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  Double s34Avg = 0.5 * (m32 + m42) - 0.25 * (m32 - m42)*(m32 - m42) / s;
  Double tHQ    = -0.5 * (s - t + u);
  Double uHQ    = -0.5 * (s + t - u);
  Double tHQ2   = tHQ * tHQ;
  Double uHQ2   = uHQ * uHQ;

  // Calculate kinematics dependence.
  Double tumHQ = tHQ * uHQ - s34Avg * s;
  sigTS = ( uHQ / tHQ - 2.25 * uHQ2 / s2 + 4.5 * s34Avg * tumHQ
    / ( s * tHQ2) + 0.5 * s34Avg * (tHQ + s34Avg) / tHQ2
    - s34Avg*s34Avg / (s * tHQ) ) / 6.;
  sigUS = ( tHQ / uHQ - 2.25 * tHQ2 / s2 + 4.5 * s34Avg * tumHQ
    / ( s * uHQ2) + 0.5 * s34Avg * (uHQ + s34Avg) / uHQ2
    - s34Avg*s34Avg / (s * uHQ) ) / 6.;
  sigSum = sigTS + sigUS;

  // Answer.
  return  g2S*g2S * sigSum ;

}


ColArr Sigma2gg2QQbar::PyColorFlow()
{
  ColArr cl;

  // Two colour flow topologies.
  Double sigRand = sigSum * Uniform(0,1);
  if (sigRand < sigTS) cl.setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 cl.setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);

  return cl;
}


void Sigma2gg2QQbar::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{

  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
	pArr[0] = sigTS / sigSum;

	clArr[1].setColAcol( 1, 2, 3, 1, 3, 0, 0, 2);
	pArr[1] = sigUS / sigSum;

}






//=========================================================================

//Derived class q qbar -> Q Qbar (c, b, t)



void Sigma2qqbar2QQbar::ColorMatrix()
{

  Double Dl1l1c, Dl1r1c;
  Double Dl1r1; 

  DataQQ( Dl1l1c, Dl1r1c, Dl1r1);

  colMatrix[0][0] =  3*Dl1r1;

  colMatrix[0][1] =  Dl1r1;

  colMatrix[1][1] =  Dl1l1c;

  colMatrix[1][0] =  Dl1r1;


}





void Sigma2qqbar2QQbar::Amplitudes()
{

  complex<Double> Exp2P = exp( Double(2)*i_*phi );
  complex<Double> Exp1P = exp( Double(1)*i_*phi );

  Double M = mP3;

  Double beta = sqrt(1.0-4.0*M*M/s);

  Double E = sqrt(s)/2;


  Double  pz = E - (M*M-t)/2./E;
  Double  pT = sqrt( beta*E*beta*E - pz*pz);

  Double Sin = pT/(E*beta);
  Double Cos = pz/(E*beta);


  //+- to ++
  JAMP[pm][pp][0] =  - 1./3*g2S* M/sqrt(s) * Sin;
  JAMP[pm][pp][1] =  + 1.  *g2S* M/sqrt(s) * Sin;

  NonZero[pm][pp] = true;
  CopyConj(pm,pp, -1);

  //+- to --
  JAMP[pm][mm][0] =  + 1./3*g2S* M/sqrt(s) * Sin * Exp2P;
  JAMP[pm][mm][1] =  - 1.  *g2S* M/sqrt(s) * Sin * Exp2P;

  NonZero[pm][mm] = true;
  CopyConj(pm,mm, -1);


  //+- to +-
  JAMP[pm][pm][0] = -g2S/6. * ( 1 + Cos )*Exp1P;
  JAMP[pm][pm][1] = +g2S/2. * ( 1 + Cos )*Exp1P;

  NonZero[pm][pm] = true;
  CopyConj(pm,pm);


  //+- to -+
  JAMP[pm][mp][0] =  g2S/6. * ( 1 - Cos )*Exp1P;
  JAMP[pm][mp][1] = -g2S/2. * ( 1 - Cos )*Exp1P;

  NonZero[pm][mp] = true;
  CopyConj(pm,mp);



}


Double Sigma2qqbar2QQbar::PyCrossSection()
{
  Double s2 = s*s;

  Double m32 = mP3*mP3; 
  Double m42 = mP3*mP3; 

  // Modified Mandelstam variables for massive kinematics with m3 = m4.
  Double s34Avg = 0.5 * (m32 + m42) - 0.25 * (m32 - m42)*(m32 - m42) / s;
  Double tHQ    = -0.5 * (s - t + u);
  Double uHQ    = -0.5 * (s + t - u);
  Double tHQ2   = tHQ * tHQ;
  Double uHQ2   = uHQ * uHQ;

  // Calculate kinematics dependence.
  sigS = (4./9.) * ((tHQ2 + uHQ2) / s2 + 2. * s34Avg / s);

  // Answer.
  return   g2S*g2S * sigS ;

}


ColArr Sigma2qqbar2QQbar::PyColorFlow()
{
  ColArr cl;

  // Colour flow topologies. Swap when antiquarks.
  cl.setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (idLeft < 0) cl.swapColAcol();

  return cl;
}


void Sigma2qqbar2QQbar::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{

  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
	pArr[0] = 1;

  if (idLeft < 0) clArr[0].swapColAcol();

}
