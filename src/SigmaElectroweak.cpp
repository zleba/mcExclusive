#include "SigmaElectroweak.h"
#include "ExclusiveHooks.h"
#include <iomanip>

void Sigma2gg2gammagamma::ColorMatrix()
{


  Double A12, A13, A14;
  Double B234, B243, B342, B432;

  DataGG( A12, A13, A14,
          B234, B243, B342, B432);

  colMatrix[0][0] = 4*A12;

}


void Sigma2gg2gammagamma::Amplitudes()
{

  complex<Double> Exp = exp(-Double(2)*i_*phi);

  Double qf;

  Pythia8::ParticleData *pd = HEPinfo::ParticleData();
  Double mb = (pd != 0) ? pd->m0(4) : 4.3;

  if(s > 4.0*mb*mb)
     qf=11.0/9.0;
  else
     qf=10.0/9.0;


  Double norm   = qf*  alphaEm * g2S /M_PI/sqrt(Nsym); //is normalization OK?


//++ to ++
  JAMP[pp][pp][0] = -norm*( 0.5*(t*t+u*u)/s/s * (pow(log(t/u),2)+M_PI*M_PI) + (t-u)/s*log(t/u) + 1.0 );
  NonZero[pp][pp] = true;
  CopyConj(pp,pp);

//++ to +-
  JAMP[pp][pm][0] = norm;
  NonZero[pp][pm] = true;
  CopyConj(pp,pm);

//++ to -+
  JAMP[pp][mp][0] = norm;
  NonZero[pp][mp] = true;
  CopyConj(pp,mp);

//++ to --
  JAMP[pp][mm][0] = norm;
  NonZero[pp][mm] = true;
  CopyConj(pp,mm);




//+- to ++
  JAMP[pm][pp][0] = Exp*norm;
  NonZero[pm][pp] = true;
  CopyConj(pm,pp);

//+- to +-
  JAMP[pm][pm][0] = -norm*( 0.5*(t*t+s*s)/u/u*(pow(log(-t/s),2) + 2.*i_*M_PI*log(-t/s))+(t-s)/u*(log(-t/s)+i_*M_PI)+1.0 )*Exp;
  NonZero[pm][pm] = true;
  CopyConj(pm,pm);

//+- to -+
  JAMP[pm][mp][0] = -norm*(0.5*(u*u+s*s)/t/t*(pow(log(-s/u),2) + 2.0*i_*M_PI*log(-s/u))+(s-u)/t*(log(-s/u)+i_*M_PI)+1.0 )*Exp;
  NonZero[pm][mp] = true;
  CopyConj(pm,mp);

//+- to --
  JAMP[pm][mm][0] = Exp*norm;
  NonZero[pm][mm] = true;
  CopyConj(pm,mm);

}




Double Sigma2gg2gammagamma::PyCrossSection()
{

  int nQuarkLoop = 5;
  // Calculate charge factor from the allowed quarks in the box.
  Double charge2Sum     = 1./9. + 4./9. + 1./9.;

  if (nQuarkLoop >= 4) charge2Sum += 4./9.;
  if (nQuarkLoop >= 5) charge2Sum += 1./9.;
  if (nQuarkLoop >= 6) charge2Sum += 4./9.;



  using Pythia8::pow2;

  // Logarithms of Mandelstam variable ratios.
  Double logST = log( -s / t );
  Double logSU = log( -s / u );
  Double logTU = log(  t / u );

  // Real and imaginary parts of separate amplitudes.
  Double b0stuRe = 1. + (t - u) / s * logTU
    + 0.5 * (t*t + u*u) / (s*s) * (pow2(logTU) + pow2(M_PI));
  Double b0stuIm = 0.;
  Double b0tsuRe = 1. + (s - u) / t * logSU
    + 0.5 * (s*s + u*u) / (t*t) * pow2(logSU);
  Double b0tsuIm = -M_PI * ( (s - u) / t + (s*s + u*u) / (t*t) * logSU);
  Double b0utsRe = 1. + (s - t) / u * logST
    + 0.5 * (s*s + t*t) / (u*u) * pow2(logST);
  Double b0utsIm = -M_PI * ( (s - t) / u + (s*s + t*t) / (u*u) * logST);
  Double b1stuRe = -1.;
  Double b1stuIm = 0.;
  Double b2stuRe = -1.;
  Double b2stuIm = 0.;

  // Calculate kinematics dependence.
  Double sigBox = pow2(b0stuRe) + pow2(b0stuIm) + pow2(b0tsuRe)
    + pow2(b0tsuIm) + pow2(b0utsRe) + pow2(b0utsIm) + 4. * pow2(b1stuRe)
    + 4. * pow2(b1stuIm) + pow2(b2stuRe) + pow2(b2stuIm);


  // Answer contains factor 1/2 from identical photons.
  Double sigma = 0.5*1./M_PI* (1. / (16. * M_PI )) * pow2(charge2Sum)
    * g2S*g2S * pow2(alphaEm) * sigBox;


  //Formula from paper
  /*
  Double res = 0.5*pow2(alphaEm)*pow2(charge2Sum) * pow2(g2S/4./M_PI) * 1./(8.*M_PI*s*s);


  Double TermUU = (s*s + t*t)/u/u;
  Double TermTT = (s*s + u*u)/t/t;
  Double TermSS = (t*t + u*u)/s/s;

  Double TermU = (s-t)/u;
  Double TermT = (s-u)/t;
  Double TermS = (t-u)/s;

  Double LogST = log(-s/t);
  Double LogSU = log(-s/u);
  Double LogTU = log( t/u);

  cout << "Before " <<setprecision(15)<< res << endl;
  res *=
  1./8 *( pow2( TermUU*pow2(LogST) + 2*TermU*LogST ) +
          pow2( TermTT*pow2(LogSU) + 2*TermT*LogSU ) +
          pow2( TermSS*(pow2(LogTU)+M_PI*M_PI) +2*TermS*LogTU) )

 +1./2 *(TermUU*LogST + 2*TermU*LogST + TermTT*pow2(LogSU)
       + 2* TermT*LogSU + TermSS*(pow2(LogTU) + M_PI*M_PI) - 2*TermS*LogTU )

 +M_PI*M_PI/2. *( pow2(TermUU*LogST + TermU) + pow2(TermTT*LogSU + TermT) ) + 4; 


  cout <<"HOPE "<< res/( sigma/(16*M_PI*s*s) ) << endl; 
  */

  return sigma;
}

ColArr Sigma2gg2gammagamma::PyColorFlow()
{

  ColArr cl;
  cl.setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

  return cl;

}


void Sigma2gg2gammagamma::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);
	pArr[0] = 1;

}













void Sigma2qqbar2gammagamma::ColorMatrix()
{

  Double Dl1l1c, Dl1r1c;

  Double Dl1r1; 

  DataQQ( Dl1l1c, Dl1r1c, Dl1r1);

  colMatrix[0][0] =  Dl1r1;

}


void Sigma2qqbar2gammagamma::Amplitudes()
{

  //Double e2 = 4*M_PI * alphaEm;

  Double cGamma;

  if(abs(idLeft) % 2 == 1)
    cGamma = -e2* 1./3 *1./3;
  else
    cGamma = -e2* 2./3 *2./3;


  complex<double> norm = exp(i_*phi)*cGamma/sqrt(Nsym);


  //+- to +-
  JAMP[pm][pm][0] = 2.*norm* sqrt(u/t);
  NonZero[pm][pm] = true;
  CopyConj(pm,pm);


  //+- to -+
  JAMP[pm][mp][0] =-2.*norm* sqrt(t/u);
  NonZero[pm][mp] = true;
  CopyConj(pm,mp);


}

Double Sigma2qqbar2gammagamma::PyCrossSection() {

  // Calculate kinematics dependence.
  sigTU  = 2. * (t*t + u*u) / (t * u);

  Double charge;
  //d, s, b
  if(abs(idLeft) % 2 == 1) 
    charge = -1./3;
  //u, c, t
  else
    charge =  2./3;

  // Answer contains factor 1/2 from identical photons.
  Double sigma0 =  e2*e2 *pow(charge,4)  * 0.5 * sigTU;
  
  //From color averaging (only 3 of 9 color combinations possible)
  return sigma0 * 1./3;
}


ColArr Sigma2qqbar2gammagamma::PyColorFlow()
{
  ColArr cl;


  // No colours at all or one flow topology. Swap if first is antiquark.
  cl.setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);

  if (idLeft < 0) cl.swapColAcol();

  return cl;
}


void Sigma2qqbar2gammagamma::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  if (idLeft < 0) clArr[0].swapColAcol();
	pArr[0] = 1;

}













void Sigma2ffbar2ffbarsgmZ::ColorMatrix()
{


  Double Dl1l1c, Dl1r1c;

  Double Dl1r1; 

  DataQQ( Dl1l1c, Dl1r1c, Dl1r1);

  colMatrix[0][0] =  Dl1r1;

}


void Sigma2ffbar2ffbarsgmZ::Amplitudes()
{

  using Pythia8::pow2;

  //MDL_MW = sqrt(MZ*MZ/2.+sqrt(MZ*MZ*MZ*MZ /4. -(MDL_AEW*M_PI*MZ*MZ)/(GF*sqrt(2))))

  //Double alphaEl = 1/132.507;


  Pythia8::ParticleData *pd = HEPinfo::ParticleData();
  Pythia8::CoupSM *co = HEPinfo::Couplings();

  Double MW, MZ, WZ, CW, SW;

  if(pd != 0 && co != 0) {
		//cout << "I am here" << endl;
    MW = pd->m0(24);
    MZ = pd->m0(23);
    WZ = pd->mWidth(23);

    SW = sqrt( co->sin2thetaW() );
    CW = sqrt( co->cos2thetaW() );
    

  } else {
    MW = 80.419002445756163;
    MZ = 91.18800;
    WZ = 2.441404;

    CW = MW/MZ;
    SW = sqrt(1-CW*CW );
  }

  //Double e2 = 4*M_PI * alphaEm;
  Double EE = sqrt(e2);


  complex<double> norm = exp(i_*phi);


  Double cGamma;


  complex<Double> GC50 = -( EE* CW*i_) /(2.*SW);
  complex<Double> GC51 =  ( EE* CW*i_) /(2.*SW);
  complex<Double> GC58 = -( EE* SW*i_) /(6.*CW);
  complex<Double> GC59 =  ( EE* SW*i_) /(2.*CW);

  //cout <<"GC50 - my " <<  GC50 << endl;
  //cout <<"GC51 - my " <<  GC51 << endl;
  //cout <<"GC58 - my " <<  GC58 << endl;
  //cout <<"GC59 - my " <<  GC59 << endl;

  complex<Double> cZpmpm, cZpmmp, cZmppm, cZmpmp;

  //d, s, b
  if(abs(idLeft) % 2 == 1) {
    cZpmpm = -(GC58*GC59)/2.;
    cZpmmp = -((GC58*GC59) + GC50*GC58)/4.;
    cZmppm = +((GC58*GC59) + GC50*GC59)/4.;
    cZmpmp = (GC58*GC59 + GC50*GC59 + GC50*GC58 + GC50*GC50)/8.;



    cGamma = -e2* 1./3;
  }
  //u c, t
  else {
    cZpmpm = +(GC58*GC59);
    cZpmmp = +((GC58*GC59) + GC50*GC58)/2.;
    cZmppm = +((GC58*GC59) + GC51*GC59)/4.;
    cZmpmp = +(GC58*GC59 + GC51*GC59+ GC50*GC58 +GC50 *GC51)/8.;
    cGamma = e2* 2./3;
  }


	//Double A = pow2(abs(cZpmpm)) + pow2(abs(cZmpmp));
	//Double B = pow2(abs(cZpmmp)) + pow2(abs(cZmppm));
	//cout << "RATIO(MY) " << (A+B)/(A-B) << endl;


	Double resM = s/MZ;
	//Double resM = MZ;

  //+- to +-
  JAMP[pm][pm][0] = norm*u*( 2.*cGamma/s  +  cZpmpm* 16./( - MZ*MZ + s + i_*resM*WZ) );
  NonZero[pm][pm] = true;


  //+- to -+
  JAMP[pm][mp][0] =-norm*t*( 2.*cGamma/s +  cZpmmp* 16./( - MZ*MZ + s + i_*resM*WZ) );
  NonZero[pm][mp] = true;


  //-+ to +-
  JAMP[mp][pm][0] =-conj(norm)*t*( 2.*cGamma/s  +  cZmppm* 16./( - MZ*MZ + s + i_*resM*WZ) );
  NonZero[mp][pm] = true;


  //-+ to -+
  JAMP[mp][mp][0] = conj(norm)*u*( 2.*cGamma/s  +  cZmpmp* 16./( - MZ*MZ + s + i_*resM*WZ) );
  NonZero[mp][mp] = true;



//   - 2*exp(i_*phi)*GC2*s^-1*t*GC3 - 8/( - MDLMZ^2 + s + i_*MDLMZ*MDLWZ)*
 //        exp(i_*phi)*t*GC58*GC59 - 8/( - MDLMZ^2 + s + i_*MDLMZ*MDLWZ)*exp(i_*phi
  //             )*t*GC50*GC58;


}


Double Sigma2ffbar2ffbarsgmZ::PyCrossSection()
{
  using Pythia8::pow2;
  
  int abId = abs(idOut1);
  if(abId != 11 && abId != 13 && abId != 15) {
    cout<<"Decay channel to pdgID "<<abId<<" not implemented."<<endl;
    exit(1);
  }


  Pythia8::ParticleData *pd = HEPinfo::ParticleData();
  Pythia8::CoupSM *co = HEPinfo::Couplings();


  const Double MW = 80.419002445756163;
  const Double MZ = 91.18800;
  const Double WZ = 2.441404;

  const Double CW = MW/MZ;
  const Double SW = sqrt(1-CW*CW );



  //initialisation
  static Double mRes      = (pd != 0) ? pd->m0(23) : MZ;
  static Double GammaRes  = (pd != 0) ? pd->mWidth(23) : WZ;
  static Double m2Res     = mRes*mRes;
  static Double GamMRat   = GammaRes / mRes;
  static Double thetaWRat = (pd != 0) ? 1./(16. * co->sin2thetaW() * co->cos2thetaW()) : 1./(16.*CW*CW*SW*SW);


  Double colQ = 3.;// * (1. +  alphaS / M_PI);

	////////////////////////////////////
  //outgoing part init
	////////////////////////////////////
 
        int idAbs = abs(idOut1);
        Double mf = (pd != 0) ? pd->m0(idAbs) : 0.;
				mf = 0;
        Double mH = sqrt(s);

        Double mr    = pow2(mf / mH);
        Double betaf = Pythia8::sqrtpos(1. - 4. * mr);

        // Combine couplings (including colour) with phase space.
				Double ef, vf, af;
				if( co != 0) {
					ef    = co->ef(idAbs);
					vf    = co->vf(idAbs);
					af    = co->af(idAbs);
					//cout << "vf " << vf<<" "<< (-0.5+2*co->sin2thetaW())*2   << endl;
				}
				else {
					ef    = -1;
					vf    = (-0.5+2*SW*SW)*2;
					af    = -1;
				}

        Double colf  = (idAbs < 6) ? colQ : 1.;
        Double gamTf = colf * ef * ef * betaf;
        Double gamLf = gamTf * 4. * mr;
        Double intTf = colf * ef * vf * betaf;
        Double intLf = intTf * 4. * mr;
        Double intAf = colf * ef * af * betaf;
        Double resTf = colf * (vf * vf * betaf + af * af * Pythia8::pow3(betaf));
        Double resLf = colf * vf * vf * betaf * 4. * mr;
        Double resAf = colf * vf * af * betaf * 4.;

				//vf = (-0.5+2*co->sin2thetaW())*2


  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  Double gamProp = M_PI * pow2(alphaEm) / s/s;
  Double intProp = gamProp * 2. * thetaWRat * s * (s - m2Res)
          / ( pow2(s - m2Res) + /*pow2(mRes*GammaRes) */ pow2(s * GamMRat)   );
  Double resProp = gamProp * pow2(thetaWRat * s)
          / ( pow2(s - m2Res) + /*pow2(mRes*GammaRes) */ pow2(s * GamMRat)   );

	//intProp = gamProp = 0; //only gamma part



  // Scattering angle in subsystem rest frame.
  Double cThe = (t - u) / s;


	////////////////////////////////////
	//ingoing part init
	////////////////////////////////////
	Double coefT, coefL, coefA;
	int id1Abs;
	{
		// Couplings for current in-flavour.
		id1Abs = abs(idLeft);
		Double ei, vi, ai;
		if(co != 0) {
			ei  = co->ef(id1Abs);
			vi  = co->vf(id1Abs);
			ai  = co->af(id1Abs);
			//cout << "ei ai "<< ei <<" "<< ai << endl;
		}
		else {

			//d, s, b
			if(abs(idLeft) % 2 == 1) {
				ei = -1./3;
				ai = -1;
				vi = (-0.5+2./3*SW*SW)*2;
			}
			//u, c, t
			else {
				ei =  2./3;
				ai = +1;
			  vi = (+0.5-4./3*SW*SW)*2;
			}

		//vi = (-0.5+2./3*co->sin2thetaW())*2; //d,s,b
		//vi = (+0.5-4./3*co->sin2thetaW())*2; //u,c,t

		}

		Double gamSumT = gamTf;
		Double gamSumL = gamLf;
		Double intSumT = intTf;
		Double intSumL = intLf;
		Double intSumA = intAf;
		Double resSumT = resTf;
		Double resSumL = resLf;
		Double resSumA = resAf;




		// Coefficients of angular expression.
		coefT = ei*ei * gamProp * gamSumT + ei*vi * intProp * intSumT
								 + (vi*vi + ai*ai) * resProp * resSumT;
		coefL = ei*ei * gamProp * gamSumL + ei*vi * intProp * intSumL
								 + (vi*vi + ai*ai) * resProp * resSumL;
		coefA = ei*ai * intProp * intSumA + vi*ai * resProp * resSumA;

	}

  // Colour factor. Answer.
  Double sigma = coefT * (1. + pow2(cThe)) + coefL * (1. - pow2(cThe))
               + 2. * coefA * cThe;


  if (id1Abs < 9) sigma /= 3.;

	
	// to ME
	sigma *= s*s / M_PI * 16*M_PI*M_PI;



  return sigma;

}

ColArr Sigma2ffbar2ffbarsgmZ::PyColorFlow()
{
  ColArr cl;


  // No colours at all or one flow topology. Swap if first is antiquark.
  cl.setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);

  if (idLeft < 0) cl.swapColAcol();

  return cl;

}



void Sigma2ffbar2ffbarsgmZ::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  if (idLeft < 0) clArr[0].swapColAcol();
	pArr[0] = 1;

}



void Sigma1gg2H::ColorMatrix()
{


  Double A12, A13, A14;
  Double B234, B243, B342, B432;

  DataGG( A12, A13, A14,
          B234, B243, B342, B432);

  colMatrix[0][0] = 4*A12;

}


void Sigma1gg2H::Amplitudes()
{

  using Pythia8::pow2;

  //Load resonance data
  Pythia8::ParticleData *pd = HEPinfo::ParticleData();
  Double mH, m2Res;
  Double widthIn, widthOut, width, BR;
  if(pd != 0) {
    int idRes = 25;//Higgs PDG
    Pythia8::ParticleDataEntry* HResPtr=pd->particleDataEntryPtr(idRes);

    mH=  HResPtr->m0();
    widthIn  = HResPtr->resWidthChan( mH, 21, 21) / 64.;

    width    = HResPtr->resWidth(idRes, mH);
    BR = HResPtr->resOpenFrac(idRes);
  }
  else {
    mH= 125;
    widthIn  = 5.426e-6;
    widthOut = 4.08e-3;
    BR = 0.5768872994; //for b bbar
  }
  m2Res = mH*mH;

  widthOut = width * BR;

  Double fact = widthIn * widthOut;
  
  complex<Double> sigma = sqrt(8)*sqrt(8*M_PI)/ ( s - m2Res + mH * width*i_ );


  sigma *= sqrt(fact);
	// to ME
	sigma *= s *sqrt(16*M_PI);/// M_PI * 16*M_PI*M_PI;


	JAMP[pp][pp][0] = sigma;
	JAMP[pp][mm][0] = sigma;
  NonZero[pp][pp] = true;
  NonZero[pp][mm] = true;

	JAMP[mm][pp][0] = sigma;
	JAMP[mm][mm][0] = sigma;
  NonZero[mm][pp] = true;
  NonZero[mm][mm] = true;

}


Double Sigma1gg2H::PyCrossSection()
{

  using Pythia8::pow2;

  if(abs(idOut1) != 5) {
    cout<<"Decay channel to pdgID "<<idOut1<<" not implemented."<<endl;
    exit(1);
  }

  Pythia8::ParticleData *pd = HEPinfo::ParticleData();
  Double mH, m2Res;
  Double widthIn, widthOut, width, BR;
  if(pd != 0) {
    int idRes = 25;//Higgs PDG
    Pythia8::ParticleDataEntry* HResPtr=pd->particleDataEntryPtr(idRes);

    mH=  HResPtr->m0();
    widthIn  = HResPtr->resWidthChan( mH, 21, 21) / 64.;

    width    = HResPtr->resWidth(idRes, mH);
    BR = HResPtr->resOpenFrac(idRes);
  }
  else {
    mH= 125;
    widthIn  = 5.426e-6;
    widthOut = 4.08e-3;
    BR = 0.5768872994; //for b bbar
  }
  m2Res = mH*mH;

  widthOut = width * BR;
  // Set up Breit-Wigner.


  Double sigBW    = 8. * M_PI/ ( pow2(s - m2Res) + pow2(mH * width) );

  // Done.
  Double sigma = widthIn * sigBW * widthOut;
  //cout <<"RADEK wIn,BW,wOut "<< widthIn<<" "<<sigBW<<" "<<widthOut<<endl;

  //cout << "sigma(RADEK)   " <<setprecision(15)<< sigma << endl;
  //cout << "radek(distance) "<< (s - m2Res) / (mH * width) <<" "<<s<<" "<<m2Res<< endl;

	// to ME
	sigma *= s*s * 16*M_PI;

	return sigma;

}






ColArr Sigma1gg2H::PyColorFlow()
{
  ColArr cl;


  // No colours at all or one flow topology. Swap if first is antiquark.
  cl.setColAcol( 1, 2, 2, 1, 3, 4, 4, 3);

  return cl;

}


void Sigma1gg2H::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 2, 2, 1, 3, 4, 4, 3);
	pArr[0] = 1;

}




void Sigma1ffbar2H::ColorMatrix()
{

  Double Dl1l1c, Dl1r1c;

  Double Dl1r1; 

  DataQQ( Dl1l1c, Dl1r1c, Dl1r1);

  colMatrix[0][0] =  Dl1r1;
	//cout << "RADECEK "<<  Dl1l1c <<" "<< Dl1r1c << " "<< Dl1r1 << endl;

}


void Sigma1ffbar2H::Amplitudes()
{

  using Pythia8::pow2;
	complex<Double> sigma;

  //Load resonance data
  Pythia8::ParticleData *pd = HEPinfo::ParticleData();
	//pd = 0;
  Double mH, m2Res;
  Double widthIn, widthOut, width, BR;
  if(pd != 0) {
    int idRes = 25;//Higgs PDG
    Pythia8::ParticleDataEntry* HResPtr=pd->particleDataEntryPtr(idRes);

    mH=  HResPtr->m0();
    widthIn  = HResPtr->resWidthChan( mH, abs(idLeft), -abs(idLeft)) / 9.;

    width    = HResPtr->resWidth(idRes, mH);
    BR = HResPtr->resOpenFrac(idRes);

		widthOut = width * BR;

  }
  else {
		cout << "I am here KOCKA" << endl;

    mH= 125;
    widthIn  = abs(idLeft)==5 ? 0.0002617450   :  0.000013089348662 ;
    widthOut = 0.0023557954; //to bBbar
		width    = 0.0040827883; //total



    //widthIn  = 0.000261741997279;
    //widthOut = 4.08e-3;
    //BR = 0.5768872994; //for b bbar
	  //width =  6.382339E-003;
  }
  m2Res = mH*mH;


  Double fact = widthIn * widthOut;
  
  complex<Double> sigBWsqrt = sqrt(4*M_PI)/ ( s - m2Res + mH * width*i_ );

	complex<Double> sigmaSqrt = sqrt(fact) * sigBWsqrt;

	complex<Double> MEtot = sigmaSqrt * sqrt(s*16*M_PI); //from sigma to ME sqrt
	//sigma *= s * 16*M_PI;

	complex<Double> MEfact = MEtot * sqrt( 9.*4. / (4.*9.) ); //Colors*Helicity / (#same elements)

	/*
  sigma *= sqrt(fact);
  sigma *= 0.003673190430791598 * (1+2.55331e-06)*(1-0.00214653)*(1+3.91122e-09);


	// to ME
	sigma *= s *sqrt(16*M_PI);/// M_PI * 16*M_PI*M_PI;

	width    = 0.0040827883; //pythia value
  sigma = 2*(s - 4*mP3*mP3) / ( s - mH*mH + mH * width*i_ );
  //sigma *= 0.000182195 * (1+-1.85414e-07);

  Double g = 0.65323293034757990;
  Double mW= 80.419002445756163;
  sigma *= 0.5*pow(g/2. *  mP3/mW,2);  
	*/

  //sigma *= (1+2.7328e-05)*(1+4.45601e-11);

  complex<Double> Exp = exp(-Double(1)*i_*phi);

	cout << "CMS energy " <<" "<< sqrt(s) << endl;

	JAMP[pp][pp][0] = +MEfact * Exp;
  NonZero[pp][pp] = true;

	JAMP[pp][mm][0] = -MEfact*conj(Exp);
  NonZero[pp][mm] = true;

	JAMP[mm][pp][0] = -MEfact* Exp;
  NonZero[mm][pp] = true;

	JAMP[mm][mm][0] = +MEfact* conj(Exp);
  NonZero[mm][mm] = true;

}


Double Sigma1ffbar2H::PyCrossSection()
{

  using Pythia8::pow2;

  if(abs(idOut1) != 5) {
    cout<<"Decay channel to pdgID "<<idOut1<<" not implemented."<<endl;
    exit(1);
  }

  Pythia8::ParticleData *pd = HEPinfo::ParticleData();
	//pd = 0;
  Double mH, m2Res;
  Double widthIn, widthOut, width;
  if(pd != 0) {
	  cout << "I am here " << endl;
    int idRes = 25;//Higgs PDG
    Pythia8::ParticleDataEntry* HResPtr=pd->particleDataEntryPtr(idRes);

    mH=  HResPtr->m0();
    widthIn  = HResPtr->resWidthChan( mH, abs(idLeft), -abs(idLeft) ) / 9.;//  64.;

    width    = HResPtr->resWidth(idRes, mH);
    Double BR = HResPtr->resOpenFrac(idRes);
		widthOut = width * BR;
  }
  else {
    mH= 125;
    widthIn  = abs(idLeft)==5 ? 0.0002617450   :  0.000013089348662 ;
    widthOut = 0.0023557954; //to bBbar
		width    = 0.0040827883; //total
		//widthIn *= 3;

  }
  m2Res = mH*mH;

	//cout << "RADEK width " << width << endl;
  // Set up Breit-Wigner.


  Double sigBW    = 4. * M_PI/ ( pow2(s - m2Res) + pow2(mH * width) );

  // Done.
  Double sigma = widthIn * sigBW * widthOut;
	cout << idLeft << endl;
  //cout <<"RADEK wIn,BW,wOut "<<setprecision(15)<< widthIn<<" "<<sigBW<<" "<<widthOut<<endl;
  //cout << "sigma(RADEK)   " <<setprecision(15)<< sigma << endl;


  //cout << "radek(distance) "<< (s - m2Res) / (mH * width) <<" "<<s<<" "<<m2Res<< endl;

	// to ME
	sigma *= s * 16*M_PI;
	//sigma *= s*s * 16*M_PI;

	return sigma;

}


ColArr Sigma1ffbar2H::PyColorFlow()
{
  ColArr cl;


  // No colours at all or one flow topology. Swap if first is antiquark.
  cl.setColAcol( 1, 0, 0, 1, 2, 0, 0, 2 );
  if (idLeft < 0) cl.swapColAcol();

  return cl;

}


void Sigma1ffbar2H::PyColorFlows(vector<ColArr> &clArr, vector<Double> &pArr)
{
  clArr.resize(nColFlow);
  pArr.resize(nColFlow);

	clArr[0].setColAcol( 1, 0, 0, 1, 2, 0, 0, 2 );
	pArr[0] = 1;
	if (idLeft < 0) clArr[0].swapColAcol();

}

