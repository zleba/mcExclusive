#include "ExclusiveHooks.h"
#include "SigmaQCD.h"
#include "SigmaElectroweak.h"
#include "ColorReconection.h"
#include "TopoSelector.h"

#include <set>

Pythia8::ParticleData* HEPinfo::particleDataPtr = 0;
Pythia8::CoupSM* HEPinfo::couplingsPtr = 0;

void ExclusiveHooks::ModifyParameters(int procId, int pRest, int _p0, int _p1, int _p2)
{
	if(procs.count(procId) == 0) {
		cout << "Process "<<procId <<" does not exist!" << endl;
		exit(1);
	}
  procs[procId]->topo.SetPs(pRest, _p0, _p1, _p2 );

}

//In case of faild hadronisation of given parton state
void ExclusiveHooks::ResetShowers(const Pythia8::Event& event)
{
	int nEm = 0;
	for(int i = 0; i < event.size(); ++i)
		if( abs(event[i].status()) == 43) ++nEm;

	if(nEm != 1) return;
	if(proc->topo.nEm == 0) return;

	if(proc->topo.nEm != 0) {
		cout << "I am reseting event " <<pyth->info.getCounter(3)<<" "<<weightEx<< " "<<isExclusive<< endl;
		if(isExclusive && weightEx != 1)
			proc->topo.RemoveEvent(weightEx);
	}

  proc->topo.nEm = 0;

  proc->RemoveShowers();

  lastEmScale = pyth->info.QRen();
  NormalSplittings = 1;

  weightEx = 1;
  nMPI = 0;

  if( pyth->info.QRen() > ExclusiveCut) {
    isExclusive = true;
    isSmallPt   = false;
  }
  else {
    isExclusive = false;
    isSmallPt   = true;
  }

  if( proc->topo.Type() == 0)
    isSmallPt = true;

}

bool ExclusiveHooks::doVetoISREmission( int sizeOld, const Pythia8::Event& event, int )
{
	//reset shower information if needed
	ResetShowers(event);

	if(isSmallPt)
		return true;


	//cout << "counter " << pyth->info.getCounter(3) << endl;
	if( isExclusive && weightEx != 1 ) {
    cout <<"isExclusive,weight "<< isExclusive <<" "<< weightEx << endl;
    cout <<"topoType "<< proc->topo.Type() <<" "<<isSmallPt << endl; 
		cout <<"Event count " << pyth->info.getCounter(3) << endl;
		event.list();
  }



	//event.list();	
  //cout << "event+size "<<pyth->info.getCounter(3)<<" "<<event.size() << endl;
	//if(4055 == pyth->info.getCounter(3)) {
		//cout <<"Event count " << pyth->info.getCounter(3) << endl;
		//event.list();
	//}

  if(proc->topo.Type() >= 0 && proc->topo.nEm >= proc->topo.Type() ) {
    isSmallPt = true;
    return true;
  }
	//cout << "OldSize, iSys " << sizeOld <<" "<< iSys <<  endl;

	//event.list();

	//Find branch
	int idOrg = -1;
	if( event[sizeOld].status() ==-41 ) {
		idOrg  = sizeOld;
	}
	else if( event[sizeOld+1].status() ==-41 ) {
		idOrg  = sizeOld+1;
	}
	else {
		cout << "Something wrong " << endl;
		exit(1);
	}

	int idEm   = event[idOrg].daughter1();
	int idDaug = event[idOrg].daughter2();



	if( event[idEm].scale() < ExclusiveCut ) {
		lastScale = event[idEm].scale();
		isSmallPt = true;
		return true;
	}
  else 
    lastEmScale = event[idEm].scale();



	//Calculate z and phi of branching
	Double z   = event[idDaug].pz() / event[idOrg].pz();
	Double phi = event[idEm].phi();

	if(z < 0 || z > 1) {
		cout << "Strange z : "<< z << endl;
		exit(1);
	}

	//Add emission
	if( event[idOrg].isAncestor(1) )
		proc->AddEmissionL(event[idOrg].id(), event[idDaug].id(), z, phi);
	else if( event[idOrg].isAncestor(2) )
		proc->AddEmissionR(event[idOrg].id(), event[idDaug].id(), z, phi);
	else {
		cout << "Something wrong " << endl;
		exit(1);
	}

	NormalSplittings *= SpinTensor::ClassicalSplitting( event[idOrg].id(), event[idDaug].id(), z );

  ++proc->topo.nEm;

	return false;
}


ExclusiveHooks::ExclusiveHooks(Pythia8::Pythia *py, Double _bSlope) : bSlope(_bSlope)
{
  pyth = py;

  pyth->settings.forceParm("SpaceShower:pT0Ref", 0.0);//change


  pyth->readString("SpaceShower:alphaSorder = 1");
  pyth->readString("SpaceShower:alphaSuseCMW = false");



  pyth->readString("SpaceShower:rapidityOrder = off");
  pyth->readString("SpaceShower:phiPolAsym = off");
  pyth->readString("SpaceShower:MEcorrections = off");
  pyth->readString("SpaceShower:phiIntAsym = off");


  pyth->readString("SpaceShower:weakShower = off");
  pyth->readString("SpaceShower:QEDshowerByQ = off");
  pyth->readString("SpaceShower:QEDshowerByL = off");


  pyth->readString("TimeShower:interleave = off");
  pyth->readString("TimeShower:allowBeamRecoil = off");

  pyth->readString("BeamRemnants:primordialKT = off");


	pyth->readString("StandardModel:sin2thetaWbar = 0.2312");


	ExclusiveCut = 5;

	WithoutEmissions = true; //default value


	procs[111] = new Sigma2gg2gg;
	procs[112] = new Sigma2gg2qqbar;
	procs[113] = new Sigma2qg2qg;
	procs[114] = new Sigma2qq2qq;
	procs[115] = new Sigma2qqbar2gg;
	procs[116] = new Sigma2qqbar2qqbarNew;

	//charm production
	procs[121] = new Sigma2gg2QQbar;
	procs[122] = new Sigma2qqbar2QQbar;

	//bottom production
	procs[123] = new Sigma2gg2QQbar;
	procs[124] = new Sigma2qqbar2QQbar;

	//top production
	procs[601] = new Sigma2gg2QQbar;
	procs[602] = new Sigma2qqbar2QQbar;


	procs[204] = new Sigma2qqbar2gammagamma;
	procs[205] = new Sigma2gg2gammagamma;
	procs[224] = new Sigma2ffbar2ffbarsgmZ;
  
  //Higgs to bb
	procs[901] = new Sigma1ffbar2H;
	procs[902] = new Sigma1gg2H;


	for(map<int,Sigma2Process *>::iterator  it = procs.begin();  it != procs.end(); ++it) {
		it->second->SetColors(&colL,  &colR );
		it->second->SetSpins (&spinL, &spinR);
	}


}

bool ExclusiveHooks::initAfterBeams()
{
    Double alphaSMZ = settingsPtr->parm("SpaceShower:alphaSvalue");

    int alphaSorder = settingsPtr->mode("SpaceShower:alphaSorder ");

    bool isISR = settingsPtr->flag("PartonLevel:ISR");
    //cout << "ISR is " << isISR << endl;
    WithoutEmissions = !isISR;

    //exit(0);

    if(alphaSorder != 1) {
        cout << "Please set order of alphaS to 1" << endl;
        exit(1);
    }

    const Double MCMIN  = 1.2;
    const Double MBMIN  = 4.0;
    //const Double MTMIN  = 167.0;
    Double  mc = max( MCMIN, particleDataPtr->m0(4));
    Double  mb = max( MBMIN, particleDataPtr->m0(5));
    //Double  mt = max( MTMIN, particleDataPtr->m0(6));

    HEPinfo::Init(particleDataPtr, coupSMPtr );


    Double sqrtS;
    int BeamMode = settingsPtr->mode("Beams:frameType");
    if(BeamMode == 1) {
        sqrtS = settingsPtr->parm("Beams:eCM");
    }
    else if(BeamMode == 2) {
        sqrtS = 2*sqrt(settingsPtr->parm("Beams:eA") * settingsPtr->parm("Beams:eB") );
    }
    ExclusiveCut = settingsPtr->parm("SpaceShower:pTmin");


    kmr = new KMRlumi(sqrtS, 0.8, alphaSMZ, mc, mb, bSlope, beamAPtr, rndmPtr);

    return true;
}



bool ExclusiveHooks::doVetoProcessLevel(Pythia8::Event& process)
{

  int code = pyth->info.code();

  if ( procs.find(code) == procs.end() ) {
    cout << "Process "<< pyth->info.name()<< " ("<<  code << ") is not implemented. " << endl;
    exit(1);
  } else {
    proc = procs[code];
  }

  //cout << "RADEK output begin" << endl;
  //process.list();
  //cout << "RADEK output end" << endl;

  int idHard1 = process.size()-2;
  int idHard2 = process.size()-1;

  if(code == 601 || code == 602) { //top production
    for(idHard1  = 1; idHard1 < process.size(); ++idHard1)
      if(process[idHard1].status() == -22) break;
    for(idHard2  = idHard1+1; idHard2 < process.size(); ++idHard2)
      if(process[idHard2].status() == -22) break;
    if(idHard1 >= process.size() || idHard2 >= process.size()) {
      cout << "Something wrong with partons identification" << endl;
      exit(1);
    }
  }
  //cout << "idHARD12 " << idHard1 << " " << idHard2 << endl;


  /*
  if( process[idHard1].status() != 23 || process[idHard2].status() != 23) {
    process.list();
    cout << "The hard process badly identified" << endl;
    exit(1);
  }
  */
  
	//process.list();

  ColRec.hardOrg.setColAcol( process[3].col(),process[3].acol(), process[4].col(),process[4].acol(),
                      process[idHard1].col(),process[idHard1].acol(), process[idHard2].col(),process[idHard2].acol() );


  proc->SetIds( pyth->info.id1(), pyth->info.id2() );
  proc->SetIdsOut(process[idHard1].id(), process[idHard2].id() );

  proc->SetSTUPhi(pyth->info.sHat(), pyth->info.tHat(), pyth->info.uHat(), pyth->info.phiHat() );
  proc->SetMasses( 0, 0, pyth->info.m3Hat(), pyth->info.m4Hat() );
  proc->SetScaleAlphas(pyth->info.QRen(), pyth->info.alphaS(), pyth->info.alphaEM() );


  lastEmScale = pyth->info.QRen();
  NormalSplittings = 1;

	weightEx = 1;
	nMPI = 0;

  if( pyth->info.QRen() > ExclusiveCut) {
    isExclusive = true;
    isSmallPt   = false;
  }
  else {
    isExclusive = false;
    isSmallPt   = true;
  }


  Double prob = kmr->NoEmProbability(pyth->info.x1(), pyth->info.x2(), 2*log(ExclusiveCut), 2*log(lastEmScale), true );
  //TopoSel.SetParam(-log(prob), 0, 0.5, 15 );
  //TopoSel.SetParam(-log(prob), 1);

	if(WithoutEmissions)
		proc->topo.SetParam(-log(prob), 0, 1, 0, 0);
	else
		proc->topo.SetA(-log(prob));

  //TopoSel.SetParam(-log(prob), 0.4, 16, 1, 2 );
  //TopoSel.SetParam(-log(prob), 1.000e+00,7.449e+01,3.365e+00,5.753e+00);

  //cout <<"TopoType:a  " <<  TopoSel.Type() <<" "<< -log(prob) <<  endl;
  if(proc->topo.Type() == 0)
    isSmallPt = true;


  return false;

}




void ExclusiveHooks::modifyEventBeforeFSR( Pythia8::Event &event )
{

  proc->PyCrossSection();

  int idStart1, idStart2;
  idStart1 = event[1].daughter1();
  idStart2 = event[2].daughter1();


  //if incoming partons can't be color singlet
	if( !(  
           ( event[idStart1].id() == 21 && event[idStart2].id() == 21 ) //both gluons
        //|| ( abs(event[idStart1].id()) <= 6 && event[idStart1].id() == -event[idStart2].id() )  //q + qbar
	) ) {
    isExclusive = false;
    return;
	}	


  Double x1, x2;
  x1 = event[idStart1].pz() / infoPtr->pzA();
  x2 = event[idStart2].pz() / infoPtr->pzB();


  Double M = sqrt(x1*x2 * infoPtr->s() );
  Double y = 1./2*log( x1/x2 );



  //cout << "Event number -- before FSR : "<< pyth->info.getCounter(3) << endl;
  //event.list();


  vector<Branch>  leftList = ColorReconection::CreateBackLists(event, 3);
  vector<Branch> rightList = ColorReconection::CreateBackLists(event, 4);


  ColArr colarr = ColRec.hardOrg;


  int isProb=ColorReconection::IsProblem(leftList, rightList, colarr);


	#if 0
  if(isProb) {
    isExclusive = false;  
    return;
  }
  else {
    int Ntry = 400000;
    int Nok = 0;
    bool isOK;
    Double pExc;
    for(int i=0; i < Ntry; ++i) {
      isOK = ColorReconection::RandomEmissions(proc, leftList, rightList, colarr,false);
      Nok += isOK;
    }
    if(Nok ==0) {
      int i = Ntry;
      Ntry += 1000;
      for(; i < Ntry; ++i) {
        isOK = ColorReconection::RandomEmissions(proc, leftList, rightList, colarr,false);
        Nok += isOK;
      }
    }
    if(Nok == 0) {
      isExclusive = false;
      cout << "No other suitaible color flow discovered" << endl;
      return;
    }
    pExc = Double(Nok)/(Ntry);

    weightEx = 1/pExc;

		Double sigma = sqrt(pExc*(1-pExc)/Ntry);
		cout << "weightEx(stat) " << setprecision(10) << 1./pExc <<" "<< endl;//pAnal << " "<< abs(pExc-pAnal)/sigma <<" "<< glEms.size() <<  endl;
  }
	#endif
	
	//weightEx = 1;
	//ColorReconection::printSpaceLikeShower(leftList,rightList,colarr);

	if(isProb) {
    isExclusive = false;  
    return;
	}
	else {
		weightEx *= ColorReconection::CountEmissions(proc, leftList, rightList);
		//cout << "weightEx " << weightEx  << endl;
	}

  //cout << "colWeight " <<TopoSel.nEm <<" "<< weightEx << endl;

  //if(WithoutEmissions)
    //weightEx *= GetExclusiveWeight(M, y, infoPtr->QRen(), infoPtr->QRen()  );
  {
    //weightEx *= GetExclusiveWeight(M, y, ExclusiveCut, lastEmScale ) * TopoSel.Weight();
    Double topoWeight = proc->topo.Weight();
    Double sudScale = (proc->topo.Type() == -1) ? ExclusiveCut : lastEmScale;
    if(topoWeight > 1e-10 && weightEx > 1e-6) {
      weightEx *= GetExclusiveWeight(M, y, sudScale, lastEmScale ) * topoWeight;
			proc->topo.FillEvent(weightEx);
		}
    else
      weightEx *= 0;
  }

  isExclusive = (weightEx > 1e-12 ) ? true : false; //minimal allowed weight

  if(!isExclusive)
    return;

  //cout << "Radek" << isExclusive << " "<<setprecision(10)<< weightEx <<  endl;


  //if problem change hard subprocess color flow
  if(isProb) {
    bool isOK = ColorReconection::ChangeColors(proc, leftList, rightList, colarr);
    isProb = !isOK;
  }

  //if still problem swap color emission
  if(isProb) {
    bool isOK = ColorReconection::SwapEmission(leftList, rightList, colarr);
    isProb = !isOK;
  }

  //if still problem randomly change shower and hard process
  if(isProb) {
    bool isOK;
    for(int i=0; i < 5000; ++i) {
      isOK = ColorReconection::RandomEmissions(proc, leftList, rightList, colarr);
      if(isOK)
        break;
    }
    isProb = !isOK;

  }


  if(isProb) {
    bool isOK;
    isOK = ColRec.ModifyColorsSingQQbarG(event);

    if(isOK)
      return;
    else {
      event.list();
      cout << "New suitaible color flow not found" << endl;
      cout << "IsExclusive " << int(isExclusive) << endl;

      isExclusive = false;
      ColorReconection::printSpaceLikeShower(leftList,rightList,colarr);
      return;
      //exit(1);
    }
  }


  ColorReconection::MakeColorSinglet(leftList, rightList, colarr);

  
  ColorReconection::CopyColorsToEvent(event, leftList, rightList, colarr);


  ColRec.StoreInColors(event);

  //cout << "Event number -- after FSR : "<< pyth->info.getCounter(3) << endl;
  //event.list();

	//cout <<"exWeight "<<setprecision(10)<< weightEx <<" "<<nMPI << endl;
}



void TransformMomenta(Pythia8::Vec4 pIn1, Pythia8::Vec4 pIn2, Pythia8::Vec4 &pOut1, Pythia8::Vec4 &pOut2);


bool ExclusiveHooks::doReconnectResonanceSystems(int, Pythia8::Event& event)
{
		/*
    for(int i=0; i < event.size(); ++i) {
			event[i].p(event[i].px(), event[i].py(), event[i].pz(), sqrt(event[i].pAbs2()+event[i].m2()));
      double massEr = abs(event[i].mCalc() - event[i].m()) / max(1., event[i].e() );
      //if(  > 0.004)
      //cout << "particleBefore " << i<<" "<<setprecision(10)<< massEr << endl;
    }
		*/

    if(!isExclusive)
      return true;

    //cout << "BeforeBefore : " <<pyth->info.getCounter(3)<< endl;
    //event.list();

		int nRemoved = 0;
	  int idLast;
    //count renmant's momenta and remove renmants
    Pythia8::Vec4 p1In, p2In;
    for(int index = event.size()-1 ; index >=  0; --index)
      if(event[index].status() == 63) {
        if(event[index].isAncestor(1)) {
          p1In += event[index].p();
        }
        else if(event[index].isAncestor(2)) {
          p2In += event[index].p();
        }
        else {
          cout << "Somthing wrong with renmant" << endl;
          exit(1);
        }
        event.remove(index, index);
				++nRemoved;
				idLast = index;
      }

		//Shift daughter info due to renmant removal
		for(int i=0; i < event.size(); ++i) {
			int daug1 = event[i].daughter1();
			int daug2 = event[i].daughter2();
      int moth1 = event[i].mother1();
      int moth2 = event[i].mother2();

			if(daug1 > idLast)
				daug1 -= nRemoved;
			if(daug2 > idLast)
				daug2 -= nRemoved;

			if(moth1 > idLast)
				moth1 -= nRemoved;
			if(moth2 > idLast)
				moth2 -= nRemoved;

			event[i].daughters(daug1, daug2);
			event[i].mothers(moth1, moth2);
		}

    Pythia8::Vec4 p1Out, p2Out;

    const Double mp = 0.938272046;
    TransformMomenta(p1In, p2In, p1Out, p2Out, mp);


    Double pz = infoPtr->eCM()/2;
    
    Pythia8::Vec4 pAMod,  pBMod;

    Pythia8::RotBstMatrix mat = CreateTransformation(4*pz*pz, 1-p1Out.pz()/pz, 1- abs(p2Out.pz())/pz, -p1, -p2, pAMod,  pBMod, 0);
      
    //cout << "Before : " <<pyth->info.getCounter(3)<< endl;
    //event.list();
    //cout <<"Rem1Rem2 "<< p1In<<" "<< p2In << endl;
    //cout <<"pAmod "<< pAMod << endl;
    //cout <<"pBmod "<< pBMod << endl;
    //cout <<"pABmod "<< pAMod+pBMod  << endl;

    //for(int i = event.size() -1; event[i].isFinal() && event[i].status() != 15  /*event[i].status() == 62 || event[i].status() == 23*/; --i)
      //event[i].rotbst(mat);
    for(int i = 0; i < event.size(); ++i)
      if(event[i].isFinal() && event[i].status() != 15)
        event[i].rotbst(mat);

    //cout << "After" << endl;
    //event.list();

    Pythia8::Vec4 pAfin,  pBfin;
    TransformMomenta(pAMod, pBMod, pAfin, pBfin, mp);
    //pAfin = pAMod;
    //pBfin = pBMod;

    //add proton
    event.append(event[1].id(), 15, 1, 0, 0, 0, 0, 0, pAfin, pAfin.mCalc() );
    event.append(event[2].id(), 15, 2, 0, 0, 0, 0, 0, pBfin, pBfin.mCalc() );

    //change colors back
    ColRec.RestoreInColors(event);



    //11      2203   uu_1                63     1     0     0     0     0   101     -0.966      0.302   1247.985   1247.985      0.771
    //12         1   d                   63     1     0     0     0   102     0     -0.832      0.701   4795.555   4795.555      0.330
    //13         2   u                   63     2     0     0     0   101     0      0.119     -0.710     -1.579      1.767      0.330
    //14      2101   ud_0                63     2     0     0     0     0   102      0.534     -1.543  -6996.910   6996.910      0.579
    return true;
}

Pythia8::RotBstMatrix ExclusiveHooks::CreateTransformation(Double s, Double x1, Double x2, Vec2 pA, Vec2 pB,
                                         Pythia8::Vec4 &pAMod, Pythia8::Vec4 &pBMod, Double mass)
{

  #define pow2(x)   (x)*(x)


  double sHat = x1*x2 * s;


  Pythia8::Vec4 pAbeam(0,0, sqrt(s/4-mass*mass),sqrt(s)/2);
  Pythia8::Vec4 pBbeam(0,0,-sqrt(s/4-mass*mass),sqrt(s)/2);


  //Vec2 pA(0.3,0.4);
  //Vec2 pB(0.2,0.3);
  //double m =1;

  // Begin kinematics of partons after primordial kT has been added.
  double sHatTAft = sHat + pow2( pA.px() + pB.px() )
                         + pow2( pA.py() + pB.py() );
  double w2A      = pA.norm2() - pA.norm2()/(1-x1);// + mass*mass*(x1*x1)/(1-x1) ;//mass*mass;
  double w2B      = pB.norm2() - pB.norm2()/(1-x2);// + mass*mass*(x2*x2)/(1-x2) ;//mass*mass;
  double w2Diff   = sHatTAft - w2A - w2B;
  double lambda   = pow2(w2Diff) - 4. * w2A * w2B;

  // Too large transverse momenta means that kinematics will not work.
  if (lambda <= 0.) {
    cout << "Unphysical" << endl;
    exit(1);
  }
  double lamRoot  = sqrt( lambda );


  double wPosBef       = x1 * sqrt(s);
  double wNegBef       = x2 * sqrt(s);


  double rescale  = sqrt(sHatTAft / sHat);
  double wPosAft  = rescale * wPosBef;
  double wNegAft  = rescale * wNegBef;

  double wPosA    = 0.5 * (sHatTAft + w2A - w2B + lamRoot) / wNegAft;
  double wNegB    = 0.5 * (sHatTAft + w2B - w2A + lamRoot) / wPosAft;


  //Pythia8::Vec4 pAMod, pBMod;

  // Store modified beam parton momenta.
  pAMod.e(  0.5 * (wPosA + w2A / wPosA) );
  pAMod.pz( 0.5 * (wPosA - w2A / wPosA) );
  pBMod.e(  0.5 * (w2B / wNegB + wNegB) );
  pBMod.pz( 0.5 * (w2B / wNegB - wNegB) );


  pAMod.px( pA.px() );
  pAMod.py( pA.py() );

  pBMod.px( pB.px() );
  pBMod.py( pB.py() );

  //cout <<  ( (pAMod+pBMod).m2Calc()-sHat )/ sHat  << endl;
  //cout << ( pAMod+pBMod).rap()  << " "<<1./2* log(x1/x2)<< endl;


  Pythia8::RotBstMatrix Msys;

  Pythia8::Vec4 pAorg(0,0, wPosBef/2, wPosBef/2 );
  Pythia8::Vec4 pBorg(0,0,-wNegBef/2, wNegBef/2 );

  // Construct system rotation and boost caused by primordial kT.
  Msys.reset();
  Msys.toCMframe(   pAorg, pBorg );
  Msys.fromCMframe( pAMod, pBMod );

  //Msys.toCMframe(   pAMod, pBMod );
  //Msys.fromCMframe( pAorg, pBorg );

  //cout << pAorg <<" "<<pBorg << endl;
  //cout << pAMod <<" "<<pBMod << endl;


  pAMod = pAbeam - pAMod;
  pBMod = pBbeam - pBMod;
  

  return Msys;
  //event[iBcopy].rotbst(Msys);
        
}



void ExclusiveHooks::TransformMomenta(Pythia8::Vec4 pIn1, Pythia8::Vec4 pIn2, Pythia8::Vec4 &pOut1, Pythia8::Vec4 &pOut2, Double mass)
{
  const Double m1 = mass;
  const Double m2 = mass;

  Pythia8::Vec4 pTotal = pIn1 + pIn2;
  pTotal.flip3();

  Pythia8::Vec4 &pIn1CMS = pIn1;
  Pythia8::Vec4 &pIn2CMS = pIn2;

  pIn1CMS.bst(pTotal);
  pIn2CMS.bst(pTotal);


  Double E = pIn1CMS.e() + pIn2CMS.e();
  Double pOld = pIn1CMS.pAbs();
  Double s= E*E;
  Double m12 = m1*m1;
  Double m22 = m2*m2;
  Double p2 = (s*s + m12*m12 + m22*m22 -2*m12*m22 - 2*s*(m12+m22) )/(4*s); 
  if(p2 < 0) {
    cout << "Troubles with momenta transormation, negative p2 = "<< p2 << endl;
    exit(1);
  }
  Double p = sqrt(p2);

  pOut1 = pIn1CMS * (p/pOld);
  pOut2 = pIn2CMS * (p/pOld);
  pOut1.e( sqrt(p2+m12) );
  pOut2.e( sqrt(p2+m22) );


  pTotal.flip3();

  pOut1.bst(pTotal);
  pOut2.bst(pTotal);


}



Double ExclusiveHooks::GetExclusiveWeight(Double M, Double y, Double Mju, Double MjuLoop )
{
  complex<Double> ExLum[3], ExLumC[3];

  id1 = pyth->event[pyth->event[1].daughter1()].id();
  //int id2 = pyth->event[pyth->event[2].daughter1()].id();

	int colAvg = abs(id1) <= 6 ? 9 : 64;

  kmr->LuminosityRand(id1, M, y, Mju, MjuLoop, p1, p2, ExLum, ExLumC );


  //ExLum[1] =0; ExLum[2] =0;
  //ExLumC[1]=0; ExLumC[2]=0;

  Double myInc  = proc->CrossSection(false, false) /9./4.;



	//cout << "Start Event" << endl;
	//cout <<"Pythia ME " <<  proc->PyCrossSection() << endl;
	//proc->MadGraphResult();
	//proc->PrintAmplitudes();
	//cout << "End Event" << endl;


	//cout <<"HHHHH Pythia X section " <<  proc->PyCrossSection() << endl;
	//cout <<"HHHHH My X section " <<  myInc << endl;

	//cout << endl << id1 <<" "<< id2 <<endl;
  //cout << "Start EventA" <<" "<< ExLum[0] <<" "<<ExLumC[0] << " "<< ExLum[1] <<" "<<ExLumC[1] << " "<< ExLum[2] <<" "<<ExLumC[2]<< endl;

  //cout <<setprecision(10)<< "Start Event" <<" "<< abs(ExLumi[0]-ExLum[0]) <<" "<< abs(ExLumiC[0]-ExLumC[0])<<" "
                             //<< abs(ExLumi[1]-ExLum[1]) <<" "<< abs(ExLumiC[2]-ExLumC[2])<<endl;

  //Exclusive cross section (multiplied by splittings and exclusive lumi)
  Double ExclusiveXsec = proc->CrossSectionFullFast(true, ExLum,  ExLumC  ) /Double(colAvg)/4.;

  //cout << "ExLum  " << ExLum[0]  <<" "<<ExLum[1]  <<" "<<ExLum[2] << endl;
  //cout << "ExLumC " << ExLumC[0] <<" "<<ExLumC[1] <<" "<<ExLumC[2] << endl;

  //cout << "RADEK normal " << NormalSplittings << endl;

  Double IncLum = kmr->LuminosityInc(id1, M, y, Mju);

  //Inclusive cross section (multiplied by splittings and inclusive lumi)
  Double InclusiveXsec = proc->PyCrossSection() * NormalSplittings * IncLum;

  //Double mySing = proc->CrossSection(true, true) /64./4.;
  //Double myInc  = proc->CrossSection(false, false) /64./4.;
	//cout << "Cross sections Ratio " << proc->PyCrossSection() / myInc <<" "<< endl;//mySing<< endl;
	//cout << "Cross sections Ratio " << (InclusiveXsec/IncLum) / myInc <<" "<< endl;//mySing<< endl;


  lumiFrac = (ExLum[0]*ExLumC[0]).real()/IncLum ;
  xSecFrac = IncLum/(ExLum[0]*ExLumC[0]).real() *   ExclusiveXsec/InclusiveXsec;

  //cout << "Cross sections " << myInc << " "<< proc->PyCrossSection() << endl;

	//cout << "myXsection "<< mySing <<" "<< ExclusiveXsec/ (ExLum[0]*ExLumC[0]).real()  << endl;
	//if( ExclusiveXsec/InclusiveXsec *IncLum/(ExLum[0]*ExLumC[0]).real()  > 1.18)
  /*
	if(proc->topo.nEm == 5) {
		cout << "xSecFrac "<< proc->topo.nEm<<" "<<M<<" "<<Mju<<" "<<ExclusiveXsec/InclusiveXsec *IncLum/(ExLum[0]*ExLumC[0]).real()  << endl;
		Double ratio = ExclusiveXsec/InclusiveXsec *IncLum/(ExLum[0]*ExLumC[0]).real();
		static Double sum = 0;
		static int    counter = 0;
		sum += ratio;
		++counter;
	  //cout << "InterMean " << sum / counter << endl;
	}
  */

  //cout << "singletXsec " << ExclusiveXsec/InclusiveXsec <<" "<< (ExLum[0]*ExLumC[0]).real()/IncLum << " "<< ExclusiveXsec/InclusiveXsec *IncLum/(ExLum[0]*ExLumC[0]).real() <<endl;



	/*
	//Create table
	for(y = -3; y <= 3; y += 0.1) {
		//y = 0;
		M = 200;
		complex<Double> ExLumM[3],ExLumCM[3];
		complex<Double> ExLumL[3],ExLumCL[3];
		complex<Double> ExLumH[3],ExLumCH[3];

		kmr->LuminosityRand(M,  y, Mju, p1, p2, ExLumM, ExLumCM );
		//kmr->LuminosityRand(M, -1., Mju, p1, p2, ExLumL, ExLumCL );
		//kmr->LuminosityRand(M, +1., Mju, p1, p2, ExLumH, ExLumCH );
		Double IncLum = kmr->LuminosityInc(M, y, Mju);

		cout <<setprecision(7)<< M<<" "<<y <<" "<< (ExLumM[0]*ExLumCM[0]).real()/IncLum << " "<< endl;
		//																						 (ExLumL[0]*ExLumCL[0]).real()/IncLum << " "<<
		//																						 (ExLumH[0]*ExLumCH[0]).real()/IncLum << endl;

	}
	exit(0);
	*/

  
  //cout << "Singlet x-sec : "<<setprecision(15)<< mySing*64*4  <<endl;//" "<<ExLum[0]*ExLumC[0]<< endl;
  //cout << endl;



  #if 0

  Double pyInc  = proc->PyCrossSection() * NormalSplittings;

  if( 1) {
		cout << "Look "<<setprecision(10)<<proc->GetCrossSectionBase()/9./4. /proc->PyCrossSection() << endl;
    cout << "Process : "<<setprecision(10)<<myInc/pyInc<<" : " << proc->getName() << endl;
    cout << "MyResult " <<  myInc << endl;
    cout << "Pythia   " <<  pyInc   << endl;
    cout << endl;
		pyth->process.list();
    /*
    cout << "MySing   " <<  mySing   << endl;
    cout << "ExX*lum   " <<  ExclusiveXsec   << endl;
    cout << "sigm*lum   " <<  mySing*ExLum[0]*ExLumC[0]   << endl;
    cout << endl;
    */
  }
  #endif

  pExc = ExclusiveXsec / InclusiveXsec;

	//cout << "RADEK pExc="<<setprecision(10) << pExc << endl;
	/*
  if((pExc > 0.1 && proc->topo.nEm == 0) || (pExc > 0.006 && proc->topo.nEm > 2) ||(pExc > 0.010 &&(proc->topo.nEm==2||proc->topo.nEm==1) )) {
    cout <<"nEm,pExc,sigExc,sigInc "<< proc->topo.nEm<<" "<<pExc <<" "<< ExclusiveXsec << " "<< InclusiveXsec << endl;

		kmr->LuminosityRand(M, y, Mju, MjuLoop, p1, p2, ExLum, ExLumC );
		Double ExclusiveXsec = proc->CrossSectionFullFast(true, ExLum,  ExLumC  ) /64./4.;
		cout << "pExc new " << ExclusiveXsec / InclusiveXsec << endl;

    cout <<"M,y,Mju "<< M<<" "<< y<<" "<< Mju << endl;


    pExc = 0;
  }
	*/

	//if(pExc > 0.08)
		//pExc = 0;


	/*
  if( abs(NormalSplittings-1)<0.00001 ) {
    Double myInc  = proc->CrossSection(false, false) /64./4.;
    //cout << "RATIO "<< myInc /proc->PyCrossSection() <<" "<< pExc << endl;
		static Double Max1 = 0;
		Max1 = max(Max1, pExc);
		cout << "Helenka1 " << Max1 << endl;
  }
	else {
		static Double MaxM = 0;
		MaxM = max(MaxM, pExc);
		cout << "HelenkaM " << MaxM << endl;
	}
  //cout << "Exclusive*Lumi : "<<setprecision(15)<< ExclusiveXsec  << endl;

  Double LumRat  = ExLum[0].real()*ExLumC[0].real()/IncLum;
  Double crossRat= mySing/proc->PyCrossSection();
  if(LumRat>0.01 &&  abs(NormalSplittings-1) <0.00001) {
    cout << "pExc"<<setprecision(15)<< pExc  <<" "<< LumRat <<" "<< crossRat<<  endl;
    cout << "M y " << M << " "<< y <<" "<< LumRat<< endl;
    cout << "RADEK"<< ExLum[0].real() <<" "<<ExLumC[0].real() << endl;

  }
	*/

  return pExc;

}


double ExclusiveHooks::biasSelectionBy( const Pythia8::SigmaProcess* ,
  const Pythia8::PhaseSpace* phaseSpacePtr, bool )
{
   //Double P0  =  1.61678e+03;
   //Double P1  = -3.10733e+00;
   //Double P2  =  5.19220e-05;

	Double P0,P1,P2;

   //P0 =  2.11859e+00;
   //P1 = -1.90577e+00;
   //P2 =  2.49797e-05;

	P0 =     0.227717;
	P1 =     0.760618;
	P2 = 0;


   Double M = sqrt( phaseSpacePtr->sHat() );

   Double &x = M;

   selBias = (P0*pow(x,P1) + P2);

   //cout << "Weight "<< selBias << endl;

   return selBias;

}



