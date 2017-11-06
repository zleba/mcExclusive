#ifndef _ExclusiveHooks_
#define _ExclusiveHooks_


#include "Pythia8/Pythia.h"
#include "SpinTensor.h"
#include "ColorTensor.h"
#include "Sigma2Process.h"
#include "ColorReconection.h"
#include "TopoSelector.h"

#include "KMRlum.h"
#include "KMRlumi.h"

#include <map>

struct Branch;

class ExclusiveHooks : public Pythia8::UserHooks{

  public:
  ExclusiveHooks(Pythia8::Pythia *py, bool _withoutEmissions);

	void ModifyParameters(int procId, int pRest, int p0, int p1, int p2);

  virtual bool initAfterBeams();

  virtual bool canVetoISREmission() {return true;}
  virtual bool doVetoISREmission( int sizeOld, const Pythia8::Event& event, int iSys);

  virtual void modifyEventBeforFSR( Pythia8::Event &event );

  virtual bool canReconnectResonanceSystems() { return  true; }
  virtual bool doReconnectResonanceSystems(int oldSizeEvent, Pythia8::Event& event);

	virtual bool canVetoStep() { return true; }
	virtual int numberVetoStep() { return 1; }
	virtual bool doVetoStep(int iPos, int nISR, int nFSR, const Pythia8::Event& event) {
        (void)iPos; (void)nFSR; (void)event; //Unused
		//cout <<"RADEK "<< iPos <<" "<< nISR << " " << nFSR << endl;
		//cout << "radecek " << event.size() <<" "<<isSmallPt<<" "<< nISR << endl;
		//return false;
		const Double prescale = 1;
		if(nISR == 1) {
			weightEx *= prescale;
			if(KMRlumi::Uniform(0,1) <1./prescale)
				return false;
			else {
				//cout << "Exiting" << endl;
				return true;
			}
		}
		return false;
	}
	
  virtual bool canVetoMPIEmission() { return true; }

  virtual bool doVetoMPIEmission( int sizeOld, const Pythia8::Event& event)
	{
    (void)sizeOld; (void)event;//Unused
    /*if(isSmallPt&& isExclusive) cout<<"ev Removed"<<endl; */
    ++nMPI;
    //weightEx=0; isSmallPt=true; isExclusive = false;
	  return true;
  }


  virtual bool canVetoProcessLevel() {return true;} 
  virtual bool doVetoProcessLevel(Pythia8::Event& event);

  virtual bool canBiasSelection() {return false; }
  virtual double biasSelectionBy( const Pythia8::SigmaProcess* sigmaProcessPtr,
  const Pythia8::PhaseSpace* phaseSpacePtr, bool inEvent);



  inline bool IsExclusive() const { return isExclusive; }
  inline Double ExclusiveWeight() const { return weightEx; }
  inline int GetnMPI() const { return nMPI;}
	inline int GetId1() const { return id1;}

	inline void EventParams(int &_nEm, Double &_a) { _nEm=TopoSel.nEm; _a=TopoSel.a; }
	void PrintBestParams() {
		for(std::map<int,Sigma2Process *>::iterator it = procs.begin(); it != procs.end(); ++it) {
			if(it->second->topo.GetEntries() > 0) {
				cout << "Process " << it->second->getName() <<" (id="<<it->first<<")"<< endl;
				it->second->topo.PrintResult();
				cout << endl;
			}
		}
	}
	Double LastScale() const { return lastEmScale; }
  Double LumiFrac() const { return lumiFrac;}
  Double XsecFrac() const { return xSecFrac;}
  Double GetPexc() const { return pExc;}
  Sigma2Process* GetProc() const {return proc;}

  Double GetExclusiveWeight(Double M, Double y, Double Mju, Double MjuLoop );


private:
  Pythia8::Pythia *pyth;

  //map containing all processes
  std::map<int,Sigma2Process *> procs;

  //pointer to current process
  Sigma2Process *proc;

  //spin and color tensors
  SpinTensor   spinL, spinR;
  ColorTensor  colL, colR;

  //product of inclusive splitting functions
  Double NormalSplittings;

  //pointer to tool calculating KMR luminosity
  KMRlumi *kmr;

  //class for reconect color lines
  ColorReconection ColRec;

  //class to sellect the generation method
  TopoSelector TopoSel;

  //tag for exluclusive event
  bool isExclusive;
  
  //event with space-like shower pt below set cut-off
  bool isSmallPt;

  //without ISR (equivalent to standard KMR approach)
  bool WithoutEmissions;

	//number of MPI
	int nMPI;

  //weight for exclusive events
  Double weightEx;

  //Treashold for exclusive space-like shower
  Double ExclusiveCut;

  //Transverse momenta of scattered protons
  Vec2 p1, p2;

	//pdgID of "left" parton incoming to the shower
	int id1;

  //Last Scale in Space-like shower
  Double lastScale, lastEmScale;

  //LumiFraction XsecFraction
  Double lumiFrac, xSecFrac, pExc;

  static void TransformMomenta(Pythia8::Vec4 pIn1, Pythia8::Vec4 pIn2, Pythia8::Vec4 &pOut1, Pythia8::Vec4 &pOut2, Double mass);
	void ResetShowers(const Pythia8::Event& event);

  static Pythia8::RotBstMatrix CreateTransformation(Double s, Double x1, Double x2, Vec2 pA, Vec2 pB,
                                         Pythia8::Vec4 &pAMod, Pythia8::Vec4 &pBMod, Double mass);

};


class HEPinfo {
  public:
  static void Init(Pythia8::ParticleData* _particleDataPtr, Pythia8::CoupSM* _couplingsPtr)
         { particleDataPtr = _particleDataPtr; couplingsPtr = _couplingsPtr; }
  static Pythia8::ParticleData* ParticleData() {return particleDataPtr; }
  static Pythia8::CoupSM*    Couplings() {return couplingsPtr; }

  static Pythia8::ParticleData*  particleDataPtr;
  static Pythia8::CoupSM*     couplingsPtr;
};


#endif
