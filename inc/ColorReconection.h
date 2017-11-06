#ifndef _ColorReconection_
#define _ColorReconection_
 
  #include "Sigma2Process.h"
  #include "Pythia8/Pythia.h"

  bool operator==(const ColArr &a, const ColArr &b);

  struct Branch;

class ColorReconection {

  int inId1, inId2;

  int inCol1, inaCol1;
  int inCol2, inaCol2;




public:
  ColArr hardOrg;

  void StoreInColors(Pythia8::Event &event);
  void RestoreInColors(Pythia8::Event &event) const;

  static vector<Branch> CreateBackLists(const Pythia8::Event &event, int start);

  static int  IsProblem(const vector<Branch > &leftList, const vector<Branch> &rightList, const ColArr &colarr);

  static void MakeColorSinglet(vector<Branch> &leftList, vector<Branch> &rightList, ColArr &hard);

  static void CopyColorsToEvent(Pythia8::Event &event, const vector<Branch> &leftList, const vector<Branch> &rightList, const ColArr &hard);

  static bool ChangeColors(Sigma2Process *proc, vector<Branch > &leftList, vector<Branch> &rightList, ColArr &hardOld);


  static bool SwapEmission(vector<Branch > &leftList, vector<Branch> &rightList, const ColArr &hardOld);

  static bool RandomEmissions(Sigma2Process *proc, vector<Branch > &leftList, vector<Branch> &rightList, ColArr &hardOld, bool replace=true);

  bool ModifyColorsSingQQbarG( Pythia8::Event &event);

  static void printSpaceLikeShower(const vector<Branch> &leftList, const vector<Branch> &rightList, const ColArr &colarr);


	static Double CountEmissions(Sigma2Process *proc, const vector<Branch > &leftList, const vector<Branch> &rightList);
	static bool TestEmission( const vector<Branch > &leftList, const vector<Branch> &rightList, const ColArr &colHard,
																	 int InverseEmId, int emCol);

	static vector<pair<int,int> > GetGluonEmissions(const vector<Branch > &leftList, const vector<Branch> &rightList);


};

  struct Branch {
    Branch(int _idIn, int _colIn, int _acolIn) {
      idIn =_idIn;
      colIn=_colIn;
      acolIn = _acolIn;
      idOut =_idIn;
      colOut=_colIn;
      acolOut = _acolIn;

      idEm = -1;
    }

    Branch(int col1, int acol1, int col2, int acol2) {
      colIn   = col1;
      acolIn  = acol1;
      colOut  = col2;
      acolOut = acol2;

      idEm = -1;
    }

    Branch() {}


    int idIn, idOut, idEm;
    int colIn, colOut, colEm;
    int acolIn, acolOut, acolEm;

    void print() const {
      if(idEm == -1)
        cout << idIn <<" ("<<colIn<< " "<<acolIn<<") -> "<< idOut <<" ("<<colOut<< " "<<acolOut<<")" << endl;
      else {
        cout << idIn <<" ("<<colIn<< " "<<acolIn<<") -> "<< idOut <<" ("<<colOut<< " "<<acolOut<<") : ";
        cout << idEm <<" ("<<colEm<< " "<<acolEm<<")" << endl;
      }

    }

    void substitute(map<int,int> rep) {

      if(idEm != -1) {
        if(rep.count(colEm) > 0)
          colEm = rep[colEm];
        if(rep.count(acolEm) > 0)
          acolEm = rep[acolEm];
      }

      if(rep.count(colIn) > 0)
        colIn = rep[colIn];
      if(rep.count(acolIn) > 0)
        acolIn = rep[acolIn];

      if(rep.count(colOut) > 0)
        colOut = rep[colOut];
      if(rep.count(acolOut) > 0)
        acolOut = rep[acolOut];

    }

    bool check(int Lc, int Lac, int Rc, int Rac) const {
      Branch b = *this;

      if(b.colIn > b.acolIn) swap(b.colIn, b.acolIn);
      if(b.colOut > b.acolOut) swap(b.colOut, b.acolOut);
      if(b.colEm > b.acolEm) swap(b.colEm, b.acolEm);

      if(Lc > Rac) swap(Lc, Rac);
      if(Rc > Lac) swap(Rc, Lac);

      bool OK = true;

			//gluon+gluon
			if( Lc > 0 && Rc > 0) {
				if(b.idEm != -1 &&  ( (Lc == b.colEm  && Rac == b.acolEm)  ||  (Rc == b.colEm  && Lac == b.acolEm)   ) )
					OK = false;

				else if( (Lc == b.colIn  && Rac == b.acolIn)  ||  (Rc == b.colIn  && Lac == b.acolIn   ) )
					OK = false;

				else if( (Lc == b.colOut  && Rac == b.acolOut) || (Rc == b.colOut  && Lac == b.acolOut   ) )
					OK = false;
			}
			//qbar+q
			else if( Rac == 0 ) {
				if(b.idEm != -1 &&   (Rc == b.colEm  && Lac == b.acolEm)  )
					OK = false;

				else if(  Rc == b.colIn  && Lac == b.acolIn  )
					OK = false;

				else if( Rc == b.colOut  && Lac == b.acolOut  )
					OK = false;
			}
			//q+qbar
			else if( Lac == 0 ) {
				if(b.idEm != -1 &&  ( Lc == b.colEm  && Rac == b.acolEm ) )
					OK = false;

				else if( Lc == b.colIn  && Rac == b.acolIn  )
					OK = false;

				else if( Lc == b.colOut  && Rac == b.acolOut  )
					OK = false;
			}
			else {
				cout << "Problem" << endl;
				cout << Lc <<" "<< Lac<<" "<< Rc <<" "<< Rac << endl;
				exit(1);
			}


      return OK;
    }

  };

#endif
