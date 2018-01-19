#include "ColorReconection.h"
#include "SigmaQCD.h"
#include "KMRlumi.h"

#include <typeinfo>



bool operator==(const ColArr &a, const ColArr &b)
{
  bool isSame = true;
  for(int i = 0; i < 8; ++i)
    isSame = isSame && (a.cols[i] == b.cols[i] );
  return isSame;
}

void ColorReconection::StoreInColors(Pythia8::Event &event)
{
  //int idStart1, idStart2;
  inId1 = event[1].daughter1();
  inId2 = event[2].daughter1();

   
  inCol1 = event[inId1].col(); inaCol1 = event[inId1].acol();
  inCol2 = event[inId2].col(); inaCol2 = event[inId2].acol();

	int  Col = (inCol1  > 0) ? 10000 : 0;
	int aCol = (inaCol1  > 0) ? 10001 : 0;

  event[inId1].cols(Col, aCol);
  event[inId2].cols(aCol, Col);

}

void ColorReconection::RestoreInColors(Pythia8::Event &event) const
{
  int m1In1, m2In1;
  int m1In2, m2In2;

  event[inId1].cols(inCol1, inaCol1);
  event[inId2].cols(inCol2, inaCol2);

  m1In1 = event[inId1].mother1(); m2In1 = event[inId1].mother2();
  m1In2 = event[inId2].mother1(); m2In2 = event[inId2].mother2();

  if( m1In1 == m2In1 && m1In1 > inId1)
    event[m1In1].cols(inCol1, inaCol1);

  if( m1In2 == m2In2 && m1In2 > inId2)
    event[m1In2].cols(inCol2, inaCol2);

}

vector<Branch> ColorReconection::CreateBackLists(const Pythia8::Event &event, int start)
{
	vector<Branch> List;

  List.push_back( Branch( start,  event[start].col(),  event[start].acol() ) ); 

  for(int i = 0; i< 100; ++i) {
    Branch br = List.back();

    int mot1 = event[br.idIn].mother1();
    int mot2 = event[br.idIn].mother2();

    if(mot1 < 3)
      break;

    if(mot1 == mot2) {

      Branch brTemp;

      brTemp.idEm = -1;

      brTemp.idOut   = br.idIn;
      brTemp.colOut  = event[br.idIn].col();
      brTemp.acolOut = event[br.idIn].acol();

      brTemp.idIn   = mot1;
      brTemp.colIn  = event[mot1].col();
      brTemp.acolIn = event[mot1].acol();



      List.push_back( brTemp );
    }

    if(mot1 > 0 && mot2 == 0) {
      int d1 = event[mot1].daughter1();
      int d2 = event[mot1].daughter2();

      Branch brTemp;
      brTemp.idEm = d1;
      brTemp.colEm = event[d1].col();
      brTemp.acolEm = event[d1].acol();

      brTemp.idOut   = d2;
      brTemp.colOut  = event[d2].col();
      brTemp.acolOut = event[d2].acol();

      brTemp.idIn   = mot1;
      brTemp.colIn  = event[mot1].col();
      brTemp.acolIn = event[mot1].acol();

      List.push_back(brTemp);

    }

  }
	return List;

}




int ColorReconection::IsProblem(const vector<Branch > &leftList, const vector<Branch> &rightList, const ColArr &colarr)
{

  int  cL = leftList.back().colIn;
  int acL = leftList.back().acolIn;

  int  cR = rightList.back().colIn;
  int acR = rightList.back().acolIn;

  int sameIndex = 0;

  for(unsigned int i = 0; i < leftList.size(); ++i) 
    if(!leftList[i].check(cL, acL, cR, acR)) {
      sameIndex = 1000 + i;
      return sameIndex;
    }

  for(unsigned int i = 0; i < rightList.size(); ++i) 
    if(!rightList[i].check(cL, acL, cR, acR)) {
      sameIndex = 2000 + i;
      return sameIndex;
    }

  Branch hardBr(colarr.c5(), colarr.c6(), colarr.c7(), colarr.c8() );

  if(!hardBr.check(cL, acL, cR, acR) ) {
    sameIndex = 3000;
    return sameIndex;
  }

  return sameIndex;

}




void ColorReconection::MakeColorSinglet(vector<Branch> &leftList, vector<Branch> &rightList, ColArr &hard)
{
  int  cL = leftList.back().colIn;
  int acL = leftList.back().acolIn;

  int  cR = rightList.back().colIn;
  int acR = rightList.back().acolIn;

  int c1 = min(cL, acR);
  int c2 = min(cR, acL);

  map<int,int> rule;

  rule[cL]  = c1;
  rule[acR] = c1;

  rule[cR]  = c2;
  rule[acL] = c2;
  
  
  for(unsigned int i = 0; i < leftList.size(); ++i) {
    leftList[i].substitute(rule);
  }
  for(unsigned int i = 0; i < rightList.size(); ++i) {
    rightList[i].substitute(rule);
  }
  hard.substitute(rule);

}



void ColorReconection::CopyColorsToEvent(Pythia8::Event &event,
     const vector<Branch> &leftList, const vector<Branch> &rightList, const ColArr &hard)
{
  int d1,d2;
  int idLoop;

  for(unsigned int i=0; i < leftList.size(); ++i) {
    const Branch & b = leftList[i];

    event[b.idIn].cols( b.colIn, b.acolIn );

    event[b.idOut].cols( b.colOut, b.acolOut );

    if(b.idEm != -1) {
      idLoop = b.idEm;
      do {
        d1 = event[idLoop].daughter1();
        d2 = event[idLoop].daughter2();

        event[idLoop].cols(b.colEm, b.acolEm);

        idLoop = d1;
      } while(d1 == d2 && d1 > 0);
    }
  }

  for(unsigned int i=0; i < rightList.size(); ++i) {
    const Branch & b = rightList[i];

    event[b.idIn].cols( b.colIn, b.acolIn );

    event[b.idOut].cols( b.colOut, b.acolOut );

    if(b.idEm != -1) {
      idLoop = b.idEm;
      do {
        d1 = event[idLoop].daughter1();
        d2 = event[idLoop].daughter2();

        event[idLoop].cols(b.colEm, b.acolEm);

        idLoop = d1;
      } while(d1 == d2 && d1 > 0);
    }
  }

  int idOut1;
  for(idOut1=0; idOut1<event.size()&&abs(event[idOut1].status())!=23 &&  (abs(event[idOut1].status())!=22  || abs(event[idOut1].id())!=6); ++idOut1)
    ;
  if(idOut1 == event.size() )
    return;

  //outgoing from hard process
    idLoop = idOut1;
    do {
      d1 = event[idLoop].daughter1();
      d2 = event[idLoop].daughter2();

      event[idLoop].cols(hard.c5(), hard.c6());

      idLoop = d1;
    } while(d1 == d2 && d1 > 0);

  //outgoing from hard process
    idLoop = idOut1+1;
    do {
      d1 = event[idLoop].daughter1();
      d2 = event[idLoop].daughter2();

      event[idLoop].cols(hard.c7(), hard.c8());

      idLoop = d1;
    } while(d1 == d2 && d1 > 0);



}


bool ColorReconection::ChangeColors(Sigma2Process *proc, vector<Branch > &leftList, vector<Branch> &rightList, ColArr &hardOld)
{
  ColArr hardNew;

  //database of used collors
  set<int> usedCols;
  ColArr  hardOldRed = hardOld;
  hardOldRed.addoffset(-100);
  usedCols.insert( hardOldRed.calcHash() );
  
  for(int nn = 0; usedCols.size() < proc->getNcolFlows() ; ++nn) {
    unsigned int i = 0;
    bool status;
    //cout << "col flow index " <<nn << endl;
    do {
      hardNew = proc->PyColorFlow();
      status = ( usedCols.insert(hardNew.calcHash() )  ).second;
      ++i;
    } while((status == false) && i < 10000);
    if(i > 9990) {
      cout << "Problem with color choice " <<i <<" "<< usedCols.size()<< endl;
      cout << proc->getName() << ", nFlows " << proc->getNcolFlows() << endl;
      if( typeid(*proc) == typeid(Sigma2gg2qqbar) ) {
        cout <<"Hura " << endl;
        (dynamic_cast<Sigma2gg2qqbar*> (proc))->printInv();
      }
      return false;
      //proc->Inv()
      //exit(1);
    }


    hardNew.addoffset(100);


    map<int,int> leftRule, rightRule;

    leftRule[hardOld.c1()] = hardNew.c1();
    leftRule[hardOld.c2()] = hardNew.c2();

    rightRule[hardOld.c3()] = hardNew.c3();
    rightRule[hardOld.c4()] = hardNew.c4();

    vector<Branch> leftListNew  = leftList;
    vector<Branch> rightListNew = rightList;

    for(i = 0; i < leftListNew.size(); ++i) {
      leftListNew[i].substitute(leftRule);
    }
    for(i = 0; i < rightListNew.size(); ++i) {
      rightListNew[i].substitute(rightRule);
    }
    
    bool notOK = IsProblem(leftListNew, rightListNew, hardNew) > 0;

    if(notOK == false) {
      //cout << "Solution found" << endl;
      leftList = leftListNew;
      rightList = rightListNew;
      hardOld = hardNew;

      return true;
    }

  }
  return false;

}


int MaxIndex(vector<Branch > &List, unsigned int ind) 
{
  if(List.size() <=ind) {
    cout <<"Index higher than list size!" << endl;
    exit(1);
  }

  int Max = 0;
  for(unsigned int i = 0; i < ind; ++i) {
    Branch &b = List[i];
    Max = max(Max, b.colIn);
    Max = max(Max, b.colOut);

    Max = max(Max, b.acolIn);
    Max = max(Max, b.acolOut);

    if(b.idEm != -1) {
      Max = max(Max, b.colEm);
      Max = max(Max, b.acolEm);
    }
  }
  Max = max(Max, List[ind].colOut);
  Max = max(Max, List[ind].acolOut);

  return Max;
}


bool ColorReconection::RandomEmissions(Sigma2Process *proc, vector<Branch > &leftList, vector<Branch> &rightList, ColArr &hardOld, bool replace)
{

  vector<Branch> leftListNew ( leftList.size() );
  vector<Branch> rightListNew( rightList.size() );

  ColArr hardNew;

  hardNew = proc->PyColorFlow();
  
  int lastCol;
  lastCol = hardNew.maximum();


  for(unsigned int i = 0; i < leftList.size(); ++i) {
    //Evolution

    Branch *bE;
    Branch *bEnew, *bCnew;

    if( i > 0) {
      leftListNew[i].colOut  = leftListNew[i-1].colIn;
      leftListNew[i].acolOut = leftListNew[i-1].acolIn;

      rightListNew[i].colOut  = rightListNew[i-1].colIn;
      rightListNew[i].acolOut = rightListNew[i-1].acolIn;
    }
    else {
      leftListNew[i].colOut  = hardNew.c1();
      leftListNew[i].acolOut = hardNew.c2();

      rightListNew[i].colOut  = hardNew.c3();
      rightListNew[i].acolOut = hardNew.c4();

    }

    if(leftList[i].idEm != -1 || rightList[i].idEm !=-1) {
      if(leftList[i].idEm != -1 && rightList[i].idEm ==-1) {
        bE = &leftList[i];
        bEnew = &leftListNew[i]; bCnew = &rightListNew[i];
      }
      else if(leftList[i].idEm == -1 && rightList[i].idEm != -1) {
        bE = &rightList[i];
        bEnew = &rightListNew[i]; bCnew = &leftListNew[i];
      }
      else{
        cout << "Emissions from both sides"<<endl;
        exit(1);
      }

      //g -> g g
      if( bE->colIn  >0 && bE->acolIn>0 && bE->colOut>0 &&
          bE->acolOut>0 && bE->colEm >0 && bE->acolEm>0 ) {
        if ( KMRlumi::Uniform(0,1) > 0.5 ) {
          bEnew->acolIn = bEnew->acolOut;
          bEnew->acolEm = bEnew->colOut;
          bEnew->colEm = ++lastCol;
          bEnew->colIn = bEnew->colEm;
        }
        else {
          bEnew->colIn = bEnew->colOut;
          bEnew->colEm = bEnew->acolOut;
          bEnew->acolEm = ++lastCol;
          bEnew->acolIn = bEnew->acolEm;
        }

      }
      else {

        bool emCol = bE->colIn  == bE->colEm && bE->acolEm == bE->colOut
                 && bE->acolIn == bE->acolOut;

        bool emaCol = bE->acolIn == bE->acolEm && bE->colEm == bE->acolOut
                 && bE->colIn  == bE->colOut;

        //
        if ( emCol ) {
          bEnew->acolIn = bEnew->acolOut;
          bEnew->acolEm = bEnew->colOut;
          if(bE->colEm != 0 ) {
            bEnew->colEm = ++lastCol;
            bEnew->colIn = bEnew->colEm;
          }
          else {
            bEnew->colEm = 0;
            bEnew->colIn = 0;
          }
        }
        //
        if ( emaCol ) {
          bEnew->colIn = bEnew->colOut;
          bEnew->colEm = bEnew->acolOut;
          if(bE->acolEm != 0 ) {
            bEnew->acolEm = ++lastCol;
            bEnew->acolIn = bEnew->acolEm;
          }
          else {
            bEnew->acolEm = 0;
            bEnew->acolIn = 0;
          }
        }

      }

      //copy the second branch
      bCnew->colIn  = bCnew->colOut;
      bCnew->acolIn = bCnew->acolOut;
    }
    //copy of both branches
    else {
      leftListNew[i].colIn = leftListNew[i].colOut;
      leftListNew[i].acolIn = leftListNew[i].acolOut;

      rightListNew[i].colIn = rightListNew[i].colOut;
      rightListNew[i].acolIn = rightListNew[i].acolOut;

    }

  }
  //cout << "Trying new emission configuration" << endl;
  //printSpaceLikeShower(leftListNew,rightListNew,hardOld);;


  int newSing = IsProblem(leftListNew, rightListNew, hardNew);
  if(newSing == 0) {
    if(replace) {
      leftList= leftListNew;
      rightList= rightListNew;
      hardOld = hardNew;
    }
    return true;
  }
  else
    return false;

}


vector<pair<int,int> > ColorReconection::GetGluonEmissions(const vector<Branch > &leftList, const vector<Branch> &rightList)
{
	vector<pair<int,int> > glEms;
  for(unsigned int i = 0; i < leftList.size(); ++i) {

		if( leftList[i].idEm  != -1 && leftList[i].colEm  > 0 && leftList[i].acolEm  > 0 ) 
			if( leftList[i].colIn  > 0 && leftList[i].acolIn  > 0 )
         glEms.push_back( make_pair(i, -1) );

		if( rightList[i].idEm != -1 && rightList[i].colEm > 0 && rightList[i].acolEm > 0 ) 
			if( rightList[i].colIn  > 0 && rightList[i].acolIn  > 0 )
				glEms.push_back( make_pair(i, +1) );

		if( leftList[i].idEm  != -1  && rightList[i].idEm != -1 ) {
			cout << "Gluon emissions from both sides!!" << endl;
			exit(1);
		}
			
	}
	return glEms;

}

Double ColorReconection::CountEmissions(Sigma2Process *proc, const vector<Branch > &leftList, const vector<Branch> &rightList)
{

  vector<Branch> leftListNew ( leftList.size() );
  vector<Branch> rightListNew( rightList.size() );

  vector<ColArr> colHardArr;
  vector<Double> probArr;

  proc->PyColorFlows(colHardArr, probArr);
  

	vector<pair<int,int> > glEms = ColorReconection::GetGluonEmissions(leftList, rightList);

	Double fracBad = 0;
	//One emission for different color line than others

	if(glEms.size() >= 3) {
		for(unsigned int i = 0; i < glEms.size(); ++i) {
			int InverseEmId = glEms[i].first;
			//int EmSide = glEms[i].second;

			for(unsigned int j=0; j < colHardArr.size(); ++j) {
				bool isBad1 = !ColorReconection::TestEmission(leftList,rightList,colHardArr[j],InverseEmId, +1);
				bool isBad2 = !ColorReconection::TestEmission(leftList,rightList,colHardArr[j],InverseEmId, -1);

				fracBad += isBad1 * probArr[j];
				fracBad += isBad2 * probArr[j];
			}
		}
	}

	if(glEms.size() <= 2 && glEms.size() > 0) {
			int InverseEmId = glEms[0].first;
			//int EmSide = glEms[i].second;

			for(unsigned int j=0; j < colHardArr.size(); ++j) {
				bool isBad1 = !ColorReconection::TestEmission(leftList,rightList,colHardArr[j],InverseEmId, +1);
				bool isBad2 = !ColorReconection::TestEmission(leftList,rightList,colHardArr[j],InverseEmId, -1);

				fracBad += isBad1 * probArr[j];
				fracBad += isBad2 * probArr[j];
			}
	}



	//All emissions from one color line
	if(glEms.size() >= 2) {
		for(unsigned int j=0; j < colHardArr.size(); ++j) {
			bool isBad1 = !ColorReconection::TestEmission(leftList,rightList,colHardArr[j],-1, +1);
			bool isBad2 = !ColorReconection::TestEmission(leftList,rightList,colHardArr[j],-1, -1);

			fracBad += isBad1 * probArr[j];
			fracBad += isBad2 * probArr[j];
		}
	}
	if(glEms.size() == 0) {
		for(unsigned int j=0; j < colHardArr.size(); ++j) {
			bool isBad1 = !ColorReconection::TestEmission(leftList,rightList,colHardArr[j],-1, +1);
			fracBad += isBad1 * probArr[j];
		}
	}

	Double totComb = pow(2, glEms.size() );

	Double corrWeight =  totComb/(totComb-fracBad);

	return corrWeight;
}

bool ColorReconection::TestEmission( const vector<Branch > &leftList, const vector<Branch> &rightList, const ColArr &colHard,
																		 int InverseEmId, int emColDir)
{

  vector<Branch> leftListNew ( leftList.size() );
  vector<Branch> rightListNew( rightList.size() );


  int lastCol;
  lastCol = colHard.maximum();


  for(unsigned int i = 0; i < leftList.size(); ++i) {
    //Evolution

    const Branch *bE;
    Branch *bEnew, *bCnew;


    if( i > 0) {
      leftListNew[i].colOut  = leftListNew[i-1].colIn;
      leftListNew[i].acolOut = leftListNew[i-1].acolIn;

      rightListNew[i].colOut  = rightListNew[i-1].colIn;
      rightListNew[i].acolOut = rightListNew[i-1].acolIn;

    }
    else {
      leftListNew[i].colOut  = colHard.c1();
      leftListNew[i].acolOut = colHard.c2();

      rightListNew[i].colOut  = colHard.c3();
      rightListNew[i].acolOut = colHard.c4();

    }

		int side;

    if(leftList[i].idEm != -1 || rightList[i].idEm !=-1) {
      if(leftList[i].idEm != -1 && rightList[i].idEm ==-1) {
        bE = &leftList[i];
        bEnew = &leftListNew[i]; bCnew = &rightListNew[i];
				side = -1;
      }
      else if(leftList[i].idEm == -1 && rightList[i].idEm != -1) {
        bE = &rightList[i];
        bEnew = &rightListNew[i]; bCnew = &leftListNew[i];
				side = +1;
      }
      else{
        cout << "Emissions from both sides"<<endl;
        exit(1);
      }

      //g -> g g
      if( bE->colIn  >0 && bE->acolIn>0 && bE->colOut>0 &&
          bE->acolOut>0 && bE->colEm >0 && bE->acolEm>0 ) {

				int EmColNow = side *   emColDir  * (InverseEmId == int(i) ? -1 : 1);

        if ( EmColNow == 1 ) {
          bEnew->acolIn = bEnew->acolOut;
          bEnew->acolEm = bEnew->colOut;
          bEnew->colEm = ++lastCol;
          bEnew->colIn = bEnew->colEm;
        }
        else {
          bEnew->colIn = bEnew->colOut;
          bEnew->colEm = bEnew->acolOut;
          bEnew->acolEm = ++lastCol;
          bEnew->acolIn = bEnew->acolEm;
        }

      }
      else {

        bool emCol = bE->colIn  == bE->colEm && bE->acolEm == bE->colOut
                 && bE->acolIn == bE->acolOut;

        bool emaCol = bE->acolIn == bE->acolEm && bE->colEm == bE->acolOut
                 && bE->colIn  == bE->colOut;

        //
        if ( emCol ) {
          bEnew->acolIn = bEnew->acolOut;
          bEnew->acolEm = bEnew->colOut;
          if(bE->colEm != 0 ) {
            bEnew->colEm = ++lastCol;
            bEnew->colIn = bEnew->colEm;
          }
          else {
            bEnew->colEm = 0;
            bEnew->colIn = 0;
          }
        }
        //
        if ( emaCol ) {
          bEnew->colIn = bEnew->colOut;
          bEnew->colEm = bEnew->acolOut;
          if(bE->acolEm != 0 ) {
            bEnew->acolEm = ++lastCol;
            bEnew->acolIn = bEnew->acolEm;
          }
          else {
            bEnew->acolEm = 0;
            bEnew->acolIn = 0;
          }
        }

      }

      //copy the second branch
      bCnew->colIn  = bCnew->colOut;
      bCnew->acolIn = bCnew->acolOut;
    }
    //copy of both branches
    else {
      leftListNew[i].colIn = leftListNew[i].colOut;
      leftListNew[i].acolIn = leftListNew[i].acolOut;
			

      rightListNew[i].colIn = rightListNew[i].colOut;
      rightListNew[i].acolIn = rightListNew[i].acolOut;
			

    }
		leftListNew[i].idIn  = leftList[i].idIn;
		leftListNew[i].idOut = leftList[i].idOut;
		leftListNew[i].idEm  = leftList[i].idEm;

		rightListNew[i].idIn  = rightList[i].idIn;
		rightListNew[i].idOut = rightList[i].idOut;
		rightListNew[i].idEm  = rightList[i].idEm;


  }
  //cout << "Trying new emission configuration" << endl;
  //printSpaceLikeShower(leftListNew,rightListNew,hardOld);;


  int newSing = IsProblem(leftListNew, rightListNew, colHard);
	//cout << "print posibilities " << InverseEmId << " "<< emColDir << endl;
  //printSpaceLikeShower(leftListNew,rightListNew,colHard);
	//cout << "print org" << endl;
  //printSpaceLikeShower(leftList,rightList,colHard);
	

  if(newSing == 0) {
    return true;
  }
  else
    return false;

}








bool ColorReconection::SwapEmission(vector<Branch > &leftList, vector<Branch> &rightList, const ColArr &hardOld)
{

  vector<Branch> leftListNew  = leftList;
  vector<Branch> rightListNew = rightList;
  

    int singlGluon = IsProblem(leftList, rightList, hardOld);

    unsigned int ind;
    if( singlGluon >= 1000 && singlGluon < 2000 ) 
      ind = singlGluon - 1000;
    else if( singlGluon >= 2000 && singlGluon < 3000 ) 
      ind = singlGluon - 2000;
    else {
      cout << "Singlet gluon in hard process" << endl;
    }
    //cout << "singGluon " << singlGluon << endl;
  
    int lastCol;
    lastCol = max( MaxIndex(leftList,ind), MaxIndex(rightList,ind) );
    lastCol = max( lastCol, hardOld.maximum() );

    //cout << "Begin lastCol " << lastCol << endl;

    for(unsigned int i = ind; i < leftList.size(); ++i) {
      //Evolution

      Branch *bE;
      Branch *bEnew, *bCnew;

      leftListNew[i].colOut  = leftListNew[i-1].colIn;
      leftListNew[i].acolOut = leftListNew[i-1].acolIn;

      rightListNew[i].colOut  = rightListNew[i-1].colIn;
      rightListNew[i].acolOut = rightListNew[i-1].acolIn;

      if(leftList[i].idEm != -1 || rightList[i].idEm !=-1) {
        if(leftList[i].idEm != -1 && rightList[i].idEm ==-1) {
          bE = &leftList[i];
          bEnew = &leftListNew[i]; bCnew = &rightListNew[i];
        }
        else if(leftList[i].idEm == -1 && rightList[i].idEm != -1) {
          bE = &rightList[i];
          bEnew = &rightListNew[i]; bCnew = &leftListNew[i];
        }
        else{
          cout << "Emissions from both sides"<<endl;
          exit(1);
        }
        bool emCol = bE->colIn  == bE->colEm && bE->acolEm == bE->colOut
                   && bE->acolIn == bE->acolOut;

        bool emaCol = bE->acolIn == bE->acolEm && bE->colEm == bE->acolOut
                    && bE->colIn  == bE->colOut;

        if(i == ind && bE->colIn *bE->acolIn == 0 ) {
          //cout << "Swapping not possible (q -> g q)" << endl;
          return false;
        }


        //test of branching type (g -> g g - acol stays or q -> q g)
        if ((emCol && i > ind) || (emaCol && i == ind) ) {
          bEnew->acolIn = bEnew->acolOut;
          bEnew->acolEm = bEnew->colOut;
          if(bE->colEm != 0 ) {
            bEnew->colEm = ++lastCol;
            bEnew->colIn = bEnew->colEm;
          }
          else {
            bEnew->colEm = 0;
            bEnew->colIn = 0;
          }
        }
        //test of branching type (g -> g g - col stays or q -> q g)
        if ( (emaCol && i > ind) || (emCol && i == ind) ) {
          bEnew->colIn = bEnew->colOut;
          bEnew->colEm = bEnew->acolOut;
          if(bE->acolEm != 0 ) {
            bEnew->acolEm = ++lastCol;
            bEnew->acolIn = bEnew->acolEm;
          }
          else {
            bEnew->acolEm = 0;
            bEnew->acolIn = 0;
          }
        }

        //copy the second branch
        bCnew->colIn  = bCnew->colOut;
        bCnew->acolIn = bCnew->acolOut;
      }
      //copy of both branches
      else {
        leftListNew[i].colIn = leftListNew[i].colOut;
        leftListNew[i].acolIn = leftListNew[i].acolOut;

        rightListNew[i].colIn = rightListNew[i].colOut;
        rightListNew[i].acolIn = rightListNew[i].acolOut;

      }

    }

  int newSing = IsProblem(leftListNew, rightListNew, hardOld);
  if(newSing == 0) {
    leftList= leftListNew;
    rightList= rightListNew;
    return true;
  }
  else
    return false;

}


bool ColorReconection::ModifyColorsSingQQbarG( Pythia8::Event &event)
{

  int idStart1, idStart2;
  for(idStart2 = event.size() -1; event[idStart2].status() > 0; --idStart2)
    ;
  idStart1 = idStart2 -1;

  int idHard1 = max(idStart1, idStart2) + 1;
  int idHard2 = idHard1 + 1;
  
  if( event[idHard1].col() == event[idHard2].acol() &&
      event[idHard1].acol() == event[idHard2].col() &&
      event.size() - idHard1 == 5  ) {
    int idQ = -1, idQbar = -1, idG = -1; 
    for(int i = idHard2+1; i < event.size(); ++i) {
      if( event[i].id() >= 1 && event[i].id() <= 6 ) idQ    = i;
      if(-event[i].id() >= 1 &&-event[i].id() <= 6 ) idQbar = i;
      if( event[i].id() == 21 ) idG = i;
    }
    if( idQ != -1 && idQbar != -1 && idG != -1 ) {
      int colMax = max(
      max( event[idHard1].col(),event[idHard1].acol() ),
      max( event[idHard2].col(),event[idHard2].acol() ) );

      event[idQ].   cols(colMax+1, 0);
      event[idQbar].cols(0, colMax+2);
      event[idG].   cols(colMax+2,colMax+1);

      event[idStart1].cols(colMax+3,colMax+4);
      event[idStart2].cols(colMax+4,colMax+3);
      StoreInColors(event);
      return true;

    }
  }
  return false;

}











void ColorReconection::printSpaceLikeShower(const vector<Branch> &leftList, const vector<Branch> &rightList, const ColArr &colarr)
{
  cout << "Left start" << endl;
  for(int i = leftList.size()-1; i >= 0; --i) {
    leftList[i].print(); 
  }
  cout << "Right start" << endl;
  for(int i = leftList.size()-1; i >= 0; --i) {
    rightList[i].print(); 
  }
  cout << "Col Arr" << endl;
  colarr.print(); 

}
